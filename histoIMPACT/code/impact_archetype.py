"""
histoIMPACT spatial archetyping module

Steps
  1. build_superpixels   
  2. signature_maps 
  3. moran_features
  4. predict_impact
"""

from __future__ import annotations

import json
from typing import Any

import numpy as np
from scipy import ndimage

LINEAGES = ['Bronchus', 'Fibrous_tissue_enriched', 'Lymph_enriched',
            'Myeloid_enriched', 'Alveoli', 'SMC & Vessel', 'Tumor']
LAYERS = ['Tumor_region', 'Interface_tumor', 'Interface_normal',
          'Normal_peri', 'Normal_distal']
LAYER_SHORT = ['TR', 'IT', 'IN', 'NP', 'ND']
LAYER_COLORS = {'TR': '#d73027', 'IT': '#fee090', 'IN': '#abd9e9',
                'NP': '#74add1', 'ND': '#4575b4'}
SIGNATURES = ['CD4_T', 'CD8_T', 'B', 'Mac', 'Fib']
THREE = {'IMPACT1a': 'IMPACT1a/2', 'IMPACT2': 'IMPACT1a/2'}

# frozen settings (see Methods)
SETTINGS = dict(
    target_um=50.0,            # super-pixel side (realised 48 um at 4 um/token)
    min_tissue_frac=0.05,      # super-pixels with < 5% tissue are dropped
    tumor_value=0,             # tumor-region model: 0 = tumor_region, 1 = non_tumor_region
    tumor_smooth_um=150.0,     # majority vote window (3 x 3 super-pixels)
    region_min_area_mm2=0.05,  # tumor islands / non-tumor holes removed
    region_min_rel=0.1,        # islands < 10% of the largest tumor removed ...
    region_keep_mm2=0.5,       # ... unless >= 0.5 mm2
    ring_edges_um=(-100.0, 0.0, 100.0, 300.0),   # five Visium layers
    airspace_close_um=100.0,
    spot_um=100.0,             # pseudo-spot (2 x 2 super-pixels, 96 um)
    window_mm=11.0,            # Visium capture-area scale
    window_step_mm=5.5,
    window_min_tissue_mm2=20.0,
    min_tumor_frac=0.05,       # tumor layers < 5% of tissue -> not assessable
)


# ---------------------------------------------------------------------
# 1. super-pixels, tumor region, five layers
# ---------------------------------------------------------------------

def block_sum(a: np.ndarray, f: int) -> np.ndarray:
    h, w = a.shape
    H, W = -(-h // f) * f, -(-w // f) * f
    a = np.pad(a, ((0, H - h), (0, W - w)))
    return a.reshape(H // f, f, W // f, f).sum((1, 3))


def _remove_small(region, valid, min_size):
    """Flip connected components (of region and of its complement within
    valid) smaller than min_size super-pixels."""
    region = region.copy()
    for target in (True, False):
        comp = (region == target) & valid
        lab, n = ndimage.label(comp)
        if n == 0:
            continue
        sizes = np.bincount(lab.ravel(), minlength=n + 1)
        small = (sizes < min_size) & (np.arange(n + 1) > 0)
        region[small[lab]] = not target
    return region & valid


def tumor_region(tumor_frac, n_tis, valid, um_sp, smooth_um, min_area_mm2,
                 min_rel, keep_mm2):
    """Tissue-weighted majority vote of token tumor calls (smooth_um window),
    small islands / holes removed, small distant islands removed."""
    size = max(1, int(round(smooth_um / um_sp)))
    num = ndimage.uniform_filter(tumor_frac * n_tis, size)
    den = ndimage.uniform_filter(n_tis, size)
    region = (num > 0.5 * den) & valid
    min_size = int(np.ceil(min_area_mm2 * 1e6 / um_sp ** 2))
    if min_size > 1:
        region = _remove_small(region, valid, min_size)
    if min_rel > 0:
        lab, n = ndimage.label(region)
        if n > 1:
            sizes = np.bincount(lab.ravel(), minlength=n + 1)
            sizes[0] = 0
            small = sizes < min_rel * sizes.max()
            if keep_mm2 > 0:
                small &= sizes < keep_mm2 * 1e6 / um_sp ** 2
            small[0] = False
            region = region & ~small[lab]
    return region


def signed_distance_um(t, valid, um):
    """+ distance to the tumor outside it, - distance to non-tumor tissue
    inside it (um)."""
    if not t.any():
        return np.full(t.shape, np.inf, dtype=np.float32)
    d_out = ndimage.distance_transform_edt(~t)
    nontumor = valid & ~t
    d_in = (ndimage.distance_transform_edt(~nontumor) if nontumor.any()
            else np.full(t.shape, np.inf))
    return (np.where(t, -d_in, d_out) * um).astype(np.float32)


def airspace_fraction(n_tis, n_pix, um_sp, close_um):
    r = max(1, int(round(close_um / um_sp)))
    occ = np.pad(n_tis > 0, r)
    env = ndimage.binary_fill_holes(
        ndimage.binary_closing(occ, iterations=r))[r:-r, r:-r]
    return np.where(env, 1 - n_tis / np.maximum(n_pix, 1), 0).astype(
        np.float32)


def build_superpixels(lineage: np.ndarray, tumor: np.ndarray,
                      um_per_token: float, s: dict = SETTINGS) -> dict:
    """lineage / tumor: token grids (-1 = background). Returns the
    super-pixel tensor: lineage fractions (7) | layer one-hot (5) | raw tumor
    call fraction | air-space fraction, plus mask, tissue counts, layer
    index and geometry."""
    lineage, tumor = np.asarray(lineage), np.asarray(tumor)
    if lineage.shape != tumor.shape:
        raise ValueError(f'lineage {lineage.shape} vs tumor {tumor.shape}')
    tissue = lineage >= 0
    f = max(1, int(round(s['target_um'] / um_per_token)))
    um_sp = f * um_per_token
    n_tis = block_sum(tissue.astype(np.float32), f)
    n_pix = block_sum(np.ones(tissue.shape, np.float32), f)
    denom = np.maximum(n_tis, 1.0)
    mask = n_tis / n_pix >= s['min_tissue_frac']
    lin = [block_sum((lineage == k).astype(np.float32), f) / denom
           for k in range(len(LINEAGES))]
    tfrac = block_sum(((tumor == s['tumor_value']) & tissue)
                      .astype(np.float32), f) / denom
    t = tumor_region(tfrac, n_tis * mask, mask, um_sp, s['tumor_smooth_um'],
                     s['region_min_area_mm2'], s['region_min_rel'],
                     s['region_keep_mm2'])
    ring = np.digitize(signed_distance_um(t, mask, um_sp), s['ring_edges_um'])
    rings = [(ring == r).astype(np.float32) for r in range(len(LAYERS))]
    air = airspace_fraction(n_tis, n_pix, um_sp, s['airspace_close_um'])
    # float16 round trip: identical to the cache the model was built with
    x = np.stack(lin + rings + [tfrac, air]).astype(np.float16)
    x[:, ~mask] = 0
    x = x.astype(np.float32)
    layer = np.where(mask, ring, -1).astype(np.int8)
    return dict(x=x, mask=mask, tissue_count=n_tis.astype(np.float32),
                pool_factor=f, um_per_superpixel=um_sp, layer=layer)


def token_layer_map(sp: dict, lineage: np.ndarray) -> np.ndarray:
    """layer index per token (each tissue token takes its super-pixel's
    layer; -1 = background / outside the tissue mask)"""
    f = int(sp['pool_factor'])
    H, W = lineage.shape
    up = np.repeat(np.repeat(sp['layer'], f, 0), f, 1)[:H, :W]
    return np.where(lineage >= 0, up, -1).astype(np.int8)


# ---------------------------------------------------------------------
# 2. Visium-supervised signature maps
# ---------------------------------------------------------------------

def load_signature_model(path) -> dict:
    z = np.load(path, allow_pickle=True)
    m = {k: z[k] for k in z.files}
    m['opts'] = json.loads(str(z['opts']))
    if m['opts'].get('context_blocks'):
        raise NotImplementedError('context features are not supported here')
    m['names'] = [str(s) for s in z['scores']]
    m['keep'] = ([int(k) for k in z['main']] if 'main' in z.files
                 else list(range(len(m['names']))))
    return m


def _gelu(x):
    return 0.5 * x * (1 + np.tanh(0.7978845608 * (x + 0.044715 * x ** 3)))


def _components(m):
    return [{key.split('__', 1)[1]: v for key, v in m.items()
             if key.startswith(f'c{k}__')} for k in range(int(m['n_comp']))]


def _apply(m, X):
    kind = str(m['kind'])
    if kind == 'ensemble':
        return np.mean([_apply(c, X) for c in _components(m)], 0)
    if kind == 'ridge':
        return X @ m['B'] + m['b0']
    h = _gelu((X - m['mu']) / m['sd'] @ m['W0'] + m['b0_'])
    h = _gelu(h @ m['W1'] + m['b1_'])
    return h @ m['W2'] + m['b2_']


def block_embeddings(feat, lineage, f, want_sq=False, band_blocks=8):
    """Per f x f block: sum (and squared sum) of the UNI embeddings of tissue
    tokens, cell-type counts and tissue counts. feat: (D, H, W) array or
    np.load(..., mmap_mode='r'); read in row bands to bound memory."""
    D, H, W = feat.shape
    if lineage.shape != (H, W):
        raise ValueError(f'tokens {feat.shape} vs lineage {lineage.shape}')
    tissue = lineage >= 0
    h, w = -(-H // f), -(-W // f)
    s = np.zeros((h, w, D), np.float32)
    sq = np.zeros((h, w, D), np.float32) if want_sq else None
    lin = np.zeros((h, w, len(LINEAGES)), np.float32)
    cnt = block_sum(tissue.astype(np.float32), f)
    onehot = np.eye(len(LINEAGES), dtype=np.float32)
    band = band_blocks * f
    for r0 in range(0, H, band):
        r1 = min(r0 + band, H)
        T = tissue[r0:r1]
        if not T.any():
            continue
        cols = np.where(T.any(0))[0]
        c0, c1 = cols[0] // f * f, min(W, -(-(cols[-1] + 1) // f) * f)
        X = np.asarray(feat[:, r0:r1, c0:c1], dtype=np.float32)
        Tw = T[:, c0:c1].astype(np.float32)
        nb_r, nb_c = -(-(r1 - r0) // f), -(-(c1 - c0) // f)
        pad = ((0, nb_r * f - (r1 - r0)), (0, nb_c * f - (c1 - c0)))

        def bsum(A):          # (C, rows, cols) -> (nb_r, nb_c, C)
            A = np.pad(A, ((0, 0),) + pad)
            return A.reshape(A.shape[0], nb_r, f, nb_c, f).sum((2, 4)) \
                .transpose(1, 2, 0)
        k0, j0 = r0 // f, c0 // f
        s[k0:k0 + nb_r, j0:j0 + nb_c] = bsum(X * Tw[None])
        if want_sq:
            sq[k0:k0 + nb_r, j0:j0 + nb_c] = bsum(X * X * Tw[None])
        L = np.clip(lineage[r0:r1, c0:c1], 0, len(LINEAGES) - 1)
        oh = onehot[L].transpose(2, 0, 1) * Tw[None]
        lin[k0:k0 + nb_r, j0:j0 + nb_c] = bsum(oh)
    return dict(sum=s, sq=sq, lin=lin, cnt=cnt)


def signature_maps(feat, lineage, sp: dict, model: dict):
    """-> maps (n_signatures, h, w) on the super-pixel grid, names"""
    f = int(sp['pool_factor'])
    opts = model['opts']
    bl = block_embeddings(feat, lineage, f, want_sq=opts.get('std', False))
    cnt = np.maximum(bl['cnt'], 1)[..., None]
    mean = bl['sum'] / cnt
    parts = [mean]
    if opts.get('std'):
        parts.append(np.sqrt(np.maximum(bl['sq'] / cnt - mean ** 2, 0)))
    if opts.get('lineage'):
        parts.append(bl['lin'] / cnt)
    F = np.concatenate(parts, -1)
    wgt = bl['cnt'] / max(bl['cnt'].sum(), 1)
    center = (F * wgt[..., None]).sum((0, 1))       # slide centring
    h, w, D = F.shape
    if (h, w) != sp['mask'].shape:
        raise ValueError(f'block grid {(h, w)} vs super-pixels '
                         f'{sp["mask"].shape}')
    X = (F - center).reshape(-1, D)
    out = np.concatenate([_apply(model, X[i:i + 65536])
                          for i in range(0, len(X), 65536)])
    maps = out.reshape(h, w, -1)[..., model['keep']].transpose(2, 0, 1)
    maps[:, bl['cnt'] == 0] = 0
    # stored as float16 in the original pipeline; keep the same rounding
    return (maps.astype(np.float16).astype(np.float32),
            [model['names'][k] for k in model['keep']])


# ---------------------------------------------------------------------
# 3. Moran's I features
# ---------------------------------------------------------------------

QUEEN = np.array([[1, 1, 1], [1, 0, 1], [1, 1, 1]], np.float32)


def _lag(z, valid):
    num = ndimage.convolve(np.where(valid, z, 0.0), QUEEN, mode='constant')
    cnt = ndimage.convolve(valid.astype(np.float32), QUEEN, mode='constant')
    return num / np.maximum(cnt, 1), valid & (cnt > 0)


def bivariate_moran(a, b, valid):
    """Mean of I(a, lag b) and I(b, lag a), row-standardised queen weights
    over valid pseudo-spots; NaN if either variable is constant."""
    za, zb = a[valid], b[valid]
    if za.std() == 0 or zb.std() == 0:
        return np.nan
    A = np.zeros_like(a, dtype=np.float64)
    B = np.zeros_like(b, dtype=np.float64)
    A[valid] = (za - za.mean()) / za.std()
    B[valid] = (zb - zb.mean()) / zb.std()
    lag_b, use = _lag(B, valid)
    lag_a, __ = _lag(A, valid)
    return 0.5 * (np.mean(A[use] * lag_b[use]) + np.mean(B[use] * lag_a[use]))


def _pseudo_spots(sp, spot_um):
    f = max(1, int(round(spot_um / sp['um_per_superpixel'])))
    w = sp['tissue_count'] * sp['mask']
    wsum = block_sum(w, f)
    valid = wsum > 0
    lay = np.stack([block_sum(sp['x'][len(LINEAGES) + r] * w, f)
                    for r in range(len(LAYERS))])
    layer = (lay / np.maximum(wsum, 1e-6)).argmax(0)
    return layer, valid, wsum, f, w


def _features(sp, maps, names, spot_um):
    layer, valid, wsum, f, w = _pseudo_spots(sp, spot_um)
    tot = wsum[valid].sum()
    feats = {f'area|{ln}': (wsum[valid & (layer == r)].sum() / tot
                            if tot > 0 else np.nan)
             for r, ln in enumerate(LAYERS)}
    for k, name in enumerate(names):
        v = block_sum(maps[k].astype(np.float64) * w, f) / np.maximum(wsum,
                                                                     1e-6)
        for r, ln in enumerate(LAYERS):
            feats[f'moran|x{name}|{ln}'] = bivariate_moran(
                v, (layer == r).astype(np.float64), valid)
    return feats


def _crop(sp, maps, i0, i1, j0, j1):
    sub = dict(sp, x=sp['x'][:, i0:i1, j0:j1], mask=sp['mask'][i0:i1, j0:j1],
               tissue_count=sp['tissue_count'][i0:i1, j0:j1])
    return sub, maps[:, i0:i1, j0:j1]


def moran_features(sp: dict, maps, names, s: dict = SETTINGS) -> dict:
    """Features on window_mm windows (tissue >= window_min_tissue_mm2),
    tissue-weighted mean over windows; a slide no larger than one window is
    one window (= whole-slide features, as for the Visium samples)."""
    um = sp['um_per_superpixel']
    wn = max(1, int(round(s['window_mm'] * 1000 / um)))
    st = max(1, int(round(s['window_step_mm'] * 1000 / um)))
    H, W = sp['mask'].shape
    px_mm2 = (um / 1000) ** 2
    starts = lambda n: sorted(set(list(range(0, max(n - wn, 0) + 1, st))
                                  + [max(n - wn, 0)]))
    rows, wts = [], []
    for i0 in starts(H):
        for j0 in starts(W):
            sub, sm = _crop(sp, maps, i0, i0 + wn, j0, j0 + wn)
            tissue = float(sub['mask'].sum()) * px_mm2
            if tissue < s['window_min_tissue_mm2']:
                continue
            rows.append(_features(sub, sm, names, s['spot_um']))
            wts.append(tissue)
    if not rows:
        rows = [_features(sp, maps, names, s['spot_um'])]
        wts = [float(sp['mask'].sum()) * px_mm2]
    wts = np.array(wts)
    out = {}
    for k in rows[0]:
        v = np.array([r[k] for r in rows], float)
        ok = np.isfinite(v)
        out[k] = float((v[ok] * wts[ok]).sum() / wts[ok].sum()) if ok.any() \
            else np.nan
    out['n_windows'] = len(rows)
    out['tissue_mm2'] = float(sp['mask'].sum()) * px_mm2
    return out


# ---------------------------------------------------------------------
# 4. IMPACT classification
# ---------------------------------------------------------------------

def impact_rule(til, mac, fib):
    if til == 'III' and mac == 'I' and fib == 'I':
        return 'IMPACT1a'
    if til == 'II' and mac == 'I' and fib == 'II':
        return 'IMPACT1b'
    if til == 'IV':
        return 'IMPACT3'
    return 'IMPACT2'


def _zscore(x, sc):
    med, mu, sd = (np.array(sc[k], float) for k in ('med', 'mu', 'sd'))
    return (np.where(np.isnan(x), med, x) - mu) / sd


def predict_impact(feats: dict, model: dict,
                   min_tumor_frac: float = SETTINGS['min_tumor_frac']) -> dict:
    """feats: one slide's features; model: impact_classifier.json.
    H&E features are z-scored with the discovery H&E statistics and read by
    classifiers trained on the Visium Moran's I (z-score calibration)."""
    cal = model.get('calibration')
    Vh = None
    if cal:
        x = np.array([feats.get(c, np.nan) for c in cal['he_columns']], float)
        Hz = np.nan_to_num(_zscore(x, cal['he_scaler']))
        Vh = dict(zip(cal['visium_columns'],
                      Hz @ np.array(cal['W']) + np.array(cal['c'])))
    probs, n_missing = {}, 0
    for t, task in model['tasks'].items():
        hx = np.array([feats.get(c, np.nan) for c in task['he_columns']],
                      float)
        n_missing += int(np.isnan(hx).sum())
        if Vh is None:
            z = _zscore(hx, task['he_scaler'])
        else:
            z = _zscore(np.array([Vh[c] for c in task['visium_columns']]),
                        task['visium_scaler'])
        W, b = np.array(task['coef']), np.array(task['intercept'])
        lg = W @ z + b
        if W.shape[0] == 1:
            p1 = 1 / (1 + np.exp(-lg[0]))
            p = np.array([1 - p1, p1])
        else:
            e = np.exp(lg - lg.max())
            p = e / e.sum()
        probs[t] = dict(zip(task['classes'], p.tolist()))
    hard = {t: max(p, key=p.get) for t, p in probs.items()}
    imp = impact_rule(hard['til_class'], hard['mac_class'], hard['fib_class'])
    til_sorted = sorted(probs['til_class'].values(), reverse=True)
    tumor = feats.get('area|Tumor_region', 0) + feats.get(
        'area|Interface_tumor', 0)
    tumor = float(np.nan_to_num(tumor))
    assessable = tumor >= min_tumor_frac
    out = dict(
        IMPACT_3class=THREE.get(imp, imp) if assessable else 'not_assessable',
        til_class=hard['til_class'], mac_class=hard['mac_class'],
        fib_class=hard['fib_class'],
        til_margin=til_sorted[0] - til_sorted[1],
        tumor_layer_frac=tumor, n_missing_features=n_missing,
        **{f'p_{t}_{c}': v for t, p in probs.items() for c, v in p.items()})
    return out


def load_impact_classifier(path) -> dict[str, Any]:
    with open(path) as fh:
        m = json.load(fh)
    if m.get('decision', 'hard') != 'hard':
        raise ValueError('only the hard decision (argmax + rules) is '
                         'implemented')
    return m
