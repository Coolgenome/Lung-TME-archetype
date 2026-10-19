"""
histoIMPACT end to end: H&E -> model 1 (cell lineage + tumor region) ->
five tissue layers -> IMPACT archetype (3-class: IMPACT1a/2, IMPACT1b,
IMPACT3).

Per sample (out_root/<sample_id>/):
  cell_lineage_grid.npy / cell_lineage_map.png   spatial tissue map (tokens)
  tumor_region_grid.npy / tumor_region_map.png   token-level tumor calls
  tissue_layers_token.npy / tissue_layers.png    5-layer segmentation
  tissue_layers_superpixel.npy                   layers on the 48 um grid
  signature_maps.npz (+ signature_maps.png)      CD4_T/CD8_T/B/Mac/Fib maps
  features.json, impact.json                     Moran's I features, call
Cohort: out_root/histoimpact_predictions.csv (one row per sample).

Inputs
  A. images (PNG / TIFF / JPEG, 40x preferred):
       python run_histoimpact.py --image ../HE_demo/ADC1.png --mpp 0.25
     Images at another resolution are resampled to 0.25 um/px first
     (the models were trained at 0.25 um/px; 2x upsampling of 20x images
     reproduced 40x predictions closely, native 20x tokens did not).
  B. a manifest CSV with column sample_id and either
       image [+ mpp]                       -> full pipeline, or
       tokens + lineage_grid + tumor_grid [+ um_per_token]
                                           -> archetyping only (re-use
                                              existing UNI tokens and
                                              model-1 grids, e.g. whole
                                              slides processed elsewhere)
       python run_histoimpact.py --manifest samples.csv

Model files (default: ../model/):
  general_uni_cls_cell_lineage.pt, general_uni_cls_tumor_region.pt (model 1)
  signature_model.npz      Visium-supervised signature MLP (expr/model_mt.npz)
  impact_classifier.json   compartment classifiers (final_model.json)
UNI weights: ../UNI/pytorch_model.bin
"""

from __future__ import annotations

import argparse
import json
import os
import pickle
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

import impact_archetype as ia

# the model-1 checkpoints are plain dicts saved by us; torch.load's
# weights_only notice is not relevant here
warnings.filterwarnings('ignore', category=FutureWarning,
                        message='You are using `torch.load`')

HERE = Path(__file__).resolve().parent
MODEL_DIR = HERE.parent / 'model'
UNI_WEIGHTS = HERE.parent / 'UNI' / 'pytorch_model.bin'
TARGET_MPP = 0.25
TOKEN_PX = 16

LINEAGE_PALETTE = {
    'Tumor': '#ef3b2c', 'Alveoli': '#e7e1ee',
    'Lymph_enriched': '#23b9dc', 'Myeloid_enriched': '#a8d36a',
    'Fibrous_tissue_enriched': '#ffd700', 'SMC & Vessel': '#76cbd6',
    'Bronchus': '#e7298a'}
TUMOR_PALETTE = {'tumor_region': '#d73027', 'non_tumor_region': '#4393c3'}


# ---------------------------------------------------------------------
# model 1
# ---------------------------------------------------------------------

def read_image(path: Path, mpp: float | None, resample: bool):
    from PIL import Image
    Image.MAX_IMAGE_PIXELS = None
    im = Image.open(path).convert('RGB')
    if mpp is None:
        print(f'  [{path.name}] no mpp given: assuming {TARGET_MPP} um/px (40x)')
        mpp = TARGET_MPP
    if resample and abs(mpp / TARGET_MPP - 1) > 0.1:
        s = mpp / TARGET_MPP
        size = (int(round(im.width * s)), int(round(im.height * s)))
        print(f'  [{path.name}] resampling {mpp:g} -> {TARGET_MPP} um/px '
              f'({im.width}x{im.height} -> {size[0]}x{size[1]})')
        im = im.resize(size, Image.LANCZOS)
        mpp = TARGET_MPP
    return np.asarray(im), mpp


def run_model1(img, uni, device, out_dir: Path, save_tokens: bool):
    """UNI tokens, tissue mask, lineage and tumor-region grids (he_pred)"""
    from he_pred import (encode_he_with_uni_highres_from_img,
                         he_background_mask, load_ckpt_classifier,
                         predict_on_feat_array_class, resize_bool_pil)
    print('  UNI features')
    feat, meta = encode_he_with_uni_highres_from_img(
        img=img, uni=uni, patch_size=224, stride=112, batch_size=64,
        device=device, normalize=True, roi_mask=None, blend_alpha=0.5,
        out_dtype='float16')
    if save_tokens:
        np.save(out_dir / 'HE_token.npy', feat)
        with open(out_dir / 'HE_token_meta.pkl', 'wb') as fh:
            pickle.dump(meta, fh, protocol=pickle.HIGHEST_PROTOCOL)
    print('  tissue mask')
    bg = he_background_mask(img, patch_size=16, density_thresh=100,
                            min_size=10, glcm_levels=64,
                            clean_background=True, texture_output_dir=None,
                            mask_output_dir=None)
    bg_full = bg['mask_fullres'][:img.shape[0], :img.shape[1]]
    grids = {}
    for name, ckpt in (('lineage', 'general_uni_cls_cell_lineage.pt'),
                       ('tumor', 'general_uni_cls_tumor_region.pt')):
        print(f'  model 1: {name}')
        clf, m = load_ckpt_classifier(MODEL_DIR / ckpt, map_location='cpu')
        labels = [str(x) for x in m['label_order']]
        g = predict_on_feat_array_class(
            feat, clf, layout=str(m.get('feat_layout', 'DHW')),
            device=device).astype(np.int16)
        g[~resize_bool_pil(bg_full, g.shape)] = -1
        grids[name] = (g, labels)
    return feat, grids


def check_model1():
    """Class names in the model-1 checkpoints must match the archetyping
    code (checked before the slow UNI step)."""
    from he_pred import load_ckpt_classifier
    __, m = load_ckpt_classifier(MODEL_DIR / 'general_uni_cls_cell_lineage.pt')
    lin = [str(x) for x in m['label_order']]
    if lin != ia.LINEAGES:
        raise SystemExit(f'cell-lineage checkpoint classes {lin}\n'
                         f'do not match impact_archetype.LINEAGES {ia.LINEAGES}:'
                         ' update code/ and model/ together')
    __, m = load_ckpt_classifier(MODEL_DIR / 'general_uni_cls_tumor_region.pt')
    tum = [str(x) for x in m['label_order']]
    if tum[0] != 'tumor_region':
        raise SystemExit(f'tumor-region checkpoint classes {tum}: expected '
                         'class 0 = tumor_region')


def load_uni(device):
    import timm
    import torch
    model = timm.create_model('vit_large_patch16_224', img_size=224,
                              patch_size=16, init_values=1e-5, num_classes=0,
                              dynamic_img_size=True)
    model.load_state_dict(torch.load(str(UNI_WEIGHTS), map_location='cpu',
                                     weights_only=True), strict=True)
    return model.eval()


def choose_device(req):
    if req != 'auto':
        return req
    import torch
    if torch.cuda.is_available():
        return 'cuda:0'
    if getattr(torch.backends, 'mps', None) and torch.backends.mps.is_available():
        return 'mps'
    return 'cpu'


# ---------------------------------------------------------------------
# figures
# ---------------------------------------------------------------------

def _style():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.family': 'sans-serif',
                         'font.sans-serif': ['Arial', 'Helvetica',
                                             'DejaVu Sans'],
                         'font.size': 8, 'pdf.fonttype': 42,
                         'svg.fonttype': 'none'})
    return plt


def _rgb(labels, colors, bg, max_px=3000):
    """categorical grid -> RGB (block-averaged if larger than max_px)"""
    from matplotlib.colors import to_rgb
    lut = np.array([to_rgb(c) for c in colors] + [to_rgb(bg)], np.float32)
    rgb = lut[np.where(labels >= 0, labels, len(colors))]
    k = int(np.ceil(max(labels.shape) / max_px))
    if k > 1:
        H, W = (labels.shape[0] // k) * k, (labels.shape[1] // k) * k
        rgb = rgb[:H, :W].reshape(H // k, k, W // k, k, 3).mean((1, 3))
    return rgb, max(k, 1)


def plot_categorical(grid, names, colors, path, title, um, bg='#000000',
                     legend_cols=1):
    plt = _style()
    from matplotlib.patches import Patch
    import matplotlib.transforms as mtransforms
    rgb, k = _rgb(grid, colors, bg)
    h, w = rgb.shape[:2]
    fig, ax = plt.subplots(figsize=(6.5, min(6.5 * h / w, 9) + 0.8))
    ax.imshow(rgb, interpolation='nearest')
    ax.set_axis_off()
    ax.set_title(title, loc='left', fontsize=8)
    # scale bar below the image
    wmm = w * k * um / 1000
    bar = 1 if wmm < 8 else (2 if wmm < 25 else 5)
    n = bar * 1000 / (um * k)
    tr = mtransforms.blended_transform_factory(ax.transData, ax.transAxes)
    ax.plot([0, n], [-0.03, -0.03], color='#0b0b0b', lw=1.2, transform=tr,
            clip_on=False, solid_capstyle='butt')
    ax.text(n + 0.01 * w, -0.03, f'{bar} mm', va='center', fontsize=7,
            transform=tr, clip_on=False)
    ax.legend(handles=[Patch(facecolor=c, edgecolor='none', label=nm)
                       for nm, c in zip(names, colors)],
              loc='center left', bbox_to_anchor=(1.01, 0.5), frameon=False,
              ncol=legend_cols)
    fig.savefig(path, dpi=int(np.clip(max(h, w) / 5, 200, 600)),
                bbox_inches='tight', facecolor='white')
    plt.close(fig)


def plot_maps(maps, names, sp, path, title):
    plt = _style()
    n = len(names)
    fig, axes = plt.subplots(1, n, figsize=(2.6 * n, 3.0))
    mask = sp['mask']
    for ax, m, nm in zip(np.atleast_1d(axes), maps, names):
        v = np.where(mask, m, np.nan)
        lo, hi = np.nanpercentile(v, [2, 98]) if mask.any() else (0, 1)
        ax.imshow(v, cmap='turbo', vmin=lo, vmax=hi if hi > lo else lo + 1e-6,
                  interpolation='nearest')
        ax.set_title(nm.replace('_', ' '), fontsize=8)
        ax.set_axis_off()
    fig.suptitle(title + '  (relative within slide; 2nd-98th percentile)',
                 fontsize=8, x=0.01, ha='left')
    fig.tight_layout()
    fig.savefig(path, dpi=200, facecolor='white')
    plt.close(fig)


# ---------------------------------------------------------------------
# one sample
# ---------------------------------------------------------------------

def process(sample: dict, models: dict, out_root: Path, args, uni=None,
            device='cpu') -> dict:
    sid = str(sample['sample_id'])
    out = out_root / sid
    out.mkdir(parents=True, exist_ok=True)
    print(f'[{sid}]')
    image = sample.get('image')
    if isinstance(image, str) and image:
        mpp = sample.get('mpp')
        mpp = float(mpp) if mpp is not None and pd.notna(mpp) else args.mpp
        img, mpp = read_image(Path(image), mpp, not args.no_resample)
        feat, grids = run_model1(img, uni, device, out, args.save_tokens)
        del img
        lineage, lin_labels = grids['lineage']
        tumor, tum_labels = grids['tumor']
        um = mpp * TOKEN_PX
    else:
        feat = np.load(sample['tokens'], mmap_mode='r')
        lineage = np.load(sample['lineage_grid'])
        tumor = np.load(sample['tumor_grid'])
        lin_labels, tum_labels = ia.LINEAGES, ['tumor_region', 'non_tumor_region']
        u = sample.get('um_per_token')
        um = float(u) if u is not None and pd.notna(u) else args.mpp * TOKEN_PX
        if feat.ndim == 3 and feat.shape[1:] != lineage.shape \
                and feat.shape[:2] == lineage.shape:
            feat = np.moveaxis(feat, 2, 0)          # (H, W, D) -> (D, H, W)
    if abs(um / (TARGET_MPP * TOKEN_PX) - 1) > 0.1:
        print(f'  WARNING: {um:g} um/token; the models were trained at '
              f'{TARGET_MPP * TOKEN_PX:g} um/token (40x)')
    np.save(out / 'cell_lineage_grid.npy', lineage.astype(np.int8))
    np.save(out / 'tumor_region_grid.npy', tumor.astype(np.int8))

    # archetyping
    sp = ia.build_superpixels(lineage, tumor, um)
    lay_tok = ia.token_layer_map(sp, lineage)
    np.save(out / 'tissue_layers_superpixel.npy', sp['layer'])
    np.save(out / 'tissue_layers_token.npy', lay_tok)
    maps, names = ia.signature_maps(feat, lineage, sp, models['signature'])
    np.savez_compressed(out / 'signature_maps.npz',
                        maps=maps.astype(np.float16), names=np.array(names),
                        um_per_superpixel=sp['um_per_superpixel'])
    feats = ia.moran_features(sp, maps, names)
    res = ia.predict_impact(feats, models['classifier'])
    with open(out / 'features.json', 'w') as fh:
        json.dump(feats, fh, indent=1)
    with open(out / 'impact.json', 'w') as fh:
        json.dump(dict(sample_id=sid, um_per_token=um, **res), fh, indent=1)

    if not args.no_plots:
        tit = (f'{sid} · {res["IMPACT_3class"]} '
               f'(tumor layers {res["tumor_layer_frac"]:.0%})')
        plot_categorical(lineage, lin_labels,
                         [LINEAGE_PALETTE.get(n, '#999999') for n in lin_labels],
                         out / 'cell_lineage_map.png', f'{sid} · cell lineage', um)
        plot_categorical(tumor, tum_labels,
                         [TUMOR_PALETTE.get(n, '#999999') for n in tum_labels],
                         out / 'tumor_region_map.png',
                         f'{sid} · tumor region (token calls)', um)
        plot_categorical(lay_tok, [f'{l} ({s})' for l, s in
                                   zip(ia.LAYERS, ia.LAYER_SHORT)],
                         [ia.LAYER_COLORS[s] for s in ia.LAYER_SHORT],
                         out / 'tissue_layers.png', tit, um, bg='#ffffff')
        if args.plot_maps:
            plot_maps(maps, names, sp, out / 'signature_maps.png', sid)
    print(f'  -> {res["IMPACT_3class"]} (TIL {res["til_class"]}, Mac {res["mac_class"]}, '
          f'Fib {res["fib_class"]}; tumor layers '
          f'{res["tumor_layer_frac"]:.1%}, {feats["n_windows"]} window(s))')
    return dict(sample_id=sid, um_per_token=um, n_windows=feats['n_windows'],
                tissue_mm2=round(feats['tissue_mm2'], 2), **res)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--image', action='append', type=Path, default=None)
    p.add_argument('--manifest', type=Path, default=None)
    p.add_argument('--mpp', type=float, default=None,
                   help='um per pixel of --image inputs (default 0.25); '
                   'manifest column mpp overrides')
    p.add_argument('--no-resample', action='store_true',
                   help='do not resample images to 0.25 um/px')
    p.add_argument('--out-root', type=Path, default=HERE.parent / 'results_histoimpact')
    p.add_argument('--signature-model', type=Path,
                   default=MODEL_DIR / 'signature_model.npz')
    p.add_argument('--classifier', type=Path,
                   default=MODEL_DIR / 'impact_classifier.json')
    p.add_argument('--device', default='auto')
    p.add_argument('--save-tokens', action='store_true',
                   help='keep the UNI token array (large)')
    p.add_argument('--plot-maps', action='store_true',
                   help='also draw the five signature maps')
    p.add_argument('--no-plots', action='store_true')
    args = p.parse_args()

    if args.manifest:
        man = pd.read_csv(args.manifest, dtype={'sample_id': str})
        samples = man.to_dict('records')
    elif args.image:
        samples = [dict(sample_id=im.stem, image=str(im), mpp=args.mpp)
                   for im in args.image]
    else:
        p.error('give --image or --manifest')
    if args.mpp is None:
        args.mpp = TARGET_MPP
    models = dict(signature=ia.load_signature_model(args.signature_model),
                  classifier=ia.load_impact_classifier(args.classifier))
    need_uni = any(isinstance(s.get('image'), str) and s.get('image')
                   for s in samples)
    uni, device = None, 'cpu'
    if need_uni:
        check_model1()
        device = choose_device(args.device)
        print('device:', device, '| loading UNI:', UNI_WEIGHTS)
        uni = load_uni(device)
    args.out_root.mkdir(parents=True, exist_ok=True)
    rows = []
    for s in samples:
        try:
            rows.append(process(s, models, args.out_root, args, uni, device))
        except Exception as e:  # keep going with the next sample
            print(f'[{s.get("sample_id")}] FAILED: {e}')
            rows.append(dict(sample_id=s.get('sample_id'),
                             IMPACT_3class='failed', error=str(e)))
    df = pd.DataFrame(rows)
    first = ['sample_id', 'IMPACT_3class', 'til_class', 'mac_class',
             'fib_class', 'til_margin', 'tumor_layer_frac', 'tissue_mm2',
             'n_windows', 'um_per_token', 'n_missing_features']
    df = df[[c for c in first if c in df] + [c for c in df if c not in first]]
    path = args.out_root / 'histoimpact_predictions.csv'
    df.round(4).to_csv(path, index=False)
    print(f'\n{len(df)} samples -> {path}')
    print(df['IMPACT_3class'].value_counts().to_string())


if __name__ == '__main__':
    main()
