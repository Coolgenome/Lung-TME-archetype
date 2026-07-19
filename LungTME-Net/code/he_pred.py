"""
Include H&E UNI feature extraction, tissue/background masking, prediction, and categorical rendering

For running the prediction, you need to provide:
  - `uni`: your loaded UNI / UNI-like torch model for feature extraction.
  - `clf`: your trained torch classifier for prediction.

"""

from __future__ import annotations

import os
import pickle
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence, Tuple

import numpy as np
import torch
from PIL import Image


# ---------------------------------------------------------------------
# H&E -> dense UNI token grid
# ---------------------------------------------------------------------


def _as_rgb_uint8(img: np.ndarray) -> np.ndarray:
    img = np.asarray(img)
    if img.ndim != 3:
        raise ValueError(f"img must be HxWxC, got {img.shape}")
    if img.shape[2] == 1:
        img = np.repeat(img, 3, axis=2)
    if img.shape[2] == 4:
        img = img[..., :3]
    if img.shape[2] != 3:
        raise ValueError(f"img must have 3 RGB channels, got {img.shape}")

    if img.dtype == np.uint8:
        return img
    if img.dtype == np.uint16:
        return np.clip(img.astype(np.float32) / 65535.0 * 255.0, 0, 255).astype(np.uint8)

    x = img.astype(np.float32)
    if np.nanmax(x) <= 1.5:
        x *= 255.0
    return np.clip(x, 0, 255).astype(np.uint8)


def _imagenet_norm_params(device: torch.device):
    mean = torch.tensor([0.485, 0.456, 0.406], device=device).view(1, 3, 1, 1)
    std = torch.tensor([0.229, 0.224, 0.225], device=device).view(1, 3, 1, 1)
    return mean, std


def _center_weights_14(sigma: float = 0.5) -> np.ndarray:
    x = np.linspace(-1, 1, 14, dtype=np.float32)
    y = np.linspace(-1, 1, 14, dtype=np.float32)
    xx, yy = np.meshgrid(x, y)
    w = np.exp(-(xx**2 + yy**2) / (2.0 * sigma**2))
    return (w / (w.max() + 1e-12)).astype(np.float32)


def _uni_forward_tokens(
    uni: torch.nn.Module, x: torch.Tensor
) -> Tuple[Optional[torch.Tensor], Optional[torch.Tensor]]:
    """Best-effort extraction of `(cls, patch_tokens)` from common UNI/timm wrappers."""

    def unwrap(out):
        if isinstance(out, (tuple, list)) and out:
            out = out[0]
        if isinstance(out, dict):
            for key in ("tokens", "last_hidden_state", "x", "feat", "features"):
                if key in out:
                    return out[key]
        return out

    def parse(out_t: torch.Tensor):
        if out_t.ndim == 2:
            return out_t, None
        if out_t.ndim != 3:
            return None, None
        _, n_tokens, _ = out_t.shape
        if n_tokens == 197:
            return out_t[:, 0], out_t[:, 1:]
        if n_tokens == 196:
            return None, out_t
        if n_tokens > 197:
            maybe = out_t[:, :197, :]
            return maybe[:, 0], maybe[:, 1:]
        return None, None

    out = unwrap(uni(x))
    if torch.is_tensor(out):
        cls, tokens = parse(out)
        if tokens is not None:
            return cls, tokens

    candidates = [uni]
    for attr in ("model", "backbone", "net", "encoder"):
        module = getattr(uni, attr, None)
        if isinstance(module, torch.nn.Module):
            candidates.append(module)

    for module in candidates:
        for method_name in ("forward_features", "forward_intermediates"):
            if not hasattr(module, method_name):
                continue
            try:
                out2 = unwrap(getattr(module, method_name)(x))
                if torch.is_tensor(out2):
                    cls, tokens = parse(out2)
                    if tokens is not None:
                        return cls, tokens
            except Exception:
                pass

    if torch.is_tensor(out) and out.ndim == 2:
        return out, None
    return None, None


@torch.inference_mode()
def encode_he_with_uni_highres_from_img(
    img: np.ndarray,
    uni: torch.nn.Module,
    *,
    patch_size: int = 224,
    stride: int = 112,
    token_size: int = 16,
    batch_size: int = 32,
    device: str = "cuda:0",
    normalize: bool = True,
    pad: Optional[int] = None,
    pad_mode: str = "constant",
    pad_value: int = 0,
    roi_mask: Optional[np.ndarray] = None,
    blend_alpha: Optional[float] = None,
    center_sigma: float = 0.5,
    out_dtype: str = "float16",
) -> Tuple[np.ndarray, dict[str, Any]]:
    """Return dense UNI features as `(D, H_tokens, W_tokens)`."""
    img_u8 = _as_rgb_uint8(img)
    h, w = img_u8.shape[:2]

    if pad is None:
        pad = patch_size // 2
    if roi_mask is not None:
        if roi_mask.shape[:2] != (h, w):
            raise ValueError(f"roi_mask must be {(h, w)}, got {roi_mask.shape}")
        roi_mask = roi_mask.astype(bool)

    if pad > 0:
        if pad_mode == "constant":
            img_pad = np.pad(
                img_u8,
                ((pad, pad), (pad, pad), (0, 0)),
                mode="constant",
                constant_values=pad_value,
            )
            mask_pad = None if roi_mask is None else np.pad(
                roi_mask, ((pad, pad), (pad, pad)), mode="constant", constant_values=0
            )
        else:
            img_pad = np.pad(img_u8, ((pad, pad), (pad, pad), (0, 0)), mode=pad_mode)
            mask_pad = None if roi_mask is None else np.pad(roi_mask, ((pad, pad), (pad, pad)), mode=pad_mode)
    else:
        img_pad = img_u8
        mask_pad = roi_mask

    hp, wp = img_pad.shape[:2]
    hp2 = int(np.ceil(hp / token_size) * token_size)
    wp2 = int(np.ceil(wp / token_size) * token_size)
    if hp2 != hp or wp2 != wp:
        img_pad = np.pad(
            img_pad,
            ((0, hp2 - hp), (0, wp2 - wp), (0, 0)),
            mode="constant",
            constant_values=pad_value,
        )
        if mask_pad is not None:
            mask_pad = np.pad(mask_pad, ((0, hp2 - hp), (0, wp2 - wp)), mode="constant", constant_values=0)
        hp, wp = hp2, wp2

    h_tokens_pad, w_tokens_pad = hp // token_size, wp // token_size
    tops = list(range(0, hp - patch_size + 1, stride))
    lefts = list(range(0, wp - patch_size + 1, stride))

    dev = torch.device(device)
    uni = uni.to(dev).eval()
    mean, std = _imagenet_norm_params(dev) if normalize else (None, None)
    weights14 = _center_weights_14(center_sigma)

    feat_sum = None
    count_sum = np.zeros((h_tokens_pad, w_tokens_pad), dtype=np.float32)
    batch_patches: list[np.ndarray] = []
    batch_pos: list[tuple[int, int]] = []

    def patch_has_roi(top_px: int, left_px: int) -> bool:
        if mask_pad is None:
            return True
        sub = mask_pad[top_px : top_px + patch_size : token_size, left_px : left_px + patch_size : token_size]
        return bool(np.any(sub))

    def flush_batch():
        nonlocal feat_sum, batch_patches, batch_pos, count_sum
        if not batch_pos:
            return

        arr = np.stack(batch_patches, axis=0)
        x = torch.from_numpy(arr).permute(0, 3, 1, 2).float().to(dev) / 255.0
        if normalize:
            x = (x - mean) / std

        cls, tokens = _uni_forward_tokens(uni, x)
        if tokens is None:
            raise RuntimeError("UNI did not expose 196 patch tokens; high-res extraction requires tokens.")
        if tokens.ndim == 3 and tokens.shape[1] == 197:
            tokens = tokens[:, 1:, :]
        if tokens.shape[1] != 196:
            raise RuntimeError(f"Expected 196 UNI tokens per patch, got {tokens.shape[1]}")

        batch_n, _, dim = tokens.shape
        tok = tokens.view(batch_n, 14, 14, dim)
        if blend_alpha is not None and cls is not None:
            alpha = float(blend_alpha)
            tok = alpha * tok + (1.0 - alpha) * cls[:, None, None, :].expand(-1, 14, 14, -1)

        tok_np = tok.detach().cpu().numpy().astype(np.float32, copy=False)
        if feat_sum is None:
            feat_sum = np.zeros((h_tokens_pad, w_tokens_pad, dim), dtype=np.float32)

        for i, (top_px, left_px) in enumerate(batch_pos):
            top_t = top_px // token_size
            left_t = left_px // token_size
            feat_sum[top_t : top_t + 14, left_t : left_t + 14, :] += tok_np[i] * weights14[..., None]
            count_sum[top_t : top_t + 14, left_t : left_t + 14] += weights14

        batch_patches = []
        batch_pos = []

    for top in tops:
        for left in lefts:
            if not patch_has_roi(top, left):
                continue
            batch_patches.append(img_pad[top : top + patch_size, left : left + patch_size])
            batch_pos.append((top, left))
            if len(batch_pos) >= batch_size:
                flush_batch()
    flush_batch()

    if feat_sum is None:
        return np.zeros((0, 0, 0), dtype=np.float32), {"reason": "no patches selected"}

    feat_map = feat_sum / np.maximum(count_sum, 1e-8)[..., None]
    start_y = pad // token_size
    start_x = pad // token_size
    h_out = h // token_size
    w_out = w // token_size
    feat_out = feat_map[start_y : start_y + h_out, start_x : start_x + w_out, :]
    feat_out = np.transpose(feat_out, (2, 0, 1))
    feat_out = feat_out.astype(np.float16 if out_dtype == "float16" else np.float32, copy=False)

    meta = {
        "H": h,
        "W": w,
        "patch_size": patch_size,
        "stride": stride,
        "token_size": token_size,
        "pad": pad,
        "grid_tokens_out": (h_out, w_out),
        "dim": int(feat_out.shape[0]),
        "blend_alpha": blend_alpha,
        "out_dtype": out_dtype,
    }
    return feat_out, meta


# ---------------------------------------------------------------------
# H&E background/tissue mask
# ---------------------------------------------------------------------


def ensure_uint8_rgb(he: np.ndarray) -> np.ndarray:
    he = np.asarray(he)
    if he.ndim == 2:
        he = np.stack([he, he, he], axis=-1)
    if he.ndim != 3 or he.shape[-1] not in (3, 4):
        raise ValueError(f"Expected image of shape (H, W, 3/4), got {he.shape}")
    if he.shape[-1] == 4:
        he = he[..., :3]
    if he.dtype != np.uint8:
        if np.issubdtype(he.dtype, np.floating):
            he_min, he_max = np.nanmin(he), np.nanmax(he)
            if he_max > he_min:
                he = (he - he_min) / (he_max - he_min)
            he = (255.0 * np.clip(he, 0, 1)).astype(np.uint8)
        else:
            he = he.astype(np.uint8)
    return he


def pad_to_multiple(he: np.ndarray, patch_size: int, pad_value: int = 255) -> np.ndarray:
    h, w = he.shape[:2]
    new_h = int(np.ceil(h / patch_size) * patch_size)
    new_w = int(np.ceil(w / patch_size) * patch_size)
    pad_width = ((0, new_h - h), (0, new_w - w), (0, 0))
    return np.pad(he, pad_width, mode="constant", constant_values=pad_value) if (new_h != h or new_w != w) else he


def patchify(x: np.ndarray, patch_size: int):
    if x.ndim != 3:
        raise ValueError("patchify expects image of shape (H, W, C)")
    shape_ori = np.array(x.shape[:2])
    shape_ext = (shape_ori + patch_size - 1) // patch_size * patch_size
    x_padded = np.pad(
        x,
        ((0, shape_ext[0] - x.shape[0]), (0, shape_ext[1] - x.shape[1]), (0, 0)),
        mode="edge",
    )
    tiles_shape = shape_ext // patch_size
    strides = (
        x_padded.strides[0] * patch_size,
        x_padded.strides[1] * patch_size,
        x_padded.strides[0],
        x_padded.strides[1],
        x_padded.strides[2],
    )
    patches = np.lib.stride_tricks.as_strided(
        x_padded,
        shape=(tiles_shape[0], tiles_shape[1], patch_size, patch_size, x.shape[2]),
        strides=strides,
    ).reshape(-1, patch_size, patch_size, x.shape[2])
    return patches, {"original": shape_ori, "padded": shape_ext, "tiles": tiles_shape}


def compute_metrics(he: np.ndarray, patch_size: int = 16):
    he_tiles, shapes = patchify(he, patch_size=patch_size)
    he_std_image = np.std(he_tiles, axis=(1, 2, 3)).reshape(shapes["tiles"])

    mean_rgb = np.mean(he_tiles, axis=(1, 2))
    channel_var = np.var(mean_rgb, axis=0)
    z_v = (mean_rgb * channel_var).sum(axis=1) / (channel_var.sum() + 1e-8)
    z_v_image = z_v.reshape(shapes["tiles"])

    flat_std = he_std_image.flatten()
    flat_std_norm = (flat_std - flat_std.min()) / (flat_std.max() - flat_std.min() + 1e-8) + 1.0
    he_std_norm_image = flat_std_norm.reshape(z_v_image.shape)

    flat_mean = z_v_image.flatten()
    flat_mean_norm = (flat_mean - flat_mean.min()) / (flat_mean.max() - flat_mean.min() + 1e-8) + 1.0
    z_v_norm_image = flat_mean_norm.reshape(z_v_image.shape)

    ratio_norm = flat_std_norm / (flat_mean_norm + 1e-8)
    ratio_norm = (ratio_norm - ratio_norm.min()) / (ratio_norm.max() - ratio_norm.min() + 1e-8)
    ratio_norm_image = ratio_norm.reshape(z_v_image.shape)
    return he_std_norm_image, he_std_image, z_v_norm_image, z_v_image, ratio_norm, ratio_norm_image


def compute_low_density_mask(
    z_v_image: np.ndarray,
    he_std_image: np.ndarray,
    ratio_norm: np.ndarray,
    density_thresh: float = 100.0,
) -> np.ndarray:
    from scipy.ndimage import gaussian_filter

    mean = z_v_image.flatten()
    std = he_std_image.flatten()
    bin_n = max(int(len(ratio_norm) * 0.00005), 32)
    hist, xedges, yedges = np.histogram2d(mean, std, bins=[bin_n, bin_n])
    smoothed_density = gaussian_filter(hist, sigma=max(1, int(bin_n * 0.02)))

    xidx = np.searchsorted(xedges, mean) - 1
    yidx = np.searchsorted(yedges, std) - 1
    valid = (
        (xidx >= 0)
        & (xidx < smoothed_density.shape[0])
        & (yidx >= 0)
        & (yidx < smoothed_density.shape[1])
    )
    density = np.zeros_like(std)
    density[valid] = smoothed_density[xidx[valid], yidx[valid]]
    return (density < density_thresh).reshape(z_v_image.shape)


def run_texture_analysis(
    image: np.ndarray,
    tissue_mask: np.ndarray,
    patch_size: int = 16,
    glcm_levels: int = 64,
    output_dir: str | None = None,
) -> np.ndarray:
    from skimage.color import rgb2gray
    from skimage.feature import graycomatrix, graycoprops
    from sklearn.mixture import GaussianMixture

    if output_dir is not None:
        os.makedirs(output_dir, exist_ok=True)

    gray_image = (rgb2gray(image) * 255).astype(np.uint8)
    gray_image = (gray_image / 255.0 * (glcm_levels - 1)).astype(np.uint8)
    mask = tissue_mask.astype(bool)
    h_mask, w_mask = mask.shape
    assert gray_image.shape[0] // patch_size == h_mask and gray_image.shape[1] // patch_size == w_mask

    energy_map = np.full(mask.shape, np.nan)
    homogeneity_map = np.full(mask.shape, np.nan)
    entropy_map = np.full(mask.shape, np.nan)
    color_filter_mask = np.zeros(mask.shape, dtype=bool)
    color_filtered = 0

    for i in range(h_mask):
        for j in range(w_mask):
            if not mask[i, j]:
                continue
            y0, y1 = i * patch_size, (i + 1) * patch_size
            x0, x1 = j * patch_size, (j + 1) * patch_size
            patch_gray = gray_image[y0:y1, x0:x1]
            patch_rgb = image[y0:y1, x0:x1, :]

            mean_r, mean_g, mean_b = np.mean(patch_rgb, axis=(0, 1))
            is_green = (mean_g > mean_r + 20) and (mean_g > mean_b + 20)
            is_gray = np.std([mean_r, mean_g, mean_b]) < 10
            is_too_bright = (mean_r > 230) and (mean_g > 230) and (mean_b > 230)
            if is_green or is_gray or is_too_bright:
                color_filter_mask[i, j] = True
                color_filtered += 1
                continue

            glcm = graycomatrix(
                patch_gray,
                distances=[1],
                angles=[0],
                levels=glcm_levels,
                symmetric=True,
                normed=True,
            )
            energy_map[i, j] = graycoprops(glcm, "energy")[0, 0]
            homogeneity_map[i, j] = graycoprops(glcm, "homogeneity")[0, 0]
            vals = glcm[:, :, 0, 0]
            vals = vals[vals > 0]
            entropy_map[i, j] = -np.sum(vals * np.log2(vals))

    if mask.sum() == 0 or mask.sum() - color_filtered < 2:
        return mask.copy()

    def normalize_map(m: np.ndarray) -> np.ndarray:
        valid = ~np.isnan(m)
        if not valid.any():
            return m
        m_min, m_max = np.nanmin(m), np.nanmax(m)
        return np.zeros_like(m) if m_max <= m_min else (m - m_min) / (m_max - m_min)

    energy = normalize_map(energy_map).flatten()
    homog = normalize_map(homogeneity_map).flatten()
    entropy = normalize_map(entropy_map).flatten()
    valid = ~(np.isnan(energy) | np.isnan(homog) | np.isnan(entropy))
    if not valid.any():
        return mask.copy()

    features = np.stack([homog[valid], energy[valid], entropy[valid]], axis=1)
    labels = GaussianMixture(n_components=4, random_state=45).fit_predict(features)

    cluster_map = np.full(mask.shape, np.nan)
    flat_indices = np.where(valid)[0]
    for idx, label in zip(flat_indices, labels):
        cluster_map[idx // w_mask, idx % w_mask] = label

    cluster_scores = {}
    energy2 = normalize_map(energy_map)
    homog2 = normalize_map(homogeneity_map)
    entropy2 = normalize_map(entropy_map)
    for k in range(4):
        coords = np.where(cluster_map == k)
        if len(coords[0]) == 0:
            cluster_scores[k] = np.inf
        else:
            cluster_scores[k] = np.nanmean(energy2[coords]) + np.nanmean(homog2[coords]) - np.nanmean(entropy2[coords])

    keep_labels = [k for k, score in sorted(cluster_scores.items(), key=lambda kv: kv[1])[:2] if np.isfinite(score)]
    updated_mask = mask.copy()
    if keep_labels:
        updated_mask[np.where(np.isin(cluster_map, keep_labels))] = False

    if output_dir is not None:
        cluster_rgb = np.zeros((h_mask, w_mask, 3), dtype=np.uint8)
        for i in range(h_mask):
            for j in range(w_mask):
                if color_filter_mask[i, j]:
                    cluster_rgb[i, j] = [30, 144, 255]
                elif np.isnan(cluster_map[i, j]):
                    cluster_rgb[i, j] = [0, 0, 0]
                elif cluster_map[i, j] in keep_labels:
                    cluster_rgb[i, j] = [0, 255, 0]
                else:
                    cluster_rgb[i, j] = [255, 0, 0]
        Image.fromarray(cluster_rgb).save(Path(output_dir) / "cluster_labels_colored.png")

    return updated_mask


def run_ratio_filtering(ratio_norm: np.ndarray, mask1_lowdensity_updated: np.ndarray):
    from skimage.filters import threshold_otsu

    kept = ~mask1_lowdensity_updated.flatten()
    data = ratio_norm[kept]
    if len(data) == 0:
        return np.zeros_like(ratio_norm, dtype=bool), 0.0
    data = data[data < np.percentile(data, 90)]
    if len(data) == 0:
        return np.zeros_like(ratio_norm, dtype=bool), 0.0
    otsu_thresh = threshold_otsu(data)
    return ratio_norm < otsu_thresh, otsu_thresh


def generate_final_mask(
    he: np.ndarray,
    mask1_updated: np.ndarray,
    mask2: np.ndarray,
    super_pixel_size: int = 16,
    min_size: int = 10,
    clean_background: bool = True,
    output_dir: str | None = None,
):
    image_height, image_width = he.shape[:2]
    ny, nx = image_height // super_pixel_size, image_width // super_pixel_size
    mask1_updated = np.asarray(mask1_updated).reshape(ny, nx)
    mask2 = np.asarray(mask2).reshape(ny, nx)
    mask_superpixels = ~(mask1_updated | mask2)

    if clean_background:
        from skimage import morphology

        mask_superpixels = morphology.remove_small_objects(mask_superpixels.astype(bool), min_size=min_size, connectivity=2)

    mask_fullres = np.kron(mask_superpixels, np.ones((super_pixel_size, super_pixel_size), dtype=bool))
    mask_fullres = mask_fullres[:image_height, :image_width]

    if output_dir is not None:
        os.makedirs(output_dir, exist_ok=True)
        Image.fromarray(mask_superpixels.astype(np.uint8) * 255).save(Path(output_dir) / "mask-small.png")
        Image.fromarray(mask_fullres.astype(np.uint8) * 255).save(Path(output_dir) / "mask.png")

    return mask_superpixels, mask_fullres


def he_background_mask(
    he: np.ndarray,
    patch_size: int = 16,
    density_thresh: float = 100.0,
    min_size: int = 10,
    glcm_levels: int = 64,
    clean_background: bool = True,
    texture_output_dir: str | None = None,
    mask_output_dir: str | None = None,
):
    he_padded = pad_to_multiple(ensure_uint8_rgb(he), patch_size=patch_size, pad_value=255)
    (
        he_std_norm_image,
        he_std_image,
        z_v_norm_image,
        z_v_image,
        ratio_norm,
        ratio_norm_image,
    ) = compute_metrics(he_padded, patch_size=patch_size)

    mask1_lowdensity = compute_low_density_mask(z_v_image, he_std_image, ratio_norm, density_thresh=density_thresh)
    mask1_lowdensity_updated = run_texture_analysis(
        he_padded,
        mask1_lowdensity,
        patch_size=patch_size,
        glcm_levels=glcm_levels,
        output_dir=texture_output_dir,
    )
    mask2_lowratio, otsu_thresh = run_ratio_filtering(ratio_norm, mask1_lowdensity_updated)
    mask_superpixels, mask_fullres = generate_final_mask(
        he_padded,
        mask1_lowdensity_updated,
        mask2_lowratio,
        super_pixel_size=patch_size,
        min_size=min_size,
        clean_background=clean_background,
        output_dir=mask_output_dir,
    )

    return {
        "he_processed": he_padded,
        "he_std_norm_image": he_std_norm_image,
        "he_std_image": he_std_image,
        "z_v_norm_image": z_v_norm_image,
        "z_v_image": z_v_image,
        "ratio_norm": ratio_norm,
        "ratio_norm_image": ratio_norm_image,
        "mask1_lowdensity": mask1_lowdensity,
        "mask1_lowdensity_updated": mask1_lowdensity_updated,
        "mask2_lowratio": mask2_lowratio,
        "mask_superpixels": mask_superpixels,
        "mask_fullres": mask_fullres,
        "otsu_thresh": otsu_thresh,
    }


# ---------------------------------------------------------------------
# Classifier prediction on feature grid
# ---------------------------------------------------------------------


class UNIMLP(torch.nn.Module):
    def __init__(self, in_dim=1024, width=512, depth=3, dropout=0.2):
        super().__init__()
        layers = [torch.nn.Linear(in_dim, width), torch.nn.GELU(), torch.nn.Dropout(dropout)]
        for _ in range(depth - 1):
            layers += [torch.nn.Linear(width, width), torch.nn.GELU(), torch.nn.Dropout(dropout)]
        self.trunk = torch.nn.Sequential(*layers)

    def forward(self, x):
        return self.trunk(x)


class UNIMLPRes(torch.nn.Module):
    def __init__(self, in_dim=1024, width=512, depth=3, dropout=0.3):
        super().__init__()
        self.in_ln = torch.nn.LayerNorm(in_dim)
        self.fc_in = torch.nn.Linear(in_dim, width)
        self.blocks = torch.nn.ModuleList()
        for _ in range(max(0, depth - 1)):
            self.blocks.append(
                torch.nn.Sequential(
                    torch.nn.LayerNorm(width),
                    torch.nn.Linear(width, width * 2),
                    torch.nn.GELU(),
                    torch.nn.Dropout(dropout),
                    torch.nn.Linear(width * 2, width),
                    torch.nn.Dropout(dropout),
                )
            )
        self.out_ln = torch.nn.LayerNorm(width)

    def forward(self, x):
        x = self.in_ln(x)
        x = self.fc_in(x)
        for block in self.blocks:
            x = x + block(x)
        return self.out_ln(x)


class UNIClassifier(torch.nn.Module):
    def __init__(self, in_dim=1024, n_classes=7, width=512, depth=3, dropout=0.3, trunk="resmlp"):
        super().__init__()
        if trunk == "resmlp":
            self.backbone = UNIMLPRes(in_dim=in_dim, width=width, depth=depth, dropout=dropout)
        elif trunk == "mlp":
            self.backbone = UNIMLP(in_dim=in_dim, width=width, depth=depth, dropout=dropout)
        elif trunk == "linear_probe":
            self.backbone = torch.nn.Sequential(torch.nn.LayerNorm(in_dim), torch.nn.Identity())
            width = in_dim
        else:
            raise ValueError("trunk must be one of ['resmlp', 'mlp', 'linear_probe']")
        self.trunk = trunk
        self.width = width
        self.head = torch.nn.Linear(width, n_classes)

    def forward(self, x):
        return self.head(self.backbone(x))


def load_ckpt_classifier(path: str | Path, map_location="cpu") -> Tuple[UNIClassifier, dict[str, Any]]:
    ckpt = torch.load(str(path), map_location=map_location)
    meta = ckpt.get("meta", {}) or {}
    model = UNIClassifier(
        in_dim=int(meta.get("in_dim", 1024)),
        n_classes=int(meta.get("n_classes", meta.get("num_classes", 2))),
        width=int(meta.get("width", 512)),
        depth=int(meta.get("depth", 3)),
        dropout=float(meta.get("dropout", 0.2)),
        trunk=str(meta.get("trunk", "mlp")),
    )
    model.load_state_dict(ckpt["state_dict"], strict=True)
    if "label_order" in meta and "id_to_name" not in meta:
        label_order = [str(x) for x in meta["label_order"]]
        meta["id_to_name"] = {i: n for i, n in enumerate(label_order)}
        meta["name_to_id"] = {n: i for i, n in enumerate(label_order)}
    return model, meta


def _get_axes(layout: str):
    layout = layout.upper()
    if set(layout) != {"D", "H", "W"} or len(layout) != 3:
        raise ValueError(f"layout must be a permutation of 'DHW', got {layout}")
    return layout.index("D"), layout.index("H"), layout.index("W")


def _get_block_x(feat: np.ndarray, layout: str, y0: int, y1: int, x0: int, x1: int) -> np.ndarray:
    ax_d, ax_h, ax_w = _get_axes(layout)
    sl = [slice(None)] * 3
    sl[ax_h] = slice(y0, y1)
    sl[ax_w] = slice(x0, x1)
    blk = feat[tuple(sl)]
    blk_dhw = np.moveaxis(blk, (ax_d, ax_h, ax_w), (0, 1, 2))
    d, h, w = blk_dhw.shape
    return np.asarray(blk_dhw.reshape(d, h * w).T, dtype=np.float32)


@torch.no_grad()
def predict_on_feat_array_class(
    feat: np.ndarray,
    model: torch.nn.Module,
    *,
    layout: str = "DHW",
    device: str = "cuda:0",
    block_hw: tuple[int, int] = (256, 256),
    batch_tokens: int = 8192,
    out_dtype=np.uint8,
) -> np.ndarray:
    """Predict a class id for each token-grid position; returns `(Ht, Wt)`."""
    if feat.ndim != 3:
        raise ValueError(f"feat must be 3D, got {feat.shape}")
    _, ax_h, ax_w = _get_axes(layout)
    ht, wt = int(feat.shape[ax_h]), int(feat.shape[ax_w])

    model = model.to(device).eval()
    out = np.zeros((ht, wt), dtype=out_dtype)
    bh, bw = map(int, block_hw)
    for y0 in range(0, ht, bh):
        y1 = min(ht, y0 + bh)
        for x0 in range(0, wt, bw):
            x1 = min(wt, x0 + bw)
            x = _get_block_x(feat, layout, y0, y1, x0, x1)
            preds = []
            for start in range(0, x.shape[0], batch_tokens):
                xb = torch.from_numpy(x[start : start + batch_tokens]).to(device, non_blocking=True)
                logits = model(xb).float()
                preds.append(logits.argmax(1).cpu().numpy().astype(out_dtype))
            out[y0:y1, x0:x1] = np.concatenate(preds, axis=0).reshape(y1 - y0, x1 - x0)
    return out


def resize_bool_pil(mask: np.ndarray, out_hw: Tuple[int, int]) -> np.ndarray:
    m = (mask > 0).astype(np.uint8) * 255
    im = Image.fromarray(m, mode="L")
    im2 = im.resize((out_hw[1], out_hw[0]), resample=Image.NEAREST)
    return np.asarray(im2) > 0


# ---------------------------------------------------------------------
# Minimal categorical pixel plot replacement
# ---------------------------------------------------------------------


def _build_id2name(id2cat: Sequence[str] | Mapping[int, str] | None) -> dict[int, str] | None:
    if id2cat is None:
        return None
    if isinstance(id2cat, Mapping):
        return {int(k): str(v) for k, v in id2cat.items()}
    return {i: str(v) for i, v in enumerate(id2cat)}


def _categorical_rgba(
    labels: np.ndarray,
    *,
    palette: Mapping[Any, str] | None,
    id2cat: Sequence[str] | Mapping[int, str] | None,
    na_value: int | None,
    na_color: str,
) -> tuple[np.ndarray, list[tuple[str, Any]]]:
    import matplotlib.pyplot as plt
    from matplotlib.colors import to_rgba

    labels = np.asarray(labels)
    if labels.ndim != 2:
        raise ValueError("labels must be a 2D array")

    id2name = _build_id2name(id2cat)
    na_mask = (labels == int(na_value)) if na_value is not None and np.issubdtype(labels.dtype, np.integer) else np.zeros(labels.shape, bool)
    valid_ids = np.unique(labels[~na_mask]).astype(int) if np.any(~na_mask) else np.array([], dtype=int)

    max_id = int(valid_ids.max()) if valid_ids.size else 0
    lut = np.zeros((max_id + 1, 4), dtype=np.float32)
    tab20 = plt.get_cmap("tab20")
    legend_items = []

    for class_id in valid_ids:
        name = str(id2name[class_id]) if id2name is not None and class_id in id2name else str(class_id)
        color = None
        if palette is not None:
            color = palette.get(class_id) or palette.get(str(class_id)) or palette.get(name)
        rgba = to_rgba(color) if color is not None else tab20(class_id % 20)
        lut[class_id] = rgba
        legend_items.append((name, rgba))

    clipped = labels.copy()
    clipped[na_mask] = 0
    clipped = np.clip(clipped, 0, max_id)
    rgba_img = lut[clipped]
    rgba_img[na_mask] = to_rgba(na_color)
    return rgba_img, legend_items


def pixel_plot_categorical(
    color: np.ndarray,
    *,
    palette: Mapping[Any, str] | None = None,
    id2cat: Sequence[str] | Mapping[int, str] | None = None,
    na_value: int | None = -1,
    na_color: str = "#000000",
    legend: bool = True,
    output_path: str | Path | None = None,
    figsize: tuple[float, float] | None = None,
    dpi: int = 300,
):
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    rgba, legend_items = _categorical_rgba(color, palette=palette, id2cat=id2cat, na_value=na_value, na_color=na_color)
    fig, ax = plt.subplots(figsize=figsize if figsize is not None else (8, 8))
    ax.imshow(rgba, origin="upper", interpolation="nearest")
    ax.set_aspect("equal")
    ax.set_axis_off()

    if legend and legend_items:
        handles = [Patch(facecolor=rgba_val, edgecolor="none", label=label) for label, rgba_val in legend_items]
        ax.legend(handles=handles, loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False)

    if output_path is not None:
        out = Path(output_path)
        out.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out, dpi=dpi, bbox_inches="tight")
    return fig, ax


def predict_slide(
    *,
    feat_path: str | Path,
    background_path: str | Path,
    clf: torch.nn.Module,
    id2cat: Sequence[str] | Mapping[int, str],
    palette: Mapping[Any, str] | None,
    output_path: str | Path,
    layout: str = "DHW",
    device: str = "cuda:0",
) -> np.ndarray:
    """Load feature and background files, run model prediction, mask background, save PNG."""
    feat = np.load(feat_path)
    with open(background_path, "rb") as f:
        back_res = pickle.load(f)

    pred_grid = predict_on_feat_array_class(feat, clf, layout=layout, device=device).astype(np.int16)
    mask_tok = resize_bool_pil(back_res["mask_fullres"], (pred_grid.shape[0], pred_grid.shape[1]))
    pred_grid[~mask_tok] = -1
    pixel_plot_categorical(
        pred_grid,
        palette=palette,
        id2cat=id2cat,
        na_value=-1,
        na_color="#000000",
        legend=True,
        output_path=output_path,
    )
    return pred_grid


if __name__ == "__main__":
    print(
        "Import this module and pass your loaded `uni` and `clf` objects. "
        "See the docstring at the top for example calls."
    )
