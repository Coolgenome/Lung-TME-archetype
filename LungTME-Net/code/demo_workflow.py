"""
Run the demo pipeline for high-resolution H&E images.

Expected inputs for demo:
  ../HE_demo/ADC1.png
  ../HE_demo/ADC2.png
  ../HE_demo/ADC3.png
"""

from __future__ import annotations

import argparse
import os
import pickle
from pathlib import Path

import numpy as np
import torch
from PIL import Image

from he_pred import (
    encode_he_with_uni_highres_from_img,
    he_background_mask,
    load_ckpt_classifier,
    pixel_plot_categorical,
    predict_on_feat_array_class,
    resize_bool_pil,
)


DEMO_DIR = Path("../HE_demo")
UNI_WEIGHTS = "../UNI/pytorch_model.bin" # download the UNI weights and put it here
TUMOR_CKPT = "../model/general_uni_cls_tumor_region.pt"
LINEAGE_CKPT = "../model/general_uni_cls_cell_lineage.pt"

DEFAULT_CROPS = [
    DEMO_DIR / "ADC1.png",
    DEMO_DIR / "ADC2.png",
    DEMO_DIR / "ADC3.png",
]


def choose_device(requested: str) -> str:
    if requested != "auto":
        return requested
    if torch.cuda.is_available():
        return "cuda:0"
    if hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
        return "mps"
    return "cpu"


def load_uni_model(weights_path: Path = UNI_WEIGHTS):
    import timm

    model = timm.create_model(
        "vit_large_patch16_224",
        img_size=224,
        patch_size=16,
        init_values=1e-5,
        num_classes=0,
        dynamic_img_size=True,
    )
    state = torch.load(str(weights_path), map_location="cpu")
    model.load_state_dict(state, strict=True)
    return model.eval()


def label_order_from_meta(meta: dict) -> list[str]:
    if "label_order" in meta:
        return [str(x) for x in meta["label_order"]]
    if "id_to_name" in meta:
        id_to_name = {int(k): str(v) for k, v in meta["id_to_name"].items()}
        return [id_to_name[i] for i in range(max(id_to_name) + 1)]
    n_classes = int(meta.get("n_classes", meta.get("num_classes", 0)))
    return [str(i) for i in range(n_classes)]


def default_palette(label_order: list[str]) -> dict[str, str]:
    original_like = {
        "Tumor_epi": "#ef3b2c",
        "tumor": "#d73027",
        "non_tumor": "#4393c3",
        "Normal Alveolar": "#e7e1ee",
        "Lymphocytes": "#23b9dc",
        "Macrophages": "#a8d36a",
        "Fibrous tissue": "#ffd700",
        "SMC & Vessel": "#76cbd6",
        "Bronchus": "#e7298a",
    }
    fallback = [
        "#d73027", "#4393c3", "#6aacd5", "#ef689e", "#4e9242", "#fb8072",
        "#8b6bae", "#ffc452", "#ca4f27", "#fcdfd1", "#80b1d3", "#b3de69",
    ]
    return {name: original_like.get(name, fallback[i % len(fallback)]) for i, name in enumerate(label_order)}


def predict_one_model(
    *,
    feat: np.ndarray,
    bg_mask_fullres: np.ndarray,
    ckpt_path: Path,
    out_prefix: Path,
    device: str,
) -> tuple[np.ndarray, list[str]]:
    clf, meta = load_ckpt_classifier(ckpt_path, map_location="cpu")
    label_order = label_order_from_meta(meta)
    palette = default_palette(label_order)
    layout = str(meta.get("feat_layout", "DHW"))

    pred_grid = predict_on_feat_array_class(feat, clf, layout=layout, device=device).astype(np.int16)

    # Match the original notebook helper: PIL nearest-neighbor resize of full-res bg mask.
    bg_mask_tok = resize_bool_pil(bg_mask_fullres, (pred_grid.shape[0], pred_grid.shape[1]))
    pred_grid[~bg_mask_tok] = -1

    np.save(str(out_prefix) + ".npy", pred_grid)
    pixel_plot_categorical(
        pred_grid,
        palette=palette,
        id2cat=label_order,
        na_value=-1,
        na_color="#000000",
        legend=True,
        output_path=str(out_prefix) + ".png",
    )
    return pred_grid, label_order


def process_crop(
    *,
    image_path: Path,
    out_dir: Path,
    uni,
    device: str,
    force: bool,
) -> None:
    os.environ.setdefault("MPLCONFIGDIR", str(out_dir / "matplotlib_cache"))
    Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)
    out_dir.mkdir(parents=True, exist_ok=True)

    img = np.asarray(Image.open(image_path).convert("RGB"))

    feat_path = out_dir / "HE_token.npy"
    feat_meta_path = out_dir / "HE_token_meta.pkl"
    if force or not feat_path.exists():
        print(f"[{image_path.name}] extracting UNI features")
        feat, feat_meta = encode_he_with_uni_highres_from_img(
            img=img,
            uni=uni,
            patch_size=224,
            stride=112,
            batch_size=64,
            device=device,
            normalize=True,
            roi_mask=None,
            blend_alpha=0.5,
            out_dtype="float16",
        )
        np.save(feat_path, feat)
        with open(feat_meta_path, "wb") as f:
            pickle.dump(feat_meta, f, protocol=pickle.HIGHEST_PROTOCOL)
    else:
        feat = np.load(feat_path)

    bg_path = out_dir / "HE_background.pkl"
    if force or not bg_path.exists():
        print(f"[{image_path.name}] computing Spyrrow-style bgmask")
        bg_res = he_background_mask(
            img,
            patch_size=16,
            density_thresh=100,
            min_size=10,
            glcm_levels=64,
            clean_background=True,
            texture_output_dir=None,
            mask_output_dir=str(out_dir / "mask_images"),
        )
        with open(bg_path, "wb") as f:
            pickle.dump(bg_res, f, protocol=pickle.HIGHEST_PROTOCOL)
    else:
        with open(bg_path, "rb") as f:
            bg_res = pickle.load(f)

    # he_background_mask pads to a multiple of patch_size; crop back before resizing to token grid.
    bg_mask_fullres = bg_res["mask_fullres"][: img.shape[0], : img.shape[1]].copy()

    print(f"[{image_path.name}] predicting tumor region")
    tumor_pred, tumor_labels = predict_one_model(
        feat=feat,
        bg_mask_fullres=bg_mask_fullres,
        ckpt_path=TUMOR_CKPT,
        out_prefix=out_dir / "tumor_region_prediction_bgmask",
        device=device,
    )

    print(f"[{image_path.name}] predicting cell lineage")
    lineage_pred, lineage_labels = predict_one_model(
        feat=feat,
        bg_mask_fullres=bg_mask_fullres,
        ckpt_path=LINEAGE_CKPT,
        out_prefix=out_dir / "cell_lineage_prediction_bgmask",
        device=device,
    )

    def summarize(name: str, arr: np.ndarray, labels: list[str]) -> None:
        vals, counts = np.unique(arr, return_counts=True)
        print(name)
        for value, count in zip(vals, counts):
            label = "background" if int(value) == -1 else labels[int(value)]
            print(f"  {int(value):>2} {label:<18} {count:>7} {count / arr.size:>7.1%}")

    summarize("tumor summary", tumor_pred, tumor_labels)
    summarize("lineage summary", lineage_pred, lineage_labels)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--image", action="append", type=Path, default=None)
    parser.add_argument("--out-root", type=Path, default=DEMO_DIR / "bgmask_only_crop_outputs")
    parser.add_argument("--device", default="auto")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    device = choose_device(args.device)
    image_paths = args.image or DEFAULT_CROPS
    print("device:", device)
    print("loading UNI:", UNI_WEIGHTS)
    uni = load_uni_model(UNI_WEIGHTS)

    for image_path in image_paths:
        crop_name = image_path.parent.name.replace("_outputs", "")
        out_dir = args.out_root / crop_name
        process_crop(image_path=image_path, out_dir=out_dir, uni=uni, device=device, force=args.force)
        print("[done]", image_path, "->", out_dir)


if __name__ == "__main__":
    main()
