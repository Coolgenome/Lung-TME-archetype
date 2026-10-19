# histoIMPACT

H&E → spatial tissue map → five tissue layers → IMPACT archetype (IMPACT1a/2, IMPACT1b, IMPACT3).

`code/run_histoimpact.py` runs the full pipeline on H&E images:

1. **Model 1** (`he_pred.py`): UNI token features (16 px = 4 µm at 0.25 µm/px), tissue mask, and per-token predictions of
   - the cell lineage (7 classes) and
   - the tumor region (tumor / non-tumor).
2. **Archetyping** (`impact_archetype.py`):
   - 48 µm super-pixels.
   - Tumor-region cleanup: majority vote over 3 × 3 super-pixels; islands and holes < 0.05 mm² removed; islands < 10 % of the largest tumor (and < 0.5 mm²) removed.
   - Five tissue layers by signed distance to the tumor border: TR < −100 µm, IT −100–0, IN 0–100, NP 100–300, ND ≥ 300 µm.
   - Visium-supervised CD4 T / CD8 T / B / macrophage / fibroblast maps.
   - 25 bivariate Moran's I features (map × layer) on 96 µm pseudo-spots. Slides larger than 11 mm use 11 mm windows (5.5 mm step).
   - TIL / Mac / Fib classifiers and the IMPACT rules.

## Outputs

Per sample (`results_histoimpact/<sample>/`):

| File | Content |
|---|---|
| `cell_lineage_grid.npy`, `cell_lineage_map.png` | Spatial tissue map of cell lineages (token grid, −1 = background) |
| `tumor_region_grid.npy`, `tumor_region_map.png` | Token-level tumor-region calls (0 = tumor_region, 1 = non_tumor_region) |
| `tissue_layers_token.npy`, `tissue_layers.png` | Five-layer segmentation on the token grid (0 TR, 1 IT, 2 IN, 3 NP, 4 ND, −1 background) |
| `tissue_layers_superpixel.npy` | The same on the 48 µm super-pixel grid, where the layers are defined |
| `signature_maps.npz` (+ `.png` with `--plot-maps`) | CD4_T, CD8_T, B, Mac, Fib maps (relative within slide) |
| `features.json` | 25 Moran's I features, layer areas, number of windows |
| `impact.json` | 3-class IMPACT call, compartment classes and probabilities |

Cohort summary: `results_histoimpact/histoimpact_predictions.csv`. Its key columns:

- `IMPACT_3class`: one of IMPACT1a/2, IMPACT1b, IMPACT3 or not_assessable. Slides whose tumor layers (TR + IT) cover less than 5 % of the tissue are not assessable.
- `til_margin`: gap between the two most likely TIL classes; below 0.05 means a borderline call.
- `tumor_layer_frac`

## Model files (`model/`)

| File | Source |
|---|---|
| `general_uni_cls_cell_lineage.pt`, `general_uni_cls_tumor_region.pt` | Model 1 |
| `signature_model.npz` | Visium-supervised multi-task MLP (copy of `expr/model_mt.npz`) |
| `impact_classifier.json` | Compartment classifiers (copy of `runs/cleanup/final_r0.1k0.5/final_model.json`) |

UNI weights: `UNI/pytorch_model.bin` (see `UNI/README.md`).

## Installation

```bash
bash setup_env.sh
conda activate histoimpact
```

`setup_env.sh`:

- creates the conda env `histoimpact` from `environment.yml` (Python 3.10, numpy, scipy, pandas, pillow, matplotlib, scikit-image, scikit-learn);
- installs PyTorch and timm with pip. The PyTorch build follows the hardware: CUDA 12.1 when an NVIDIA GPU is found, the CPU build otherwise, and the standard wheel with Apple-silicon (MPS) support on macOS;
- checks that all packages import, reports the compute device, and confirms the UNI weights and the four model files are present.

Options:

- `TORCH_CUDA=cu118 bash setup_env.sh`: use another CUDA build.
- `ENV_NAME=myenv bash setup_env.sh`: use another env name.
- `RECREATE=1 bash setup_env.sh`: rebuild the env.

PyTorch and timm are needed for model 1 (UNI features and token classification). The archetyping step (`impact_archetype.py`) needs only numpy, scipy and pandas. A GPU is strongly recommended for UNI feature extraction.

## Usage

```bash
cd code
# images at 40x (0.25 um/px)
python run_histoimpact.py --image ../HE_demo/ADC1.png --image ../HE_demo/ADC2.png --mpp 0.25
# images at another resolution are resampled to 0.25 um/px first
python run_histoimpact.py --image slide_20x.tif --mpp 0.5
# a cohort: manifest with sample_id, image, mpp
python run_histoimpact.py --manifest samples.csv --out-root ../results_cohort
```

**Whole-slide images.** UNI extraction from an image holds the full image and its token array in memory, so whole slides are better processed with a tiled token extractor. To run only the archetyping on existing outputs, use a manifest with these columns:

- `sample_id`
- `tokens`: (1024, H, W) `.npy`, read memory-mapped
- `lineage_grid` and `tumor_grid`: (H, W) `.npy`, −1 = background
- `um_per_token`

**Notes**

- The models were trained at 0.25 µm/px. Native 20× tokens shift the TIL features toward IMPACT3; 2× upsampling before tokenisation reproduces 40× predictions closely.
- Signature maps are relative within a slide (each slide is centred on its own mean embedding).
- Single-sample calls near class boundaries are uncertain. For cohort analyses, consider flagging slides with `til_margin` < 0.05 (the probability of TIL IV, `p_til_class_IV`, is available as a continuous measure).
