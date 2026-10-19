# histoIMPACT

H&E → spatial tissue map → five tissue layers → IMPACT archetype (IMPACT1a/2, IMPACT1b, IMPACT3).

`code/run_histoimpact.py` runs the full pipeline on H&E images:

1. **Spatial tissue map** (`he_pred.py`): UNI token features (16 px = 4 µm at 0.25 µm/px), tissue mask, and per-token predictions of
   - the cell lineage (7 classes) and
   - the tumor region (tumor / non-tumor).
2. **Archetyping** (`impact_archetype.py`):
   - Five tissue layers by signed distance to the tumor border: TR < −100 µm, IT −100–0, IN 0–100, NP 100–300, ND ≥ 300 µm.
   - Visium-supervised CD4 T / CD8 T / B / macrophage / fibroblast maps.
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



## Usage

```bash
cd code
# images at 40x (0.25 um/px)
python run_histoimpact.py --image /path/to/image.tif --mpp 0.25
# images at another resolution are resampled to 0.25 um/px first
python run_histoimpact.py --image slide_20x.tif --mpp 0.5
# a cohort: manifest with sample_id, image, mpp
python run_histoimpact.py --manifest samples.csv --out-root ../results_cohort
```
