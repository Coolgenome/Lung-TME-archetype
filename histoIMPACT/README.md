# histoIMPACT

Please download the demo H&E image here: https://zenodo.org/records/23483997

`code/run_histoimpact.py` runs the full pipeline on H&E images:

1. **Spatial tissue map** (`he_pred.py`): UNI token features, tissue mask, and per-token predictions of
   - spatial pathological features (7 classes) and
   - tumor region-associated tokens (tumor / non-tumor).
2. **Archetyping** (`impact_archetype.py`):
   - Five tissue layers relative to the tumor boundary: TR, IT, IN, NP, ND.
   - Major TME cell type signatures: CD4 T, CD8 T, B, macrophage, fibroblast.
   - Spatial patterns of TIL / Mac / Fib
   - IMPACT archetype prediction.

## Outputs

| File | Content |
|---|---|
| `cell_lineage_grid.npy`, `cell_lineage_map.png` | Spatial tissue map of pathological features |
| `tumor_region_grid.npy`, `tumor_region_map.png` | Tumor-region-associated tokens |
| `tissue_layers_token.npy`, `tissue_layers.png` | Five-layer segmentation on the token grid |
| `tissue_layers_superpixel.npy` | Tissue layer segmentation on the 48 µm super-pixel grid |
| `signature_maps.npz` (+ `.png` with `--plot-maps`) | CD4_T, CD8_T, B, macrophage, fibroblast signature maps |
| `features.json` | Moran's I of major TME cell types |
| `impact.json` | 3-class IMPACT prediction, TME spatial pattern classes |

Cohort summary: `results_histoimpact/histoimpact_predictions.csv`. Its key columns:

- `IMPACT_3class`: one of IMPACT1a/2, IMPACT1b, IMPACT3 or not_assessable.
- `til_margin`: gap between the two most likely TIL classes; below 0.05 means a borderline call.
- `tumor_layer_frac`


## Installation

```bash
bash setup_env.sh
conda activate histoimpact
```

`setup_env.sh`:

- creates the conda env `histoimpact` from `environment.yml`



## Usage

```bash
cd code
# single image
python run_histoimpact.py --image /path/to/image.tif --mpp 0.25
# a cohort: manifest with sample_id, image, mpp
python run_histoimpact.py --manifest samples.csv --out-root ../results_cohort
```
