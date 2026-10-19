#!/usr/bin/env bash
# Set up the conda environment for histoIMPACT and check that the full
# pipeline (UNI + model 1 + archetyping) can run.
#
#   bash setup_env.sh                 create / update env "histoimpact"
#   ENV_NAME=myenv bash setup_env.sh  another env name
#   TORCH_CUDA=cu118 bash setup_env.sh  PyTorch CUDA build (default: cu121
#                                       when an NVIDIA GPU is found; cpu
#                                       otherwise; macOS uses the default
#                                       wheel with Apple-silicon MPS)
#   RECREATE=1 bash setup_env.sh      delete and rebuild the env
#
# Then:  conda activate histoimpact
#        cd code && python run_histoimpact.py --image ../HE_demo/ADC1.png --mpp 0.25
set -euo pipefail
cd "$(dirname "$0")"
ROOT=$(pwd)
ENV_NAME=${ENV_NAME:-histoimpact}

if ! command -v conda >/dev/null 2>&1; then
  echo "conda not found: install Miniconda / Miniforge first" >&2
  exit 1
fi
eval "$(conda shell.bash hook)"
SOLVER=conda
command -v mamba >/dev/null 2>&1 && SOLVER=mamba

# 1. conda env with the scientific stack
if conda env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
  if [ "${RECREATE:-0}" = 1 ]; then
    echo ">> removing env $ENV_NAME"
    conda env remove -y -n "$ENV_NAME"
    $SOLVER env create -n "$ENV_NAME" -f environment.yml
  else
    echo ">> updating env $ENV_NAME"
    $SOLVER env update -n "$ENV_NAME" -f environment.yml
  fi
else
  echo ">> creating env $ENV_NAME"
  $SOLVER env create -n "$ENV_NAME" -f environment.yml
fi
conda activate "$ENV_NAME"

# 2. PyTorch (build matching the hardware) and timm (UNI backbone)
if [ "$(uname -s)" = Darwin ]; then
  TORCH_INDEX=""
elif [ -n "${TORCH_CUDA:-}" ]; then
  TORCH_INDEX="--index-url https://download.pytorch.org/whl/${TORCH_CUDA}"
elif command -v nvidia-smi >/dev/null 2>&1 && nvidia-smi >/dev/null 2>&1; then
  TORCH_INDEX="--index-url https://download.pytorch.org/whl/cu121"
else
  TORCH_INDEX="--index-url https://download.pytorch.org/whl/cpu"
fi
echo ">> installing torch + torchvision ${TORCH_INDEX:-(default wheel)}"
# torch and torchvision must come from the same index in one command, or pip
# pairs them with mismatched versions ("operator torchvision::nms does not
# exist"). Remove any torchvision from an earlier run first: a newer PyPI
# build would otherwise count as already satisfied.
python -m pip uninstall -y torchvision >/dev/null 2>&1 || true
# shellcheck disable=SC2086
python -m pip install --upgrade "torch>=2.1" torchvision $TORCH_INDEX
# timm last; torchvision is already satisfied, so pip leaves it alone
python -m pip install --upgrade "timm>=0.9.16"

# 3. check: imports, device, model files
python - "$ROOT" <<'EOF'
import sys
from pathlib import Path
root = Path(sys.argv[1])
import numpy, scipy, pandas, PIL, matplotlib, skimage, sklearn, torch, timm
print(f'python {sys.version.split()[0]} | numpy {numpy.__version__} | '
      f'torch {torch.__version__} | timm {timm.__version__}')
dev = ('cuda' if torch.cuda.is_available() else
       'mps' if getattr(torch.backends, 'mps', None)
       and torch.backends.mps.is_available() else 'cpu')
print('device available:', dev + (f' ({torch.cuda.get_device_name(0)})'
                                  if dev == 'cuda' else ''))
need = {'UNI weights': root / 'UNI' / 'pytorch_model.bin',
        'model 1, cell lineage': root / 'model' / 'general_uni_cls_cell_lineage.pt',
        'model 1, tumor region': root / 'model' / 'general_uni_cls_tumor_region.pt',
        'signature model': root / 'model' / 'signature_model.npz',
        'IMPACT classifier': root / 'model' / 'impact_classifier.json'}
missing = [k for k, p in need.items() if not p.exists()]
for k, p in need.items():
    print(f'  {"ok" if p.exists() else "MISSING":8s} {k}: {p.relative_to(root)}')
sys.path.insert(0, str(root / 'code'))
import impact_archetype as ia
ia.load_signature_model(need['signature model']) if need['signature model'].exists() else None
ia.load_impact_classifier(need['IMPACT classifier']) if need['IMPACT classifier'].exists() else None
if dev == 'cpu':
    print('note: no GPU found; UNI feature extraction will be slow on CPU')
if missing:
    print('missing files:', ', '.join(missing))
    if 'UNI weights' in missing:
        print('  UNI: request access at https://huggingface.co/MahmoodLab/UNI and '
              'save pytorch_model.bin in UNI/')
    sys.exit(1)
print('histoIMPACT environment ready')
EOF

echo
echo "Done. Run:"
echo "  conda activate $ENV_NAME"
echo "  cd code && python run_histoimpact.py --image ../HE_demo/ADC1.png --mpp 0.25"
