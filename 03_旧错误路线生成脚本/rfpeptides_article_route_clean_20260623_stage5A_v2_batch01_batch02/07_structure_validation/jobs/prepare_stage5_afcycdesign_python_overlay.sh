#!/bin/bash
set -euo pipefail

PYTHON_BIN="${AFCYCDESIGN_PYTHON:-$HOME/fga_model_envs/colabdesign-py310/.pixi/envs/default/bin/python}"
OVERLAY_DIR="${AFCYCDESIGN_PYTHON_OVERLAY:-$HOME/fga_model_envs/stage5_afcycdesign_python_overlay}"

if [[ ! -x "$PYTHON_BIN" ]]; then
  echo "FAIL: AfCycDesign Python is missing or not executable." >&2
  exit 2
fi

mkdir -p "$OVERLAY_DIR"
if PYTHONPATH="$OVERLAY_DIR" "$PYTHON_BIN" -c 'import IPython' >/dev/null 2>&1; then
  echo "PASS: Stage 5 Python overlay already provides IPython."
  exit 0
fi

"$PYTHON_BIN" -m pip install --upgrade --target "$OVERLAY_DIR" "IPython==8.37.0"
PYTHONPATH="$OVERLAY_DIR" "$PYTHON_BIN" -c 'import IPython; print("Prepared Stage 5 IPython overlay:", IPython.__version__)'
