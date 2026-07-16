#!/bin/bash
set -euo pipefail

PYTHON_BIN="${AFCYCDESIGN_PYTHON:-$HOME/fga_model_envs/colabdesign-py310/.pixi/envs/default/bin/python}"
SOURCE_DIR="${COLABDESIGN_GAMMA_SOURCE:-$HOME/fga_model_envs/sources/ColabDesign-gamma-stage5}"
OVERLAY_DIR="${AFCYCDESIGN_PYTHON_OVERLAY:-$HOME/fga_model_envs/stage5_afcycdesign_python_overlay}"
AF_PARAMS_DIR="${AF_PARAMS:-$HOME/fga_model_envs/af_params}"

echo "[Stage5 preflight] python: $PYTHON_BIN"
echo "[Stage5 preflight] source: $SOURCE_DIR"
echo "[Stage5 preflight] overlay: $OVERLAY_DIR"
echo "[Stage5 preflight] params: $AF_PARAMS_DIR"

if [[ ! -x "$PYTHON_BIN" ]]; then
  echo "FAIL: AfCycDesign Python is missing or not executable." >&2
  exit 2
fi
if [[ ! -d "$AF_PARAMS_DIR" ]]; then
  echo "FAIL: AlphaFold parameter directory is missing." >&2
  exit 2
fi
if [[ ! -d "$OVERLAY_DIR" ]]; then
  echo "FAIL: Stage 5 Python overlay is missing. Run prepare_stage5_afcycdesign_python_overlay.sh first." >&2
  exit 2
fi
if [[ ! -f "$SOURCE_DIR/colabdesign/af/contrib/predict.py" ]] || [[ ! -f "$SOURCE_DIR/colabdesign/af/contrib/cyclic.py" ]]; then
  echo "FAIL: pinned ColabDesign gamma prediction source is missing. Run prepare_pinned_colabdesign_gamma_source.sh first." >&2
  exit 2
fi

COMMIT_MARKER="$SOURCE_DIR/.stage5_colabdesign_commit"
if [[ ! -f "$COMMIT_MARKER" ]]; then
  echo "FAIL: Stage 5 ColabDesign commit marker is missing." >&2
  exit 2
fi
current_commit="$(tr -d '[:space:]' < "$COMMIT_MARKER")"
if [[ "$current_commit" != "5ab4efaba2321a6c3c314b82d2fff8e0241f5c2d" ]]; then
  echo "FAIL: ColabDesign source commit is $current_commit, expected 5ab4efaba2321a6c3c314b82d2fff8e0241f5c2d." >&2
  exit 2
fi

export PYTHONPATH="$SOURCE_DIR:$OVERLAY_DIR${PYTHONPATH:+:$PYTHONPATH}"

"$PYTHON_BIN" - <<'PY'
import importlib.util
import sys

required = [
    "colabdesign",
    "colabdesign.af.contrib.predict",
    "colabdesign.af.contrib.cyclic",
]
missing = [name for name in required if importlib.util.find_spec(name) is None]
if missing:
    print("FAIL: missing modules: " + ", ".join(missing), file=sys.stderr)
    raise SystemExit(3)

from colabdesign.af.contrib.cyclic import add_cyclic_offset
print("PASS: gamma prediction modules and add_cyclic_offset are importable.")
print("Protocol requirement: template_mode=none; use_initial_guess=false; cyclic chain index=1.")
PY

echo "PASS: Stage 5 AfCycDesign protocol preflight completed."
