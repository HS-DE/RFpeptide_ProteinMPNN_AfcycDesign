#!/bin/bash
set -euo pipefail

if [[ "${RUN_STAGE5_PREDICTIONS:-NO}" != "YES" ]]; then
  echo "Stage 5 prediction is review-gated. Set RUN_STAGE5_PREDICTIONS=YES only after inspecting the manifest and protocol audit." >&2
  exit 3
fi

PYTHON_BIN="${AFCYCDESIGN_PYTHON:-$HOME/fga_model_envs/colabdesign-py310/.pixi/envs/default/bin/python}"
SOURCE_DIR="${COLABDESIGN_GAMMA_SOURCE:-$HOME/fga_model_envs/sources/ColabDesign-gamma-stage5}"
OVERLAY_DIR="${AFCYCDESIGN_PYTHON_OVERLAY:-$HOME/fga_model_envs/stage5_afcycdesign_python_overlay}"
AF_PARAMS_DIR="${AF_PARAMS:-$HOME/fga_model_envs/af_params}"
export PYTHONPATH="$SOURCE_DIR:$OVERLAY_DIR${PYTHONPATH:+:$PYTHONPATH}"

cd '/mnt/c/SH/fga_cyclic_peptide_design'
"$PYTHON_BIN" '/mnt/c/SH/fga_cyclic_peptide_design/scripts/external/run_afcycdesign_independent_recovery.py' \
  --job-spec '/mnt/c/SH/fga_cyclic_peptide_design/results/rfpeptides_article_route_clean_20260623_stage5_batch01_batch02/07_structure_validation/inputs/job_specs/S5_02_batch01_RFpep_Site_2_2193_seed01.json' \
  --af-params "$AF_PARAMS_DIR"
