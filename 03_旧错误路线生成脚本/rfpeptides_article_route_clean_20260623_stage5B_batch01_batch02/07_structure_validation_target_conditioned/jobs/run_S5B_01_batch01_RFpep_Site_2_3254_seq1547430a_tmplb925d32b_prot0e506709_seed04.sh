#!/bin/bash
set -euo pipefail

if [[ "${RUN_STAGE5B_PREDICTIONS:-NO}" != "YES" ]]; then
  echo "Stage 5B prediction is review-gated. Set RUN_STAGE5B_PREDICTIONS=YES after reviewing the manifest and preflight." >&2
  exit 3
fi

PYTHON_BIN="${AFCYCDESIGN_PYTHON:-$HOME/fga_model_envs/colabdesign-py310/.pixi/envs/default/bin/python}"
SOURCE_DIR="${COLABDESIGN_GAMMA_SOURCE:-$HOME/fga_model_envs/sources/ColabDesign-gamma-stage5}"
OVERLAY_DIR="${AFCYCDESIGN_PYTHON_OVERLAY:-$HOME/fga_model_envs/stage5_afcycdesign_python_overlay}"
AF_PARAMS_DIR="${AF_PARAMS:-$HOME/fga_model_envs/af_params}"
export COLABDESIGN_GAMMA_SOURCE="$SOURCE_DIR"
export PYTHONPATH="$SOURCE_DIR:$OVERLAY_DIR${PYTHONPATH:+:$PYTHONPATH}"

cd '/mnt/c/SH/fga_cyclic_peptide_design'
"$PYTHON_BIN" '/mnt/c/SH/fga_cyclic_peptide_design/scripts/external/run_afcycdesign_target_conditioned_recovery.py'   --job-spec '/mnt/c/SH/fga_cyclic_peptide_design/results/rfpeptides_article_route_clean_20260623_stage5B_batch01_batch02/07_structure_validation_target_conditioned/inputs/job_specs/S5B_01_batch01_RFpep_Site_2_3254_seq1547430a_tmplb925d32b_prot0e506709_seed04.json'   --af-params "$AF_PARAMS_DIR"
