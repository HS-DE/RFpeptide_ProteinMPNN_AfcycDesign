#!/bin/bash
set -euo pipefail

PYTHON_BIN="${AFCYCDESIGN_PYTHON:-$HOME/fga_model_envs/colabdesign-py310/.pixi/envs/default/bin/python}"
SOURCE_DIR="${COLABDESIGN_GAMMA_SOURCE:-$HOME/fga_model_envs/sources/ColabDesign-gamma-stage5}"
OVERLAY_DIR="${AFCYCDESIGN_PYTHON_OVERLAY:-$HOME/fga_model_envs/stage5_afcycdesign_python_overlay}"
AF_PARAMS_DIR="${AF_PARAMS:-$HOME/fga_model_envs/af_params}"
echo "[Stage5B preflight] python: $PYTHON_BIN"
echo "[Stage5B preflight] source: $SOURCE_DIR"
echo "[Stage5B preflight] params: $AF_PARAMS_DIR"
[[ -x "$PYTHON_BIN" ]] || { echo "FAIL: AfCycDesign Python is missing." >&2; exit 2; }
[[ -d "$SOURCE_DIR" ]] || { echo "FAIL: pinned gamma source is missing." >&2; exit 2; }
[[ -d "$OVERLAY_DIR" ]] || { echo "FAIL: Stage 5 Python overlay is missing." >&2; exit 2; }
[[ -d "$AF_PARAMS_DIR" ]] || { echo "FAIL: AlphaFold parameter directory is missing." >&2; exit 2; }
[[ -f "$SOURCE_DIR/.stage5_colabdesign_commit" ]] || { echo "FAIL: commit marker is missing." >&2; exit 2; }
[[ "$(tr -d '[:space:]' < "$SOURCE_DIR/.stage5_colabdesign_commit")" == "5ab4efaba2321a6c3c314b82d2fff8e0241f5c2d" ]] || { echo "FAIL: gamma commit mismatch." >&2; exit 2; }
export COLABDESIGN_GAMMA_SOURCE="$SOURCE_DIR"
export PYTHONPATH="$SOURCE_DIR:$OVERLAY_DIR${PYTHONPATH:+:$PYTHONPATH}"
cd '/mnt/c/SH/fga_cyclic_peptide_design'
"$PYTHON_BIN" '/mnt/c/SH/fga_cyclic_peptide_design/scripts/external/run_afcycdesign_target_conditioned_recovery.py' --audit-imports

"$PYTHON_BIN" '/mnt/c/SH/fga_cyclic_peptide_design/scripts/external/run_afcycdesign_target_conditioned_recovery.py' --preflight-only --job-spec '/mnt/c/SH/fga_cyclic_peptide_design/results/rfpeptides_article_route_clean_20260623_stage5B_batch01_batch02/07_structure_validation_target_conditioned/inputs/job_specs/S5B_01_batch01_RFpep_Site_2_3254_seq1547430a_tmplb925d32b_prot0e506709_seed00.json' --af-params "$AF_PARAMS_DIR"
"$PYTHON_BIN" '/mnt/c/SH/fga_cyclic_peptide_design/scripts/external/run_afcycdesign_target_conditioned_recovery.py' --preflight-only --job-spec '/mnt/c/SH/fga_cyclic_peptide_design/results/rfpeptides_article_route_clean_20260623_stage5B_batch01_batch02/07_structure_validation_target_conditioned/inputs/job_specs/S5B_02_batch01_RFpep_Site_2_2193_seq4e9c5740_tmplb925d32b_prot0e506709_seed00.json' --af-params "$AF_PARAMS_DIR"
"$PYTHON_BIN" '/mnt/c/SH/fga_cyclic_peptide_design/scripts/external/run_afcycdesign_target_conditioned_recovery.py' --preflight-only --job-spec '/mnt/c/SH/fga_cyclic_peptide_design/results/rfpeptides_article_route_clean_20260623_stage5B_batch01_batch02/07_structure_validation_target_conditioned/inputs/job_specs/S5B_03_batch01_RFpep_Site_2_9132_seqff47d13b_tmplb925d32b_prot0e506709_seed00.json' --af-params "$AF_PARAMS_DIR"
"$PYTHON_BIN" '/mnt/c/SH/fga_cyclic_peptide_design/scripts/external/run_afcycdesign_target_conditioned_recovery.py' --preflight-only --job-spec '/mnt/c/SH/fga_cyclic_peptide_design/results/rfpeptides_article_route_clean_20260623_stage5B_batch01_batch02/07_structure_validation_target_conditioned/inputs/job_specs/S5B_04_batch02_RFpep_Site_2_1518_seqb7139021_tmplb925d32b_prot0e506709_seed00.json' --af-params "$AF_PARAMS_DIR"
"$PYTHON_BIN" '/mnt/c/SH/fga_cyclic_peptide_design/scripts/external/run_afcycdesign_target_conditioned_recovery.py' --preflight-only --job-spec '/mnt/c/SH/fga_cyclic_peptide_design/results/rfpeptides_article_route_clean_20260623_stage5B_batch01_batch02/07_structure_validation_target_conditioned/inputs/job_specs/S5B_05_batch02_RFpep_Site_2_5529_seq5c351973_tmplb925d32b_prot0e506709_seed00.json' --af-params "$AF_PARAMS_DIR"

echo "PASS: Stage 5B static and padded-template tensor preflight completed for all candidates."
