#!/usr/bin/env bash
set -u

PROJECT_DIR="${PROJECT_DIR:-/mnt/c/SH/fga_cyclic_peptide_design}"
JOBS_CSV="${JOBS_CSV:-results/colabfold_predictions_top30_seed1/colabfold_jobs.csv}"
GPU_LIST="${GPU_LIST:-0}"
MAX_JOBS="${MAX_JOBS:-0}"
MAX_RETRIES="${MAX_RETRIES:-1}"
COLABFOLD_TIMEOUT_SECONDS="${COLABFOLD_TIMEOUT_SECONDS:-7200}"
COLABFOLD_OVERRIDE="${COLABFOLD_OVERRIDE:-false}"
CF_BIN="${CF_BIN:-$HOME/fga_model_envs/localcolabfold/.pixi/envs/default/bin/colabfold_batch}"
CF_PY="${CF_PY:-$(dirname "$CF_BIN")/python}"
AF_PARAMS="${AF_PARAMS:-$HOME/fga_model_envs/af_params}"
COLABFOLD_PREDOWNLOAD_PARAMS="${COLABFOLD_PREDOWNLOAD_PARAMS:-true}"
XLA_PYTHON_CLIENT_PREALLOCATE="${XLA_PYTHON_CLIENT_PREALLOCATE:-false}"
XLA_PYTHON_CLIENT_MEM_FRACTION="${XLA_PYTHON_CLIENT_MEM_FRACTION:-0.45}"
XLA_FLAGS="${XLA_FLAGS:---xla_gpu_enable_triton_gemm=false}"

cd "$PROJECT_DIR" || exit 2

if [ ! -x "$CF_BIN" ]; then
  echo "Missing executable ColabFold binary: $CF_BIN" >&2
  exit 2
fi
if [ ! -x "$CF_PY" ]; then
  echo "Missing executable ColabFold Python: $CF_PY" >&2
  exit 2
fi

if [ ! -f "$JOBS_CSV" ]; then
  echo "Missing JOBS_CSV: $JOBS_CSV" >&2
  exit 2
fi

LOG_ROOT="$(dirname "$JOBS_CSV")/logs"
mkdir -p "$LOG_ROOT"
summary="$LOG_ROOT/colabfold_batch_summary.tsv"
failures="$LOG_ROOT/colabfold_batch_failures.tsv"
touch "$summary" "$failures"

IFS=',' read -r -a GPUS <<< "$GPU_LIST"
if [ "${#GPUS[@]}" -eq 0 ]; then
  GPUS=(0)
fi

first_model_type="$("$CF_PY" - "$JOBS_CSV" <<'PY'
import csv
import sys

with open(sys.argv[1], newline="", encoding="utf-8-sig") as handle:
    for row in csv.DictReader(handle):
        print(row.get("model_type", "auto") or "auto")
        break
PY
)"

if [ "$COLABFOLD_PREDOWNLOAD_PARAMS" = "true" ] && [ -n "$first_model_type" ]; then
  echo "Checking ColabFold params for model_type=${first_model_type} in ${AF_PARAMS}"
  if ! "$CF_PY" - "$first_model_type" "$AF_PARAMS" <<'PY'
import sys
from pathlib import Path
from colabfold.download import download_alphafold_params

model_type = sys.argv[1]
data_dir = Path(sys.argv[2])
download_alphafold_params(model_type, data_dir)
print(f"ColabFold params ready: {model_type} -> {data_dir}")
PY
  then
    echo "ColabFold parameter check/download failed; stop before launching GPU jobs." >&2
    exit 2
  fi
fi

output_complete() {
  local out_dir="$1"
  find "$out_dir" -maxdepth 1 -type f \( -name '*scores_rank_001*.json' -o -name '*scores*.json' \) | grep -q . || return 1
  find "$out_dir" -maxdepth 1 -type f \( -name '*rank_001*.pdb' -o -name '*.pdb' \) | grep -q . || return 1
  return 0
}

job_stream() {
  "$CF_PY" - "$JOBS_CSV" "$MAX_JOBS" <<'PY'
import csv
import sys

path = sys.argv[1]
max_jobs = int(sys.argv[2])
emitted = 0
with open(path, newline="", encoding="utf-8-sig") as handle:
    for row in csv.DictReader(handle):
        if max_jobs and emitted >= max_jobs:
            break
        print("\t".join([
            row.get("colabfold_job_id", ""),
            row.get("input_fasta", ""),
            row.get("output_dir", ""),
            row.get("seed", ""),
            row.get("msa_mode", ""),
            row.get("model_type", ""),
            row.get("num_models", ""),
            row.get("num_recycle", ""),
        ]))
        emitted += 1
PY
}

run_job() {
  local job_id="$1"
  local input_fasta="$2"
  local out_dir="$3"
  local seed="$4"
  local msa_mode="$5"
  local model_type="$6"
  local num_models="$7"
  local num_recycle="$8"
  local gpu="$9"
  local attempt="${10}"

  input_fasta="${input_fasta//\\//}"
  out_dir="${out_dir//\\//}"
  mkdir -p "$out_dir"
  if [ "$COLABFOLD_OVERRIDE" != "true" ] && output_complete "$out_dir"; then
    echo -e "${job_id}\tseed=${seed}\tgpu=${gpu}\tskipped_existing" | tee -a "$summary"
    return 0
  fi

  local log="$LOG_ROOT/${job_id}.attempt${attempt}.gpu${gpu}.log"
  echo "=== Running ${job_id} attempt=${attempt} gpu=${gpu} ==="
  (
    export CUDA_VISIBLE_DEVICES="$gpu"
    export XLA_PYTHON_CLIENT_PREALLOCATE
    export XLA_PYTHON_CLIENT_MEM_FRACTION
    export XLA_FLAGS
    timeout "$COLABFOLD_TIMEOUT_SECONDS" "$CF_BIN" \
      --msa-mode "$msa_mode" \
      --model-type "$model_type" \
      --num-models "$num_models" \
      --num-recycle "$num_recycle" \
      --num-seeds 1 \
      --random-seed "$seed" \
      --data "$AF_PARAMS" \
      --overwrite-existing-results \
      "$input_fasta" "$out_dir"
  ) 2>&1 | tee "$log"
  local exit_code=${PIPESTATUS[0]}

  if [ "$exit_code" -eq 0 ] && output_complete "$out_dir"; then
    echo -e "${job_id}\tseed=${seed}\tgpu=${gpu}\tpass\tattempt=${attempt}" | tee -a "$summary"
    return 0
  fi
  echo -e "${job_id}\tseed=${seed}\tgpu=${gpu}\tfail\texit=${exit_code}\tattempt=${attempt}\tlog=${log}" | tee -a "$failures"
  return 1
}

slot=0
while IFS=$'\t' read -r job_id input_fasta out_dir seed msa_mode model_type num_models num_recycle; do
  if [ -z "$job_id" ]; then
    continue
  fi
  gpu="${GPUS[$((slot % ${#GPUS[@]}))]}"
  (
    ok=0
    for attempt in $(seq 1 "$MAX_RETRIES"); do
      if run_job "$job_id" "$input_fasta" "$out_dir" "$seed" "$msa_mode" "$model_type" "$num_models" "$num_recycle" "$gpu" "$attempt"; then
        ok=1
        break
      fi
      sleep 5
    done
    if [ "$ok" -ne 1 ]; then
      echo "WARNING: ${job_id} failed after ${MAX_RETRIES} attempt(s); continuing." >&2
    fi
  ) &
  slot=$((slot + 1))
  if [ $((slot % ${#GPUS[@]})) -eq 0 ]; then
    wait
  fi
done < <(job_stream)
wait

echo "ColabFold batch complete."
echo "Summary: $summary"
echo "Failures: $failures"
