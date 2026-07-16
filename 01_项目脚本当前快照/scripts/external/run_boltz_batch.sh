#!/usr/bin/env bash
set -u

# Boltz batch runner for the FGA cyclic peptide Boltz branch.
# Reads results/boltz_predictions/boltz_jobs.csv and dispatches jobs across GPUs.

PROJECT_DIR="${PROJECT_DIR:-/mnt/c/SH/fga_cyclic_peptide_design}"
BOLTZ_ENV="${BOLTZ_ENV:-$HOME/fga_model_envs/boltz2}"
BOLTZ_CACHE="${BOLTZ_CACHE:-$HOME/fga_model_envs/boltz_cache}"
JOBS_CSV="${JOBS_CSV:-results/boltz_predictions/boltz_jobs.csv}"
GPU_LIST="${GPU_LIST:-0}"
MAX_JOBS="${MAX_JOBS:-0}"
MAX_RETRIES="${MAX_RETRIES:-1}"
OUTPUT_FORMAT="${OUTPUT_FORMAT:-pdb}"
BOLTZ_MODEL="${BOLTZ_MODEL:-boltz2}"
BOLTZ_USE_MSA_SERVER="${BOLTZ_USE_MSA_SERVER:-auto}"
BOLTZ_USE_POTENTIALS="${BOLTZ_USE_POTENTIALS:-true}"
BOLTZ_OVERRIDE="${BOLTZ_OVERRIDE:-false}"
BOLTZ_EXTRA_ARGS="${BOLTZ_EXTRA_ARGS:-}"
LOG_ROOT="${LOG_ROOT:-$(dirname "$JOBS_CSV")/logs}"

cd "$PROJECT_DIR" || exit 2
mkdir -p "$LOG_ROOT" "$BOLTZ_CACHE"

if [ ! -f "$JOBS_CSV" ]; then
  echo "Missing jobs csv: $JOBS_CSV" >&2
  exit 2
fi

if [ ! -f "$BOLTZ_ENV/bin/activate" ]; then
  echo "Missing Boltz environment: $BOLTZ_ENV" >&2
  exit 2
fi

# shellcheck disable=SC1090
source "$BOLTZ_ENV/bin/activate"
export BOLTZ_CACHE
export PYTHONUNBUFFERED="${PYTHONUNBUFFERED:-1}"

IFS=',' read -r -a GPUS <<< "$GPU_LIST"
if [ "${#GPUS[@]}" -eq 0 ]; then
  echo "GPU_LIST is empty" >&2
  exit 2
fi

summary="$LOG_ROOT/boltz_batch_summary.tsv"
failures="$LOG_ROOT/boltz_batch_failures.tsv"
touch "$summary" "$failures"

boltz_output_complete() {
  local out_dir="$1"
  local confidence
  local structure

  [ -d "$out_dir" ] || return 1
  confidence=$(find "$out_dir" -type f \( -name 'confidence_*.json' -o -name '*confidence*.json' \) | head -n 1)
  structure=$(find "$out_dir" -type f \( -name '*_model_0.pdb' -o -name '*_model_0.cif' -o -name '*.pdb' -o -name '*.cif' \) | head -n 1)
  [ -n "$confidence" ] && [ -n "$structure" ]
}

job_lines() {
  python - "$JOBS_CSV" "$MAX_JOBS" "$@" <<'PY'
import csv
import sys

csv_path = sys.argv[1]
max_jobs = int(sys.argv[2])
selected_ids = set(sys.argv[3:])
emitted = 0

with open(csv_path, "r", encoding="utf-8-sig", newline="") as handle:
    for row in csv.DictReader(handle):
        job_id = row.get("boltz_job_id", "")
        if selected_ids and job_id not in selected_ids:
            continue
        if not job_id:
            continue
        print("\t".join([
            job_id,
            row.get("input_yaml", ""),
            row.get("output_dir", ""),
            row.get("seed", ""),
            row.get("msa_mode", ""),
        ]))
        emitted += 1
        if max_jobs > 0 and emitted >= max_jobs:
            break
PY
}

run_job() {
  local job_id="$1"
  local input_yaml="$2"
  local out_dir="$3"
  local seed="$4"
  local msa_mode="$5"
  local gpu="$6"
  local attempt="$7"
  local log
  local cmd
  local exit_code

  log="${LOG_ROOT}/${job_id}.attempt${attempt}.gpu${gpu}.log"

  if boltz_output_complete "$out_dir"; then
    echo -e "${job_id}\tseed=${seed}\tgpu=${gpu}\tskipped_existing" | tee -a "$summary"
    return 0
  fi

  cmd=(boltz predict "$input_yaml" --out_dir "$out_dir" --model "$BOLTZ_MODEL" --seed "$seed" --output_format "$OUTPUT_FORMAT" --write_full_pae)

  if [ "$BOLTZ_USE_MSA_SERVER" = "true" ] || { [ "$BOLTZ_USE_MSA_SERVER" = "auto" ] && [ "$msa_mode" = "server" ]; }; then
    cmd+=(--use_msa_server)
  fi
  if [ "$BOLTZ_USE_POTENTIALS" = "true" ]; then
    cmd+=(--use_potentials)
  fi
  if [ "$BOLTZ_OVERRIDE" = "true" ]; then
    cmd+=(--override)
  fi
  if [ -n "$BOLTZ_EXTRA_ARGS" ]; then
    # Intentional word splitting for user-provided extra CLI args.
    # shellcheck disable=SC2206
    extra=($BOLTZ_EXTRA_ARGS)
    cmd+=("${extra[@]}")
  fi

  echo "=== Running ${job_id} attempt=${attempt} gpu=${gpu} ==="
  CUDA_VISIBLE_DEVICES="$gpu" "${cmd[@]}" 2>&1 | tee "$log"
  exit_code=${PIPESTATUS[0]}

  if [ "$exit_code" -eq 0 ] && boltz_output_complete "$out_dir"; then
    echo -e "${job_id}\tseed=${seed}\tgpu=${gpu}\tpass\tattempt=${attempt}" | tee -a "$summary"
    return 0
  fi

  echo -e "${job_id}\tseed=${seed}\tgpu=${gpu}\tfail\texit=${exit_code}\tattempt=${attempt}\tlog=${log}" | tee -a "$failures"
  return "$exit_code"
}

running=0
slot=0

while IFS=$'\t' read -r job_id input_yaml out_dir seed msa_mode; do
  if [ -z "$job_id" ]; then
    continue
  fi

  gpu="${GPUS[$((slot % ${#GPUS[@]}))]}"
  slot=$((slot + 1))

  (
    ok=0
    for attempt in $(seq 1 "$MAX_RETRIES"); do
      if run_job "$job_id" "$input_yaml" "$out_dir" "$seed" "$msa_mode" "$gpu" "$attempt"; then
        ok=1
        break
      fi
      sleep 3
    done
    if [ "$ok" -ne 1 ]; then
      echo "WARNING: ${job_id} failed after ${MAX_RETRIES} attempt(s); continuing." >&2
    fi
  ) &

  running=$((running + 1))
  if [ "$running" -ge "${#GPUS[@]}" ]; then
    wait -n
    running=$((running - 1))
  fi
done < <(job_lines "$@")

wait

echo "Boltz batch complete."
echo "Summary: $summary"
echo "Failures: $failures"
