#!/usr/bin/env bash
set -u

# Crash-tolerant ColabDesign batch runner.
# Runs one design per Python/JAX process, so a native JAX/CUDA segfault only
# loses one seed instead of the whole production batch.

PROJECT_DIR="${PROJECT_DIR:-/mnt/c/SH/fga_cyclic_peptide_design}"
DESIGNS_PER_JOB="${DESIGNS_PER_JOB:-100}"
MAX_RETRIES="${MAX_RETRIES:-1}"
PSSM_ITERS="${PSSM_ITERS:-80}"
GREEDY_ITERS="${GREEDY_ITERS:-32}"
START_SEED_BASE="${START_SEED_BASE:-0}"

export LCF_PY="${LCF_PY:-$HOME/fga_model_envs/colabdesign-py310/.pixi/envs/default/bin/python}"
export AF_PARAMS="${AF_PARAMS:-$HOME/fga_model_envs/af_params}"
export CUDA_VISIBLE_DEVICES="${CUDA_VISIBLE_DEVICES:-0}"
export CUDA_HOME="${CUDA_HOME:-/usr/local/cuda}"
export LD_LIBRARY_PATH="${CUDA_HOME}/lib64:${CUDA_HOME}/targets/x86_64-linux/lib:/usr/lib/x86_64-linux-gnu:/usr/lib/wsl/lib:${LD_LIBRARY_PATH:-}"
export XLA_PYTHON_CLIENT_PREALLOCATE="${XLA_PYTHON_CLIENT_PREALLOCATE:-false}"
export XLA_PYTHON_CLIENT_MEM_FRACTION="${XLA_PYTHON_CLIENT_MEM_FRACTION:-0.30}"
export XLA_PYTHON_CLIENT_ALLOCATOR="${XLA_PYTHON_CLIENT_ALLOCATOR:-platform}"
export XLA_FLAGS="${XLA_FLAGS:---xla_gpu_enable_triton_gemm=false --xla_gpu_autotune_level=0}"
export JAX_PLATFORMS="${JAX_PLATFORMS:-cuda}"
export PYTHONUNBUFFERED="${PYTHONUNBUFFERED:-1}"
export PYTHONFAULTHANDLER="${PYTHONFAULTHANDLER:-1}"
export TF_FORCE_UNIFIED_MEMORY="${TF_FORCE_UNIFIED_MEMORY:-0}"

cd "$PROJECT_DIR" || exit 2
mkdir -p logs

if [ "$#" -eq 0 ]; then
  set -- Patch_B_L12_colabdesign Patch_B_L14_colabdesign Patch_B_L16_colabdesign
fi

summary="logs/colabdesign_safe_batch_summary.tsv"
failures="logs/colabdesign_safe_batch_failures.tsv"
touch "$summary" "$failures"

for JOB in "$@"; do
  job_script="results/raw_designs/colabdesign_jobs/${JOB}.sh"
  if [ ! -f "$job_script" ]; then
    echo "Missing job script: $job_script" | tee -a "$failures"
    continue
  fi

  for IDX in $(seq 0 $((DESIGNS_PER_JOB - 1))); do
    SEED=$((START_SEED_BASE + IDX))
    TAG=$(printf "%s_seed%04d" "$JOB" "$SEED")
    OUT_DIR="results/raw_designs/colabdesign_outputs/${TAG}"
    CANDIDATES="${OUT_DIR}/candidates.csv"

    if [ -s "$CANDIDATES" ] && [ "$(wc -l < "$CANDIDATES")" -gt 1 ]; then
      echo -e "${TAG}\tseed=${SEED}\tskipped_existing" | tee -a "$summary"
      continue
    fi

    ok=0
    for ATTEMPT in $(seq 1 "$MAX_RETRIES"); do
      LOG="logs/${TAG}.attempt${ATTEMPT}.log"
      echo "=== Running ${TAG} attempt=${ATTEMPT} seed=${SEED} ==="

      NUM_DESIGNS=1 \
      START_SEED="$SEED" \
      JOB_ID="$TAG" \
      OUT_DIR="$OUT_DIR" \
      PSSM_ITERS="$PSSM_ITERS" \
      GREEDY_ITERS="$GREEDY_ITERS" \
      bash "$job_script" 2>&1 | tee "$LOG"
      exit_code=${PIPESTATUS[0]}

      if [ "$exit_code" -eq 0 ] && [ -s "$CANDIDATES" ] && [ "$(wc -l < "$CANDIDATES")" -gt 1 ]; then
        echo -e "${TAG}\tseed=${SEED}\tpass\tattempt=${ATTEMPT}" | tee -a "$summary"
        ok=1
        break
      fi

      echo -e "${TAG}\tseed=${SEED}\tfail\texit=${exit_code}\tattempt=${ATTEMPT}\tlog=${LOG}" | tee -a "$failures"
      sleep 3
    done

    if [ "$ok" -ne 1 ]; then
      echo "WARNING: ${TAG} failed after ${MAX_RETRIES} attempt(s); continuing." >&2
    fi
  done
done

echo "Safe batch complete."
echo "Summary: $summary"
echo "Failures: $failures"
