#!/usr/bin/env bash
set -u

# Faster ColabDesign batch runner.
# Runs several designs per Python/JAX process and can dispatch chunks across
# multiple GPUs. This is faster than run_colabdesign_safe_batch.sh, but a native
# JAX/CUDA crash loses the whole chunk instead of only one seed.

PROJECT_DIR="${PROJECT_DIR:-/mnt/c/SH/fga_cyclic_peptide_design}"
DESIGNS_TOTAL_PER_JOB="${DESIGNS_TOTAL_PER_JOB:-100}"
DESIGNS_PER_PROCESS="${DESIGNS_PER_PROCESS:-5}"
MAX_RETRIES="${MAX_RETRIES:-1}"
PSSM_ITERS="${PSSM_ITERS:-80}"
GREEDY_ITERS="${GREEDY_ITERS:-32}"
START_SEED_BASE="${START_SEED_BASE:-0}"
GPU_LIST="${GPU_LIST:-0}"

export LCF_PY="${LCF_PY:-$HOME/fga_model_envs/colabdesign-py310/.pixi/envs/default/bin/python}"
export AF_PARAMS="${AF_PARAMS:-$HOME/fga_model_envs/af_params}"
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
  set -- Patch_A_L12_colabdesign Patch_A_L14_colabdesign Patch_A_L16_colabdesign
fi

IFS=',' read -r -a GPUS <<< "$GPU_LIST"
if [ "${#GPUS[@]}" -eq 0 ]; then
  echo "GPU_LIST is empty" >&2
  exit 2
fi

summary="logs/colabdesign_chunk_batch_summary.tsv"
failures="logs/colabdesign_chunk_batch_failures.tsv"
touch "$summary" "$failures"

running=0
slot=0

run_chunk() {
  local job="$1"
  local seed="$2"
  local n_designs="$3"
  local gpu="$4"
  local attempt="$5"
  local tag
  local out_dir
  local candidates
  local log
  local exit_code

  tag=$(printf "%s_chunk_seed%04d_n%03d" "$job" "$seed" "$n_designs")
  out_dir="results/raw_designs/colabdesign_outputs/${tag}"
  candidates="${out_dir}/candidates.csv"
  log="logs/${tag}.attempt${attempt}.gpu${gpu}.log"

  if [ -s "$candidates" ] && [ "$(wc -l < "$candidates")" -gt 1 ]; then
    echo -e "${tag}\tseed=${seed}\tn=${n_designs}\tgpu=${gpu}\tskipped_existing" | tee -a "$summary"
    return 0
  fi

  echo "=== Running ${tag} attempt=${attempt} gpu=${gpu} ==="

  CUDA_VISIBLE_DEVICES="$gpu" \
  NUM_DESIGNS="$n_designs" \
  START_SEED="$seed" \
  JOB_ID="$tag" \
  OUT_DIR="$out_dir" \
  PSSM_ITERS="$PSSM_ITERS" \
  GREEDY_ITERS="$GREEDY_ITERS" \
  bash "results/raw_designs/colabdesign_jobs/${job}.sh" 2>&1 | tee "$log"
  exit_code=${PIPESTATUS[0]}

  if [ "$exit_code" -eq 0 ] && [ -s "$candidates" ] && [ "$(wc -l < "$candidates")" -gt 1 ]; then
    echo -e "${tag}\tseed=${seed}\tn=${n_designs}\tgpu=${gpu}\tpass\tattempt=${attempt}" | tee -a "$summary"
    return 0
  fi

  echo -e "${tag}\tseed=${seed}\tn=${n_designs}\tgpu=${gpu}\tfail\texit=${exit_code}\tattempt=${attempt}\tlog=${log}" | tee -a "$failures"
  return "$exit_code"
}

for JOB in "$@"; do
  job_script="results/raw_designs/colabdesign_jobs/${JOB}.sh"
  if [ ! -f "$job_script" ]; then
    echo "Missing job script: $job_script" | tee -a "$failures"
    continue
  fi

  offset=0
  while [ "$offset" -lt "$DESIGNS_TOTAL_PER_JOB" ]; do
    remaining=$((DESIGNS_TOTAL_PER_JOB - offset))
    n_designs="$DESIGNS_PER_PROCESS"
    if [ "$remaining" -lt "$n_designs" ]; then
      n_designs="$remaining"
    fi

    seed=$((START_SEED_BASE + offset))
    gpu="${GPUS[$((slot % ${#GPUS[@]}))]}"
    slot=$((slot + 1))

    (
      ok=0
      for ATTEMPT in $(seq 1 "$MAX_RETRIES"); do
        if run_chunk "$JOB" "$seed" "$n_designs" "$gpu" "$ATTEMPT"; then
          ok=1
          break
        fi
        sleep 3
      done
      if [ "$ok" -ne 1 ]; then
        tag=$(printf "%s_chunk_seed%04d_n%03d" "$JOB" "$seed" "$n_designs")
        echo "WARNING: ${tag} failed after ${MAX_RETRIES} attempt(s); continuing." >&2
      fi
    ) &

    running=$((running + 1))
    if [ "$running" -ge "${#GPUS[@]}" ]; then
      wait -n
      running=$((running - 1))
    fi

    offset=$((offset + n_designs))
  done
done

wait
echo "Chunk batch complete."
echo "Summary: $summary"
echo "Failures: $failures"
