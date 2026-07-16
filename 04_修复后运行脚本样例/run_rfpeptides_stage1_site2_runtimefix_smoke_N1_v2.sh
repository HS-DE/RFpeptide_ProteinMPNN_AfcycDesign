#!/bin/bash
set -eo pipefail

cd '/home/luomi/fga_model_envs/rfpeptides/RFdiffusion'
source ~/fga_model_envs/miniforge3/etc/profile.d/conda.sh
conda activate 'SE3nv'
set -u
export RFPEPTIDES_RUNTIME_ROOT='/home/luomi/fga_model_envs/rfpeptides/RFdiffusion'
export PYTHONPATH="${RFPEPTIDES_RUNTIME_ROOT}${PYTHONPATH:+:${PYTHONPATH}}"
export EXPECTED_RFPEPTIDES_COMMIT='2d0c003df46b9db41d119321f15403dec3716cd9'
export EXPECTED_INFERENCE_UTILS_SHA256='7a3f4deeffd679fb95331806c26b4e4d476932a5c309ffbbd53316c044c1addf'
export EXPECTED_MODEL_RUNNERS_SHA256='47bee8cab513f7b49159de3e819ce41b05c255584f7edbfdb79cede02054ae71'
export EXPECTED_RUN_INFERENCE_SHA256='7d91212d71558bb8cc35a89e5dab4f3282a4b27a04b4cae434c3f334bf254680'
export EXPECTED_UTIL_SHA256='688dca4c8f07ef480f117b8cae1b3fa692a1d2f3b299ef1217918db5a3dc6908'

python - <<'PY'
import hashlib
import os
from pathlib import Path
import subprocess
from types import SimpleNamespace
import rfdiffusion.util as rfd_util
import rfdiffusion.inference.model_runners as model_runners
import rfdiffusion.inference.utils as inference_utils
from rfdiffusion.inference.utils import get_idx0_hotspots
runtime_root = Path(os.environ['RFPEPTIDES_RUNTIME_ROOT']).resolve()
expected_paths = {
    'inference_utils': runtime_root / 'rfdiffusion/inference/utils.py',
    'model_runners': runtime_root / 'rfdiffusion/inference/model_runners.py',
    'run_inference': runtime_root / 'scripts/run_inference.py',
    'util': runtime_root / 'rfdiffusion/util.py',
}
loaded_paths = {
    'inference_utils': Path(inference_utils.__file__).resolve(),
    'model_runners': Path(model_runners.__file__).resolve(),
    'util': Path(rfd_util.__file__).resolve(),
}
for key, loaded_path in loaded_paths.items():
    if loaded_path != expected_paths[key]:
        raise RuntimeError(f'RFpeptides runtime split-brain: {key} loaded from {loaded_path}, expected {expected_paths[key]}')
expected_hashes = {
    'inference_utils': os.environ['EXPECTED_INFERENCE_UTILS_SHA256'],
    'model_runners': os.environ['EXPECTED_MODEL_RUNNERS_SHA256'],
    'run_inference': os.environ['EXPECTED_RUN_INFERENCE_SHA256'],
    'util': os.environ['EXPECTED_UTIL_SHA256'],
}
for key, path in expected_paths.items():
    observed_hash = hashlib.sha256(path.read_bytes()).hexdigest()
    if observed_hash != expected_hashes[key]:
        raise RuntimeError(f'RFpeptides runtime file changed: {path}; observed={observed_hash}, expected={expected_hashes[key]}')
observed_commit = subprocess.check_output(['git', '-C', str(runtime_root), 'rev-parse', 'HEAD'], text=True).strip()
if observed_commit != os.environ['EXPECTED_RFPEPTIDES_COMMIT']:
    raise RuntimeError(f'RFpeptides runtime commit changed: observed={observed_commit}, expected={os.environ["EXPECTED_RFPEPTIDES_COMMIT"]}')
mappings = {
    'receptor_con_ref_pdb_idx': [('A', 82), ('A', 84)],
    'receptor_con_hal_idx0': [81, 83],
}
observed = get_idx0_hotspots(mappings, SimpleNamespace(hotspot_res=['A82', 'A84']), binderlen=17)
expected = [98, 100]
if list(observed or []) != expected:
    raise RuntimeError(f'RFpeptides hotspot indexing preflight failed: observed={observed}, expected={expected}')
print(f'[RFpeptides preflight] runtime={runtime_root}')
print(f'[RFpeptides preflight] commit={observed_commit}')
print('[RFpeptides preflight] module identity, hashes, and hotspot binder-length offset: PASS')
PY

mkdir -p '/mnt/c/SH/fga_cyclic_peptide_design/results/rfpeptides_article_route_clean_20260716_runtimefix_smoke_N1_v2/02_rfpeptides_backbones/RFpep_Site_2'
mkdir -p '/mnt/c/SH/fga_cyclic_peptide_design/results/rfpeptides_article_route_clean_20260716_runtimefix_smoke_N1_v2/logs'

echo '[RFpeptides] starting RFpep_Site_2_L17_17_N1'
{
python ./scripts/run_inference.py --config-name base \
inference.output_prefix=/mnt/c/SH/fga_cyclic_peptide_design/results/rfpeptides_article_route_clean_20260716_runtimefix_smoke_N1_v2/02_rfpeptides_backbones/RFpep_Site_2/RFpep_Site_2_L17_17 \
inference.num_designs=1 \
'contigmap.contigs=[17-17 A1-86/0]' \
inference.input_pdb=/mnt/c/SH/fga_cyclic_peptide_design/results/rfpeptides_article_route_clean_20260615_fpocket/00_target_inputs/RFpep_Site_2_target.pdb \
inference.cyclic=True \
diffuser.T=50 \
inference.cyc_chains='a' \
'ppi.hotspot_res=[A85,A82,A86,A84]' \
inference.write_trajectory=False
} 2>&1 | tee '/mnt/c/SH/fga_cyclic_peptide_design/results/rfpeptides_article_route_clean_20260716_runtimefix_smoke_N1_v2/logs/rfpeptides_stage1_RFpep_Site_2_L17_17_N1.log'
echo '[RFpeptides] finished RFpep_Site_2_L17_17_N1'
