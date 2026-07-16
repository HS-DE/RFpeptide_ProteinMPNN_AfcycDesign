#!/bin/bash
set -eo pipefail

cd '/home/luomi/fga_model_envs/rfpeptides/RFdiffusion'
source ~/fga_model_envs/miniforge3/etc/profile.d/conda.sh
conda activate 'SE3nv'
set -u
mkdir -p '/mnt/c/SH/fga_cyclic_peptide_design/results/rfpeptides_article_route_clean_20260716_runtimefix_smoke_N1_v3/02_rfpeptides_backbones/RFpep_Site_2'
mkdir -p '/mnt/c/SH/fga_cyclic_peptide_design/results/rfpeptides_article_route_clean_20260716_runtimefix_smoke_N1_v3/logs'
exec > >(tee '/mnt/c/SH/fga_cyclic_peptide_design/results/rfpeptides_article_route_clean_20260716_runtimefix_smoke_N1_v3/logs/run_rfpeptides_stage1_site2_runtimefix_smoke_N1_v3.complete.log') 2>&1
export RFPEPTIDES_RUNTIME_ROOT='/home/luomi/fga_model_envs/rfpeptides/RFdiffusion'
export RFPEPTIDES_REQUIRE_PROVENANCE_CLOSURE='1'
export PYTHONPATH="${RFPEPTIDES_RUNTIME_ROOT}${PYTHONPATH:+:${PYTHONPATH}}"
export EXPECTED_RFPEPTIDES_COMMIT='2d0c003df46b9db41d119321f15403dec3716cd9'
export EXPECTED_INFERENCE_UTILS_SHA256='7a3f4deeffd679fb95331806c26b4e4d476932a5c309ffbbd53316c044c1addf'
export EXPECTED_MODEL_RUNNERS_SHA256='ffebaf0ffedc2710ed9ef983279a6da97f4f5465591cd9adf8d79e86bdca8e08'
export EXPECTED_RUN_INFERENCE_SHA256='d1432bbba91e646531213f1aa39152230cba230aba649616fb5718a829c59123'
export EXPECTED_UTIL_SHA256='688dca4c8f07ef480f117b8cae1b3fa692a1d2f3b299ef1217918db5a3dc6908'

python - <<'PY'
import csv
import hashlib
import json
import os
from pathlib import Path
import subprocess
from types import SimpleNamespace
import rfdiffusion.util as rfd_util
from rfdiffusion.contigs import ContigMap
import rfdiffusion.inference.model_runners as model_runners
import rfdiffusion.inference.utils as inference_utils
from rfdiffusion.inference.utils import get_idx0_hotspots, process_target
job_specs = json.loads('[{"contig": "17-17 A1-86/0", "hotspots": ["A82", "A84", "A85", "A86"], "hotspots_normalized": "A82,A84,A85,A86", "hotspots_sha256": "85ad077878621b8709a23cb8d0c86fe201d8889ae3645ccb9909c5912f05f3dc", "job_id": "RFpep_Site_2_L17_17_N1", "mapping_csv": "/mnt/c/SH/fga_cyclic_peptide_design/results/rfpeptides_article_route_clean_20260615_fpocket/00_target_inputs/RFpep_Site_2_crop_renumbering_mapping.csv", "mapping_csv_sha256": "b407e30325721318f18d1bd9dbbbc4ebe65f225916ab7c7eb436c71aa9bc5cd3", "target_pdb": "/mnt/c/SH/fga_cyclic_peptide_design/results/rfpeptides_article_route_clean_20260615_fpocket/00_target_inputs/RFpep_Site_2_target.pdb", "target_pdb_sha256": "4a2253962e5b380215180c2f463af9d13ffaf2bd9470f092a19e36b6fce87ccd"}]')
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
for job in job_specs:
    target_pdb = Path(job['target_pdb']).resolve()
    mapping_csv = Path(job['mapping_csv']).resolve()
    if hashlib.sha256(target_pdb.read_bytes()).hexdigest() != job['target_pdb_sha256']:
        raise RuntimeError(f"{job['job_id']}: Stage 0 target PDB hash mismatch")
    if hashlib.sha256(mapping_csv.read_bytes()).hexdigest() != job['mapping_csv_sha256']:
        raise RuntimeError(f"{job['job_id']}: Stage 0 mapping CSV hash mismatch")
    if hashlib.sha256(job['hotspots_normalized'].encode('utf-8')).hexdigest() != job['hotspots_sha256']:
        raise RuntimeError(f"{job['job_id']}: normalized hotspot hash mismatch")
    with mapping_csv.open('r', encoding='utf-8-sig', newline='') as handle:
        mapping_rows = list(csv.DictReader(handle))
    mapping_rows.sort(key=lambda row: int(row['rfpeptides_residue_number']))
    if not mapping_rows:
        raise RuntimeError(f"{job['job_id']}: Stage 0 mapping CSV is empty")
    expected_refs = [
        (str(row['rfpeptides_chain']), int(row['rfpeptides_residue_number']))
        for row in mapping_rows
    ]
    if len(set(expected_refs)) != len(expected_refs):
        raise RuntimeError(f"{job['job_id']}: duplicate Stage 0 mapping residues")
    mapping_hotspots = sorted(
        f"{row['rfpeptides_chain']}{int(row['rfpeptides_residue_number'])}"
        for row in mapping_rows
        if str(row.get('is_selected_hotspot', '')).strip().lower() == 'true'
    )
    if mapping_hotspots != sorted(job['hotspots']):
        raise RuntimeError(
            f"{job['job_id']}: Stage 0 mapping hotspots differ from the job: "
            f"mapping={mapping_hotspots} job={sorted(job['hotspots'])}"
        )
    target_feats = process_target(str(target_pdb), parse_hetatom=True, center=False)
    target_pdb_idx = [(str(chain), int(number)) for chain, number in target_feats['pdb_idx']]
    if target_pdb_idx != expected_refs:
        raise RuntimeError(
            f"{job['job_id']}: real Stage 0 target PDB residue order differs from mapping CSV"
        )
    contig_map = ContigMap(target_feats, contigs=[job['contig']])
    mappings = contig_map.get_mappings()
    binderlen = len(contig_map.inpaint)
    receptor_refs = [(str(chain), int(number)) for chain, number in mappings['receptor_con_ref_pdb_idx']]
    receptor_hal = [int(value) for value in mappings['receptor_con_hal_idx0']]
    complex_refs = [(str(chain), int(number)) for chain, number in mappings['complex_con_ref_pdb_idx']]
    complex_hal = [int(value) for value in mappings['complex_con_hal_idx0']]
    expected_receptor_hal = list(range(len(expected_refs)))
    expected_complex_hal = list(range(binderlen, binderlen + len(expected_refs)))
    if receptor_refs != expected_refs or complex_refs != expected_refs:
        raise RuntimeError(f"{job['job_id']}: real ContigMap target references differ from Stage 0 mapping")
    if receptor_hal != expected_receptor_hal:
        raise RuntimeError(f"{job['job_id']}: receptor-local ContigMap indices are not contiguous from zero")
    if complex_hal != expected_complex_hal:
        raise RuntimeError(f"{job['job_id']}: complex-global ContigMap indices do not include binderlen")
    helper_indices = sorted(
        int(value)
        for value in (get_idx0_hotspots(mappings, SimpleNamespace(hotspot_res=job['hotspots']), binderlen) or [])
    )
    hotspot_set = set(job['hotspots'])
    expected_indices = sorted(
        binderlen + index
        for index, (chain, number) in enumerate(expected_refs)
        if f'{chain}{number}' in hotspot_set
    )
    contig_indices = sorted(
        complex_hal[index]
        for index, (chain, number) in enumerate(complex_refs)
        if f'{chain}{number}' in hotspot_set
    )
    if len(expected_indices) != len(hotspot_set):
        raise RuntimeError(f"{job['job_id']}: not all hotspots occur in the real Stage 0 mapping")
    if helper_indices != expected_indices or contig_indices != expected_indices:
        raise RuntimeError(
            f"{job['job_id']}: real ContigMap hotspot provenance mismatch: "
            f"helper={helper_indices} contig={contig_indices} expected={expected_indices}"
        )
    print('[RFpeptides preflight] real_contigmap=' + json.dumps({
        'job_id': job['job_id'],
        'binderlen': binderlen,
        'receptor_con_ref_pdb_idx': receptor_refs,
        'receptor_con_hal_idx0': receptor_hal,
        'complex_con_ref_pdb_idx': complex_refs,
        'complex_con_hal_idx0': complex_hal,
        'expected_hotspot_complex_global_indices': expected_indices,
    }, separators=(',', ':')))
print(f'[RFpeptides preflight] runtime={runtime_root}')
print(f'[RFpeptides preflight] commit={observed_commit}')
print('[RFpeptides preflight] module identity, hashes, Stage 0 inputs, and real ContigMap: PASS')
PY

echo '[RFpeptides] starting RFpep_Site_2_L17_17_N1'
export RFPEPTIDES_JOB_ID='RFpep_Site_2_L17_17_N1'
export STAGE0_TARGET_PDB_SHA256='4a2253962e5b380215180c2f463af9d13ffaf2bd9470f092a19e36b6fce87ccd'
export STAGE0_MAPPING_CSV='/mnt/c/SH/fga_cyclic_peptide_design/results/rfpeptides_article_route_clean_20260615_fpocket/00_target_inputs/RFpep_Site_2_crop_renumbering_mapping.csv'
export STAGE0_MAPPING_CSV_SHA256='b407e30325721318f18d1bd9dbbbc4ebe65f225916ab7c7eb436c71aa9bc5cd3'
export STAGE0_HOTSPOTS_NORMALIZED='A82,A84,A85,A86'
export STAGE0_HOTSPOTS_SHA256='85ad077878621b8709a23cb8d0c86fe201d8889ae3645ccb9909c5912f05f3dc'
{
python ./scripts/run_inference.py --config-name base \
inference.output_prefix=/mnt/c/SH/fga_cyclic_peptide_design/results/rfpeptides_article_route_clean_20260716_runtimefix_smoke_N1_v3/02_rfpeptides_backbones/RFpep_Site_2/RFpep_Site_2_L17_17 \
inference.num_designs=1 \
'contigmap.contigs=[17-17 A1-86/0]' \
inference.input_pdb=/mnt/c/SH/fga_cyclic_peptide_design/results/rfpeptides_article_route_clean_20260615_fpocket/00_target_inputs/RFpep_Site_2_target.pdb \
inference.cyclic=True \
diffuser.T=50 \
inference.cyc_chains='a' \
'ppi.hotspot_res=[A82,A84,A85,A86]' \
inference.write_trajectory=False
} 2>&1 | tee '/mnt/c/SH/fga_cyclic_peptide_design/results/rfpeptides_article_route_clean_20260716_runtimefix_smoke_N1_v3/logs/rfpeptides_stage1_RFpep_Site_2_L17_17_N1.log'
echo '[RFpeptides] finished RFpep_Site_2_L17_17_N1'
