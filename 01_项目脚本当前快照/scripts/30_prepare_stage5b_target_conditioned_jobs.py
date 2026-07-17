from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
from pathlib import Path
from typing import Any, Mapping

from common import assert_active_route_path, append_run_header, read_csv, resolve_path, rows_to_markdown, setup_logger, write_csv, write_markdown
from pdb_utils import parse_residues, residue_sequence


COLABDESIGN_GAMMA_COMMIT = "5ab4efaba2321a6c3c314b82d2fff8e0241f5c2d"
PROTOCOL_VERSION = "stage5B_v1_target_only_template_single_sequence_no_mlm_dropout"
LEGAL_AA = set("ACDEFGHIKLMNPQRSTVWY")

MANIFEST_FIELDS = [
    "stage5B_candidate_id",
    "source_stage5A_candidate_id",
    "peptide_sequence_hash",
    "target_sequence_hash",
    "target_template_sha1",
    "protocol_hash",
    "protocol_version",
    "selection_order",
    "batch",
    "backbone_id",
    "source_stage4_design_id",
    "peptide_sequence",
    "peptide_length",
    "site_label",
    "site_id",
    "reference_design_pdb",
    "staged_reference_design_pdb",
    "target_template_pdb",
    "target_chain",
    "peptide_chain",
    "ddg_proxy_no_repack",
    "ddg_proxy_per_peptide_residue",
    "ddg_proxy_per_target_contact",
    "ddg_proxy_per_site_contact",
    "num_target_contacts",
    "num_target_site_contacts",
    "num_hotspot_contacts",
    "macrocycle_terminal_cn_distance",
    "clash_status",
    "selection_reason",
    "validation_test_type",
    "target_msa_mode",
    "peptide_msa_mode",
    "use_mlm",
    "use_dropout",
    "template_mode",
    "template_sequence_masked",
    "template_sidechains_masked",
    "template_interchain_features_masked",
    "use_initial_guess",
    "cyclic_chain_index",
    "cyclic_topology_encoding",
    "requested_recycles",
    "forward_passes",
    "seeds_per_candidate",
    "models_per_seed",
    "prediction_job_count",
    "status",
]

JOB_FIELDS = [
    "stage5B_job_id",
    "stage5B_candidate_id",
    "peptide_sequence_hash",
    "target_sequence_hash",
    "target_template_sha1",
    "protocol_hash",
    "batch",
    "backbone_id",
    "seed",
    "target_sequence_length",
    "peptide_sequence",
    "peptide_length",
    "target_template_pdb",
    "reference_design_pdb_for_posthoc",
    "job_spec_json",
    "prediction_output_dir",
    "run_script",
    "validation_test_type",
    "template_mode",
    "use_initial_guess",
    "target_msa_mode",
    "peptide_msa_mode",
    "use_mlm",
    "use_dropout",
    "requested_recycles",
    "forward_passes",
    "status",
]


def _safe_token(value: Any) -> str:
    return "".join(ch if ch.isalnum() or ch in {"_", "-", "."} else "_" for ch in str(value)).strip("_") or "item"


def _sha1_text(value: str) -> str:
    return hashlib.sha1(value.encode("utf-8")).hexdigest()


def _sha1_file(path: Path) -> str:
    digest = hashlib.sha1()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _protocol_hash(payload: Mapping[str, Any]) -> str:
    canonical = json.dumps(dict(payload), sort_keys=True, separators=(",", ":"))
    return _sha1_text(canonical)[:8]


def _resolve_mixed_path(value: str | Path) -> Path:
    text = str(value).strip().replace("\\", "/")
    if os.name == "nt" and text.startswith("/mnt/") and len(text) > 7 and text[6] == "/":
        return Path(f"{text[5].upper()}:/{text[7:]}")
    if os.name != "nt" and len(text) >= 3 and text[1] == ":" and text[2] == "/":
        return Path(f"/mnt/{text[0].lower()}{text[2:]}")
    path = Path(text)
    return path if path.is_absolute() else resolve_path(path)


def _to_wsl_path(path: str | Path) -> str:
    text = str(path).replace("\\", "/")
    if len(text) >= 3 and text[1] == ":" and text[2] == "/":
        return f"/mnt/{text[0].lower()}{text[2:]}"
    return text


def _quote_bash(value: str | Path) -> str:
    return "'" + str(value).replace("'", "'\"'\"'") + "'"


def _write_text(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(content.rstrip() + "\n")


def _validate_sequence(value: Any, label: str, expected_length: int | None = None) -> str:
    sequence = str(value).strip().upper()
    if not sequence:
        raise RuntimeError(f"{label} is empty")
    invalid = sorted(set(sequence) - LEGAL_AA)
    if invalid:
        raise RuntimeError(f"{label} contains unsupported amino acids: {','.join(invalid)}")
    if expected_length is not None and len(sequence) != expected_length:
        raise RuntimeError(f"{label} length is {len(sequence)}, expected {expected_length}")
    return sequence


def _pdb_context(path: Path, expected_target_length: int | None = None) -> tuple[str, str]:
    chains = parse_residues(path)
    if "A" not in chains:
        raise RuntimeError(f"Target chain A is missing: {path}")
    target = _validate_sequence(residue_sequence(chains["A"]), f"target chain A in {path}", expected_target_length)
    peptide = _validate_sequence(residue_sequence(chains["B"]), f"peptide chain B in {path}") if "B" in chains else ""
    return target, peptide


def _target_template_sequence(path: Path, expected_length: int) -> str:
    chains = parse_residues(path)
    if set(chains) != {"A"}:
        raise RuntimeError(f"Target template must contain only chain A; observed chains={sorted(chains)} in {path}")
    residues = chains["A"]
    if len(residues) != expected_length:
        raise RuntimeError(f"Target template chain A has {len(residues)} residues, expected {expected_length}")
    missing_ca = [str(row.get("pdb_residue_number", "")) for row in residues if "CA" not in row.get("atoms", {})]
    if missing_ca:
        raise RuntimeError(f"Target template has residues without CA atoms: {','.join(missing_ca)}")
    return _validate_sequence(residue_sequence(residues), "target template sequence", expected_length)


def _load_mapping(path: Path, expected_length: int) -> tuple[list[int], list[int], str]:
    rows = read_csv(path)
    if len(rows) != expected_length:
        raise RuntimeError(f"Stage 0 mapping has {len(rows)} rows, expected {expected_length}: {path}")
    ordered = sorted(rows, key=lambda row: int(row["rfpeptides_residue_number"]))
    numbers = [int(row["rfpeptides_residue_number"]) for row in ordered]
    if numbers != list(range(1, expected_length + 1)):
        raise RuntimeError("Stage 0 target mapping is not contiguous from 1 to 86")
    three_to_one = {
        "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
        "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
        "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    }
    sequence = "".join(three_to_one.get(str(row.get("rfpeptides_residue_name", "")).upper(), "X") for row in ordered)
    sequence = _validate_sequence(sequence, "Stage 0 mapping sequence", expected_length)
    site = [int(row["rfpeptides_residue_number"]) for row in ordered if str(row.get("is_target_site_residue", "")).lower() == "true"]
    hotspots = [int(row["rfpeptides_residue_number"]) for row in ordered if str(row.get("is_selected_hotspot", "")).lower() == "true"]
    if not site or not hotspots:
        raise RuntimeError("Stage 0 mapping lacks Site_2 or hotspot residues")
    return site, hotspots, sequence


def _job_script(
    project_root: Path,
    runner: Path,
    spec: Path,
    python_bin: str,
    source_dir: str,
    overlay_dir: str,
    af_params: str,
) -> str:
    return f"""#!/bin/bash
set -euo pipefail

if [[ "${{RUN_STAGE5B_PREDICTIONS:-NO}}" != "YES" ]]; then
  echo "Stage 5B prediction is review-gated. Set RUN_STAGE5B_PREDICTIONS=YES after reviewing the manifest and preflight." >&2
  exit 3
fi

PYTHON_BIN="${{AFCYCDESIGN_PYTHON:-{python_bin}}}"
SOURCE_DIR="${{COLABDESIGN_GAMMA_SOURCE:-{source_dir}}}"
OVERLAY_DIR="${{AFCYCDESIGN_PYTHON_OVERLAY:-{overlay_dir}}}"
AF_PARAMS_DIR="${{AF_PARAMS:-{af_params}}}"
export COLABDESIGN_GAMMA_SOURCE="$SOURCE_DIR"
export PYTHONPATH="$SOURCE_DIR:$OVERLAY_DIR${{PYTHONPATH:+:$PYTHONPATH}}"

cd {_quote_bash(_to_wsl_path(project_root))}
"$PYTHON_BIN" {_quote_bash(_to_wsl_path(runner))} \
  --job-spec {_quote_bash(_to_wsl_path(spec))} \
  --af-params "$AF_PARAMS_DIR"
"""


def _preflight_script(
    project_root: Path,
    runner: Path,
    representative_specs: list[Path],
    python_bin: str,
    source_dir: str,
    overlay_dir: str,
    af_params: str,
) -> str:
    lines = [
        "#!/bin/bash",
        "set -euo pipefail",
        "",
        f'PYTHON_BIN="${{AFCYCDESIGN_PYTHON:-{python_bin}}}"',
        f'SOURCE_DIR="${{COLABDESIGN_GAMMA_SOURCE:-{source_dir}}}"',
        f'OVERLAY_DIR="${{AFCYCDESIGN_PYTHON_OVERLAY:-{overlay_dir}}}"',
        f'AF_PARAMS_DIR="${{AF_PARAMS:-{af_params}}}"',
        'echo "[Stage5B preflight] python: $PYTHON_BIN"',
        'echo "[Stage5B preflight] source: $SOURCE_DIR"',
        'echo "[Stage5B preflight] params: $AF_PARAMS_DIR"',
        '[[ -x "$PYTHON_BIN" ]] || { echo "FAIL: AfCycDesign Python is missing." >&2; exit 2; }',
        '[[ -d "$SOURCE_DIR" ]] || { echo "FAIL: pinned gamma source is missing." >&2; exit 2; }',
        '[[ -d "$OVERLAY_DIR" ]] || { echo "FAIL: Stage 5 Python overlay is missing." >&2; exit 2; }',
        '[[ -d "$AF_PARAMS_DIR" ]] || { echo "FAIL: AlphaFold parameter directory is missing." >&2; exit 2; }',
        '[[ -f "$SOURCE_DIR/.stage5_colabdesign_commit" ]] || { echo "FAIL: commit marker is missing." >&2; exit 2; }',
        f'[[ "$(tr -d \'[:space:]\' < "$SOURCE_DIR/.stage5_colabdesign_commit")" == "{COLABDESIGN_GAMMA_COMMIT}" ]] || {{ echo "FAIL: gamma commit mismatch." >&2; exit 2; }}',
        'export COLABDESIGN_GAMMA_SOURCE="$SOURCE_DIR"',
        'export PYTHONPATH="$SOURCE_DIR:$OVERLAY_DIR${PYTHONPATH:+:$PYTHONPATH}"',
        f"cd {_quote_bash(_to_wsl_path(project_root))}",
        f'"$PYTHON_BIN" {_quote_bash(_to_wsl_path(runner))} --audit-imports',
        "",
    ]
    for spec in representative_specs:
        lines.append(
            f'"$PYTHON_BIN" {_quote_bash(_to_wsl_path(runner))} --preflight-only '
            f'--job-spec {_quote_bash(_to_wsl_path(spec))} --af-params "$AF_PARAMS_DIR"'
        )
    lines.extend(["", 'echo "PASS: Stage 5B static and padded-template tensor preflight completed for all candidates."'])
    return "\n".join(lines) + "\n"


def _master_script(preflight: Path, job_scripts: list[Path]) -> str:
    lines = [
        "#!/bin/bash",
        "set -euo pipefail",
        "",
        'if [[ "${RUN_STAGE5B_PREDICTIONS:-NO}" != "YES" ]]; then',
        '  echo "Stage 5B full prediction is review-gated. No model will run." >&2',
        '  echo "After review, set RUN_STAGE5B_PREDICTIONS=YES for this command only." >&2',
        "  exit 3",
        "fi",
        "",
        f"bash {_quote_bash(_to_wsl_path(preflight))}",
        "",
    ]
    lines.extend(f"bash {_quote_bash(_to_wsl_path(path))}" for path in job_scripts)
    return "\n".join(lines) + "\n"


def _report(rows: list[Mapping[str, Any]], output_dir: Path, protocol_hash: str) -> str:
    columns = [
        "stage5B_candidate_id", "batch", "backbone_id", "peptide_sequence", "peptide_length",
        "ddg_proxy_no_repack", "num_target_site_contacts", "num_hotspot_contacts", "selection_reason",
    ]
    return f"""# FGA RFpeptides Stage 5B Target-Structure-Conditioned Recovery

Status: inputs and 25 review-gated seed jobs prepared; predictions have not been run.

This is `target_structure_conditioned_recovery`, not sequence-only independent
recovery and not experimental binding validation. The 86-aa Stage 0 target
structure is supplied as a target-only backbone template. No peptide template,
peptide coordinates, complete design complex, or initial guess is supplied to
the prediction runner. Stage 4 design PDBs are staged only for post-prediction
alignment and pose-recovery analysis.

```text
protocol_hash: {protocol_hash}
target_msa_mode: single_sequence
peptide_msa_mode: single_sequence
use_mlm: false
use_dropout: true
dropout_interpretation: stochastic_prediction_ensemble
template_mode: target_only
template_sequence_masked: true
template_sidechains_masked: true
template_interchain_features_masked: true
use_initial_guess: false
cyclic_chain_index: 1
requested_recycles: 6
forward_passes: 7
candidates: {len(rows)}
seed_jobs: {len(rows) * 5}
model_predictions_planned: {len(rows) * 25}
```

The target template is padded to the full target-plus-peptide model length.
Every peptide template coordinate and atom mask is zero. The master runner
executes a static and padded-template tensor preflight before starting any job.

## Candidates

{rows_to_markdown(rows, columns, "No candidates were prepared.")}

## Output

```text
{output_dir}
```

These are validation inputs, not final peptide candidates.
"""


def main() -> int:
    parser = argparse.ArgumentParser(description="Prepare Stage 5B target-structure-conditioned AfCycDesign jobs.")
    parser.add_argument("--stage5a-root", required=True)
    parser.add_argument("--stage0-root", required=True)
    parser.add_argument("--output-root", required=True)
    parser.add_argument("--candidate-count", type=int, default=5)
    parser.add_argument("--seeds-per-candidate", type=int, default=5)
    parser.add_argument("--models-per-seed", type=int, default=5)
    parser.add_argument("--recycles", type=int, default=6)
    parser.add_argument("--expected-target-length", type=int, default=86)
    parser.add_argument("--afcycdesign-python", default="$HOME/fga_model_envs/colabdesign-py310/.pixi/envs/default/bin/python")
    parser.add_argument("--colabdesign-source", default="$HOME/fga_model_envs/sources/ColabDesign-gamma-stage5")
    parser.add_argument("--python-overlay", default="$HOME/fga_model_envs/stage5_afcycdesign_python_overlay")
    parser.add_argument("--af-params", default="$HOME/fga_model_envs/af_params")
    args = parser.parse_args()

    logger = setup_logger("30_prepare_stage5b_target_conditioned_jobs")
    append_run_header(logger, "30_prepare_stage5b_target_conditioned_jobs.py")
    if args.candidate_count != 5 or args.seeds_per_candidate != 5 or args.models_per_seed != 5:
        raise RuntimeError("Stage 5B v1 is fixed at 5 candidates, 5 seeds, and 5 model parameter sets.")
    if args.recycles != 6 or args.expected_target_length != 86:
        raise RuntimeError("Stage 5B v1 is fixed at 6 recycles and an 86-aa target template.")

    project_root = resolve_path(".")
    stage5a_root = _resolve_mixed_path(args.stage5a_root)
    assert_active_route_path(stage5a_root, "Stage 30 Stage 5A root")
    stage5a_dir = stage5a_root / "07_structure_validation"
    source_manifest = stage5a_dir / "FGA_rfpeptides_stage5_candidate_manifest.csv"
    assert_active_route_path(source_manifest, "Stage 30 Stage 5A candidate manifest CSV")
    source_rows = read_csv(source_manifest)
    if len(source_rows) < args.candidate_count:
        raise RuntimeError(f"Stage 5A-v2 manifest has {len(source_rows)} rows, expected at least 5: {source_manifest}")
    source_rows = sorted(source_rows, key=lambda row: int(row.get("selection_order", 999999)))[: args.candidate_count]

    output_root = _resolve_mixed_path(args.output_root)
    assert_active_route_path(output_root, "Stage 30 output root", must_exist=False)
    output_dir = output_root / "07_structure_validation_target_conditioned"
    inputs_dir = output_dir / "inputs"
    target_dir = inputs_dir / "target"
    reference_dir = inputs_dir / "reference_design_poses"
    specs_dir = inputs_dir / "job_specs"
    fasta_dir = inputs_dir / "fastas"
    jobs_dir = output_dir / "jobs"
    predictions_dir = output_dir / "predictions"

    stage0_root = _resolve_mixed_path(args.stage0_root)
    assert_active_route_path(stage0_root, "Stage 30 Stage 0 root")
    stage0_dir = stage0_root / "00_target_inputs"
    source_target = stage0_dir / "RFpep_Site_2_target.pdb"
    assert_active_route_path(source_target, "Stage 30 Stage 0 target template")
    if not source_target.is_file():
        raise RuntimeError(f"Missing Stage 0 target template: {source_target}")
    target_template = target_dir / source_target.name
    target_template.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source_target, target_template)
    target_sequence = _target_template_sequence(target_template, args.expected_target_length)
    mapping_csv = stage0_dir / "RFpep_Site_2_crop_renumbering_mapping.csv"
    assert_active_route_path(mapping_csv, "Stage 30 Stage 0 mapping CSV")
    site2_indices, hotspot_indices, mapping_sequence = _load_mapping(
        mapping_csv,
        args.expected_target_length,
    )
    if mapping_sequence != target_sequence:
        raise RuntimeError("Stage 0 mapping sequence does not match target template chain A")

    target_sequence_hash = _sha1_text(target_sequence)[:8]
    target_template_sha1 = _sha1_file(target_template)
    protocol_payload = {
        "protocol_version": PROTOCOL_VERSION,
        "validation_test_type": "target_structure_conditioned_recovery",
        "colabdesign_commit": COLABDESIGN_GAMMA_COMMIT,
        "model_type": "alphafold2_multimer_v3",
        "target_template_sha1": target_template_sha1,
        "template_mode": "target_only",
        "use_templates": True,
        "use_batch_as_template": False,
        "template_sequence_masked": True,
        "template_sidechains_masked": True,
        "template_interchain_features_masked": True,
        "use_initial_guess": False,
        "target_msa_mode": "single_sequence",
        "peptide_msa_mode": "single_sequence",
        "use_mlm": False,
        "mlm_replace_fraction": 0.0,
        "use_dropout": True,
        "cyclic_chain_index": 1,
        "cyclic_topology_encoding": "peptide_chain_relative_position_cyclic_offset",
        "models_per_seed": args.models_per_seed,
        "requested_recycles": args.recycles,
    }
    protocol_hash = _protocol_hash(protocol_payload)
    _write_text(target_dir / "RFpep_Site_2_target.fasta", f">RFpep_Site_2_target\n{target_sequence}")

    runner = project_root / "scripts" / "external" / "run_afcycdesign_target_conditioned_recovery.py"
    if not runner.is_file():
        raise RuntimeError(f"Missing Stage 5B runner: {runner}")

    manifest_rows: list[dict[str, Any]] = []
    job_rows: list[dict[str, Any]] = []
    job_scripts: list[Path] = []
    representative_specs: list[Path] = []
    for order, source in enumerate(source_rows, start=1):
        sequence = _validate_sequence(source.get("peptide_sequence", ""), f"candidate {order} peptide sequence")
        recorded_length = int(source.get("peptide_length", len(sequence)))
        if recorded_length != len(sequence):
            raise RuntimeError(f"Candidate {order} manifest peptide length does not match its sequence")
        source_design = _resolve_mixed_path(source.get("design_pdb", ""))
        if not source_design.is_file():
            source_design = _resolve_mixed_path(source.get("staged_design_pdb", ""))
        if not source_design.is_file():
            raise RuntimeError(f"Missing Stage 4 reference design PDB for candidate {order}")
        assert_active_route_path(source_design, f"Stage 30 Stage 4 reference PDB for candidate {order}")
        design_target, design_peptide = _pdb_context(source_design, args.expected_target_length)
        if design_target != target_sequence or design_peptide != sequence:
            raise RuntimeError(f"Candidate {order} Stage 4 PDB sequences do not match target/manifest sequences")

        peptide_hash = _sha1_text(sequence)[:8]
        candidate_id = f"S5B_{order:02d}_{source['batch']}_{_safe_token(source['backbone_id'])}_seq{peptide_hash}"
        staged_reference = reference_dir / f"{candidate_id}_stage4_reference.pdb"
        staged_reference.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source_design, staged_reference)
        fasta = fasta_dir / f"{candidate_id}.fasta"
        _write_text(
            fasta,
            f">{candidate_id}|validation=target_structure_conditioned_recovery|template=target_only|initial_guess=false\n"
            f"{target_sequence}:({sequence})",
        )

        row = {
            "stage5B_candidate_id": candidate_id,
            "source_stage5A_candidate_id": source.get("stage5_candidate_id", ""),
            "peptide_sequence_hash": peptide_hash,
            "target_sequence_hash": target_sequence_hash,
            "target_template_sha1": target_template_sha1,
            "protocol_hash": protocol_hash,
            "protocol_version": PROTOCOL_VERSION,
            "selection_order": order,
            "batch": source.get("batch", ""),
            "backbone_id": source.get("backbone_id", ""),
            "source_stage4_design_id": source.get("source_stage4_design_id", ""),
            "peptide_sequence": sequence,
            "peptide_length": len(sequence),
            "site_label": source.get("site_label", ""),
            "site_id": source.get("site_id", ""),
            "reference_design_pdb": source_design,
            "staged_reference_design_pdb": staged_reference,
            "target_template_pdb": target_template,
            "target_chain": "A",
            "peptide_chain": "B",
            "ddg_proxy_no_repack": source.get("ddg_proxy_no_repack", ""),
            "ddg_proxy_per_peptide_residue": source.get("ddg_proxy_per_peptide_residue", ""),
            "ddg_proxy_per_target_contact": source.get("ddg_proxy_per_target_contact", ""),
            "ddg_proxy_per_site_contact": source.get("ddg_proxy_per_site_contact", ""),
            "num_target_contacts": source.get("num_target_contacts", ""),
            "num_target_site_contacts": source.get("num_target_site_contacts", ""),
            "num_hotspot_contacts": source.get("num_hotspot_contacts", ""),
            "macrocycle_terminal_cn_distance": source.get("macrocycle_terminal_cn_distance", ""),
            "clash_status": source.get("clash_status", ""),
            "selection_reason": source.get("selection_reason", ""),
            "validation_test_type": "target_structure_conditioned_recovery",
            "target_msa_mode": "single_sequence",
            "peptide_msa_mode": "single_sequence",
            "use_mlm": "false",
            "use_dropout": "true",
            "template_mode": "target_only",
            "template_sequence_masked": "true",
            "template_sidechains_masked": "true",
            "template_interchain_features_masked": "true",
            "use_initial_guess": "false",
            "cyclic_chain_index": 1,
            "cyclic_topology_encoding": "peptide_chain_relative_position_cyclic_offset",
            "requested_recycles": args.recycles,
            "forward_passes": args.recycles + 1,
            "seeds_per_candidate": args.seeds_per_candidate,
            "models_per_seed": args.models_per_seed,
            "prediction_job_count": args.seeds_per_candidate,
            "status": "prepared_not_run",
        }
        manifest_rows.append(row)

        for seed in range(args.seeds_per_candidate):
            job_id = f"{candidate_id}_tmpl{target_template_sha1[:8]}_prot{protocol_hash}_seed{seed:02d}"
            prediction_output = predictions_dir / candidate_id / f"protocol_{protocol_hash}" / f"seed_{seed:02d}"
            spec_path = specs_dir / f"{job_id}.json"
            spec = {
                "stage5B_job_id": job_id,
                "stage5B_candidate_id": candidate_id,
                "protocol_hash": protocol_hash,
                "protocol_version": PROTOCOL_VERSION,
                "validation_test_type": "target_structure_conditioned_recovery",
                "batch": source.get("batch", ""),
                "backbone_id": source.get("backbone_id", ""),
                "target_sequence": target_sequence,
                "target_sequence_length": len(target_sequence),
                "target_sequence_hash": target_sequence_hash,
                "peptide_sequence": sequence,
                "peptide_sequence_length": len(sequence),
                "peptide_sequence_hash": peptide_hash,
                "target_template_path": _to_wsl_path(target_template),
                "target_template_sha1": target_template_sha1,
                "target_template_chain": "A",
                "target_template_expected_coverage": len(target_sequence),
                "peptide_template_expected_coverage": 0,
                "reference_design_pdb_for_posthoc": _to_wsl_path(staged_reference),
                "reference_design_loaded_by_prediction_runner": False,
                "site2_target_indices_1based": site2_indices,
                "hotspot_target_indices_1based": hotspot_indices,
                "seed": seed,
                "models_per_seed": args.models_per_seed,
                "model_names_expected": [f"model_{index}_multimer_v3" for index in range(1, 6)],
                "requested_recycles": args.recycles,
                "forward_passes": args.recycles + 1,
                "target_msa_mode": "single_sequence",
                "peptide_msa_mode": "single_sequence",
                "msa_rows_input": 1,
                "use_mlm": False,
                "mlm_replace_fraction": 0.0,
                "use_dropout": True,
                "dropout_interpretation": "stochastic_prediction_ensemble",
                "template_mode": "target_only",
                "use_templates": True,
                "use_batch_as_template": False,
                "template_sequence_masked": True,
                "template_sidechains_masked": True,
                "template_interchain_features_masked": True,
                "template_padding_mode": "target_atoms_then_zero_masked_peptide_positions",
                "use_initial_guess": False,
                "cyclic_chain_index": 1,
                "cyclic_topology_encoding": "peptide_chain_relative_position_cyclic_offset",
                "prediction_output_dir": _to_wsl_path(prediction_output),
                "colabdesign_commit": COLABDESIGN_GAMMA_COMMIT,
            }
            _write_text(spec_path, json.dumps(spec, indent=2, sort_keys=True))
            if seed == 0:
                representative_specs.append(spec_path)

            run_script = jobs_dir / f"run_{job_id}.sh"
            _write_text(
                run_script,
                _job_script(
                    project_root,
                    runner,
                    spec_path,
                    args.afcycdesign_python,
                    args.colabdesign_source,
                    args.python_overlay,
                    args.af_params,
                ),
            )
            job_scripts.append(run_script)
            job_rows.append(
                {
                    "stage5B_job_id": job_id,
                    "stage5B_candidate_id": candidate_id,
                    "peptide_sequence_hash": peptide_hash,
                    "target_sequence_hash": target_sequence_hash,
                    "target_template_sha1": target_template_sha1,
                    "protocol_hash": protocol_hash,
                    "batch": source.get("batch", ""),
                    "backbone_id": source.get("backbone_id", ""),
                    "seed": seed,
                    "target_sequence_length": len(target_sequence),
                    "peptide_sequence": sequence,
                    "peptide_length": len(sequence),
                    "target_template_pdb": target_template,
                    "reference_design_pdb_for_posthoc": staged_reference,
                    "job_spec_json": spec_path,
                    "prediction_output_dir": prediction_output,
                    "run_script": run_script,
                    "validation_test_type": "target_structure_conditioned_recovery",
                    "template_mode": "target_only",
                    "use_initial_guess": "false",
                    "target_msa_mode": "single_sequence",
                    "peptide_msa_mode": "single_sequence",
                    "use_mlm": "false",
                    "use_dropout": "true",
                    "requested_recycles": args.recycles,
                    "forward_passes": args.recycles + 1,
                    "status": "prepared_not_run",
                }
            )

    manifest_csv = output_dir / "FGA_rfpeptides_stage5B_candidate_manifest.csv"
    jobs_csv = output_dir / "FGA_rfpeptides_stage5B_prediction_jobs.csv"
    preflight = jobs_dir / "check_stage5B_target_conditioned_protocol.sh"
    master = jobs_dir / "run_stage5B_target_conditioned_recovery_all.sh"
    write_csv(manifest_csv, manifest_rows, MANIFEST_FIELDS)
    write_csv(jobs_csv, job_rows, JOB_FIELDS)
    write_markdown(output_dir / "FGA_rfpeptides_stage5B_candidate_manifest.md", _report(manifest_rows, output_dir, protocol_hash))
    _write_text(
        preflight,
        _preflight_script(
            project_root,
            runner,
            representative_specs,
            args.afcycdesign_python,
            args.colabdesign_source,
            args.python_overlay,
            args.af_params,
        ),
    )
    _write_text(master, _master_script(preflight, job_scripts))

    logger.info("Stage 5B candidates prepared: %s", len(manifest_rows))
    logger.info("Stage 5B seed jobs prepared: %s", len(job_rows))
    logger.info("Planned model predictions: %s", len(job_rows) * args.models_per_seed)
    logger.info("Output directory: %s", output_dir)
    logger.info("No Stage 5B prediction was run.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
