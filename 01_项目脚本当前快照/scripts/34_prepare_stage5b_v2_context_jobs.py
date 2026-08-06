from __future__ import annotations

import argparse
import hashlib
import json
import os
import shlex
import shutil
from pathlib import Path
from typing import Any, Iterable, Mapping

from common import (
    ROUTE_PROVENANCE_FIELDS,
    SOURCE_ROUTE_PROVENANCE_FIELDS,
    assert_active_route_path,
    append_run_header,
    canonical_json_sha256,
    load_route_manifest,
    read_csv,
    resolve_path,
    route_provenance_fields,
    setup_logger,
    sha256_file,
    validate_route_project_config,
    validate_row_route_provenance,
    write_csv,
    write_markdown,
    write_route_manifest,
)
from pdb_utils import parse_residues, residue_sequence
from stage5_contract import (
    STAGE4_IDENTITY_FIELDS,
    load_stage4_validation_contract,
    stage4_identity_values,
)


COLABDESIGN_GAMMA_COMMIT = "5ab4efaba2321a6c3c314b82d2fff8e0241f5c2d"
TOP5_PROTOCOL_VERSION = "stage5B_v2_C1_C3_full_target_template_top5_v1"
ALL_PASS_PROTOCOL_VERSION = "stage5B_v2_C1_C3_full_target_template_allpass_v1"
PROTOCOL_VERSION = TOP5_PROTOCOL_VERSION
VALIDATION_TEST_TYPE = "target_context_structure_conditioned_recovery"
MODEL_NAMES = [f"model_{index}_multimer_v3" for index in range(1, 6)]
LEGAL_AA = set("ACDEFGHIKLMNPQRSTVWY")
DEFAULT_CONTEXT_IDS = (
    "C1_crop86_full_template",
    "C3_native_GHI301_full_template",
)

MANIFEST_FIELDS = [
    "stage5B_v2_candidate_context_id",
    "stage5B_v2_candidate_id",
    "context_id",
    "context_description",
    "context_pdb",
    "context_pdb_sha256",
    "context_mapping_csv",
    "context_mapping_csv_sha256",
    "context_chain_ids",
    "context_chain_lengths",
    "context_total_length",
    "stage0_crop_context_indices_1based",
    "site2_fga_indices_1based",
    "hotspot_fga_indices_1based",
    "peptide_chain_id",
    "cyclic_chain_index",
    "peptide_sequence",
    "peptide_sequence_hash",
    "peptide_length",
    "reference_target_chain",
    "reference_peptide_chain",
    "staged_reference_design_pdb",
    "staged_reference_design_pdb_sha256",
    "selection_order",
    "batch",
    "backbone_id",
    *STAGE4_IDENTITY_FIELDS,
    "validation_test_type",
    "template_mode",
    "template_sequence_masked",
    "template_sidechains_masked",
    "template_interchain_features_masked",
    "peptide_template_expected_coverage",
    "use_initial_guess",
    "target_msa_mode",
    "peptide_msa_mode",
    "use_mlm",
    "use_dropout",
    "cyclic_topology_encoding",
    "requested_recycles",
    "forward_passes",
    "seeds_per_candidate_context",
    "models_per_seed",
    "protocol_hash",
    "stage5_campaign_id",
    "stage5_selection_mode",
    "status",
    "notes",
    *SOURCE_ROUTE_PROVENANCE_FIELDS,
    *ROUTE_PROVENANCE_FIELDS,
]

JOB_FIELDS = [
    "stage5B_v2_job_id",
    "stage5B_v2_candidate_context_id",
    "stage5B_v2_candidate_id",
    "context_id",
    "peptide_sequence_hash",
    "protocol_hash",
    "stage5_campaign_id",
    "stage5_selection_mode",
    "job_shard",
    "seed",
    "requested_recycles",
    "forward_passes",
    "models_per_seed",
    "job_spec_json",
    "job_spec_sha256",
    "prediction_output_dir",
    "run_script",
    "template_mode",
    "peptide_template_expected_coverage",
    "use_initial_guess",
    "cyclic_chain_index",
    "status",
    *SOURCE_ROUTE_PROVENANCE_FIELDS,
    *ROUTE_PROVENANCE_FIELDS,
]


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


def _write_text(path: Path, value: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(value.rstrip() + "\n", encoding="utf-8", newline="\n")


def _sha1_text(value: str) -> str:
    return hashlib.sha1(value.encode("utf-8")).hexdigest()


def _safe_token(value: Any) -> str:
    return "".join(ch if ch.isalnum() or ch in "_-" else "_" for ch in str(value)).strip("_")


def _truthy(value: Any) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes"}


def _split_csv(value: str) -> list[str]:
    return [item.strip() for item in str(value).split(",") if item.strip()]


def _validate_sequence(value: Any, label: str, expected_length: int | None = None) -> str:
    sequence = str(value).strip().upper()
    if not sequence:
        raise RuntimeError(f"{label} is empty")
    invalid = sorted(set(sequence) - LEGAL_AA)
    if invalid:
        raise RuntimeError(f"{label} contains unsupported amino acids: {','.join(invalid)}")
    if expected_length is not None and len(sequence) != expected_length:
        raise RuntimeError(f"{label} has length {len(sequence)}, expected {expected_length}")
    return sequence


def _reference_sequences(
    path: Path,
    target_chain: str,
    peptide_chain: str,
    expected_target: str,
    expected_peptide: str,
) -> None:
    chains = parse_residues(path)
    if set(chains) != {target_chain, peptide_chain} or target_chain == peptide_chain:
        raise RuntimeError(
            f"Stage 4 reference chains {list(chains)} do not match "
            f"target={target_chain}, peptide={peptide_chain}: {path}"
        )
    target = _validate_sequence(
        residue_sequence(chains[target_chain]), f"target chain {target_chain} in {path}"
    )
    peptide = _validate_sequence(
        residue_sequence(chains[peptide_chain]), f"peptide chain {peptide_chain} in {path}"
    )
    if target != expected_target or peptide != expected_peptide:
        raise RuntimeError(f"Stage 4 reference sequence mismatch: {path}")


def _context_contract(row: Mapping[str, str]) -> dict[str, Any]:
    context_id = str(row["context_id"]).strip()
    context_pdb = assert_active_route_path(
        _resolve_mixed_path(row["context_pdb"]), f"Stage 34 {context_id} context PDB"
    )
    mapping_csv = assert_active_route_path(
        _resolve_mixed_path(row["context_mapping_csv"]), f"Stage 34 {context_id} mapping CSV"
    )
    if sha256_file(context_pdb) != str(row["context_pdb_sha256"]):
        raise RuntimeError(f"{context_id} context PDB SHA-256 mismatch")
    if sha256_file(mapping_csv) != str(row["context_mapping_csv_sha256"]):
        raise RuntimeError(f"{context_id} mapping CSV SHA-256 mismatch")
    for field in (
        "template_sequence_masked",
        "template_sidechains_masked",
        "template_interchain_features_masked",
    ):
        if _truthy(row.get(field)):
            raise RuntimeError(f"{context_id} is not a full target template: {field}=true")
    if str(row.get("template_mode")) != "target_context_full":
        raise RuntimeError(f"{context_id} template_mode must be target_context_full")

    mapping = read_csv(mapping_csv)
    chain_ids = _split_csv(row["chain_ids"])
    chain_lengths = [int(value) for value in _split_csv(row["chain_lengths"])]
    if len(mapping) != sum(chain_lengths) or len(chain_ids) != len(chain_lengths):
        raise RuntimeError(f"{context_id} chain lengths do not match its mapping")
    sequences: dict[str, list[str]] = {chain: [] for chain in chain_ids}
    numbers: dict[str, list[int]] = {chain: [] for chain in chain_ids}
    crop_indices: list[int] = []
    site_indices: list[int] = []
    hotspot_indices: list[int] = []
    for mapping_row in mapping:
        chain = str(mapping_row["context_chain"])
        if chain not in sequences:
            raise RuntimeError(f"{context_id} mapping contains unexpected chain {chain}")
        number = int(mapping_row["context_residue_number"])
        sequences[chain].append(str(mapping_row["context_residue_aa"]).upper())
        numbers[chain].append(number)
        if (
            chain == "A"
            and str(mapping_row["original_chain_id"]) == "G"
            and 115 <= int(mapping_row["original_pdb_residue_number"]) <= 200
        ):
            crop_indices.append(number)
        if chain == "A" and _truthy(mapping_row.get("is_target_site_residue")):
            site_indices.append(number)
        if chain == "A" and _truthy(mapping_row.get("is_selected_hotspot")):
            hotspot_indices.append(number)
    chain_sequences = []
    for chain, length in zip(chain_ids, chain_lengths):
        if numbers[chain] != list(range(1, length + 1)):
            raise RuntimeError(f"{context_id} chain {chain} mapping numbering is not contiguous")
        chain_sequences.append(_validate_sequence("".join(sequences[chain]), f"{context_id} chain {chain}", length))
    pdb_chains = parse_residues(context_pdb)
    if list(pdb_chains) != chain_ids:
        raise RuntimeError(f"{context_id} PDB chain order {list(pdb_chains)} != mapping {chain_ids}")
    if [residue_sequence(pdb_chains[chain]) for chain in chain_ids] != chain_sequences:
        raise RuntimeError(f"{context_id} PDB sequences do not match mapping sequences")
    if len(crop_indices) != 86 or len(hotspot_indices) != 4 or not site_indices:
        raise RuntimeError(f"{context_id} lost the 86-aa crop, Site_2, or four hotspots")
    if context_id == DEFAULT_CONTEXT_IDS[0] and (chain_ids != ["A"] or chain_lengths != [86]):
        raise RuntimeError("C1 must be the single-chain 86-aa crop")
    if context_id == DEFAULT_CONTEXT_IDS[1] and (chain_ids != ["A", "B", "C"] or chain_lengths != [174, 67, 60]):
        raise RuntimeError("C3 must contain the 174/67/60-aa native G/H/I context")
    return {
        "context_id": context_id,
        "context_description": row["context_description"],
        "context_pdb": context_pdb,
        "context_mapping_csv": mapping_csv,
        "context_chain_ids": chain_ids,
        "context_chain_lengths": chain_lengths,
        "context_chain_sequences": chain_sequences,
        "stage0_crop_context_indices_1based": crop_indices,
        "site2_fga_indices_1based": site_indices,
        "hotspot_fga_indices_1based": hotspot_indices,
    }


def _job_script(
    *,
    project_root: Path,
    runner: Path,
    job_spec: Path,
    python_bin: str,
    source_dir: str,
    overlay_dir: str,
    af_params: str,
) -> str:
    return f'''#!/bin/bash
set -euo pipefail
if [[ "${{RUN_STAGE5B_V2_PREDICTIONS:-NO}}" != "YES" ]]; then
  echo "Stage 5B-v2 predictions are review-gated. Set RUN_STAGE5B_V2_PREDICTIONS=YES after inspection." >&2
  exit 3
fi
PYTHON_BIN="${{AFCYCDESIGN_PYTHON:-{python_bin}}}"
SOURCE_DIR="${{COLABDESIGN_GAMMA_SOURCE:-{source_dir}}}"
OVERLAY_DIR="${{AFCYCDESIGN_PYTHON_OVERLAY:-{overlay_dir}}}"
AF_PARAMS_DIR="${{AF_PARAMS:-{af_params}}}"
export COLABDESIGN_GAMMA_SOURCE="$SOURCE_DIR"
export PYTHONPATH="$SOURCE_DIR:$OVERLAY_DIR${{PYTHONPATH:+:$PYTHONPATH}}"
cd {shlex.quote(_to_wsl_path(project_root))}
"$PYTHON_BIN" {shlex.quote(_to_wsl_path(runner))} \
  --job-spec {shlex.quote(_to_wsl_path(job_spec))} \
  --af-params "$AF_PARAMS_DIR"
'''


def _preflight_script(
    *,
    project_root: Path,
    runner: Path,
    specs: Iterable[Path],
    python_bin: str,
    source_dir: str,
    overlay_dir: str,
    af_params: str,
) -> str:
    commands = "\n".join(
        f'"$PYTHON_BIN" {shlex.quote(_to_wsl_path(runner))} --job-spec '
        f'{shlex.quote(_to_wsl_path(spec))} --af-params "$AF_PARAMS_DIR" --preflight-only'
        for spec in specs
    )
    return f'''#!/bin/bash
set -euo pipefail
PYTHON_BIN="${{AFCYCDESIGN_PYTHON:-{python_bin}}}"
SOURCE_DIR="${{COLABDESIGN_GAMMA_SOURCE:-{source_dir}}}"
OVERLAY_DIR="${{AFCYCDESIGN_PYTHON_OVERLAY:-{overlay_dir}}}"
AF_PARAMS_DIR="${{AF_PARAMS:-{af_params}}}"
export COLABDESIGN_GAMMA_SOURCE="$SOURCE_DIR"
export PYTHONPATH="$SOURCE_DIR:$OVERLAY_DIR${{PYTHONPATH:+:$PYTHONPATH}}"
cd {shlex.quote(_to_wsl_path(project_root))}
{commands}
echo "PASS: all Stage 5B-v2 C1/C3 preflights completed."
'''


def _master_script(preflight: Path, scripts: Iterable[Path]) -> str:
    commands = "\n".join(f"bash {shlex.quote(_to_wsl_path(path))}" for path in scripts)
    return f'''#!/bin/bash
set -euo pipefail
bash {shlex.quote(_to_wsl_path(preflight))}
{commands}
'''


def _plan_markdown(
    rows: list[Mapping[str, Any]],
    output_dir: Path,
    selection_mode: str,
    seed_jobs: int,
    model_predictions: int,
    job_shards: int,
    representative_preflights: int,
) -> str:
    candidates = sorted({str(row["stage5B_v2_candidate_id"]) for row in rows})
    contexts = sorted({str(row["context_id"]) for row in rows})
    return f"""# Stage 5B-v2 C1/C3 Recovery Plan

This campaign compares corrected-route Stage 4 candidates under two full
target-template contexts: C1 (the 86-aa crop) and C3 (the native G/H/I 301-aa
context).

Protocol invariants:

- target template sequence, backbone, sidechains, and target-target interchain
  features are retained;
- peptide sequence is supplied, but peptide template coverage is exactly zero;
- no peptide design coordinates and no initial guess are loaded by prediction;
- the cyclic positional offset is applied only to the peptide chain;
- target and peptide MSAs are one-row single-sequence inputs;
- dropout and MLM are disabled for this controlled comparison;
- all five AlphaFold multimer-v3 parameter sets and six requested recycles are used.

Planned matrix:

```text
candidates: {len(candidates)}
selection_mode: {selection_mode}
contexts: {', '.join(contexts)}
candidate_context_pairs: {len(rows)}
seed_jobs: {seed_jobs}
model_predictions: {model_predictions}
job_shards: {job_shards}
representative_preflights: {representative_preflights}
output_directory: {output_dir}
```

C1 is primarily a protocol sanity check. C3 is the more biologically relevant
native-context site-recovery test. Neither result alone makes a peptide a final
candidate. Prediction scripts remain gated by `RUN_STAGE5B_V2_PREDICTIONS=YES`.
"""


def main() -> int:
    parser = argparse.ArgumentParser(description="Prepare Stage 5B-v2 peptide recovery in C1 and C3 contexts.")
    parser.add_argument("--source-run-root", required=True)
    parser.add_argument("--stage4-scores-csv", required=True)
    parser.add_argument("--stage4-top-candidates-csv", required=True)
    parser.add_argument("--context-control-root", required=True)
    parser.add_argument("--output-root", required=True)
    parser.add_argument("--project-config", required=True)
    parser.add_argument("--context-ids", default=",".join(DEFAULT_CONTEXT_IDS))
    parser.add_argument(
        "--selection-mode",
        choices=["top_validation", "all_stage4_pass"],
        default="top_validation",
        help="Use the exact Stage 4 top table or every row that passes the Stage 4 hard gates.",
    )
    parser.add_argument("--candidate-count", type=int, default=5)
    parser.add_argument("--seeds-per-candidate-context", type=int, default=1)
    parser.add_argument("--models-per-seed", type=int, default=5)
    parser.add_argument("--recycles", type=int, default=6)
    parser.add_argument("--job-shards", type=int, default=1)
    parser.add_argument(
        "--allow-large-campaign",
        action="store_true",
        help="Required to write an all_stage4_pass campaign. Validation-only does not require it.",
    )
    parser.add_argument("--validate-inputs-only", action="store_true")
    parser.add_argument(
        "--afcycdesign-python",
        default="/home/luomi/fga_model_envs/colabdesign-py310/.pixi/envs/default/bin/python",
    )
    parser.add_argument(
        "--colabdesign-source",
        default="/home/luomi/fga_model_envs/sources/ColabDesign-gamma-stage5",
    )
    parser.add_argument(
        "--python-overlay",
        default="/home/luomi/fga_model_envs/stage5_afcycdesign_python_overlay",
    )
    parser.add_argument("--af-params", default="/home/luomi/fga_model_envs/af_params")
    args = parser.parse_args()

    logger = setup_logger("34_prepare_stage5b_v2_context_jobs")
    append_run_header(logger, "34_prepare_stage5b_v2_context_jobs.py")
    selected_contexts = _split_csv(args.context_ids)
    if selected_contexts != list(DEFAULT_CONTEXT_IDS):
        raise RuntimeError(f"This protocol is fixed to ordered contexts {DEFAULT_CONTEXT_IDS}")
    if args.candidate_count != 5:
        raise RuntimeError("--candidate-count is the fixed five-row Stage 4 audit top table size")
    if args.seeds_per_candidate_context < 1:
        raise RuntimeError("--seeds-per-candidate-context must be positive")
    if args.models_per_seed != 5 or args.recycles != 6:
        raise RuntimeError("Stage 5B-v2 requires five model sets and six requested recycles")
    if args.job_shards < 1:
        raise RuntimeError("--job-shards must be >= 1")
    if args.selection_mode == "all_stage4_pass" and not args.validate_inputs_only and not args.allow_large_campaign:
        raise RuntimeError(
            "Writing an all_stage4_pass C1/C3 campaign requires --allow-large-campaign after reviewing its size."
        )

    stage4 = load_stage4_validation_contract(
        source_run_root=args.source_run_root,
        stage4_scores_csv=args.stage4_scores_csv,
        stage4_top_candidates_csv=args.stage4_top_candidates_csv,
        project_config=args.project_config,
        candidate_count=5,
        selection_mode=args.selection_mode,
    )
    source_manifest = stage4["source_route_manifest"]
    source_manifest_path = stage4["source_route_manifest_path"]
    source_manifest_sha256 = stage4["source_route_manifest_sha256"]
    source_rows = list(stage4["selected_rows"])
    if args.job_shards > len(source_rows) * len(selected_contexts) * args.seeds_per_candidate_context:
        raise RuntimeError("--job-shards cannot exceed the planned seed-job count")

    context_root = assert_active_route_path(
        _resolve_mixed_path(args.context_control_root), "Stage 34 context-control root"
    )
    context_manifest_path, context_route, context_route_sha256 = load_route_manifest(context_root.parent)
    validate_route_project_config(args.project_config, context_route)
    for field in (
        "route_protocol_version",
        "hotspot_mapping_version",
        "cyclization",
        "project_config_sha256",
        "effective_project_config_sha256",
    ):
        if str(context_route[field]) != str(source_manifest[field]):
            raise RuntimeError(f"Stage 4 and target-context routes disagree on {field}")
    if context_route["stage0_sites"] != source_manifest["stage0_sites"]:
        raise RuntimeError("Stage 4 and target-context routes do not share identical Stage 0 provenance")
    context_provenance = route_provenance_fields(context_manifest_path, context_route, context_route_sha256)
    context_rows = read_csv(context_root / "FGA_rfpeptides_stage5_target_context_control_manifest.csv")
    context_lookup: dict[str, dict[str, Any]] = {}
    for row in context_rows:
        validate_row_route_provenance(row, context_provenance, f"Stage 34 context row {row.get('context_id')}")
        if row.get("context_id") in selected_contexts:
            context_lookup[str(row["context_id"])] = _context_contract(row)
    if list(context_lookup) != selected_contexts:
        raise RuntimeError(f"Missing or reordered C1/C3 context rows: {list(context_lookup)}")

    stage0_site = source_manifest["stage0_sites"][0]
    stage0_target = assert_active_route_path(stage0_site["target_pdb"], "Stage 34 Stage 0 target PDB")
    stage0_target_sequence = residue_sequence(parse_residues(stage0_target)["A"])
    if len(stage0_target_sequence) != 86:
        raise RuntimeError("Stage 0 target sequence must contain 86 residues")
    for context in context_lookup.values():
        crop_sequence = "".join(
            context["context_chain_sequences"][0][index - 1]
            for index in context["stage0_crop_context_indices_1based"]
        )
        if crop_sequence != stage0_target_sequence:
            raise RuntimeError(f"{context['context_id']} crop sequence differs from Stage 0")

    validated_sources: list[dict[str, Any]] = []
    for order, source in enumerate(source_rows, start=1):
        peptide = _validate_sequence(source["peptide_sequence"], f"Stage 4 candidate {order}")
        if int(source["peptide_length"]) != len(peptide):
            raise RuntimeError(f"Stage 4 candidate {order} peptide length mismatch")
        reference = assert_active_route_path(
            _resolve_mixed_path(source["scored_pdb"]), f"Stage 34 Stage 4 reference {order}"
        )
        if sha256_file(reference) != str(source["scored_pdb_sha256"]):
            raise RuntimeError(f"Stage 4 candidate {order} reference SHA-256 mismatch")
        target_chain = str(source.get("target_chain", "")).strip()
        peptide_chain = str(source.get("peptide_chain", "")).strip()
        if not target_chain or not peptide_chain:
            raise RuntimeError(f"Stage 4 candidate {order} lacks target_chain/peptide_chain provenance")
        _reference_sequences(
            reference,
            target_chain,
            peptide_chain,
            stage0_target_sequence,
            peptide,
        )
        validated_sources.append(
            {
                **source,
                "_peptide": peptide,
                "_reference": reference,
                "_reference_target_chain": target_chain,
                "_reference_peptide_chain": peptide_chain,
            }
        )

    planned_pairs = len(validated_sources) * len(selected_contexts)
    planned_jobs = planned_pairs * args.seeds_per_candidate_context
    planned_models = planned_jobs * args.models_per_seed
    if args.validate_inputs_only:
        logger.info("Stage 4 candidates validated: %s", len(validated_sources))
        logger.info("Stage 5B-v2 selection mode: %s", args.selection_mode)
        logger.info("Full target contexts validated: %s", len(selected_contexts))
        logger.info("Candidate-context pairs planned: %s", planned_pairs)
        logger.info("Seed jobs planned: %s", planned_jobs)
        logger.info("Model predictions planned: %s", planned_models)
        logger.info("No Stage 5B-v2 files or predictions were written.")
        return 0

    output_root = assert_active_route_path(
        _resolve_mixed_path(args.output_root), "Stage 34 output root", must_exist=False
    )
    output_dir = output_root / "07_structure_validation_target_context_conditioned"
    context_dir = output_dir / "inputs" / "target_contexts"
    mapping_dir = output_dir / "inputs" / "mappings"
    reference_dir = output_dir / "inputs" / "reference_design_poses"
    spec_dir = output_dir / "inputs" / "job_specs"
    jobs_dir = output_dir / "jobs"
    predictions_dir = output_dir / "predictions"
    project_root = resolve_path(".")
    runner = project_root / "scripts" / "external" / "run_afcycdesign_stage5b_v2_context_recovery.py"
    if not runner.is_file():
        raise RuntimeError(f"Missing Stage 5B-v2 runner: {runner}")

    staged_contexts: dict[str, dict[str, Any]] = {}
    for context_id, context in context_lookup.items():
        staged_pdb = context_dir / f"{context_id}.pdb"
        staged_mapping = mapping_dir / f"{context_id}_mapping.csv"
        staged_pdb.parent.mkdir(parents=True, exist_ok=True)
        staged_mapping.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(context["context_pdb"], staged_pdb)
        shutil.copy2(context["context_mapping_csv"], staged_mapping)
        staged_contexts[context_id] = {
            **context,
            "context_pdb": staged_pdb,
            "context_mapping_csv": staged_mapping,
        }

    protocol_version = (
        TOP5_PROTOCOL_VERSION if args.selection_mode == "top_validation" else ALL_PASS_PROTOCOL_VERSION
    )
    protocol_payload = {
        "protocol_version": protocol_version,
        "validation_test_type": VALIDATION_TEST_TYPE,
        "contexts": selected_contexts,
        "context_pdb_sha256": {
            key: sha256_file(value["context_pdb"]) for key, value in staged_contexts.items()
        },
        "template_mode": "target_context_full",
        "template_sequence_masked": False,
        "template_sidechains_masked": False,
        "template_interchain_features_masked": False,
        "peptide_template_expected_coverage": 0,
        "use_initial_guess": False,
        "reference_design_loaded_by_prediction_runner": False,
        "target_msa_mode": "single_sequence",
        "peptide_msa_mode": "single_sequence",
        "use_mlm": False,
        "use_dropout": False,
        "cyclic_topology_encoding": "peptide_chain_relative_position_cyclic_offset",
        "models_per_seed": args.models_per_seed,
        "requested_recycles": args.recycles,
        "colabdesign_commit": COLABDESIGN_GAMMA_COMMIT,
        "stage4_run_id": stage4["stage4_run_id"],
        "stage4_protocol_identity_sha256": stage4["stage4_protocol_identity_sha256"],
        "stage4_top_candidates_csv_sha256": stage4["stage4_top_candidates_csv_sha256"],
        "context_control_route_manifest_sha256": context_route_sha256,
        "candidate_selection_mode": args.selection_mode,
        "candidate_count": len(validated_sources),
    }
    protocol_hash = canonical_json_sha256(protocol_payload)[:12]
    campaign_label = "top5" if args.selection_mode == "top_validation" else f"all{len(validated_sources)}"
    campaign_id = f"stage5B_v2_{campaign_label}_C1_C3_{stage4['stage4_run_id']}_{protocol_hash}"
    route_manifest_path, route_manifest, route_manifest_sha256 = write_route_manifest(
        output_root,
        {
            "batch_id": campaign_id,
            "site_labels": list(source_manifest["site_labels"]),
            "protocol_peptide_length_min": int(source_manifest["protocol_peptide_length_min"]),
            "protocol_peptide_length_max": int(source_manifest["protocol_peptide_length_max"]),
            "run_peptide_length_min": int(source_manifest["run_peptide_length_min"]),
            "run_peptide_length_max": int(source_manifest["run_peptide_length_max"]),
            "num_designs_requested": int(source_manifest["num_designs_requested"]),
            "project_config": source_manifest["project_config"],
            "project_config_sha256": source_manifest["project_config_sha256"],
            "effective_project_config_sha256": source_manifest["effective_project_config_sha256"],
            "stage0_sites": list(source_manifest["stage0_sites"]),
            "source_route_manifests": [
                {
                    "run_id": source_manifest["run_id"],
                    "batch_id": source_manifest["batch_id"],
                    "manifest_path": str(source_manifest_path),
                    "manifest_sha256": source_manifest_sha256,
                },
                {
                    "run_id": context_route["run_id"],
                    "batch_id": context_route["batch_id"],
                    "manifest_path": str(context_manifest_path),
                    "manifest_sha256": context_route_sha256,
                },
            ],
            "stage5B_v2_protocol": protocol_payload,
            "stage5_selection_mode": args.selection_mode,
            "stage5_campaign_id": campaign_id,
            "stage5_candidate_count": len(validated_sources),
            "stage5_context_count": len(selected_contexts),
            "stage5_seed_job_count": planned_jobs,
            "stage5_model_prediction_count": planned_models,
            "stage5_job_shards": args.job_shards,
            "stage4_scores_csv": str(stage4["stage4_scores_csv"]),
            "stage4_scores_csv_sha256": stage4["stage4_scores_csv_sha256"],
            "stage4_top_candidates_csv": str(stage4["stage4_top_candidates_csv"]),
            "stage4_top_candidates_csv_sha256": stage4["stage4_top_candidates_csv_sha256"],
        },
    )
    provenance = route_provenance_fields(route_manifest_path, route_manifest, route_manifest_sha256)

    manifest_rows: list[dict[str, Any]] = []
    job_rows: list[dict[str, Any]] = []
    representative_specs: list[Path] = []
    representative_spec_keys: set[tuple[str, int]] = set()
    scripts_by_context: dict[str, list[Path]] = {context_id: [] for context_id in selected_contexts}
    scripts_by_shard: list[list[Path]] = [[] for _ in range(args.job_shards)]
    context_shard_counters: dict[str, int] = {context_id: 0 for context_id in selected_contexts}
    for order, source in enumerate(validated_sources, start=1):
        peptide = source["_peptide"]
        peptide_hash = _sha1_text(peptide)[:8]
        backbone_hash = _sha1_text(str(source["global_backbone_id"]))[:10]
        candidate_prefix = "S5B2CTX" if args.selection_mode == "top_validation" else "S5B2CTXALL"
        order_width = 2 if args.selection_mode == "top_validation" else 4
        candidate_id = (
            f"{candidate_prefix}_{order:0{order_width}d}_{source['batch']}_"
            f"gb{backbone_hash}_seq{peptide_hash}"
        )
        staged_reference = reference_dir / f"{candidate_id}_stage4_reference.pdb"
        staged_reference.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source["_reference"], staged_reference)
        reference_sha256 = sha256_file(staged_reference)
        for context_id in selected_contexts:
            context = staged_contexts[context_id]
            candidate_context_id = f"{candidate_id}_{_safe_token(context_id)}"
            peptide_chain_id = chr(ord("A") + len(context["context_chain_ids"]))
            cyclic_chain_index = len(context["context_chain_ids"])
            row = {
                "stage5B_v2_candidate_context_id": candidate_context_id,
                "stage5B_v2_candidate_id": candidate_id,
                "context_id": context_id,
                "context_description": context["context_description"],
                "context_pdb": context["context_pdb"],
                "context_pdb_sha256": sha256_file(context["context_pdb"]),
                "context_mapping_csv": context["context_mapping_csv"],
                "context_mapping_csv_sha256": sha256_file(context["context_mapping_csv"]),
                "context_chain_ids": ",".join(context["context_chain_ids"]),
                "context_chain_lengths": ",".join(str(value) for value in context["context_chain_lengths"]),
                "context_total_length": sum(context["context_chain_lengths"]),
                "stage0_crop_context_indices_1based": ",".join(str(value) for value in context["stage0_crop_context_indices_1based"]),
                "site2_fga_indices_1based": ",".join(str(value) for value in context["site2_fga_indices_1based"]),
                "hotspot_fga_indices_1based": ",".join(str(value) for value in context["hotspot_fga_indices_1based"]),
                "peptide_chain_id": peptide_chain_id,
                "cyclic_chain_index": cyclic_chain_index,
                "peptide_sequence": peptide,
                "peptide_sequence_hash": peptide_hash,
                "peptide_length": len(peptide),
                "reference_target_chain": source["_reference_target_chain"],
                "reference_peptide_chain": source["_reference_peptide_chain"],
                "staged_reference_design_pdb": staged_reference,
                "staged_reference_design_pdb_sha256": reference_sha256,
                "selection_order": order,
                "batch": source["batch"],
                "backbone_id": source["backbone_id"],
                **stage4_identity_values(source),
                "validation_test_type": VALIDATION_TEST_TYPE,
                "template_mode": "target_context_full",
                "template_sequence_masked": "false",
                "template_sidechains_masked": "false",
                "template_interchain_features_masked": "false",
                "peptide_template_expected_coverage": 0,
                "use_initial_guess": "false",
                "target_msa_mode": "single_sequence",
                "peptide_msa_mode": "single_sequence",
                "use_mlm": "false",
                "use_dropout": "false",
                "cyclic_topology_encoding": "peptide_chain_relative_position_cyclic_offset",
                "requested_recycles": args.recycles,
                "forward_passes": args.recycles + 1,
                "seeds_per_candidate_context": args.seeds_per_candidate_context,
                "models_per_seed": args.models_per_seed,
                "protocol_hash": protocol_hash,
                "stage5_campaign_id": campaign_id,
                "stage5_selection_mode": args.selection_mode,
                "status": "prepared_not_run",
                "notes": "C1=protocol sanity; C3=native-context recovery; not a final peptide decision",
                **{field: source.get(field, "") for field in SOURCE_ROUTE_PROVENANCE_FIELDS},
                **provenance,
            }
            manifest_rows.append(row)
            for seed in range(args.seeds_per_candidate_context):
                job_id = f"{candidate_context_id}_prot{protocol_hash}_seed{seed:02d}"
                output = predictions_dir / context_id / candidate_id / f"protocol_{protocol_hash}" / f"seed_{seed:02d}"
                spec_path = spec_dir / f"{job_id}.json"
                spec = {
                    "stage5B_v2_job_id": job_id,
                    "stage5B_v2_candidate_context_id": candidate_context_id,
                    "stage5B_v2_candidate_id": candidate_id,
                    "protocol_hash": protocol_hash,
                    "protocol_version": protocol_version,
                    "stage5_campaign_id": campaign_id,
                    "stage5_selection_mode": args.selection_mode,
                    "validation_test_type": VALIDATION_TEST_TYPE,
                    "context_id": context_id,
                    "context_description": context["context_description"],
                    "context_pdb": _to_wsl_path(context["context_pdb"]),
                    "context_pdb_sha256": sha256_file(context["context_pdb"]),
                    "context_mapping_csv": _to_wsl_path(context["context_mapping_csv"]),
                    "context_mapping_csv_sha256": sha256_file(context["context_mapping_csv"]),
                    "context_chain_ids": context["context_chain_ids"],
                    "context_chain_lengths": context["context_chain_lengths"],
                    "context_chain_sequences": context["context_chain_sequences"],
                    "stage0_crop_context_indices_1based": context["stage0_crop_context_indices_1based"],
                    "site2_fga_indices_1based": context["site2_fga_indices_1based"],
                    "hotspot_fga_indices_1based": context["hotspot_fga_indices_1based"],
                    "fga_chain_index": 0,
                    "peptide_chain_id": peptide_chain_id,
                    "peptide_sequence": peptide,
                    "peptide_sequence_length": len(peptide),
                    "peptide_sequence_hash": peptide_hash,
                    "cyclic_chain_index": cyclic_chain_index,
                    "cyclic_topology_encoding": "peptide_chain_relative_position_cyclic_offset",
                    "template_mode": "target_context_full",
                    "use_templates": True,
                    "use_batch_as_template": False,
                    "template_sequence_masked": False,
                    "template_sidechains_masked": False,
                    "template_interchain_features_masked": False,
                    "template_padding_mode": "full_target_context_then_zero_unknown_peptide_positions",
                    "target_template_expected_coverage": sum(context["context_chain_lengths"]),
                    "peptide_template_expected_coverage": 0,
                    "use_initial_guess": False,
                    "reference_design_pdb_for_posthoc": _to_wsl_path(staged_reference),
                    "reference_design_pdb_sha256": reference_sha256,
                    "reference_target_chain": source["_reference_target_chain"],
                    "reference_peptide_chain": source["_reference_peptide_chain"],
                    "reference_design_loaded_by_prediction_runner": False,
                    "batch": source["batch"],
                    "backbone_id": source["backbone_id"],
                    **stage4_identity_values(source),
                    "target_msa_mode": "single_sequence",
                    "peptide_msa_mode": "single_sequence",
                    "msa_rows_expected": 1,
                    "use_mlm": False,
                    "use_dropout": False,
                    "seed": seed,
                    "requested_recycles": args.recycles,
                    "forward_passes": args.recycles + 1,
                    "models_per_seed": args.models_per_seed,
                    "model_names_expected": MODEL_NAMES,
                    "colabdesign_commit": COLABDESIGN_GAMMA_COMMIT,
                    "prediction_output_dir": _to_wsl_path(output),
                    **{field: source.get(field, "") for field in SOURCE_ROUTE_PROVENANCE_FIELDS},
                    **provenance,
                }
                _write_text(spec_path, json.dumps(spec, indent=2, sort_keys=True))
                run_script = jobs_dir / f"run_{job_id}.sh"
                _write_text(
                    run_script,
                    _job_script(
                        project_root=project_root,
                        runner=runner,
                        job_spec=spec_path,
                        python_bin=args.afcycdesign_python,
                        source_dir=args.colabdesign_source,
                        overlay_dir=args.python_overlay,
                        af_params=args.af_params,
                    ),
                )
                if args.selection_mode == "top_validation":
                    representative_specs.append(spec_path)
                elif seed == 0:
                    representative_key = (context_id, len(peptide))
                    if representative_key not in representative_spec_keys:
                        representative_specs.append(spec_path)
                        representative_spec_keys.add(representative_key)
                scripts_by_context[context_id].append(run_script)
                shard_index = context_shard_counters[context_id] % args.job_shards
                context_shard_counters[context_id] += 1
                scripts_by_shard[shard_index].append(run_script)
                job_rows.append(
                    {
                        "stage5B_v2_job_id": job_id,
                        "stage5B_v2_candidate_context_id": candidate_context_id,
                        "stage5B_v2_candidate_id": candidate_id,
                        "context_id": context_id,
                        "peptide_sequence_hash": peptide_hash,
                        "protocol_hash": protocol_hash,
                        "stage5_campaign_id": campaign_id,
                        "stage5_selection_mode": args.selection_mode,
                        "job_shard": shard_index + 1,
                        "seed": seed,
                        "requested_recycles": args.recycles,
                        "forward_passes": args.recycles + 1,
                        "models_per_seed": args.models_per_seed,
                        "job_spec_json": spec_path,
                        "job_spec_sha256": sha256_file(spec_path),
                        "prediction_output_dir": output,
                        "run_script": run_script,
                        "template_mode": "target_context_full",
                        "peptide_template_expected_coverage": 0,
                        "use_initial_guess": "false",
                        "cyclic_chain_index": cyclic_chain_index,
                        "status": "prepared_not_run",
                        **{field: source.get(field, "") for field in SOURCE_ROUTE_PROVENANCE_FIELDS},
                        **provenance,
                    }
                )

    preflight = jobs_dir / "check_stage5B_v2_C1_C3_protocol.sh"
    _write_text(
        preflight,
        _preflight_script(
            project_root=project_root,
            runner=runner,
            specs=representative_specs,
            python_bin=args.afcycdesign_python,
            source_dir=args.colabdesign_source,
            overlay_dir=args.python_overlay,
            af_params=args.af_params,
        ),
    )
    all_scripts = [path for context_id in selected_contexts for path in scripts_by_context[context_id]]
    _write_text(jobs_dir / "run_stage5B_v2_C1_C3_all.sh", _master_script(preflight, all_scripts))
    for context_id in selected_contexts:
        _write_text(
            jobs_dir / f"run_stage5B_v2_{_safe_token(context_id)}_all.sh",
            _master_script(preflight, scripts_by_context[context_id]),
        )
    for shard_index, shard_scripts in enumerate(scripts_by_shard, start=1):
        _write_text(
            jobs_dir / f"run_stage5B_v2_C1_C3_shard_{shard_index:02d}_of_{args.job_shards:02d}.sh",
            _master_script(preflight, shard_scripts),
        )

    write_csv(output_dir / "FGA_rfpeptides_stage5B_v2_candidate_context_manifest.csv", manifest_rows, MANIFEST_FIELDS)
    write_csv(output_dir / "FGA_rfpeptides_stage5B_v2_prediction_jobs.csv", job_rows, JOB_FIELDS)
    write_markdown(
        output_dir / "FGA_rfpeptides_stage5B_v2_plan.md",
        _plan_markdown(
            manifest_rows,
            output_dir,
            args.selection_mode,
            len(job_rows),
            len(job_rows) * args.models_per_seed,
            args.job_shards,
            len(representative_specs),
        ),
    )
    logger.info("Stage 5B-v2 candidates prepared: %s", len(validated_sources))
    logger.info("Stage 5B-v2 contexts prepared: %s", len(selected_contexts))
    logger.info("Candidate-context pairs prepared: %s", len(manifest_rows))
    logger.info("Seed jobs prepared: %s", len(job_rows))
    logger.info("Planned model predictions: %s", len(job_rows) * args.models_per_seed)
    logger.info("Stage 5B-v2 selection mode: %s", args.selection_mode)
    logger.info("Stage 5B-v2 job shards: %s", args.job_shards)
    logger.info("Representative preflight specs: %s", len(representative_specs))
    logger.info("Output directory: %s", output_dir)
    logger.info("No Stage 5B-v2 prediction was run.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
