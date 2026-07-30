from __future__ import annotations

import argparse
import hashlib
import importlib.util
import math
import os
from collections import Counter
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

from common import (
    ROUTE_PROVENANCE_FIELDS,
    SOURCE_ROUTE_PROVENANCE_FIELDS,
    assert_active_route_path,
    append_run_header,
    load_route_manifest,
    read_csv,
    resolve_path,
    route_provenance_fields,
    rows_to_markdown,
    setup_logger,
    validate_route_project_config,
    write_csv,
    write_markdown,
    write_route_manifest,
)
from pdb_utils import parse_residues, residue_sequence


SCRIPT_DIR = Path(__file__).resolve().parent
STAGE3C_PATH = SCRIPT_DIR / "23_collect_proteinmpnn_sequences.py"
BACKBONE_ATOMS = ("N", "CA", "C", "O")


STAGE3D1_FIELDS = [
    "stage3d1_design_id",
    "source_sequence_design_id",
    "stage3_job_id",
    "run_group_id",
    "protocol_identity_sha256",
    "stage3_mode",
    "global_backbone_id",
    "source_local_design_id",
    "source_batch_label",
    "backbone_family_id",
    "backbone_id",
    "site_label",
    "site_id",
    "stage3_jobs_csv",
    "stage3_jobs_csv_sha256",
    "stage2_selection_csv",
    "stage2_selection_csv_sha256",
    "stage3c_qc_csv",
    "stage3c_qc_csv_sha256",
    "stage3_input_pdb",
    "stage3_input_pdb_sha256",
    "source_backbone_pdb",
    "source_backbone_pdb_sha256",
    "input_pdb",
    "input_pdb_sha256",
    "repacked_pdb",
    "repacked_pdb_sha256",
    "target_chain",
    "peptide_chain",
    "peptide_sequence",
    "target_sequence",
    "sequence_duplicate_status",
    "repack_status",
    "repack_mode",
    "input_pose_total_score",
    "repacked_pose_total_score",
    "repack_score_delta",
    "repack_residue_count",
    "peptide_repack_residue_count",
    "target_repack_residue_count",
    "target_repack_residue_labels",
    "peptide_backbone_atom_count",
    "peptide_backbone_rmsd_A",
    "peptide_backbone_max_displacement_A",
    "target_backbone_atom_count",
    "target_backbone_rmsd_A",
    "target_backbone_max_displacement_A",
    "backbone_max_displacement_tolerance_A",
    "backbone_geometry_status",
    "repacked_peptide_sequence",
    "repacked_target_sequence",
    "pre_target_site_recovery_status",
    "pre_hotspot_recovery_status",
    "pre_macrocycle_geometry_status",
    "pre_clash_status",
    "pre_num_target_site_contacts",
    "pre_num_hotspot_contacts",
    "pre_peptide_site_min_distance",
    "pre_peptide_hotspot_min_distance",
    "pre_macrocycle_terminal_cn_distance",
    "post_target_site_recovery_status",
    "post_hotspot_recovery_status",
    "post_macrocycle_geometry_status",
    "post_clash_status",
    "post_num_target_site_contacts",
    "post_num_hotspot_contacts",
    "post_peptide_site_min_distance",
    "post_peptide_hotspot_min_distance",
    "post_macrocycle_terminal_cn_distance",
    "post_rosetta_pose_total",
    "pass_stage3d1_qc",
    "stage3d1_failure_reasons",
    "notes",
] + ROUTE_PROVENANCE_FIELDS + SOURCE_ROUTE_PROVENANCE_FIELDS


def _load_stage3c_module() -> Any:
    spec = importlib.util.spec_from_file_location("stage3c_qc", STAGE3C_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not import Stage 3C QC helper: {STAGE3C_PATH}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _split_csv(value: str) -> list[str]:
    return [item.strip() for item in str(value or "").split(",") if item.strip()]


def _safe_token(value: str) -> str:
    keep = []
    for ch in str(value):
        if ch.isalnum() or ch in {"_", "-", "."}:
            keep.append(ch)
        else:
            keep.append("_")
    return "".join(keep).strip("_") or "item"


def _resolve_mixed_path(value: str | Path) -> Path:
    text = str(value).strip().replace("\\", "/")
    if not text:
        return resolve_path(text)
    if os.name == "nt" and text.startswith("/mnt/") and len(text) > 7 and text[6] == "/":
        return Path(f"{text[5].upper()}:/{text[7:]}")
    if os.name != "nt" and len(text) >= 3 and text[1] == ":" and text[2] == "/":
        return Path(f"/mnt/{text[0].lower()}{text[2:]}")
    path = Path(text)
    if path.is_absolute():
        return path
    return resolve_path(path)


def _read_required_csv(path: Path) -> list[dict[str, str]]:
    rows = read_csv(path)
    if not rows:
        raise RuntimeError(f"Missing or empty CSV: {path}")
    return rows


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _parse_int(value: Any, default: int = 0) -> int:
    try:
        return int(str(value).strip())
    except ValueError:
        return default


def _same_path(left: str | Path, right: str | Path) -> bool:
    return _resolve_mixed_path(left).resolve() == _resolve_mixed_path(right).resolve()


def _require_equal(observed: Any, expected: Any, label: str) -> None:
    if str(observed).strip() != str(expected).strip():
        raise RuntimeError(
            f"{label} mismatch: observed={str(observed).strip()!r}, expected={str(expected).strip()!r}"
        )


def _require_recorded_file(
    *,
    row: Mapping[str, Any],
    path_field: str,
    sha256_field: str,
    expected_path: Path,
    expected_sha256: str,
    label: str,
) -> None:
    recorded_text = str(row.get(path_field, "")).strip()
    if not recorded_text:
        raise RuntimeError(f"{label} is missing {path_field}")
    recorded_path = assert_active_route_path(
        _resolve_mixed_path(recorded_text),
        f"{label} {path_field}",
    )
    if recorded_path.resolve() != expected_path.resolve():
        raise RuntimeError(
            f"{label} {path_field} mismatch: observed={recorded_path}, expected={expected_path}"
        )
    recorded_sha256 = str(row.get(sha256_field, "")).strip().lower()
    if recorded_sha256 != expected_sha256:
        raise RuntimeError(
            f"{label} {sha256_field} mismatch: observed={recorded_sha256}, expected={expected_sha256}"
        )


def _validate_stage3c_table(
    *,
    qchelper: Any,
    stage3c_rows: Sequence[Mapping[str, str]],
    stage2_selection_lookup: Mapping[str, Mapping[str, str]],
    stage3_job_lookup: Mapping[str, Mapping[str, str]],
    expected_by_backbone: Mapping[str, Mapping[Path, tuple[Path, str]]],
    stage3_route_provenance: Mapping[str, str],
    aggregate_sources: Mapping[Path, Mapping[str, str]],
    stage2_selection_csv: Path,
    stage2_selection_csv_sha256: str,
    stage3_jobs_csv: Path,
    stage3_jobs_csv_sha256: str,
    run_group_id: str,
    protocol_identity_sha256: str,
    stage3_mode: str,
) -> list[dict[str, str]]:
    expected_paths: dict[Path, tuple[str, Path]] = {}
    for backbone_id, outputs in expected_by_backbone.items():
        for output_pdb, (input_pdb, _) in outputs.items():
            resolved = output_pdb.resolve()
            if resolved in expected_paths:
                raise RuntimeError(f"Two Stage 3 jobs expect the same output PDB: {resolved}")
            expected_paths[resolved] = (backbone_id, input_pdb.resolve())

    if len(stage3c_rows) != len(expected_paths):
        raise RuntimeError(
            "Stage 3C table is not the complete authoritative Stage 3B output set: "
            f"rows={len(stage3c_rows)}, expected={len(expected_paths)}"
        )

    seen_sequence_ids: set[str] = set()
    seen_output_paths: set[Path] = set()
    validated: list[dict[str, str]] = []
    for raw_row in stage3c_rows:
        row = dict(raw_row)
        sequence_design_id = str(row.get("sequence_design_id", "")).strip()
        if not sequence_design_id:
            raise RuntimeError("Stage 3C row is missing sequence_design_id")
        if sequence_design_id in seen_sequence_ids:
            raise RuntimeError(f"Duplicate Stage 3C sequence_design_id: {sequence_design_id}")
        seen_sequence_ids.add(sequence_design_id)

        global_backbone_id = str(row.get("global_backbone_id", "")).strip()
        if not global_backbone_id or "," in global_backbone_id:
            raise RuntimeError(
                f"Stage 3C row has invalid global_backbone_id: {global_backbone_id!r}"
            )
        backbone_row = stage2_selection_lookup.get(global_backbone_id)
        job_row = stage3_job_lookup.get(global_backbone_id)
        if backbone_row is None or job_row is None:
            raise RuntimeError(
                f"Stage 3C row is not linked to the authoritative Stage 2.5/jobs tables: "
                f"{global_backbone_id}"
            )

        output_pdb = assert_active_route_path(
            _resolve_mixed_path(str(row.get("relaxed_pdb", ""))),
            f"Stage 24 Stage 3C PDB {sequence_design_id}",
        )
        output_resolved = output_pdb.resolve()
        if output_resolved in seen_output_paths:
            raise RuntimeError(f"Duplicate Stage 3C relaxed_pdb: {output_pdb}")
        seen_output_paths.add(output_resolved)
        expected_record = expected_paths.get(output_resolved)
        if expected_record is None:
            raise RuntimeError(
                f"Stage 3C row references an output not authorized by the jobs table: {output_pdb}"
            )
        expected_backbone_id, expected_input_pdb = expected_record
        if expected_backbone_id != global_backbone_id:
            raise RuntimeError(
                f"Stage 3C output belongs to {expected_backbone_id}, not {global_backbone_id}: "
                f"{output_pdb}"
            )
        if output_pdb.stem != sequence_design_id:
            raise RuntimeError(
                f"Stage 3C sequence_design_id/PDB stem mismatch: "
                f"{sequence_design_id} != {output_pdb.stem}"
            )
        output_sha256 = _sha256_file(output_pdb)
        _require_equal(
            row.get("relaxed_pdb_sha256", ""),
            output_sha256,
            f"Stage 3C relaxed PDB SHA-256 for {sequence_design_id}",
        )

        _require_equal(
            row.get("backbone_id", ""),
            global_backbone_id,
            f"Stage 3C backbone_id for {sequence_design_id}",
        )
        _require_equal(
            row.get("stage3_job_id", ""),
            job_row.get("stage3_job_id", ""),
            f"Stage 3C stage3_job_id for {sequence_design_id}",
        )
        _require_equal(
            row.get("run_group_id", ""),
            run_group_id,
            f"Stage 3C run_group_id for {sequence_design_id}",
        )
        _require_equal(
            row.get("protocol_identity_sha256", ""),
            protocol_identity_sha256,
            f"Stage 3C protocol identity for {sequence_design_id}",
        )
        _require_equal(
            row.get("stage3_mode", ""),
            stage3_mode,
            f"Stage 3C mode for {sequence_design_id}",
        )
        for field in [
            "source_local_design_id",
            "source_batch_label",
            "backbone_family_id",
            "site_label",
            "site_id",
        ]:
            _require_equal(
                row.get(field, ""),
                backbone_row.get(field, ""),
                f"Stage 3C {field} for {sequence_design_id}",
            )
        for field in ["target_chain", "peptide_chain"]:
            _require_equal(
                row.get(field, ""),
                backbone_row.get(field, ""),
                f"Stage 3C {field} for {sequence_design_id}",
            )

        _require_recorded_file(
            row=row,
            path_field="stage3_jobs_csv",
            sha256_field="stage3_jobs_csv_sha256",
            expected_path=stage3_jobs_csv,
            expected_sha256=stage3_jobs_csv_sha256,
            label=f"Stage 3C row {sequence_design_id}",
        )
        _require_recorded_file(
            row=row,
            path_field="stage2_selection_csv",
            sha256_field="stage2_selection_csv_sha256",
            expected_path=stage2_selection_csv,
            expected_sha256=stage2_selection_csv_sha256,
            label=f"Stage 3C row {sequence_design_id}",
        )

        input_pdb = assert_active_route_path(
            _resolve_mixed_path(str(row.get("stage3_input_pdb", ""))),
            f"Stage 24 Stage 3 input PDB for {sequence_design_id}",
        )
        if input_pdb.resolve() != expected_input_pdb:
            raise RuntimeError(
                f"Stage 3C input PDB does not match the jobs table for {sequence_design_id}: "
                f"observed={input_pdb}, expected={expected_input_pdb}"
            )
        _require_equal(
            row.get("stage3_input_pdb_sha256", ""),
            _sha256_file(input_pdb),
            f"Stage 3C input PDB SHA-256 for {sequence_design_id}",
        )

        source_pdb = assert_active_route_path(
            _resolve_mixed_path(str(row.get("source_backbone_pdb", ""))),
            f"Stage 24 source backbone PDB for {sequence_design_id}",
        )
        expected_source_pdb = assert_active_route_path(
            _resolve_mixed_path(str(job_row.get("source_backbone_pdb", ""))),
            f"Stage 24 jobs-table source backbone for {global_backbone_id}",
        )
        if source_pdb.resolve() != expected_source_pdb.resolve():
            raise RuntimeError(
                f"Stage 3C source backbone path mismatch for {sequence_design_id}"
            )
        source_sha256 = _sha256_file(source_pdb)
        _require_equal(
            row.get("source_backbone_pdb_sha256", ""),
            source_sha256,
            f"Stage 3C source backbone SHA-256 for {sequence_design_id}",
        )
        _require_equal(
            job_row.get("source_backbone_pdb_sha256", ""),
            source_sha256,
            f"Stage 3 job source backbone SHA-256 for {global_backbone_id}",
        )

        expected_count = len(expected_by_backbone[global_backbone_id])
        _require_equal(
            row.get("expected_stage3b_outputs_for_backbone", ""),
            expected_count,
            f"Stage 3C expected output count for {sequence_design_id}",
        )
        _require_equal(
            row.get("observed_stage3b_outputs_for_backbone", ""),
            expected_count,
            f"Stage 3C observed output count for {sequence_design_id}",
        )

        qchelper._validate_cached_route_provenance(
            row,
            stage3_route_provenance,
            f"Stage 24 Stage 3C row {sequence_design_id}",
        )
        qchelper._validate_aggregate_source_membership(
            row,
            aggregate_sources,
            f"Stage 24 Stage 3C row {sequence_design_id}",
        )
        for field in SOURCE_ROUTE_PROVENANCE_FIELDS:
            _require_equal(
                row.get(field, ""),
                job_row.get(field, ""),
                f"Stage 3C {field} for {sequence_design_id}",
            )
        validated.append(row)

    missing = sorted(set(expected_paths) - seen_output_paths, key=str)
    if missing:
        raise RuntimeError(
            f"Stage 3C table is missing {len(missing)} jobs-table outputs; first={missing[:5]}"
        )
    return validated


def _validate_authoritative_stage3_inputs(
    *,
    qchelper: Any,
    stage3_root: Path,
    stage2_selection_csv: Path,
    stage3_jobs_csv: Path,
    stage3c_qc_csv: Path,
    stage3_mode: str,
    stage3_route_provenance: Mapping[str, str],
    aggregate_sources: Mapping[Path, Mapping[str, str]],
    manifest_site_label: str,
) -> dict[str, Any]:
    stage2_selection_rows = _read_required_csv(stage2_selection_csv)
    stage2_selection_lookup = qchelper._strict_lookup_rows(
        stage2_selection_rows,
        "global_backbone_id",
        "Stage 2.5 selection",
    )
    job_rows = _read_required_csv(stage3_jobs_csv)
    stage3_job_lookup = qchelper._strict_lookup_rows(
        job_rows,
        "global_backbone_id",
        "Stage 3 job",
    )
    qchelper._strict_lookup_rows(job_rows, "stage3_job_id", "Stage 3 job")
    if set(stage2_selection_lookup) != set(stage3_job_lookup):
        missing_jobs = sorted(set(stage2_selection_lookup) - set(stage3_job_lookup))
        extra_jobs = sorted(set(stage3_job_lookup) - set(stage2_selection_lookup))
        raise RuntimeError(
            "Stage 2.5 selection and Stage 3 jobs must contain the same global backbone IDs: "
            f"missing_jobs={missing_jobs[:5]}, extra_jobs={extra_jobs[:5]}"
        )

    run_group_id = qchelper._single_value(job_rows, "run_group_id", "Stage 3 jobs")
    protocol_identity_sha256 = qchelper._single_value(
        job_rows,
        "protocol_identity_sha256",
        "Stage 3 jobs",
    )
    if len(protocol_identity_sha256) != 64:
        raise RuntimeError("Stage 3 jobs protocol_identity_sha256 is not a full SHA-256 digest")
    expected_run_group_id = f"stage3_{len(job_rows)}bp_{protocol_identity_sha256[:12]}"
    if run_group_id != expected_run_group_id:
        raise RuntimeError(
            f"Stage 3 run_group_id formula mismatch: observed={run_group_id}, "
            f"expected={expected_run_group_id}"
        )
    locked_mode = qchelper._single_value(job_rows, "stage3_mode", "Stage 3 jobs")
    if locked_mode != stage3_mode:
        raise RuntimeError(
            f"--stage3-mode={stage3_mode} does not match jobs table mode={locked_mode}"
        )
    if locked_mode != "proteinmpnn_only":
        raise RuntimeError(
            "Stage 3D-1 accepts only ProteinMPNN-only outputs; FastRelax outputs are not valid inputs"
        )

    stage2_selection_csv_sha256 = _sha256_file(stage2_selection_csv)
    stage3_jobs_csv_sha256 = _sha256_file(stage3_jobs_csv)
    expected_by_backbone: dict[str, dict[Path, tuple[Path, str]]] = {}
    family_ids: set[str] = set()
    for backbone_id, job_row in stage3_job_lookup.items():
        backbone_row = stage2_selection_lookup[backbone_id]
        if str(backbone_row.get("pass_backbone_qc", "")).strip().lower() != "true":
            raise RuntimeError(f"Stage 2.5 row is not pass_backbone_qc=true: {backbone_id}")
        if str(backbone_row.get("stage2_5_selected", "")).strip().lower() != "true":
            raise RuntimeError(f"Stage 2.5 row is not stage2_5_selected=true: {backbone_id}")
        if str(backbone_row.get("site_label", "")).strip() != manifest_site_label:
            raise RuntimeError(
                f"Stage 2.5 site label does not match the manifest-locked site: {backbone_id}"
            )
        family_id = str(backbone_row.get("backbone_family_id", "")).strip()
        if not family_id:
            raise RuntimeError(f"Stage 2.5 row has no backbone_family_id: {backbone_id}")
        if family_id in family_ids:
            raise RuntimeError(
                f"Stage 3 jobs contain more than one representative from family {family_id}"
            )
        family_ids.add(family_id)
        if str(backbone_row.get("family_representative_global_backbone_id", "")).strip() != backbone_id:
            raise RuntimeError(
                f"Stage 2.5 selected row is not its family representative: {backbone_id}"
            )
        qchelper._validate_cached_route_provenance(
            backbone_row,
            stage3_route_provenance,
            f"Stage 24 Stage 2.5 row {backbone_id}",
        )
        qchelper._validate_cached_route_provenance(
            job_row,
            stage3_route_provenance,
            f"Stage 24 Stage 3 job {backbone_id}",
        )
        expected_by_backbone[backbone_id] = qchelper._validate_stage3_job_row(
            backbone_id=backbone_id,
            backbone_row=backbone_row,
            job_row=job_row,
            expected_stage3_mode=stage3_mode,
            stage3_root=stage3_root,
            stage2_selection_csv=stage2_selection_csv,
            stage2_selection_csv_sha256=stage2_selection_csv_sha256,
            aggregate_sources=aggregate_sources,
        )

    qchelper._validate_runlist_contract(job_rows=job_rows, stage3_root=stage3_root)
    expected_output_count = sum(len(outputs) for outputs in expected_by_backbone.values())
    observed_output_count = qchelper._validate_complete_stage3_outputs(expected_by_backbone)
    if observed_output_count != expected_output_count:
        raise RuntimeError(
            "Stage 3B output count changed after completeness validation: "
            f"expected={expected_output_count}, observed={observed_output_count}"
        )

    expected_stage3c_dir = (
        stage3_root
        / "05_proteinmpnn_sequences"
        / "stage3c"
        / _safe_token(run_group_id)
    )
    expected_stage3c_csv = (
        expected_stage3c_dir
        / f"FGA_rfpeptides_{_safe_token(run_group_id)}_stage3C_sequences_qc.csv"
    )
    if stage3c_qc_csv.resolve() != expected_stage3c_csv.resolve():
        raise RuntimeError(
            "Stage 3D-1 requires the run-group-scoped authoritative Stage 3C QC table: "
            f"observed={stage3c_qc_csv}, expected={expected_stage3c_csv}"
        )
    stage3c_rows = _validate_stage3c_table(
        qchelper=qchelper,
        stage3c_rows=_read_required_csv(stage3c_qc_csv),
        stage2_selection_lookup=stage2_selection_lookup,
        stage3_job_lookup=stage3_job_lookup,
        expected_by_backbone=expected_by_backbone,
        stage3_route_provenance=stage3_route_provenance,
        aggregate_sources=aggregate_sources,
        stage2_selection_csv=stage2_selection_csv,
        stage2_selection_csv_sha256=stage2_selection_csv_sha256,
        stage3_jobs_csv=stage3_jobs_csv,
        stage3_jobs_csv_sha256=stage3_jobs_csv_sha256,
        run_group_id=run_group_id,
        protocol_identity_sha256=protocol_identity_sha256,
        stage3_mode=stage3_mode,
    )
    return {
        "stage2_selection_lookup": stage2_selection_lookup,
        "stage3_job_lookup": stage3_job_lookup,
        "stage3c_rows": stage3c_rows,
        "run_group_id": run_group_id,
        "protocol_identity_sha256": protocol_identity_sha256,
        "stage2_selection_csv_sha256": stage2_selection_csv_sha256,
        "stage3_jobs_csv_sha256": stage3_jobs_csv_sha256,
        "stage3c_qc_csv_sha256": _sha256_file(stage3c_qc_csv),
        "expected_output_count": expected_output_count,
    }


def _selected_stage3c_rows(
    rows: Sequence[Mapping[str, str]],
    selected_global_backbones: set[str],
    include_duplicate_sequences: bool,
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    selected: list[dict[str, str]] = []
    skipped: list[dict[str, str]] = []
    seen_sequences: set[tuple[str, str]] = set()
    required = {
        "file_status": "pass",
        "parse_status": "pass",
        "sequence_status": "pass_sequence",
        "target_site_recovery_status": "site_contact_pass",
        "hotspot_recovery_status": "hotspot_contact_pass",
        "macrocycle_geometry_status": "pass_head_to_tail_macrocycle",
    }
    for raw_row in rows:
        row = dict(raw_row)
        backbone_id = str(row.get("global_backbone_id", "")).strip()
        if selected_global_backbones and backbone_id not in selected_global_backbones:
            continue
        sequence = str(row.get("peptide_sequence", "")).strip()
        if not sequence:
            skipped.append({**row, "skip_reason": "missing_peptide_sequence"})
            continue
        failed = [field for field, expected in required.items() if row.get(field, "") != expected]
        if failed:
            skipped.append(
                {
                    **row,
                    "skip_reason": "failed_pre_repack_filter:" + ",".join(failed),
                }
            )
            continue
        duplicate_key = (backbone_id, sequence)
        is_duplicate = duplicate_key in seen_sequences
        if is_duplicate and not include_duplicate_sequences:
            skipped.append(
                {
                    **row,
                    "skip_reason": "duplicate_sequence_within_global_backbone",
                }
            )
            continue
        row["_sequence_duplicate_status"] = (
            "included_duplicate_within_global_backbone"
            if is_duplicate
            else "unique_within_global_backbone"
        )
        selected.append(row)
        seen_sequences.add(duplicate_key)
    return selected, skipped


def _validate_selected_input_structures(
    *,
    selected_rows: Sequence[Mapping[str, str]],
    stage2_selection_lookup: Mapping[str, Mapping[str, str]],
    stage0_target_sequence: str,
) -> dict[str, str]:
    target_sequences: dict[str, str] = {}
    for row in selected_rows:
        sequence_id = str(row.get("sequence_design_id", "")).strip()
        backbone_id = str(row.get("global_backbone_id", "")).strip()
        backbone_row = stage2_selection_lookup[backbone_id]
        input_pdb = assert_active_route_path(
            _resolve_mixed_path(str(row.get("relaxed_pdb", ""))),
            f"Stage 24 selected Stage 3C PDB {sequence_id}",
        )
        _require_equal(
            row.get("relaxed_pdb_sha256", ""),
            _sha256_file(input_pdb),
            f"Stage 24 selected PDB SHA-256 for {sequence_id}",
        )
        chains = parse_residues(input_pdb)
        target_chain = str(row.get("target_chain", "")).strip()
        peptide_chain = str(row.get("peptide_chain", "")).strip()
        target_residues = list(chains.get(target_chain, []))
        peptide_residues = list(chains.get(peptide_chain, []))
        if not target_residues or not peptide_residues:
            raise RuntimeError(
                f"Stage 24 selected input is missing target/peptide chain for {sequence_id}"
            )
        peptide_sequence = residue_sequence(peptide_residues)
        target_sequence = residue_sequence(target_residues)
        _require_equal(
            peptide_sequence,
            row.get("peptide_sequence", ""),
            f"Stage 24 peptide sequence for {sequence_id}",
        )
        expected_length = _parse_int(backbone_row.get("peptide_length", "0"))
        if len(peptide_residues) != expected_length:
            raise RuntimeError(
                f"Stage 24 peptide length mismatch for {sequence_id}: "
                f"observed={len(peptide_residues)}, expected={expected_length}"
            )
        if target_sequence != stage0_target_sequence:
            raise RuntimeError(
                f"Stage 24 target sequence differs from manifest-locked Stage 0 target for {sequence_id}"
            )
        target_sequences[sequence_id] = target_sequence
    return target_sequences


def _load_pyrosetta() -> Any:
    try:
        import pyrosetta  # type: ignore
    except ModuleNotFoundError as exc:
        raise RuntimeError(
            "PyRosetta is required for Stage 3D-1 side-chain repack. "
            "Run in the proteinmpnn_binder_design environment."
        ) from exc
    return pyrosetta


def _init_pyrosetta(pyrosetta: Any) -> None:
    pyrosetta.init(
        "-beta_nov16 -mute all -use_terminal_residues true "
        "-ex1 -ex2aro -packing:use_input_sc"
    )


def _xyz_tuple(xyz: Any) -> tuple[float, float, float]:
    try:
        return float(xyz.x), float(xyz.y), float(xyz.z)
    except TypeError:
        return float(xyz[0]), float(xyz[1]), float(xyz[2])


def _residue_atom_coords(
    residue: Any,
    heavy_only: bool = True,
) -> list[tuple[float, float, float]]:
    coords: list[tuple[float, float, float]] = []
    for atom_idx in range(1, residue.natoms() + 1):
        if heavy_only and residue.atom_is_hydrogen(atom_idx):
            continue
        coords.append(_xyz_tuple(residue.xyz(atom_idx)))
    return coords


def _sq_distance(
    left: tuple[float, float, float],
    right: tuple[float, float, float],
) -> float:
    return sum((left[idx] - right[idx]) ** 2 for idx in range(3))


def _pose_chain(pose: Any, pose_idx: int) -> str:
    pdb_info = pose.pdb_info()
    chain = pdb_info.chain(pose_idx) if pdb_info is not None else ""
    return str(chain).strip() or "_"


def _pose_residue_label(pose: Any, pose_idx: int) -> str:
    pdb_info = pose.pdb_info()
    chain = _pose_chain(pose, pose_idx)
    number = pdb_info.number(pose_idx) if pdb_info is not None else pose_idx
    icode = pdb_info.icode(pose_idx).strip() if pdb_info is not None else ""
    return f"{chain}{number}{icode}"


def _pose_number(pose: Any, pose_idx: int) -> str:
    pdb_info = pose.pdb_info()
    if pdb_info is None:
        return str(pose_idx)
    icode = pdb_info.icode(pose_idx).strip()
    return f"{pdb_info.number(pose_idx)}{icode}"


def _mapped_site_hotspots(
    *,
    qchelper: Any,
    backbone_row: Mapping[str, str],
    input_pdb: Path,
    target_chain: str,
    site_numbers: set[str],
    hotspot_numbers: set[str],
) -> tuple[set[str], set[str]]:
    input_chains = parse_residues(input_pdb)
    input_target_residues = list(input_chains.get(target_chain, []))
    source_backbone_pdb = assert_active_route_path(
        _resolve_mixed_path(str(backbone_row.get("rf_pdb", ""))),
        "Stage 24 source backbone PDB",
    )
    target_number_map, _, _ = qchelper._target_number_mapping(
        source_backbone_pdb=source_backbone_pdb,
        target_chain=target_chain,
        relaxed_target_residues=input_target_residues,
    )
    mapped_site_numbers = qchelper._apply_number_mapping(site_numbers, target_number_map)
    mapped_hotspot_numbers = qchelper._apply_number_mapping(
        hotspot_numbers,
        target_number_map,
    )
    return mapped_site_numbers, mapped_hotspot_numbers


def _repack_positions(
    *,
    pose: Any,
    peptide_chain: str,
    target_chain: str,
    target_repack_radius: float,
    mapped_hotspot_numbers: set[str],
) -> tuple[set[int], set[int], set[int]]:
    peptide_positions = {
        idx
        for idx in range(1, pose.total_residue() + 1)
        if _pose_chain(pose, idx) == peptide_chain
    }
    target_positions = {
        idx
        for idx in range(1, pose.total_residue() + 1)
        if _pose_chain(pose, idx) == target_chain
    }
    peptide_coords: list[tuple[float, float, float]] = []
    for idx in sorted(peptide_positions):
        peptide_coords.extend(
            _residue_atom_coords(pose.residue(idx), heavy_only=True)
        )

    radius_sq = target_repack_radius * target_repack_radius
    interface_target_positions: set[int] = set()
    for idx in sorted(target_positions):
        target_coords = _residue_atom_coords(
            pose.residue(idx),
            heavy_only=True,
        )
        if any(
            _sq_distance(left, right) <= radius_sq
            for left in target_coords
            for right in peptide_coords
        ):
            interface_target_positions.add(idx)

    for idx in sorted(target_positions):
        if _pose_number(pose, idx) in mapped_hotspot_numbers:
            interface_target_positions.add(idx)

    return (
        peptide_positions | interface_target_positions,
        peptide_positions,
        interface_target_positions,
    )


def _chain_backbone_records(
    pdb_path: Path,
    chain_id: str,
) -> tuple[str, list[str], list[tuple[float, float, float]]]:
    residues = list(parse_residues(pdb_path).get(chain_id, []))
    if not residues:
        raise RuntimeError(f"Chain {chain_id} is missing from {pdb_path}")
    residue_numbers: list[str] = []
    coords: list[tuple[float, float, float]] = []
    for residue in residues:
        residue_number = str(residue.get("pdb_residue_number", ""))
        atoms = residue.get("atoms", {})
        missing_atoms = [atom for atom in BACKBONE_ATOMS if atom not in atoms]
        if missing_atoms:
            raise RuntimeError(
                f"Missing backbone atoms {missing_atoms} at {chain_id}{residue_number} "
                f"in {pdb_path}"
            )
        residue_numbers.append(residue_number)
        coords.extend(atoms[atom] for atom in BACKBONE_ATOMS)
    return residue_sequence(residues), residue_numbers, coords


def _backbone_displacement_metrics(
    *,
    input_pdb: Path,
    output_pdb: Path,
    chain_id: str,
    chain_label: str,
    max_displacement_tolerance: float,
) -> dict[str, Any]:
    input_sequence, input_numbers, input_coords = _chain_backbone_records(
        input_pdb,
        chain_id,
    )
    output_sequence, output_numbers, output_coords = _chain_backbone_records(
        output_pdb,
        chain_id,
    )
    if input_numbers != output_numbers:
        raise RuntimeError(
            f"{chain_label} residue numbering changed during Stage 3D-1 repack"
        )
    if input_sequence != output_sequence:
        raise RuntimeError(
            f"{chain_label} sequence changed during Stage 3D-1 repack"
        )
    if len(input_coords) != len(output_coords) or not input_coords:
        raise RuntimeError(
            f"{chain_label} backbone atom correspondence changed during Stage 3D-1 repack"
        )
    squared = [
        _sq_distance(left, right)
        for left, right in zip(input_coords, output_coords)
    ]
    rmsd = math.sqrt(sum(squared) / len(squared))
    max_displacement = math.sqrt(max(squared))
    return {
        f"{chain_label}_backbone_atom_count": len(input_coords),
        f"{chain_label}_backbone_rmsd_A": round(rmsd, 6),
        f"{chain_label}_backbone_max_displacement_A": round(
            max_displacement,
            6,
        ),
        f"repacked_{chain_label}_sequence": output_sequence,
        f"{chain_label}_backbone_fixed": max_displacement
        <= max_displacement_tolerance,
    }


def _run_repack(
    *,
    pyrosetta: Any,
    input_pdb: Path,
    output_pdb: Path,
    peptide_chain: str,
    target_chain: str,
    target_repack_radius: float,
    mapped_hotspot_numbers: set[str],
    backbone_max_displacement_tolerance: float,
) -> dict[str, Any]:
    pose = pyrosetta.pose_from_pdb(str(input_pdb))
    scorefxn = pyrosetta.get_fa_scorefxn()
    before_score = float(scorefxn(pose))
    repack_positions, peptide_positions, target_positions = _repack_positions(
        pose=pose,
        peptide_chain=peptide_chain,
        target_chain=target_chain,
        target_repack_radius=target_repack_radius,
        mapped_hotspot_numbers=mapped_hotspot_numbers,
    )
    if not peptide_positions:
        raise RuntimeError(f"Peptide chain {peptide_chain} was not found in {input_pdb}")
    if not target_positions:
        raise RuntimeError(f"No target interface residues selected for {input_pdb}")

    task = pyrosetta.standard_packer_task(pose)
    task.restrict_to_repacking()
    task.or_include_current(True)
    for idx in range(1, pose.total_residue() + 1):
        if idx not in repack_positions:
            task.nonconst_residue_task(idx).prevent_repacking()

    mover = pyrosetta.rosetta.protocols.minimization_packing.PackRotamersMover(
        scorefxn,
        task,
    )
    mover.apply(pose)
    after_score = float(scorefxn(pose))

    output_pdb.parent.mkdir(parents=True, exist_ok=True)
    temporary_pdb = output_pdb.with_name(f"{output_pdb.stem}.tmp.pdb")
    if temporary_pdb.exists():
        temporary_pdb.unlink()
    try:
        pose.dump_pdb(str(temporary_pdb))
        peptide_metrics = _backbone_displacement_metrics(
            input_pdb=input_pdb,
            output_pdb=temporary_pdb,
            chain_id=peptide_chain,
            chain_label="peptide",
            max_displacement_tolerance=backbone_max_displacement_tolerance,
        )
        target_metrics = _backbone_displacement_metrics(
            input_pdb=input_pdb,
            output_pdb=temporary_pdb,
            chain_id=target_chain,
            chain_label="target",
            max_displacement_tolerance=backbone_max_displacement_tolerance,
        )
        temporary_pdb.replace(output_pdb)
    except Exception:
        if temporary_pdb.exists():
            temporary_pdb.unlink()
        raise

    peptide_backbone_fixed = bool(
        peptide_metrics.pop("peptide_backbone_fixed")
    )
    target_backbone_fixed = bool(
        target_metrics.pop("target_backbone_fixed")
    )
    backbone_geometry_status = (
        "pass_fixed_backbone"
        if peptide_backbone_fixed and target_backbone_fixed
        else "fail_backbone_moved"
    )
    return {
        "input_pose_total_score": round(before_score, 3),
        "repacked_pose_total_score": round(after_score, 3),
        "repack_score_delta": round(after_score - before_score, 3),
        "repack_residue_count": len(repack_positions),
        "peptide_repack_residue_count": len(peptide_positions),
        "target_repack_residue_count": len(target_positions),
        "target_repack_residue_labels": ",".join(
            _pose_residue_label(pose, idx)
            for idx in sorted(target_positions)
        ),
        **peptide_metrics,
        **target_metrics,
        "backbone_max_displacement_tolerance_A": (
            backbone_max_displacement_tolerance
        ),
        "backbone_geometry_status": backbone_geometry_status,
    }


def _qc_repacked_pdb(
    *,
    qchelper: Any,
    repacked_pdb: Path,
    backbone_row: Mapping[str, str],
    site_numbers: set[str],
    hotspot_numbers: set[str],
    args: argparse.Namespace,
) -> dict[str, Any]:
    return qchelper._qc_row_for_relaxed_pdb(
        relaxed_pdb=repacked_pdb,
        backbone_row=backbone_row,
        site_numbers=site_numbers,
        hotspot_numbers=hotspot_numbers,
        contact_cutoff=args.contact_cutoff,
        site_near_distance=args.site_near_distance,
        hotspot_near_distance=args.hotspot_near_distance,
        severe_clash_distance=args.severe_clash_distance,
        min_target_contacts=args.min_target_contacts,
        min_site_contacts=args.min_site_contacts,
        min_hotspot_contacts=args.min_hotspot_contacts,
        macrocycle_pass_distance=args.macrocycle_pass_distance,
        macrocycle_warn_distance=args.macrocycle_warn_distance,
        forbidden_aas=args.forbidden_aas,
    )


def _stage3d1_failure_reasons(
    *,
    post_qc: Mapping[str, Any],
    input_peptide_sequence: str,
    input_target_sequence: str,
    repack_metrics: Mapping[str, Any],
) -> list[str]:
    reasons: list[str] = []
    if str(post_qc.get("peptide_sequence", "")) != input_peptide_sequence:
        reasons.append("peptide_sequence_changed")
    if str(repack_metrics.get("repacked_peptide_sequence", "")) != input_peptide_sequence:
        reasons.append("repacked_peptide_sequence_changed")
    if str(repack_metrics.get("repacked_target_sequence", "")) != input_target_sequence:
        reasons.append("target_sequence_changed")
    if repack_metrics.get("backbone_geometry_status") != "pass_fixed_backbone":
        reasons.append(
            str(
                repack_metrics.get(
                    "backbone_geometry_status",
                    "backbone_geometry_not_verified",
                )
            )
        )
    if post_qc.get("sequence_status") != "pass_sequence":
        reasons.append(str(post_qc.get("sequence_status", "sequence_not_pass")))
    if post_qc.get("target_site_recovery_status") != "site_contact_pass":
        reasons.append(
            str(post_qc.get("target_site_recovery_status", "site_not_pass"))
        )
    if post_qc.get("hotspot_recovery_status") != "hotspot_contact_pass":
        reasons.append(
            str(post_qc.get("hotspot_recovery_status", "hotspot_not_pass"))
        )
    if post_qc.get("macrocycle_geometry_status") != "pass_head_to_tail_macrocycle":
        reasons.append(
            str(post_qc.get("macrocycle_geometry_status", "macrocycle_not_pass"))
        )
    if post_qc.get("clash_status") != "pass_no_severe_clash":
        reasons.append(str(post_qc.get("clash_status", "clash_not_pass")))
    return list(dict.fromkeys(reasons))


def _stage3d1_provenance(
    *,
    input_row: Mapping[str, str],
    job_row: Mapping[str, str],
    stage2_selection_csv: Path,
    stage2_selection_csv_sha256: str,
    stage3_jobs_csv: Path,
    stage3_jobs_csv_sha256: str,
    stage3c_qc_csv: Path,
    stage3c_qc_csv_sha256: str,
    output_route_provenance: Mapping[str, str],
) -> dict[str, Any]:
    return {
        "stage3_job_id": str(job_row.get("stage3_job_id", "")),
        "run_group_id": str(job_row.get("run_group_id", "")),
        "protocol_identity_sha256": str(
            job_row.get("protocol_identity_sha256", "")
        ),
        "stage3_mode": str(job_row.get("stage3_mode", "")),
        "global_backbone_id": str(job_row.get("global_backbone_id", "")),
        "source_local_design_id": str(
            job_row.get("source_local_design_id", "")
        ),
        "source_batch_label": str(job_row.get("source_batch_label", "")),
        "backbone_family_id": str(job_row.get("backbone_family_id", "")),
        "stage3_jobs_csv": str(stage3_jobs_csv),
        "stage3_jobs_csv_sha256": stage3_jobs_csv_sha256,
        "stage2_selection_csv": str(stage2_selection_csv),
        "stage2_selection_csv_sha256": stage2_selection_csv_sha256,
        "stage3c_qc_csv": str(stage3c_qc_csv),
        "stage3c_qc_csv_sha256": stage3c_qc_csv_sha256,
        "stage3_input_pdb": str(input_row.get("stage3_input_pdb", "")),
        "stage3_input_pdb_sha256": str(
            input_row.get("stage3_input_pdb_sha256", "")
        ),
        "source_backbone_pdb": str(
            input_row.get("source_backbone_pdb", "")
        ),
        "source_backbone_pdb_sha256": str(
            input_row.get("source_backbone_pdb_sha256", "")
        ),
        **output_route_provenance,
        **{
            field: str(job_row.get(field, ""))
            for field in SOURCE_ROUTE_PROVENANCE_FIELDS
        },
    }


def _summary_markdown(
    *,
    rows: list[Mapping[str, Any]],
    skipped_rows: list[Mapping[str, Any]],
    output_dir: Path,
    args: argparse.Namespace,
    run_group_id: str,
    expected_stage3c_rows: int,
) -> str:
    pass_rows = [
        row for row in rows if row.get("pass_stage3d1_qc") == "true"
    ]
    columns = [
        "stage3d1_design_id",
        "global_backbone_id",
        "peptide_sequence",
        "post_target_site_recovery_status",
        "post_hotspot_recovery_status",
        "post_macrocycle_geometry_status",
        "post_clash_status",
        "peptide_backbone_max_displacement_A",
        "target_backbone_max_displacement_A",
        "backbone_geometry_status",
        "pass_stage3d1_qc",
        "stage3d1_failure_reasons",
    ]
    skip_counts = Counter(
        str(row.get("skip_reason", ""))
        for row in skipped_rows
    )
    skip_lines = (
        "\n".join(
            f"- {key or 'blank'}: {skip_counts[key]}"
            for key in sorted(skip_counts)
        )
        or "- none: 0"
    )
    return f"""# FGA RFpeptides Stage 3D-1 Side-Chain Repack QC

Status: ProteinMPNN-only Stage 3C structures were processed by side-chain
repack only. This stage does not redesign sequence, run FastRelax, or minimize
the backbone.

Authoritative input contract:

```text
run_group_id: {run_group_id}
stage3_mode: proteinmpnn_only
complete_stage3c_rows_validated: {expected_stage3c_rows}
stage2_selection_csv: {args.stage2_selection_csv}
stage3_jobs_csv: {args.stage3_jobs_csv}
stage3c_qc_csv: {args.stage3c_qc_csv}
```

Method:

```text
fixed sequence: true
backbone minimization: false
FastRelax: false
peptide chain repack: all side chains
target chain repack: side chains within {args.target_repack_radius:g} A of peptide plus mapped hotspots
backbone max-displacement tolerance: {args.backbone_max_displacement_tolerance:g} A
```

The fixed-backbone claim is checked from the written PDB coordinates for N, CA,
C, and O atoms on both peptide and target chains. Any sequence change, backbone
movement above tolerance, lost Site_2/hotspot contact, damaged head-to-tail
geometry, or remaining severe clash fails Stage 3D-1.

Output directory:

```text
{output_dir}
```

## Counts

```text
input_rows_repacked: {len(rows)}
pass_stage3d1_qc: {len(pass_rows)}
skipped_stage3c_rows: {len(skipped_rows)}
```

Skipped input rows:

{skip_lines}

## Repacked Rows

{rows_to_markdown(rows, columns, "No Stage 3D-1 rows were generated.")}

Passing Stage 3D-1 is an intermediate structural QC result, not a final peptide
candidate claim.
"""


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Stage 3D-1 fixed-backbone side-chain repack for authoritative "
            "ProteinMPNN-only Stage 3C outputs."
        )
    )
    parser.add_argument("--stage0-root", required=True)
    parser.add_argument("--stage3-root", required=True)
    parser.add_argument("--project-config", required=True)
    parser.add_argument("--stage2-selection-csv", required=True)
    parser.add_argument("--stage3-jobs-csv", required=True)
    parser.add_argument("--stage3c-qc-csv", required=True)
    parser.add_argument(
        "--stage3-mode",
        choices=["proteinmpnn_only"],
        required=True,
    )
    parser.add_argument(
        "--output-root",
        default="",
        help="Defaults to --stage3-root.",
    )
    parser.add_argument(
        "--selected-global-backbones",
        default="",
        help=(
            "Optional comma-separated global_backbone_id subset. The complete "
            "Stage 3C table is still validated before the subset is repacked."
        ),
    )
    parser.add_argument("--target-repack-radius", type=float, default=8.0)
    parser.add_argument(
        "--backbone-max-displacement-tolerance",
        type=float,
        default=0.002,
    )
    parser.add_argument("--include-duplicate-sequences", action="store_true")
    parser.add_argument("--validate-inputs-only", action="store_true")
    parser.add_argument("--contact-cutoff", type=float, default=5.0)
    parser.add_argument("--site-near-distance", type=float, default=6.0)
    parser.add_argument("--hotspot-near-distance", type=float, default=8.0)
    parser.add_argument("--severe-clash-distance", type=float, default=1.2)
    parser.add_argument("--min-target-contacts", type=int, default=1)
    parser.add_argument("--min-site-contacts", type=int, default=1)
    parser.add_argument("--min-hotspot-contacts", type=int, default=1)
    parser.add_argument("--macrocycle-pass-distance", type=float, default=2.0)
    parser.add_argument("--macrocycle-warn-distance", type=float, default=3.0)
    parser.add_argument("--forbidden-aas", default="CX")
    parser.add_argument(
        "--pymol-path-style",
        choices=["windows", "wsl", "native"],
        default="windows",
        help="Path style for generated PyMOL review scripts.",
    )
    args = parser.parse_args()

    logger = setup_logger("24_stage3d1_sidechain_repack")
    append_run_header(logger, "24_stage3d1_sidechain_repack.py")

    if args.target_repack_radius <= 0:
        raise RuntimeError("--target-repack-radius must be > 0")
    if args.backbone_max_displacement_tolerance <= 0:
        raise RuntimeError(
            "--backbone-max-displacement-tolerance must be > 0"
        )
    if args.contact_cutoff <= 0:
        raise RuntimeError("--contact-cutoff must be > 0")
    if args.site_near_distance < args.contact_cutoff:
        raise RuntimeError("--site-near-distance must be >= --contact-cutoff")
    if args.hotspot_near_distance < args.contact_cutoff:
        raise RuntimeError(
            "--hotspot-near-distance must be >= --contact-cutoff"
        )
    if args.macrocycle_warn_distance < args.macrocycle_pass_distance:
        raise RuntimeError(
            "--macrocycle-warn-distance must be >= "
            "--macrocycle-pass-distance"
        )

    qchelper = _load_stage3c_module()
    stage0_root = assert_active_route_path(
        _resolve_mixed_path(args.stage0_root),
        "Stage 24 Stage 0 root",
    )
    stage3_root = assert_active_route_path(
        _resolve_mixed_path(args.stage3_root),
        "Stage 24 Stage 3 root",
    )
    output_root = assert_active_route_path(
        _resolve_mixed_path(args.output_root)
        if args.output_root
        else stage3_root,
        "Stage 24 output root",
        must_exist=False,
    )
    stage2_selection_csv = assert_active_route_path(
        _resolve_mixed_path(args.stage2_selection_csv),
        "Stage 24 Stage 2.5 selection CSV",
    )
    stage3_jobs_csv = assert_active_route_path(
        _resolve_mixed_path(args.stage3_jobs_csv),
        "Stage 24 Stage 3 jobs CSV",
    )
    stage3c_qc_csv = assert_active_route_path(
        _resolve_mixed_path(args.stage3c_qc_csv),
        "Stage 24 Stage 3C QC CSV",
    )
    qchelper._require_within(
        stage2_selection_csv,
        stage3_root,
        "Stage 24 Stage 2.5 selection CSV",
    )
    qchelper._require_within(
        stage3_jobs_csv,
        stage3_root,
        "Stage 24 Stage 3 jobs CSV",
    )
    qchelper._require_within(
        stage3c_qc_csv,
        stage3_root,
        "Stage 24 Stage 3C QC CSV",
    )

    route_manifest_path, route_manifest, route_manifest_sha256 = (
        load_route_manifest(stage3_root)
    )
    validate_route_project_config(args.project_config, route_manifest)
    stage3_route_provenance = route_provenance_fields(
        route_manifest_path,
        route_manifest,
        route_manifest_sha256,
    )
    aggregate_sources = qchelper._aggregate_source_manifest_lookup(
        route_manifest
    )
    manifest_site_label, stage0_target_pdb, site_mapping_csv = (
        qchelper._stage0_site_inputs(
            route_manifest=route_manifest,
            stage0_root=stage0_root,
        )
    )
    site_numbers, hotspot_numbers = qchelper._load_site_mapping(
        site_mapping_csv
    )

    contract = _validate_authoritative_stage3_inputs(
        qchelper=qchelper,
        stage3_root=stage3_root,
        stage2_selection_csv=stage2_selection_csv,
        stage3_jobs_csv=stage3_jobs_csv,
        stage3c_qc_csv=stage3c_qc_csv,
        stage3_mode=args.stage3_mode,
        stage3_route_provenance=stage3_route_provenance,
        aggregate_sources=aggregate_sources,
        manifest_site_label=manifest_site_label,
    )
    stage2_selection_lookup = contract["stage2_selection_lookup"]
    stage3_job_lookup = contract["stage3_job_lookup"]
    run_group_id = str(contract["run_group_id"])

    requested = _split_csv(args.selected_global_backbones)
    if len(requested) != len(set(requested)):
        raise RuntimeError(
            "--selected-global-backbones contains duplicate global backbone IDs"
        )
    missing_requested = sorted(set(requested) - set(stage3_job_lookup))
    if missing_requested:
        raise RuntimeError(
            "Requested global backbone IDs are absent from the authoritative "
            f"Stage 3 jobs table: {missing_requested[:10]}"
        )
    selected_rows, skipped_rows = _selected_stage3c_rows(
        contract["stage3c_rows"],
        set(requested),
        include_duplicate_sequences=args.include_duplicate_sequences,
    )
    if not selected_rows:
        raise RuntimeError(
            "No contact-preserving Stage 3C rows were selected for Stage 3D-1"
        )

    target_chains = {
        str(row.get("target_chain", "")).strip()
        for row in selected_rows
    }
    if len(target_chains) != 1:
        raise RuntimeError(
            f"Stage 3D-1 selected rows do not share one target chain: "
            f"{sorted(target_chains)}"
        )
    stage0_target_residues = list(
        parse_residues(stage0_target_pdb).get(next(iter(target_chains)), [])
    )
    if not stage0_target_residues:
        raise RuntimeError(
            "Manifest-locked Stage 0 target PDB does not contain the selected "
            "target chain"
        )
    stage0_target_sequence = residue_sequence(stage0_target_residues)
    target_sequences = _validate_selected_input_structures(
        selected_rows=selected_rows,
        stage2_selection_lookup=stage2_selection_lookup,
        stage0_target_sequence=stage0_target_sequence,
    )

    logger.info("Authoritative Stage 3 jobs validated: %s", len(stage3_job_lookup))
    logger.info(
        "Complete Stage 3C rows validated: %s",
        contract["expected_output_count"],
    )
    logger.info(
        "Contact-preserving unique Stage 3C rows selected: %s",
        len(selected_rows),
    )
    if args.validate_inputs_only:
        logger.info(
            "Input validation completed; PyRosetta was not loaded and no "
            "Stage 3D-1 outputs were written."
        )
        return 0

    if output_root.resolve() != stage3_root.resolve():
        output_manifest_path, output_manifest, output_manifest_sha256 = (
            write_route_manifest(output_root, route_manifest)
        )
        output_route_provenance = route_provenance_fields(
            output_manifest_path,
            output_manifest,
            output_manifest_sha256,
        )
    else:
        output_route_provenance = stage3_route_provenance

    output_dir = (
        output_root
        / "05_proteinmpnn_sequences"
        / "stage3d1"
        / _safe_token(run_group_id)
    )
    output_pdb_dir = output_dir / "repacked_pdbs"
    output_prefix = (
        f"FGA_rfpeptides_{_safe_token(run_group_id)}_stage3D1"
    )

    pyrosetta = _load_pyrosetta()
    _init_pyrosetta(pyrosetta)
    output_pdb_dir.mkdir(parents=True, exist_ok=True)

    result_rows: list[dict[str, Any]] = []
    post_qc_rows: list[dict[str, Any]] = []
    target_chain_values: set[str] = set()
    peptide_chain_values: set[str] = set()
    for input_row in selected_rows:
        backbone_id = str(input_row.get("global_backbone_id", "")).strip()
        backbone_row = stage2_selection_lookup[backbone_id]
        job_row = stage3_job_lookup[backbone_id]
        sequence_id = str(input_row.get("sequence_design_id", "")).strip()
        peptide_chain = str(input_row.get("peptide_chain", "")).strip()
        target_chain = str(input_row.get("target_chain", "")).strip()
        target_chain_values.add(target_chain)
        peptide_chain_values.add(peptide_chain)
        input_pdb = assert_active_route_path(
            _resolve_mixed_path(str(input_row.get("relaxed_pdb", ""))),
            f"Stage 24 input PDB {sequence_id}",
        )
        input_sequence = str(input_row.get("peptide_sequence", "")).strip()
        input_target_sequence = target_sequences[sequence_id]
        stage3d1_id = f"{sequence_id}_stage3d1_repack"
        repacked_pdb = (
            output_pdb_dir / f"{_safe_token(stage3d1_id)}.pdb"
        )

        mapped_site_numbers, mapped_hotspot_numbers = (
            _mapped_site_hotspots(
                qchelper=qchelper,
                backbone_row=backbone_row,
                input_pdb=input_pdb,
                target_chain=target_chain,
                site_numbers=site_numbers,
                hotspot_numbers=hotspot_numbers,
            )
        )
        try:
            repack_metrics = _run_repack(
                pyrosetta=pyrosetta,
                input_pdb=input_pdb,
                output_pdb=repacked_pdb,
                peptide_chain=peptide_chain,
                target_chain=target_chain,
                target_repack_radius=args.target_repack_radius,
                mapped_hotspot_numbers=mapped_hotspot_numbers,
                backbone_max_displacement_tolerance=(
                    args.backbone_max_displacement_tolerance
                ),
            )
            repack_status = "success"
            notes = ""
        except Exception as exc:
            repack_metrics = {}
            repack_status = f"failed:{exc.__class__.__name__}"
            notes = str(exc)

        if repack_status == "success":
            post_qc = _qc_repacked_pdb(
                qchelper=qchelper,
                repacked_pdb=repacked_pdb,
                backbone_row=backbone_row,
                site_numbers=site_numbers,
                hotspot_numbers=hotspot_numbers,
                args=args,
            )
            post_qc_rows.append(post_qc)
            failure_reasons = _stage3d1_failure_reasons(
                post_qc=post_qc,
                input_peptide_sequence=input_sequence,
                input_target_sequence=input_target_sequence,
                repack_metrics=repack_metrics,
            )
            repacked_pdb_sha256 = _sha256_file(repacked_pdb)
        else:
            post_qc = {}
            failure_reasons = [repack_status]
            repacked_pdb_sha256 = ""

        provenance = _stage3d1_provenance(
            input_row=input_row,
            job_row=job_row,
            stage2_selection_csv=stage2_selection_csv,
            stage2_selection_csv_sha256=contract[
                "stage2_selection_csv_sha256"
            ],
            stage3_jobs_csv=stage3_jobs_csv,
            stage3_jobs_csv_sha256=contract["stage3_jobs_csv_sha256"],
            stage3c_qc_csv=stage3c_qc_csv,
            stage3c_qc_csv_sha256=contract["stage3c_qc_csv_sha256"],
            output_route_provenance=output_route_provenance,
        )
        result_rows.append(
            {
                "stage3d1_design_id": stage3d1_id,
                "source_sequence_design_id": sequence_id,
                **provenance,
                "backbone_id": backbone_id,
                "site_label": backbone_row.get("site_label", ""),
                "site_id": backbone_row.get("site_id", ""),
                "input_pdb": input_pdb,
                "input_pdb_sha256": input_row.get(
                    "relaxed_pdb_sha256",
                    "",
                ),
                "repacked_pdb": (
                    repacked_pdb if repack_status == "success" else ""
                ),
                "repacked_pdb_sha256": repacked_pdb_sha256,
                "target_chain": target_chain,
                "peptide_chain": peptide_chain,
                "peptide_sequence": input_sequence,
                "target_sequence": input_target_sequence,
                "sequence_duplicate_status": input_row.get(
                    "_sequence_duplicate_status",
                    "",
                ),
                "repack_status": repack_status,
                "repack_mode": (
                    "side_chain_repack_only_fixed_sequence_no_backbone_"
                    "minimization_no_fastrelax"
                ),
                **repack_metrics,
                "pre_target_site_recovery_status": input_row.get(
                    "target_site_recovery_status",
                    "",
                ),
                "pre_hotspot_recovery_status": input_row.get(
                    "hotspot_recovery_status",
                    "",
                ),
                "pre_macrocycle_geometry_status": input_row.get(
                    "macrocycle_geometry_status",
                    "",
                ),
                "pre_clash_status": input_row.get("clash_status", ""),
                "pre_num_target_site_contacts": input_row.get(
                    "num_target_site_contacts",
                    "",
                ),
                "pre_num_hotspot_contacts": input_row.get(
                    "num_hotspot_contacts",
                    "",
                ),
                "pre_peptide_site_min_distance": input_row.get(
                    "peptide_site_min_distance",
                    "",
                ),
                "pre_peptide_hotspot_min_distance": input_row.get(
                    "peptide_hotspot_min_distance",
                    "",
                ),
                "pre_macrocycle_terminal_cn_distance": input_row.get(
                    "macrocycle_terminal_cn_distance",
                    "",
                ),
                "post_target_site_recovery_status": post_qc.get(
                    "target_site_recovery_status",
                    "",
                ),
                "post_hotspot_recovery_status": post_qc.get(
                    "hotspot_recovery_status",
                    "",
                ),
                "post_macrocycle_geometry_status": post_qc.get(
                    "macrocycle_geometry_status",
                    "",
                ),
                "post_clash_status": post_qc.get("clash_status", ""),
                "post_num_target_site_contacts": post_qc.get(
                    "num_target_site_contacts",
                    "",
                ),
                "post_num_hotspot_contacts": post_qc.get(
                    "num_hotspot_contacts",
                    "",
                ),
                "post_peptide_site_min_distance": post_qc.get(
                    "peptide_site_min_distance",
                    "",
                ),
                "post_peptide_hotspot_min_distance": post_qc.get(
                    "peptide_hotspot_min_distance",
                    "",
                ),
                "post_macrocycle_terminal_cn_distance": post_qc.get(
                    "macrocycle_terminal_cn_distance",
                    "",
                ),
                "post_rosetta_pose_total": post_qc.get(
                    "rosetta_pose_total",
                    "",
                ),
                "pass_stage3d1_qc": (
                    "true" if not failure_reasons else "false"
                ),
                "stage3d1_failure_reasons": ";".join(failure_reasons),
                "notes": notes,
            }
        )

    if len(target_chain_values) != 1 or len(peptide_chain_values) != 1:
        raise RuntimeError(
            "Stage 3D-1 selected rows changed chain contract during processing"
        )
    pass_rows = [
        row for row in result_rows if row.get("pass_stage3d1_qc") == "true"
    ]
    write_csv(
        output_dir / f"{output_prefix}_sidechain_repack_qc.csv",
        result_rows,
        STAGE3D1_FIELDS,
    )
    write_csv(
        output_dir / f"{output_prefix}_sidechain_repack_qc_pass.csv",
        pass_rows,
        STAGE3D1_FIELDS,
    )
    skipped_fields = list(qchelper.STAGE3C_FIELDS) + ["skip_reason"]
    write_csv(
        output_dir / f"{output_prefix}_skipped_inputs.csv",
        skipped_rows,
        skipped_fields,
    )
    write_markdown(
        output_dir / f"{output_prefix}_sidechain_repack_qc.md",
        _summary_markdown(
            rows=result_rows,
            skipped_rows=skipped_rows,
            output_dir=output_dir,
            args=args,
            run_group_id=run_group_id,
            expected_stage3c_rows=contract["expected_output_count"],
        ),
    )
    if post_qc_rows:
        qchelper._write_pymol_review(
            rows=post_qc_rows,
            output_path=(
                output_dir
                / f"{_safe_token(manifest_site_label)}_"
                f"{_safe_token(run_group_id)}_stage3D1_review.pml"
            ),
            target_chain=next(iter(target_chain_values)),
            peptide_chain=next(iter(peptide_chain_values)),
            site_numbers=site_numbers,
            hotspot_numbers=hotspot_numbers,
            top_n=20,
            pymol_path_style=args.pymol_path_style,
        )

    logger.info("Stage 3D-1 input rows repacked: %s", len(result_rows))
    logger.info("Stage 3D-1 passed QC: %s", len(pass_rows))
    logger.info("Skipped Stage 3C rows: %s", len(skipped_rows))
    logger.info("Output directory: %s", output_dir)
    if not pass_rows:
        logger.warning("No Stage 3D-1 outputs passed QC.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
