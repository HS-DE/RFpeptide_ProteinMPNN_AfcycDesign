from __future__ import annotations

import argparse
import hashlib
import importlib.util
import math
import os
import shutil
from collections import Counter
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

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
STAGE3D1_PATH = SCRIPT_DIR / "24_stage3d1_sidechain_repack.py"
HYDROPHOBIC_AAS = set("AVILMFWY")
AROMATIC_AAS = set("FWY")
BASIC_AAS = set("KR")
ACIDIC_AAS = set("DE")


STAGE4_FIELDS = [
    "stage4_design_id",
    "stage4_run_id",
    "stage4_protocol_identity_sha256",
    "source_stage3d1_design_id",
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
    "stage2_selection_csv",
    "stage2_selection_csv_sha256",
    "stage3_jobs_csv",
    "stage3_jobs_csv_sha256",
    "stage3c_qc_csv",
    "stage3c_qc_csv_sha256",
    "stage3d1_qc_csv",
    "stage3d1_qc_csv_sha256",
    "stage3d1_pass_csv",
    "stage3d1_pass_csv_sha256",
    "source_backbone_pdb",
    "source_backbone_pdb_sha256",
    "input_pdb",
    "input_pdb_sha256",
    "scored_pdb",
    "scored_pdb_sha256",
    "target_chain",
    "peptide_chain",
    "peptide_sequence",
    "peptide_length",
    "net_charge_approx",
    "hydrophobic_fraction",
    "aromatic_count",
    "pro_count",
    "gly_count",
    "sequence_liability_notes",
    "sequence_liability_count",
    "score_mode",
    "pyrosetta_score_status",
    "complex_rosetta_total_score",
    "separated_rosetta_total_score",
    "ddg_proxy_no_repack",
    "ddg_proxy_unit",
    "ddg_proxy_per_peptide_residue",
    "ddg_proxy_per_target_contact",
    "ddg_proxy_per_site_contact",
    "ddg_proxy_status",
    "CMS",
    "SAP",
    "target_contact_status",
    "target_site_recovery_status",
    "hotspot_recovery_status",
    "macrocycle_geometry_status",
    "clash_status",
    "num_target_contacts",
    "num_target_site_contacts",
    "num_hotspot_contacts",
    "peptide_target_min_distance",
    "peptide_site_min_distance",
    "peptide_hotspot_min_distance",
    "closest_target_residue",
    "closest_site_residue",
    "closest_hotspot_residue",
    "macrocycle_terminal_cn_distance",
    "peptide_radius_of_gyration",
    "detached_or_collapsed_flag",
    "stage4_energy_rank",
    "stage4_priority_rank",
    "stage4_validation_selection_rank",
    "stage4_priority_class",
    "stage4_priority_notes",
    "pass_stage4_qc",
    "stage4_failure_reasons",
    "notes",
] + ROUTE_PROVENANCE_FIELDS + SOURCE_ROUTE_PROVENANCE_FIELDS

TOP_VALIDATION_FIELDS = [
    "stage4_validation_selection_rank",
    "stage4_priority_rank",
    "stage4_priority_class",
    "stage4_design_id",
    "stage4_run_id",
    "stage4_protocol_identity_sha256",
    "run_group_id",
    "global_backbone_id",
    "backbone_family_id",
    "source_batch_label",
    "peptide_sequence",
    "peptide_length",
    "scored_pdb",
    "scored_pdb_sha256",
    "ddg_proxy_no_repack",
    "ddg_proxy_unit",
    "ddg_proxy_per_peptide_residue",
    "ddg_proxy_per_target_contact",
    "ddg_proxy_per_site_contact",
    "site_contacts",
    "hotspot_contacts",
    "macrocycle_terminal_cn_distance",
    "clash_status",
    "sequence_liability_notes",
    "reason_selected",
]


def _load_stage3c_module() -> Any:
    spec = importlib.util.spec_from_file_location("stage3c_qc", STAGE3C_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not import Stage 3C QC helper: {STAGE3C_PATH}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _load_stage3d1_module() -> Any:
    spec = importlib.util.spec_from_file_location("stage3d1_qc", STAGE3D1_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not import Stage 3D-1 helper: {STAGE3D1_PATH}")
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


def _require_equal(observed: Any, expected: Any, label: str) -> None:
    if str(observed).strip() != str(expected).strip():
        raise RuntimeError(
            f"{label} mismatch: observed={str(observed).strip()!r}, "
            f"expected={str(expected).strip()!r}"
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
    recorded_path = assert_active_route_path(
        _resolve_mixed_path(str(row.get(path_field, ""))),
        f"{label} {path_field}",
    )
    if recorded_path.resolve() != expected_path.resolve():
        raise RuntimeError(
            f"{label} {path_field} mismatch: observed={recorded_path}, "
            f"expected={expected_path}"
        )
    _require_equal(
        row.get(sha256_field, ""),
        expected_sha256,
        f"{label} {sha256_field}",
    )


def _expected_stage3d1_pass(row: Mapping[str, Any]) -> bool:
    return (
        str(row.get("repack_status", "")) == "success"
        and str(row.get("backbone_geometry_status", "")) == "pass_fixed_backbone"
        and str(row.get("post_target_site_recovery_status", ""))
        == "site_contact_pass"
        and str(row.get("post_hotspot_recovery_status", ""))
        == "hotspot_contact_pass"
        and str(row.get("post_macrocycle_geometry_status", ""))
        == "pass_head_to_tail_macrocycle"
        and str(row.get("post_clash_status", "")) == "pass_no_severe_clash"
        and str(row.get("peptide_sequence", ""))
        == str(row.get("repacked_peptide_sequence", ""))
        and str(row.get("target_sequence", ""))
        == str(row.get("repacked_target_sequence", ""))
        and not str(row.get("stage3d1_failure_reasons", "")).strip()
    )


def _validate_stage3d1_pass_subset(
    *,
    stage3d1_rows: Sequence[Mapping[str, str]],
    stage3d1_pass_rows: Sequence[Mapping[str, str]],
    stage3d1_fields: Sequence[str],
) -> list[dict[str, str]]:
    full_lookup: dict[str, dict[str, str]] = {}
    for raw_row in stage3d1_rows:
        row = dict(raw_row)
        design_id = str(row.get("stage3d1_design_id", "")).strip()
        if not design_id:
            raise RuntimeError("Stage 3D-1 QC row is missing stage3d1_design_id")
        if design_id in full_lookup:
            raise RuntimeError(f"Duplicate Stage 3D-1 design ID: {design_id}")
        observed_pass = str(row.get("pass_stage3d1_qc", "")).strip().lower()
        expected_pass = "true" if _expected_stage3d1_pass(row) else "false"
        if observed_pass != expected_pass:
            raise RuntimeError(
                f"Stage 3D-1 pass flag disagrees with hard QC fields for "
                f"{design_id}: observed={observed_pass}, expected={expected_pass}"
            )
        full_lookup[design_id] = row

    pass_lookup: dict[str, dict[str, str]] = {}
    for raw_row in stage3d1_pass_rows:
        row = dict(raw_row)
        design_id = str(row.get("stage3d1_design_id", "")).strip()
        if not design_id:
            raise RuntimeError(
                "Stage 3D-1 pass-table row is missing stage3d1_design_id"
            )
        if design_id in pass_lookup:
            raise RuntimeError(
                f"Duplicate Stage 3D-1 pass-table design ID: {design_id}"
            )
        pass_lookup[design_id] = row

    expected_pass_ids = {
        design_id
        for design_id, row in full_lookup.items()
        if str(row.get("pass_stage3d1_qc", "")).strip().lower() == "true"
    }
    if set(pass_lookup) != expected_pass_ids:
        missing = sorted(expected_pass_ids - set(pass_lookup))
        extra = sorted(set(pass_lookup) - expected_pass_ids)
        raise RuntimeError(
            "Stage 3D-1 pass table is not the exact pass subset of the full QC "
            f"table: missing={missing[:5]}, extra={extra[:5]}"
        )

    for design_id, pass_row in pass_lookup.items():
        full_row = full_lookup[design_id]
        for field in stage3d1_fields:
            _require_equal(
                pass_row.get(field, ""),
                full_row.get(field, ""),
                f"Stage 3D-1 pass/full table {field} for {design_id}",
            )
    return [pass_lookup[design_id] for design_id in sorted(pass_lookup)]


def _parse_float(value: Any) -> float | None:
    try:
        if str(value).strip() == "":
            return None
        out = float(str(value).strip())
    except ValueError:
        return None
    if math.isnan(out) or math.isinf(out):
        return None
    return out


def _parse_int(value: Any, default: int = 0) -> int:
    try:
        return int(str(value).strip())
    except ValueError:
        return default


def _validate_stage3d1_contract(
    *,
    qchelper: Any,
    stage3d1helper: Any,
    stage3_root: Path,
    stage3_contract: Mapping[str, Any],
    stage3_route_provenance: Mapping[str, str],
    aggregate_sources: Mapping[Path, Mapping[str, str]],
    stage2_selection_csv: Path,
    stage3_jobs_csv: Path,
    stage3c_qc_csv: Path,
    stage3d1_qc_csv: Path,
    stage3d1_pass_csv: Path,
    stage0_target_pdb: Path,
) -> dict[str, Any]:
    run_group_id = str(stage3_contract["run_group_id"])
    safe_run_group = _safe_token(run_group_id)
    stage3d1_dir = (
        stage3_root
        / "05_proteinmpnn_sequences"
        / "stage3d1"
        / safe_run_group
    )
    expected_qc_csv = (
        stage3d1_dir
        / f"FGA_rfpeptides_{safe_run_group}_stage3D1_sidechain_repack_qc.csv"
    )
    expected_pass_csv = (
        stage3d1_dir
        / f"FGA_rfpeptides_{safe_run_group}_stage3D1_sidechain_repack_qc_pass.csv"
    )
    if stage3d1_qc_csv.resolve() != expected_qc_csv.resolve():
        raise RuntimeError(
            "Stage 4 requires the run-group-scoped full Stage 3D-1 QC table: "
            f"observed={stage3d1_qc_csv}, expected={expected_qc_csv}"
        )
    if stage3d1_pass_csv.resolve() != expected_pass_csv.resolve():
        raise RuntimeError(
            "Stage 4 requires the run-group-scoped Stage 3D-1 pass table: "
            f"observed={stage3d1_pass_csv}, expected={expected_pass_csv}"
        )

    stage3d1_qc_sha256 = _sha256_file(stage3d1_qc_csv)
    stage3d1_pass_sha256 = _sha256_file(stage3d1_pass_csv)
    stage3c_selected_rows, _ = stage3d1helper._selected_stage3c_rows(
        stage3_contract["stage3c_rows"],
        set(),
        include_duplicate_sequences=False,
    )
    expected_by_sequence = {
        str(row.get("sequence_design_id", "")).strip(): dict(row)
        for row in stage3c_selected_rows
    }
    if len(expected_by_sequence) != len(stage3c_selected_rows):
        raise RuntimeError(
            "Stage 3D-1 expected Stage 3C sequence IDs are blank or duplicated"
        )

    stage3d1_rows = _read_required_csv(stage3d1_qc_csv)
    stage3d1_pass_rows = _validate_stage3d1_pass_subset(
        stage3d1_rows=stage3d1_rows,
        stage3d1_pass_rows=_read_required_csv(stage3d1_pass_csv),
        stage3d1_fields=stage3d1helper.STAGE3D1_FIELDS,
    )
    if len(stage3d1_rows) != len(expected_by_sequence):
        raise RuntimeError(
            "Stage 3D-1 full QC table is not the complete contact-preserving "
            "unique Stage 3C set: "
            f"rows={len(stage3d1_rows)}, expected={len(expected_by_sequence)}"
        )

    stage2_selection_lookup = stage3_contract["stage2_selection_lookup"]
    stage3_job_lookup = stage3_contract["stage3_job_lookup"]
    full_sequence_ids: set[str] = set()
    full_stage3d1_ids: set[str] = set()
    repacked_paths: set[Path] = set()
    repacked_dir = stage3d1_dir / "repacked_pdbs"
    stage2_selection_sha256 = str(
        stage3_contract["stage2_selection_csv_sha256"]
    )
    stage3_jobs_sha256 = str(stage3_contract["stage3_jobs_csv_sha256"])
    stage3c_qc_sha256 = str(stage3_contract["stage3c_qc_csv_sha256"])
    expected_repack_mode = (
        "side_chain_repack_only_fixed_sequence_no_backbone_minimization_"
        "no_fastrelax"
    )

    for raw_row in stage3d1_rows:
        row = dict(raw_row)
        stage3d1_id = str(row.get("stage3d1_design_id", "")).strip()
        sequence_id = str(row.get("source_sequence_design_id", "")).strip()
        global_backbone_id = str(row.get("global_backbone_id", "")).strip()
        if not stage3d1_id or stage3d1_id in full_stage3d1_ids:
            raise RuntimeError(
                f"Blank or duplicate Stage 3D-1 design ID: {stage3d1_id!r}"
            )
        if not sequence_id or sequence_id in full_sequence_ids:
            raise RuntimeError(
                f"Blank or duplicate Stage 3D-1 source sequence ID: "
                f"{sequence_id!r}"
            )
        full_stage3d1_ids.add(stage3d1_id)
        full_sequence_ids.add(sequence_id)
        stage3c_row = expected_by_sequence.get(sequence_id)
        if stage3c_row is None:
            raise RuntimeError(
                f"Stage 3D-1 row is not authorized by the selected Stage 3C "
                f"set: {sequence_id}"
            )
        backbone_row = stage2_selection_lookup.get(global_backbone_id)
        job_row = stage3_job_lookup.get(global_backbone_id)
        if backbone_row is None or job_row is None:
            raise RuntimeError(
                "Stage 3D-1 row is not linked to the authoritative Stage 2.5 "
                f"and Stage 3 job tables: {global_backbone_id}"
            )

        _require_equal(
            stage3d1_id,
            f"{sequence_id}_stage3d1_repack",
            f"Stage 3D-1 design ID for {sequence_id}",
        )
        _require_equal(
            row.get("backbone_id", ""),
            global_backbone_id,
            f"Stage 3D-1 backbone_id for {sequence_id}",
        )
        _require_equal(
            stage3c_row.get("global_backbone_id", ""),
            global_backbone_id,
            f"Stage 3D-1/Stage 3C global backbone for {sequence_id}",
        )
        for field in [
            "stage3_job_id",
            "run_group_id",
            "protocol_identity_sha256",
            "stage3_mode",
        ]:
            _require_equal(
                row.get(field, ""),
                job_row.get(field, ""),
                f"Stage 3D-1 {field} for {sequence_id}",
            )
        for field in [
            "source_local_design_id",
            "source_batch_label",
            "backbone_family_id",
            "site_label",
            "site_id",
            "target_chain",
            "peptide_chain",
            "peptide_sequence",
        ]:
            _require_equal(
                row.get(field, ""),
                stage3c_row.get(field, ""),
                f"Stage 3D-1 {field} for {sequence_id}",
            )
        for field in SOURCE_ROUTE_PROVENANCE_FIELDS:
            _require_equal(
                row.get(field, ""),
                job_row.get(field, ""),
                f"Stage 3D-1 {field} for {sequence_id}",
            )

        _require_recorded_file(
            row=row,
            path_field="stage2_selection_csv",
            sha256_field="stage2_selection_csv_sha256",
            expected_path=stage2_selection_csv,
            expected_sha256=stage2_selection_sha256,
            label=f"Stage 3D-1 row {sequence_id}",
        )
        _require_recorded_file(
            row=row,
            path_field="stage3_jobs_csv",
            sha256_field="stage3_jobs_csv_sha256",
            expected_path=stage3_jobs_csv,
            expected_sha256=stage3_jobs_sha256,
            label=f"Stage 3D-1 row {sequence_id}",
        )
        _require_recorded_file(
            row=row,
            path_field="stage3c_qc_csv",
            sha256_field="stage3c_qc_csv_sha256",
            expected_path=stage3c_qc_csv,
            expected_sha256=stage3c_qc_sha256,
            label=f"Stage 3D-1 row {sequence_id}",
        )

        input_pdb = assert_active_route_path(
            _resolve_mixed_path(str(row.get("input_pdb", ""))),
            f"Stage 4 Stage 3D-1 input PDB for {sequence_id}",
        )
        expected_input_pdb = assert_active_route_path(
            _resolve_mixed_path(str(stage3c_row.get("relaxed_pdb", ""))),
            f"Stage 4 Stage 3C PDB for {sequence_id}",
        )
        if input_pdb.resolve() != expected_input_pdb.resolve():
            raise RuntimeError(
                f"Stage 3D-1 input PDB does not match Stage 3C for "
                f"{sequence_id}"
            )
        input_sha256 = _sha256_file(input_pdb)
        _require_equal(
            row.get("input_pdb_sha256", ""),
            input_sha256,
            f"Stage 3D-1 input PDB SHA-256 for {sequence_id}",
        )
        _require_equal(
            stage3c_row.get("relaxed_pdb_sha256", ""),
            input_sha256,
            f"Stage 3C PDB SHA-256 for {sequence_id}",
        )
        _require_equal(
            row.get("repack_mode", ""),
            expected_repack_mode,
            f"Stage 3D-1 repack mode for {sequence_id}",
        )
        _require_equal(
            row.get("sequence_duplicate_status", ""),
            "unique_within_global_backbone",
            f"Stage 3D-1 sequence duplicate status for {sequence_id}",
        )

        if str(row.get("repack_status", "")) == "success":
            repacked_pdb = assert_active_route_path(
                _resolve_mixed_path(str(row.get("repacked_pdb", ""))),
                f"Stage 4 repacked PDB for {sequence_id}",
            )
            qchelper._require_within(
                repacked_pdb,
                repacked_dir,
                f"Stage 4 repacked PDB for {sequence_id}",
            )
            expected_repacked = repacked_dir / f"{_safe_token(stage3d1_id)}.pdb"
            if repacked_pdb.resolve() != expected_repacked.resolve():
                raise RuntimeError(
                    f"Stage 3D-1 repacked PDB path mismatch for {sequence_id}: "
                    f"observed={repacked_pdb}, expected={expected_repacked}"
                )
            if repacked_pdb.resolve() in repacked_paths:
                raise RuntimeError(
                    f"Duplicate Stage 3D-1 repacked PDB: {repacked_pdb}"
                )
            repacked_paths.add(repacked_pdb.resolve())
            _require_equal(
                row.get("repacked_pdb_sha256", ""),
                _sha256_file(repacked_pdb),
                f"Stage 3D-1 repacked PDB SHA-256 for {sequence_id}",
            )
        elif str(row.get("pass_stage3d1_qc", "")).lower() == "true":
            raise RuntimeError(
                f"Stage 3D-1 pass row has no successful repack: {sequence_id}"
            )

        qchelper._validate_cached_route_provenance(
            row,
            stage3_route_provenance,
            f"Stage 4 Stage 3D-1 row {sequence_id}",
        )
        qchelper._validate_aggregate_source_membership(
            row,
            aggregate_sources,
            f"Stage 4 Stage 3D-1 row {sequence_id}",
        )

    if full_sequence_ids != set(expected_by_sequence):
        missing = sorted(set(expected_by_sequence) - full_sequence_ids)
        extra = sorted(full_sequence_ids - set(expected_by_sequence))
        raise RuntimeError(
            "Stage 3D-1 full table does not match the Stage 3C selected set: "
            f"missing={missing[:5]}, extra={extra[:5]}"
        )

    pass_target_chains = {
        str(row.get("target_chain", "")).strip() for row in stage3d1_pass_rows
    }
    pass_peptide_chains = {
        str(row.get("peptide_chain", "")).strip() for row in stage3d1_pass_rows
    }
    if len(pass_target_chains) != 1 or len(pass_peptide_chains) != 1:
        raise RuntimeError(
            "Stage 3D-1 pass rows do not share one target chain and one peptide "
            f"chain: target={sorted(pass_target_chains)}, "
            f"peptide={sorted(pass_peptide_chains)}"
        )
    target_chain = next(iter(pass_target_chains))
    peptide_chain = next(iter(pass_peptide_chains))
    stage0_target_residues = list(
        parse_residues(stage0_target_pdb).get(target_chain, [])
    )
    if not stage0_target_residues:
        raise RuntimeError(
            "Manifest-locked Stage 0 target PDB does not contain the Stage "
            f"3D-1 target chain {target_chain}"
        )
    stage0_target_sequence = residue_sequence(stage0_target_residues)

    for row in stage3d1_pass_rows:
        design_id = str(row["stage3d1_design_id"])
        global_backbone_id = str(row["global_backbone_id"])
        repacked_pdb = assert_active_route_path(
            _resolve_mixed_path(str(row["repacked_pdb"])),
            f"Stage 4 pass PDB {design_id}",
        )
        chains = parse_residues(repacked_pdb)
        target_residues = list(chains.get(target_chain, []))
        peptide_residues = list(chains.get(peptide_chain, []))
        if not target_residues or not peptide_residues:
            raise RuntimeError(
                f"Stage 3D-1 pass PDB is missing target/peptide chain: "
                f"{design_id}"
            )
        _require_equal(
            residue_sequence(target_residues),
            stage0_target_sequence,
            f"Stage 3D-1 target sequence for {design_id}",
        )
        _require_equal(
            residue_sequence(peptide_residues),
            row.get("peptide_sequence", ""),
            f"Stage 3D-1 peptide sequence for {design_id}",
        )
        expected_length = _parse_int(
            stage2_selection_lookup[global_backbone_id].get(
                "peptide_length", "0"
            )
        )
        if len(peptide_residues) != expected_length:
            raise RuntimeError(
                f"Stage 3D-1 peptide length mismatch for {design_id}: "
                f"observed={len(peptide_residues)}, expected={expected_length}"
            )

    return {
        "stage3d1_rows": stage3d1_rows,
        "stage3d1_pass_rows": stage3d1_pass_rows,
        "stage3d1_qc_csv_sha256": stage3d1_qc_sha256,
        "stage3d1_pass_csv_sha256": stage3d1_pass_sha256,
        "target_chain": target_chain,
        "peptide_chain": peptide_chain,
    }


def _sequence_properties(sequence: str) -> dict[str, Any]:
    seq = "".join(aa for aa in str(sequence).upper() if aa.isalpha())
    length = len(seq)
    if not seq:
        return {
            "peptide_length": 0,
            "net_charge_approx": "",
            "hydrophobic_fraction": "",
            "aromatic_count": 0,
            "pro_count": 0,
            "gly_count": 0,
            "sequence_liability_notes": "missing_sequence",
            "sequence_liability_count": 1,
        }

    hydrophobic_count = sum(1 for aa in seq if aa in HYDROPHOBIC_AAS)
    aromatic_count = sum(1 for aa in seq if aa in AROMATIC_AAS)
    pro_count = seq.count("P")
    gly_count = seq.count("G")
    net_charge = sum(1.0 for aa in seq if aa in BASIC_AAS)
    net_charge += 0.1 * seq.count("H")
    net_charge -= sum(1.0 for aa in seq if aa in ACIDIC_AAS)
    hydrophobic_fraction = hydrophobic_count / length

    notes = []
    if hydrophobic_fraction >= 0.60:
        notes.append("warning_high_hydrophobic_fraction")
    if pro_count >= max(4, math.ceil(length * 0.30)):
        notes.append("warning_high_pro_fraction")
    if gly_count >= max(4, math.ceil(length * 0.30)):
        notes.append("warning_high_gly_fraction")
    most_common_count = max(Counter(seq).values())
    if most_common_count >= max(5, math.ceil(length * 0.45)):
        notes.append("warning_low_complexity_single_residue_dominant")
    if any(seq[idx] == seq[idx + 1] == seq[idx + 2] == seq[idx + 3] for idx in range(max(0, length - 3))):
        notes.append("warning_low_complexity_four_residue_run")

    return {
        "peptide_length": length,
        "net_charge_approx": round(net_charge, 3),
        "hydrophobic_fraction": round(hydrophobic_fraction, 3),
        "aromatic_count": aromatic_count,
        "pro_count": pro_count,
        "gly_count": gly_count,
        "sequence_liability_notes": ";".join(notes) if notes else "none",
        "sequence_liability_count": len(notes),
    }


def _load_pyrosetta() -> Any:
    try:
        import pyrosetta  # type: ignore
    except ModuleNotFoundError as exc:
        raise RuntimeError(
            "PyRosetta is required for Stage 4 score-only interface scoring. "
            "Run in the proteinmpnn_binder_design environment."
        ) from exc
    return pyrosetta


def _init_pyrosetta(pyrosetta: Any) -> None:
    pyrosetta.init("-beta_nov16 -mute all -use_terminal_residues true")


def _xyz_tuple(xyz: Any) -> tuple[float, float, float]:
    try:
        return float(xyz.x), float(xyz.y), float(xyz.z)
    except TypeError:
        return float(xyz[0]), float(xyz[1]), float(xyz[2])


def _pose_chain(pose: Any, pose_idx: int) -> str:
    pdb_info = pose.pdb_info()
    chain = pdb_info.chain(pose_idx) if pdb_info is not None else ""
    return str(chain).strip() or "_"


def _chain_positions(pose: Any, chain_id: str) -> list[int]:
    return [idx for idx in range(1, pose.total_residue() + 1) if _pose_chain(pose, idx) == chain_id]


def _translate_positions(pyrosetta: Any, pose: Any, positions: Sequence[int], separation_distance: float) -> None:
    atom_id_cls = pyrosetta.rosetta.core.id.AtomID
    xyz_vector_cls = pyrosetta.rosetta.numeric.xyzVector_double_t
    for res_idx in positions:
        residue = pose.residue(res_idx)
        for atom_idx in range(1, residue.natoms() + 1):
            atom_id = atom_id_cls(atom_idx, res_idx)
            x, y, z = _xyz_tuple(pose.xyz(atom_id))
            pose.set_xyz(atom_id, xyz_vector_cls(x + separation_distance, y, z))


def _score_complex_and_separated(
    *,
    pyrosetta: Any,
    input_pdb: Path,
    peptide_chain: str,
    separation_distance: float,
) -> dict[str, Any]:
    pose = pyrosetta.pose_from_pdb(str(input_pdb))
    peptide_positions = _chain_positions(pose, peptide_chain)
    if not peptide_positions:
        raise RuntimeError(f"Peptide chain {peptide_chain} was not found in {input_pdb}")

    scorefxn = pyrosetta.get_fa_scorefxn()
    complex_score = float(scorefxn(pose))
    separated_pose = pose.clone()
    _translate_positions(pyrosetta, separated_pose, peptide_positions, separation_distance)
    separated_score = float(scorefxn(separated_pose))
    return {
        "complex_rosetta_total_score": round(complex_score, 3),
        "separated_rosetta_total_score": round(separated_score, 3),
        "ddg_proxy_no_repack": round(complex_score - separated_score, 3),
    }


def _normalized_ddg_metrics(
    *,
    ddg_proxy: float | None,
    peptide_length: int,
    target_contacts: int,
    site_contacts: int,
) -> dict[str, float | str]:
    if ddg_proxy is None:
        return {
            "ddg_proxy_per_peptide_residue": "",
            "ddg_proxy_per_target_contact": "",
            "ddg_proxy_per_site_contact": "",
        }
    return {
        "ddg_proxy_per_peptide_residue": (
            round(ddg_proxy / peptide_length, 4) if peptide_length > 0 else ""
        ),
        "ddg_proxy_per_target_contact": (
            round(ddg_proxy / target_contacts, 4) if target_contacts > 0 else ""
        ),
        "ddg_proxy_per_site_contact": (
            round(ddg_proxy / site_contacts, 4) if site_contacts > 0 else ""
        ),
    }


def _qc_scored_pdb(
    *,
    qchelper: Any,
    scored_pdb: Path,
    backbone_row: Mapping[str, str],
    site_numbers: set[str],
    hotspot_numbers: set[str],
    args: argparse.Namespace,
) -> dict[str, Any]:
    return qchelper._qc_row_for_relaxed_pdb(
        relaxed_pdb=scored_pdb,
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


def _detached_or_collapsed_flag(qc: Mapping[str, Any], min_peptide_rog: float) -> str:
    flags = []
    if qc.get("target_contact_status") not in {"target_contact_pass", "target_contact_low_count"}:
        flags.append(str(qc.get("target_contact_status", "detached_from_target_crop")))
    if qc.get("target_site_recovery_status") != "site_contact_pass":
        flags.append(str(qc.get("target_site_recovery_status", "site_not_recovered")))
    if qc.get("hotspot_recovery_status") != "hotspot_contact_pass":
        flags.append(str(qc.get("hotspot_recovery_status", "hotspot_not_recovered")))
    rog = _parse_float(qc.get("peptide_radius_of_gyration", ""))
    if rog is not None and rog < min_peptide_rog:
        flags.append("peptide_collapsed_low_radius_of_gyration")
    return ";".join(flags) if flags else "pass_basic_pose_geometry"


def _ddg_proxy_status(ddg_proxy: float | None) -> str:
    if ddg_proxy is None:
        return "not_available"
    if ddg_proxy <= 0.0:
        return "favorable_or_neutral_interface_proxy"
    return "unfavorable_interface_proxy"


def _stage4_failure_reasons(
    *,
    qc: Mapping[str, Any],
    score_status: str,
    detached_flag: str,
    ddg_status: str,
    require_favorable_interface: bool,
) -> list[str]:
    reasons: list[str] = []
    if score_status != "success":
        reasons.append(score_status)
    if qc.get("sequence_status") != "pass_sequence":
        reasons.append(str(qc.get("sequence_status", "sequence_not_pass")))
    if qc.get("target_site_recovery_status") != "site_contact_pass":
        reasons.append(str(qc.get("target_site_recovery_status", "site_not_pass")))
    if qc.get("hotspot_recovery_status") != "hotspot_contact_pass":
        reasons.append(str(qc.get("hotspot_recovery_status", "hotspot_not_pass")))
    if qc.get("macrocycle_geometry_status") != "pass_head_to_tail_macrocycle":
        reasons.append(str(qc.get("macrocycle_geometry_status", "macrocycle_not_pass")))
    if qc.get("clash_status") != "pass_no_severe_clash":
        reasons.append(str(qc.get("clash_status", "clash_not_pass")))
    if detached_flag != "pass_basic_pose_geometry":
        reasons.append(detached_flag)
    if require_favorable_interface and ddg_status != "favorable_or_neutral_interface_proxy":
        reasons.append(ddg_status)
    return reasons


def _assign_energy_ranks(rows: list[dict[str, Any]]) -> None:
    ranked = sorted(
        [row for row in rows if _parse_float(row.get("ddg_proxy_no_repack", "")) is not None],
        key=lambda row: (
            float(row.get("ddg_proxy_no_repack", 999999.0)),
            float(row.get("complex_rosetta_total_score", 999999.0)),
            str(row.get("stage4_design_id", "")),
        ),
    )
    for rank, row in enumerate(ranked, start=1):
        row["stage4_energy_rank"] = rank


def _priority_hard_qc_pass(row: Mapping[str, Any]) -> bool:
    return (
        row.get("macrocycle_geometry_status") == "pass_head_to_tail_macrocycle"
        and row.get("target_site_recovery_status") == "site_contact_pass"
        and row.get("hotspot_recovery_status") == "hotspot_contact_pass"
        and row.get("clash_status") == "pass_no_severe_clash"
        and row.get("detached_or_collapsed_flag") == "pass_basic_pose_geometry"
    )


def _priority_sort_key(row: Mapping[str, Any]) -> tuple[Any, ...]:
    hard_qc_rank = 0 if _priority_hard_qc_pass(row) else 1
    ddg = _parse_float(row.get("ddg_proxy_no_repack", ""))
    return (
        hard_qc_rank,
        ddg if ddg is not None else 999999.0,
        -_parse_int(row.get("num_hotspot_contacts", 0)),
        -_parse_int(row.get("num_target_site_contacts", 0)),
        _parse_int(row.get("sequence_liability_count", 0)),
        str(row.get("stage4_design_id", "")),
    )


def _priority_notes(row: Mapping[str, Any]) -> str:
    notes = []
    if _priority_hard_qc_pass(row):
        notes.append("hard geometry/contact QC passed")
    else:
        notes.append("hard geometry/contact QC not fully passed")
    if row.get("ddg_proxy_status") == "favorable_or_neutral_interface_proxy":
        notes.append("no-repack ddG proxy favorable_or_neutral")
    elif row.get("ddg_proxy_status") == "unfavorable_interface_proxy":
        notes.append("no-repack ddG proxy weak_or_not_clearly_favorable")
    else:
        notes.append("no-repack ddG proxy not_available")
    seq_notes = str(row.get("sequence_liability_notes", ""))
    if seq_notes and seq_notes != "none":
        notes.append(f"sequence warnings: {seq_notes}")
    else:
        notes.append("no obvious sequence liability warning")
    return "; ".join(notes)


def _assign_priority_ranks(rows: list[dict[str, Any]]) -> None:
    ranked = sorted(rows, key=_priority_sort_key)
    for rank, row in enumerate(ranked, start=1):
        row["stage4_priority_rank"] = rank
        row["stage4_validation_selection_rank"] = ""
        row["stage4_priority_class"] = "priority_3_low_priority"
        row["stage4_priority_notes"] = _priority_notes(row)


def _select_validation_source_rows(
    rows: list[dict[str, Any]],
    top_validation_count: int,
    max_per_global_backbone: int,
) -> list[dict[str, Any]]:
    selected: list[dict[str, Any]] = []
    counts: Counter[str] = Counter()
    ranked = sorted(
        rows,
        key=lambda item: _parse_int(
            item.get("stage4_priority_rank", 999999),
            999999,
        ),
    )
    for row in ranked:
        if not _priority_hard_qc_pass(row):
            continue
        backbone_id = str(row.get("global_backbone_id", "")).strip()
        if not backbone_id:
            raise RuntimeError(
                "Stage 4 row is missing global_backbone_id during validation "
                "candidate selection"
            )
        if counts[backbone_id] >= max_per_global_backbone:
            continue
        selected.append(row)
        counts[backbone_id] += 1
        if len(selected) >= top_validation_count:
            break
    for selection_rank, row in enumerate(selected, start=1):
        row["stage4_validation_selection_rank"] = selection_rank
        row["stage4_priority_class"] = (
            "priority_1_validation_ready"
            if selection_rank <= min(3, top_validation_count)
            else "priority_2_backup"
        )
        row["stage4_priority_notes"] = (
            f"{row['stage4_priority_notes']}; diversity-aware validation "
            f"selection rank {selection_rank}"
        )
    return selected


def _top_validation_rows(
    selected: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    top_rows = []
    for row in selected:
        ddg = _parse_float(row.get("ddg_proxy_no_repack", ""))
        if ddg is None:
            energy_note = "interface-energy proxy is not available"
        elif ddg <= 0:
            energy_note = "interface-energy proxy is favorable or neutral"
        else:
            energy_note = "interface-energy proxy is weak or not clearly favorable"
        top_rows.append(
            {
                "stage4_validation_selection_rank": row.get(
                    "stage4_validation_selection_rank", ""
                ),
                "stage4_priority_rank": row.get("stage4_priority_rank", ""),
                "stage4_priority_class": row.get("stage4_priority_class", ""),
                "stage4_design_id": row.get("stage4_design_id", ""),
                "stage4_run_id": row.get("stage4_run_id", ""),
                "stage4_protocol_identity_sha256": row.get(
                    "stage4_protocol_identity_sha256", ""
                ),
                "run_group_id": row.get("run_group_id", ""),
                "global_backbone_id": row.get("global_backbone_id", ""),
                "backbone_family_id": row.get("backbone_family_id", ""),
                "source_batch_label": row.get("source_batch_label", ""),
                "peptide_sequence": row.get("peptide_sequence", ""),
                "peptide_length": row.get("peptide_length", ""),
                "scored_pdb": row.get("scored_pdb", ""),
                "scored_pdb_sha256": row.get("scored_pdb_sha256", ""),
                "ddg_proxy_no_repack": row.get("ddg_proxy_no_repack", ""),
                "ddg_proxy_unit": row.get("ddg_proxy_unit", ""),
                "ddg_proxy_per_peptide_residue": row.get(
                    "ddg_proxy_per_peptide_residue", ""
                ),
                "ddg_proxy_per_target_contact": row.get(
                    "ddg_proxy_per_target_contact", ""
                ),
                "ddg_proxy_per_site_contact": row.get(
                    "ddg_proxy_per_site_contact", ""
                ),
                "site_contacts": row.get("num_target_site_contacts", ""),
                "hotspot_contacts": row.get("num_hotspot_contacts", ""),
                "macrocycle_terminal_cn_distance": row.get("macrocycle_terminal_cn_distance", ""),
                "clash_status": row.get("clash_status", ""),
                "sequence_liability_notes": row.get("sequence_liability_notes", ""),
                "reason_selected": (
                    "geometry/contact/macrocycle QC passed; no severe clash; "
                    f"{energy_note}; selected for downstream structure validation."
                ),
            }
        )
    return top_rows


def _status_lines(rows: list[Mapping[str, Any]], field: str) -> str:
    counts = Counter(str(row.get(field, "")) for row in rows)
    if not counts:
        return "- none: 0"
    return "\n".join(f"- {key or 'blank'}: {counts[key]}" for key in sorted(counts))


def _summary_markdown(
    *,
    rows: list[Mapping[str, Any]],
    top_validation_rows: list[Mapping[str, Any]],
    output_dir: Path,
    args: argparse.Namespace,
    run_group_id: str,
    stage3d1_full_count: int,
) -> str:
    pass_rows = [row for row in rows if row.get("pass_stage4_qc") == "true"]
    top_rows = sorted(
        rows,
        key=lambda row: (
            row.get("pass_stage4_qc") != "true",
            float(row.get("stage4_priority_rank") or 999999.0),
            float(row.get("peptide_hotspot_min_distance") or 999.0),
            str(row.get("stage4_design_id", "")),
        ),
    )[: args.top_report]
    ddg_values = [_parse_float(row.get("ddg_proxy_no_repack", "")) for row in rows]
    ddg_values = [value for value in ddg_values if value is not None]
    all_ddg_positive = bool(ddg_values) and all(value > 0.0 for value in ddg_values)
    ddg_interpretation = (
        "All parsed ddg_proxy_no_repack values are positive. Geometry/contact/macrocycle QC passed, "
        "but the interface-energy proxy is weak or not clearly favorable; selected rows are only inputs "
        "for downstream structure validation."
        if all_ddg_positive
        else "ddg_proxy_no_repack is used only as a lightweight ranking proxy for downstream validation triage."
    )
    columns = [
        "stage4_design_id",
        "global_backbone_id",
        "backbone_family_id",
        "peptide_sequence",
        "ddg_proxy_no_repack",
        "ddg_proxy_per_peptide_residue",
        "ddg_proxy_status",
        "num_target_site_contacts",
        "num_hotspot_contacts",
        "peptide_hotspot_min_distance",
        "macrocycle_terminal_cn_distance",
        "clash_status",
        "detached_or_collapsed_flag",
        "stage4_energy_rank",
        "stage4_priority_rank",
        "stage4_priority_class",
        "sequence_liability_notes",
        "pass_stage4_qc",
        "stage4_failure_reasons",
    ]
    return f"""# FGA RFpeptides Stage 4 Rosetta Score-Only Interface Scoring

Status: Stage 3D-1 repack-only structures were scored without ordinary
FastRelax, backbone movement, minimization, sequence design, or additional
side-chain repacking.

Important rule: Stage 4 ranks and filters structures for downstream validation.
It does not make final peptide candidates. Ordinary unconstrained FastRelax is
not used here because it opened/damaged the RFpep_Site_2 macrocycle control.

`ddg_proxy_no_repack` is a no-repack separated-state proxy. It is not a formal
experimental binding energy and is not a final affinity judgment.

Method:

```text
run_group_id: {run_group_id}
input_stage3d1_pass_csv: {args.stage3d1_pass_csv}
input_stage3d1_full_qc_csv: {args.stage3d1_qc_csv}
score_mode: score_only_no_repack_no_minimization
separation_distance_A_for_ddg_proxy: {args.separation_distance:g}
require_favorable_interface_proxy: {args.require_favorable_interface}
top_validation_candidates: {args.top_validation_candidates}
max_top_candidates_per_global_backbone: {args.max_top_per_global_backbone}
CMS: not_available
SAP: not_available
```

Output directory:

```text
{output_dir}
```

## Counts

```text
full_stage3d1_rows_validated: {stage3d1_full_count}
input_rows_scored: {len(rows)}
pass_stage4_qc: {len(pass_rows)}
```

## QC Status Counts

Target-site recovery:

{_status_lines(rows, "target_site_recovery_status")}

Hotspot recovery:

{_status_lines(rows, "hotspot_recovery_status")}

Macrocycle geometry:

{_status_lines(rows, "macrocycle_geometry_status")}

Clash:

{_status_lines(rows, "clash_status")}

ddG proxy:

{_status_lines(rows, "ddg_proxy_status")}

Sequence liability notes:

{_status_lines(rows, "sequence_liability_notes")}

## ddG Proxy Interpretation

{ddg_interpretation}

The proxy is reported in Rosetta Energy Units (REU), not kcal/mol. Per-residue
and per-contact values are normalization aids only and are not affinity
estimates.

## Recommended Top Validation Candidates

{rows_to_markdown(top_validation_rows, TOP_VALIDATION_FIELDS, "No top validation candidates were selected.")}

## Ranked Stage 4 Rows

{rows_to_markdown(top_rows, columns, "No Stage 4 rows were generated.")}
"""


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Stage 4A-v2 score-only interface scoring for the authoritative "
            "run-group-scoped Stage 3D-1 pass set."
        )
    )
    parser.add_argument("--stage0-root", required=True)
    parser.add_argument("--stage3-root", required=True)
    parser.add_argument("--project-config", required=True)
    parser.add_argument("--stage2-selection-csv", required=True)
    parser.add_argument("--stage3-jobs-csv", required=True)
    parser.add_argument("--stage3c-qc-csv", required=True)
    parser.add_argument("--stage3d1-qc-csv", required=True)
    parser.add_argument("--stage3d1-pass-csv", required=True)
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
    parser.add_argument("--validate-inputs-only", action="store_true")
    parser.add_argument("--separation-distance", type=float, default=1000.0)
    parser.add_argument("--require-favorable-interface", action="store_true")
    parser.add_argument("--min-peptide-rog", type=float, default=1.0)
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
    parser.add_argument("--top-report", type=int, default=50)
    parser.add_argument("--top-pymol", type=int, default=20)
    parser.add_argument("--top-validation-candidates", type=int, default=5)
    parser.add_argument(
        "--max-top-per-global-backbone",
        type=int,
        default=1,
        help=(
            "Maximum top-validation rows selected from one global backbone. "
            "The full score table always retains every row."
        ),
    )
    parser.add_argument(
        "--pymol-path-style",
        choices=["windows", "wsl", "native"],
        default="windows",
        help="Path style for generated PyMOL review scripts.",
    )
    args = parser.parse_args()

    logger = setup_logger("25_stage4_rosetta_interface_scoring")
    append_run_header(logger, "25_stage4_rosetta_interface_scoring.py")

    if args.separation_distance <= 0:
        raise RuntimeError("--separation-distance must be > 0")
    if args.contact_cutoff <= 0:
        raise RuntimeError("--contact-cutoff must be > 0")
    if args.site_near_distance < args.contact_cutoff:
        raise RuntimeError("--site-near-distance must be >= --contact-cutoff")
    if args.hotspot_near_distance < args.contact_cutoff:
        raise RuntimeError("--hotspot-near-distance must be >= --contact-cutoff")
    if args.macrocycle_warn_distance < args.macrocycle_pass_distance:
        raise RuntimeError(
            "--macrocycle-warn-distance must be >= "
            "--macrocycle-pass-distance"
        )
    if args.top_validation_candidates <= 0:
        raise RuntimeError("--top-validation-candidates must be > 0")
    if args.max_top_per_global_backbone <= 0:
        raise RuntimeError("--max-top-per-global-backbone must be > 0")

    qchelper = _load_stage3c_module()
    stage3d1helper = _load_stage3d1_module()
    stage0_root = assert_active_route_path(
        _resolve_mixed_path(args.stage0_root),
        "Stage 25 Stage 0 root",
    )
    stage3_root = assert_active_route_path(
        _resolve_mixed_path(args.stage3_root),
        "Stage 25 Stage 3 root",
    )
    output_root = assert_active_route_path(
        _resolve_mixed_path(args.output_root)
        if args.output_root
        else stage3_root,
        "Stage 25 output root",
        must_exist=False,
    )
    stage2_selection_csv = assert_active_route_path(
        _resolve_mixed_path(args.stage2_selection_csv),
        "Stage 25 Stage 2.5 selection CSV",
    )
    stage3_jobs_csv = assert_active_route_path(
        _resolve_mixed_path(args.stage3_jobs_csv),
        "Stage 25 Stage 3 jobs CSV",
    )
    stage3c_qc_csv = assert_active_route_path(
        _resolve_mixed_path(args.stage3c_qc_csv),
        "Stage 25 Stage 3C QC CSV",
    )
    stage3d1_qc_csv = assert_active_route_path(
        _resolve_mixed_path(args.stage3d1_qc_csv),
        "Stage 25 Stage 3D-1 full QC CSV",
    )
    stage3d1_pass_csv = assert_active_route_path(
        _resolve_mixed_path(args.stage3d1_pass_csv),
        "Stage 25 Stage 3D-1 pass CSV",
    )
    for path, label in [
        (stage2_selection_csv, "Stage 2.5 selection CSV"),
        (stage3_jobs_csv, "Stage 3 jobs CSV"),
        (stage3c_qc_csv, "Stage 3C QC CSV"),
        (stage3d1_qc_csv, "Stage 3D-1 full QC CSV"),
        (stage3d1_pass_csv, "Stage 3D-1 pass CSV"),
    ]:
        qchelper._require_within(path, stage3_root, f"Stage 25 {label}")
    args.stage3d1_qc_csv = str(stage3d1_qc_csv)
    args.stage3d1_pass_csv = str(stage3d1_pass_csv)

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
    stage3_contract = stage3d1helper._validate_authoritative_stage3_inputs(
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
    stage3d1_contract = _validate_stage3d1_contract(
        qchelper=qchelper,
        stage3d1helper=stage3d1helper,
        stage3_root=stage3_root,
        stage3_contract=stage3_contract,
        stage3_route_provenance=stage3_route_provenance,
        aggregate_sources=aggregate_sources,
        stage2_selection_csv=stage2_selection_csv,
        stage3_jobs_csv=stage3_jobs_csv,
        stage3c_qc_csv=stage3c_qc_csv,
        stage3d1_qc_csv=stage3d1_qc_csv,
        stage3d1_pass_csv=stage3d1_pass_csv,
        stage0_target_pdb=stage0_target_pdb,
    )
    run_group_id = str(stage3_contract["run_group_id"])
    stage3d1_pass_rows = stage3d1_contract["stage3d1_pass_rows"]
    logger.info(
        "Authoritative Stage 3 jobs validated: %s",
        len(stage3_contract["stage3_job_lookup"]),
    )
    logger.info(
        "Complete Stage 3C rows validated: %s",
        stage3_contract["expected_output_count"],
    )
    logger.info(
        "Complete Stage 3D-1 QC rows validated: %s",
        len(stage3d1_contract["stage3d1_rows"]),
    )
    logger.info(
        "Stage 3D-1 pass rows selected for Stage 4: %s",
        len(stage3d1_pass_rows),
    )
    if args.validate_inputs_only:
        logger.info(
            "Input validation completed; PyRosetta was not loaded and no "
            "Stage 4 outputs were written."
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

    stage4_protocol_payload = {
        "stage": "stage4A_v2_score_only",
        "run_group_id": run_group_id,
        "stage3d1_qc_csv_sha256": stage3d1_contract[
            "stage3d1_qc_csv_sha256"
        ],
        "stage3d1_pass_csv_sha256": stage3d1_contract[
            "stage3d1_pass_csv_sha256"
        ],
        "score_mode": "score_only_no_repack_no_minimization",
        "separation_distance_A": args.separation_distance,
        "require_favorable_interface": args.require_favorable_interface,
        "min_peptide_rog_A": args.min_peptide_rog,
        "contact_cutoff_A": args.contact_cutoff,
        "site_near_distance_A": args.site_near_distance,
        "hotspot_near_distance_A": args.hotspot_near_distance,
        "severe_clash_distance_A": args.severe_clash_distance,
        "min_target_contacts": args.min_target_contacts,
        "min_site_contacts": args.min_site_contacts,
        "min_hotspot_contacts": args.min_hotspot_contacts,
        "macrocycle_pass_distance_A": args.macrocycle_pass_distance,
        "macrocycle_warn_distance_A": args.macrocycle_warn_distance,
        "forbidden_aas": args.forbidden_aas,
        "top_validation_candidates": args.top_validation_candidates,
        "max_top_per_global_backbone": args.max_top_per_global_backbone,
        "script_sha256": _sha256_file(Path(__file__).resolve()),
    }
    stage4_protocol_identity_sha256 = canonical_json_sha256(
        stage4_protocol_payload
    )
    stage4_run_id = (
        f"stage4A_v2_{_safe_token(run_group_id)}_"
        f"{stage4_protocol_identity_sha256[:12]}"
    )
    output_dir = (
        output_root
        / "06_rosetta_scoring"
        / _safe_token(run_group_id)
        / _safe_token(stage4_run_id)
    )
    output_pdb_dir = output_dir / "scored_pdbs"
    output_prefix = f"FGA_rfpeptides_{_safe_token(stage4_run_id)}"

    pyrosetta = _load_pyrosetta()
    _init_pyrosetta(pyrosetta)
    output_pdb_dir.mkdir(parents=True, exist_ok=True)

    result_rows: list[dict[str, Any]] = []
    pymol_rows: list[dict[str, Any]] = []
    stage2_selection_lookup = stage3_contract["stage2_selection_lookup"]
    target_chain = str(stage3d1_contract["target_chain"])
    peptide_chain = str(stage3d1_contract["peptide_chain"])

    for input_row in stage3d1_pass_rows:
        global_backbone_id = str(
            input_row.get("global_backbone_id", "")
        ).strip()
        backbone_row = stage2_selection_lookup.get(global_backbone_id)
        if backbone_row is None:
            raise RuntimeError(
                f"Missing Stage 2.5 row for global backbone: "
                f"{global_backbone_id}"
            )

        input_pdb = assert_active_route_path(
            _resolve_mixed_path(str(input_row.get("repacked_pdb", ""))),
            f"Stage 25 repacked PDB for {global_backbone_id}",
        )
        input_pdb_sha256 = _sha256_file(input_pdb)
        _require_equal(
            input_row.get("repacked_pdb_sha256", ""),
            input_pdb_sha256,
            f"Stage 25 input PDB SHA-256 for {global_backbone_id}",
        )
        source_stage3d1_id = str(
            input_row.get("stage3d1_design_id", "")
        ).strip()
        stage4_id = (
            f"{source_stage3d1_id}_stage4A_v2_"
            f"{stage4_protocol_identity_sha256[:12]}"
        )
        scored_pdb = output_pdb_dir / f"{_safe_token(stage4_id)}.pdb"
        shutil.copy2(input_pdb, scored_pdb)
        scored_pdb_sha256 = _sha256_file(scored_pdb)
        if scored_pdb_sha256 != input_pdb_sha256:
            raise RuntimeError(
                f"Stage 4 score-only PDB copy changed bytes for {stage4_id}"
            )

        score_status = "success"
        notes = (
            "Score-only; no FastRelax, repack, sequence change, backbone "
            "movement, or minimization."
        )
        try:
            score_metrics = _score_complex_and_separated(
                pyrosetta=pyrosetta,
                input_pdb=scored_pdb,
                peptide_chain=peptide_chain,
                separation_distance=args.separation_distance,
            )
        except Exception as exc:
            score_status = f"failed:{exc.__class__.__name__}"
            score_metrics = {
                "complex_rosetta_total_score": "",
                "separated_rosetta_total_score": "",
                "ddg_proxy_no_repack": "",
            }
            notes = str(exc)

        qc = _qc_scored_pdb(
            qchelper=qchelper,
            scored_pdb=scored_pdb,
            backbone_row=backbone_row,
            site_numbers=site_numbers,
            hotspot_numbers=hotspot_numbers,
            args=args,
        )
        pymol_rows.append(qc)
        ddg_proxy = _parse_float(score_metrics.get("ddg_proxy_no_repack", ""))
        ddg_status = _ddg_proxy_status(ddg_proxy)
        detached_flag = _detached_or_collapsed_flag(
            qc,
            args.min_peptide_rog,
        )
        failure_reasons = _stage4_failure_reasons(
            qc=qc,
            score_status=score_status,
            detached_flag=detached_flag,
            ddg_status=ddg_status,
            require_favorable_interface=args.require_favorable_interface,
        )
        peptide_sequence = str(
            input_row.get(
                "peptide_sequence",
                qc.get("peptide_sequence", ""),
            )
        ).strip()
        sequence_metrics = _sequence_properties(peptide_sequence)
        ddg_normalized = _normalized_ddg_metrics(
            ddg_proxy=ddg_proxy,
            peptide_length=_parse_int(sequence_metrics["peptide_length"]),
            target_contacts=_parse_int(qc.get("num_target_contacts", 0)),
            site_contacts=_parse_int(
                qc.get("num_target_site_contacts", 0)
            ),
        )

        result_rows.append(
            {
                "stage4_design_id": stage4_id,
                "stage4_run_id": stage4_run_id,
                "stage4_protocol_identity_sha256": (
                    stage4_protocol_identity_sha256
                ),
                "source_stage3d1_design_id": source_stage3d1_id,
                "source_sequence_design_id": input_row.get(
                    "source_sequence_design_id", ""
                ),
                "stage3_job_id": input_row.get("stage3_job_id", ""),
                "run_group_id": run_group_id,
                "protocol_identity_sha256": input_row.get(
                    "protocol_identity_sha256", ""
                ),
                "stage3_mode": input_row.get("stage3_mode", ""),
                "global_backbone_id": global_backbone_id,
                "source_local_design_id": input_row.get(
                    "source_local_design_id", ""
                ),
                "source_batch_label": input_row.get(
                    "source_batch_label", ""
                ),
                "backbone_family_id": input_row.get(
                    "backbone_family_id", ""
                ),
                "backbone_id": global_backbone_id,
                "site_label": backbone_row.get("site_label", ""),
                "site_id": backbone_row.get("site_id", ""),
                "stage2_selection_csv": stage2_selection_csv,
                "stage2_selection_csv_sha256": stage3_contract[
                    "stage2_selection_csv_sha256"
                ],
                "stage3_jobs_csv": stage3_jobs_csv,
                "stage3_jobs_csv_sha256": stage3_contract[
                    "stage3_jobs_csv_sha256"
                ],
                "stage3c_qc_csv": stage3c_qc_csv,
                "stage3c_qc_csv_sha256": stage3_contract[
                    "stage3c_qc_csv_sha256"
                ],
                "stage3d1_qc_csv": stage3d1_qc_csv,
                "stage3d1_qc_csv_sha256": stage3d1_contract[
                    "stage3d1_qc_csv_sha256"
                ],
                "stage3d1_pass_csv": stage3d1_pass_csv,
                "stage3d1_pass_csv_sha256": stage3d1_contract[
                    "stage3d1_pass_csv_sha256"
                ],
                "source_backbone_pdb": input_row.get(
                    "source_backbone_pdb", ""
                ),
                "source_backbone_pdb_sha256": input_row.get(
                    "source_backbone_pdb_sha256", ""
                ),
                "input_pdb": input_pdb,
                "input_pdb_sha256": input_pdb_sha256,
                "scored_pdb": scored_pdb,
                "scored_pdb_sha256": scored_pdb_sha256,
                "target_chain": target_chain,
                "peptide_chain": peptide_chain,
                "peptide_sequence": peptide_sequence,
                **sequence_metrics,
                "score_mode": "score_only_no_repack_no_minimization",
                "pyrosetta_score_status": score_status,
                **score_metrics,
                "ddg_proxy_unit": "REU",
                **ddg_normalized,
                "ddg_proxy_status": ddg_status,
                "CMS": "not_available",
                "SAP": "not_available",
                "target_contact_status": qc.get(
                    "target_contact_status", ""
                ),
                "target_site_recovery_status": qc.get(
                    "target_site_recovery_status", ""
                ),
                "hotspot_recovery_status": qc.get(
                    "hotspot_recovery_status", ""
                ),
                "macrocycle_geometry_status": qc.get(
                    "macrocycle_geometry_status", ""
                ),
                "clash_status": qc.get("clash_status", ""),
                "num_target_contacts": qc.get("num_target_contacts", ""),
                "num_target_site_contacts": qc.get(
                    "num_target_site_contacts", ""
                ),
                "num_hotspot_contacts": qc.get(
                    "num_hotspot_contacts", ""
                ),
                "peptide_target_min_distance": qc.get(
                    "peptide_target_min_distance", ""
                ),
                "peptide_site_min_distance": qc.get(
                    "peptide_site_min_distance", ""
                ),
                "peptide_hotspot_min_distance": qc.get(
                    "peptide_hotspot_min_distance", ""
                ),
                "closest_target_residue": qc.get(
                    "closest_target_residue", ""
                ),
                "closest_site_residue": qc.get(
                    "closest_site_residue", ""
                ),
                "closest_hotspot_residue": qc.get(
                    "closest_hotspot_residue", ""
                ),
                "macrocycle_terminal_cn_distance": qc.get(
                    "macrocycle_terminal_cn_distance", ""
                ),
                "peptide_radius_of_gyration": qc.get(
                    "peptide_radius_of_gyration", ""
                ),
                "detached_or_collapsed_flag": detached_flag,
                "stage4_energy_rank": "",
                "stage4_priority_rank": "",
                "stage4_validation_selection_rank": "",
                "stage4_priority_class": "",
                "stage4_priority_notes": "",
                "pass_stage4_qc": (
                    "true" if not failure_reasons else "false"
                ),
                "stage4_failure_reasons": ";".join(failure_reasons),
                "notes": notes,
                **output_route_provenance,
                **{
                    field: str(input_row.get(field, ""))
                    for field in SOURCE_ROUTE_PROVENANCE_FIELDS
                },
            }
        )

    _assign_energy_ranks(result_rows)
    _assign_priority_ranks(result_rows)
    selected_validation_rows = _select_validation_source_rows(
        result_rows,
        args.top_validation_candidates,
        args.max_top_per_global_backbone,
    )
    top_validation_rows = _top_validation_rows(selected_validation_rows)
    pass_rows = [
        row for row in result_rows if row.get("pass_stage4_qc") == "true"
    ]
    write_csv(
        output_dir / f"{output_prefix}_rosetta_interface_scores.csv",
        result_rows,
        STAGE4_FIELDS,
    )
    write_csv(
        output_dir / f"{output_prefix}_rosetta_interface_scores_pass.csv",
        pass_rows,
        STAGE4_FIELDS,
    )
    write_csv(
        output_dir / f"{output_prefix}_top_validation_candidates.csv",
        top_validation_rows,
        TOP_VALIDATION_FIELDS,
    )
    write_markdown(
        output_dir / f"{output_prefix}_rosetta_interface_scores.md",
        _summary_markdown(
            rows=result_rows,
            top_validation_rows=top_validation_rows,
            output_dir=output_dir,
            args=args,
            run_group_id=run_group_id,
            stage3d1_full_count=len(stage3d1_contract["stage3d1_rows"]),
        ),
    )
    qchelper._write_pymol_review(
        rows=pymol_rows,
        output_path=(
            output_dir
            / f"{_safe_token(manifest_site_label)}_"
            f"{_safe_token(stage4_run_id)}_score_only_review.pml"
        ),
        target_chain=target_chain,
        peptide_chain=peptide_chain,
        site_numbers=site_numbers,
        hotspot_numbers=hotspot_numbers,
        top_n=args.top_pymol,
        pymol_path_style=args.pymol_path_style,
    )

    logger.info("Stage 4 input rows scored: %s", len(result_rows))
    logger.info("Stage 4 passed QC: %s", len(pass_rows))
    logger.info(
        "Stage 4 top validation candidates: %s",
        len(top_validation_rows),
    )
    logger.info("Stage 4 run ID: %s", stage4_run_id)
    logger.info("Output directory: %s", output_dir)
    if not pass_rows:
        logger.warning("No Stage 4 outputs passed score-only QC.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
