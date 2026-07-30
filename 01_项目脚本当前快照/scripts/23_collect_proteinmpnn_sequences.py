from __future__ import annotations

import argparse
import hashlib
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
    write_fasta,
    write_markdown,
    write_route_manifest,
)
from pdb_utils import ca_coord, centroid, distance, parse_residues, residue_sequence


STAGE3C_FIELDS = [
    "sequence_design_id",
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
    "stage3_input_pdb",
    "stage3_input_pdb_sha256",
    "relaxed_pdb",
    "relaxed_pdb_sha256",
    "source_backbone_pdb",
    "source_backbone_pdb_sha256",
    "expected_stage3b_outputs_for_backbone",
    "observed_stage3b_outputs_for_backbone",
    "file_status",
    "parse_status",
    "target_chain",
    "peptide_chain",
    "peptide_sequence",
    "peptide_length",
    "expected_peptide_length",
    "forbidden_aas",
    "forbidden_aas_present",
    "sequence_status",
    "target_residue_count",
    "target_residue_numbering_status",
    "target_residue_numbering_offset",
    "mapped_target_site_pdb_residue_numbers",
    "mapped_hotspot_pdb_residue_numbers",
    "target_site_residue_count",
    "hotspot_residue_count",
    "num_target_contacts",
    "num_target_site_contacts",
    "num_hotspot_contacts",
    "peptide_target_min_distance",
    "peptide_site_min_distance",
    "peptide_hotspot_min_distance",
    "closest_target_residue",
    "closest_site_residue",
    "closest_hotspot_residue",
    "target_contact_status",
    "target_site_recovery_status",
    "hotspot_recovery_status",
    "macrocycle_terminal_cn_distance",
    "macrocycle_geometry_status",
    "clash_status",
    "peptide_ca_centroid_x",
    "peptide_ca_centroid_y",
    "peptide_ca_centroid_z",
    "peptide_radius_of_gyration",
    "rosetta_score_status",
    "rosetta_pose_total",
    "rosetta_peptide_total_sum",
    "rosetta_target_total_sum",
    "pass_stage3c_qc",
    "qc_failure_reasons",
    "qc_notes",
] + ROUTE_PROVENANCE_FIELDS + SOURCE_ROUTE_PROVENANCE_FIELDS


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
    text = str(value).strip()
    if not text:
        return resolve_path(text)
    text = text.replace("\\", "/")
    if os.name == "nt" and text.startswith("/mnt/") and len(text) > 7 and text[6] == "/":
        return Path(f"{text[5].upper()}:/{text[7:]}")
    if os.name != "nt" and len(text) >= 3 and text[1] == ":" and text[2] == "/":
        return Path(f"/mnt/{text[0].lower()}{text[2:]}")
    path = Path(text)
    if path.is_absolute():
        return path
    return resolve_path(path)


def _format_pymol_path(path: str | Path, style: str) -> str:
    text = str(path).replace("\\", "/")
    if style == "windows":
        if text.startswith("/mnt/") and len(text) > 7 and text[6] == "/":
            return f"{text[5].upper()}:/{text[7:]}"
        return text
    if style == "wsl":
        if len(text) >= 3 and text[1] == ":" and text[2] == "/":
            return f"/mnt/{text[0].lower()}{text[2:]}"
        return text
    if style == "native":
        return text
    raise RuntimeError(f"Unsupported PyMOL path style: {style}")


def _parse_int(value: Any, default: int = 0) -> int:
    try:
        return int(str(value).strip())
    except ValueError:
        return default


def _round_float(value: float | None, digits: int = 3) -> float | str:
    if value is None or math.isinf(value) or math.isnan(value):
        return ""
    return round(value, digits)


def _residue_sort_key(value: str) -> tuple[int, str]:
    digits = "".join(ch for ch in str(value) if ch.isdigit() or ch == "-")
    suffix = "".join(ch for ch in str(value) if not (ch.isdigit() or ch == "-"))
    return _parse_int(digits), suffix


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


def _strict_lookup_rows(
    rows: Iterable[Mapping[str, str]],
    key: str,
    label: str,
) -> dict[str, dict[str, str]]:
    lookup: dict[str, dict[str, str]] = {}
    for row in rows:
        item = dict(row)
        value = str(item.get(key, "")).strip()
        if not value:
            raise RuntimeError(f"{label} row is missing {key}")
        if "," in value:
            raise RuntimeError(f"{label} row contains aggregate {key}: {value}")
        if value in lookup:
            raise RuntimeError(f"Duplicate {label} {key}: {value}")
        lookup[value] = item
    return lookup


def _single_value(rows: Sequence[Mapping[str, str]], field: str, label: str) -> str:
    values = {str(row.get(field, "")).strip() for row in rows if str(row.get(field, "")).strip()}
    if len(values) != 1:
        raise RuntimeError(f"Expected exactly one non-empty {field} in {label}; found {sorted(values)}")
    return next(iter(values))


def _require_within(path: Path, root: Path, label: str) -> None:
    try:
        path.resolve().relative_to(root.resolve())
    except ValueError as exc:
        raise RuntimeError(f"{label} is outside its required root: path={path}, root={root}") from exc


def _validate_cached_route_provenance(
    row: Mapping[str, Any],
    expected: Mapping[str, str],
    label: str,
) -> None:
    observed_manifest = _resolve_mixed_path(str(row.get("route_manifest_path", ""))).resolve()
    expected_manifest = _resolve_mixed_path(str(expected["route_manifest_path"])).resolve()
    if observed_manifest != expected_manifest:
        raise RuntimeError(
            f"{label} route provenance mismatch for route_manifest_path: "
            f"observed={observed_manifest}, expected={expected_manifest}"
        )
    for field in ROUTE_PROVENANCE_FIELDS:
        if field == "route_manifest_path":
            continue
        observed = str(row.get(field, "")).strip()
        if observed != str(expected[field]):
            raise RuntimeError(
                f"{label} route provenance mismatch for {field}: "
                f"observed={observed!r}, expected={expected[field]!r}"
            )


def _aggregate_source_manifest_lookup(
    route_manifest: Mapping[str, Any],
) -> dict[Path, dict[str, str]]:
    records = route_manifest.get("source_route_manifests")
    if not isinstance(records, list) or not records:
        raise RuntimeError("Stage 23 requires a Stage 2.5 aggregate route manifest with source_route_manifests")
    expected_count = int(route_manifest.get("stage2_5_source_count", len(records)))
    if len(records) != expected_count:
        raise RuntimeError(
            f"Stage 23 aggregate source manifest count mismatch: records={len(records)}, expected={expected_count}"
        )

    lookup: dict[Path, dict[str, str]] = {}
    for index, record in enumerate(records, start=1):
        if not isinstance(record, Mapping):
            raise RuntimeError(f"Stage 23 aggregate source manifest record {index} is invalid")
        required = ["run_id", "batch_id", "manifest_path", "manifest_sha256"]
        missing = [field for field in required if not str(record.get(field, "")).strip()]
        if missing:
            raise RuntimeError(f"Stage 23 aggregate source manifest record {index} is missing {missing}")
        manifest_path = assert_active_route_path(
            _resolve_mixed_path(str(record["manifest_path"])),
            f"Stage 23 aggregate source route manifest {index}",
        ).resolve()
        if manifest_path in lookup:
            raise RuntimeError(f"Stage 23 aggregate route lists a source manifest more than once: {manifest_path}")
        actual_sha256 = _sha256_file(manifest_path)
        expected_sha256 = str(record["manifest_sha256"]).strip()
        if actual_sha256 != expected_sha256:
            raise RuntimeError(
                f"Stage 23 aggregate source route manifest SHA-256 mismatch: {manifest_path}"
            )
        _, source_manifest, loaded_sha256 = load_route_manifest(manifest_path.parent)
        if loaded_sha256 != expected_sha256:
            raise RuntimeError(f"Stage 23 source route manifest changed during validation: {manifest_path}")
        if str(source_manifest["run_id"]) != str(record["run_id"]).strip():
            raise RuntimeError(f"Stage 23 aggregate source run_id mismatch: {manifest_path}")
        if str(source_manifest["batch_id"]) != str(record["batch_id"]).strip():
            raise RuntimeError(f"Stage 23 aggregate source batch_id mismatch: {manifest_path}")
        lookup[manifest_path] = {
            "source_run_id": str(record["run_id"]).strip(),
            "source_batch_id": str(record["batch_id"]).strip(),
            "source_route_manifest_sha256": expected_sha256,
        }
    return lookup


def _validate_aggregate_source_membership(
    row: Mapping[str, Any],
    aggregate_sources: Mapping[Path, Mapping[str, str]],
    label: str,
) -> None:
    source_manifest = _resolve_mixed_path(str(row.get("source_route_manifest", ""))).resolve()
    expected = aggregate_sources.get(source_manifest)
    if expected is None:
        raise RuntimeError(f"{label} source route manifest is not listed by the Stage 2.5 aggregate route")
    for field in ["source_run_id", "source_batch_id", "source_route_manifest_sha256"]:
        observed = str(row.get(field, "")).strip()
        if observed != str(expected[field]):
            raise RuntimeError(
                f"{label} aggregate source provenance mismatch for {field}: "
                f"observed={observed!r}, expected={expected[field]!r}"
            )


def _expected_stage3_outputs(
    *,
    job_row: Mapping[str, str],
    stage3_root: Path,
    expected_stage3_mode: str,
    backbone_id: str,
) -> dict[Path, tuple[Path, str]]:
    tags = _split_csv(str(job_row.get("input_tags", "")))
    input_pdbs = [_resolve_mixed_path(value) for value in _split_csv(str(job_row.get("input_pdbs", "")))]
    input_hashes = _split_csv(str(job_row.get("input_pdb_sha256s", "")))
    if not tags or len(tags) != len(set(tags)):
        raise RuntimeError(f"Stage 3 job has blank or duplicate input tags for {backbone_id}")
    if len(input_pdbs) != len(tags) or len(input_hashes) != len(tags):
        raise RuntimeError(f"Stage 3 job input tag/PDB/hash counts disagree for {backbone_id}")

    input_dir = assert_active_route_path(
        _resolve_mixed_path(str(job_row.get("input_pdb_dir", ""))),
        f"Stage 23 input PDB directory for {backbone_id}",
    )
    _require_within(input_dir, stage3_root, f"Stage 23 input PDB directory for {backbone_id}")
    for tag, input_pdb, expected_sha256 in zip(tags, input_pdbs, input_hashes):
        input_pdb = assert_active_route_path(input_pdb, f"Stage 23 input PDB {tag}")
        _require_within(input_pdb, input_dir, f"Stage 23 input PDB {tag}")
        if input_pdb.stem != tag:
            raise RuntimeError(f"Stage 3 input tag/PDB stem mismatch for {backbone_id}: {tag} != {input_pdb.stem}")
        if _sha256_file(input_pdb) != expected_sha256:
            raise RuntimeError(f"Stage 3 input PDB SHA-256 mismatch for {backbone_id}: {input_pdb}")

    if str(job_row.get("input_tag", "")).strip() != tags[0]:
        raise RuntimeError(f"Stage 3 job input_tag does not match first input_tags entry for {backbone_id}")
    if _resolve_mixed_path(str(job_row.get("input_pdb", ""))).resolve() != input_pdbs[0].resolve():
        raise RuntimeError(f"Stage 3 job input_pdb does not match first input_pdbs entry for {backbone_id}")
    if str(job_row.get("input_pdb_sha256", "")).strip() != input_hashes[0]:
        raise RuntimeError(f"Stage 3 job input_pdb_sha256 does not match first hash for {backbone_id}")

    output_dir = _resolve_mixed_path(str(job_row.get("output_pdb_dir", "")))
    if not str(job_row.get("output_pdb_dir", "")).strip():
        raise RuntimeError(f"Stage 3 job is missing output_pdb_dir for {backbone_id}")
    assert_active_route_path(output_dir, f"Stage 23 output PDB directory for {backbone_id}", must_exist=False)
    _require_within(output_dir, stage3_root, f"Stage 23 output PDB directory for {backbone_id}")

    seqs_per_backbone = _parse_int(job_row.get("seqs_per_backbone", "0"))
    relax_cycles = _parse_int(job_row.get("relax_cycles", "-1"), -1)
    if seqs_per_backbone <= 0 or relax_cycles < 0:
        raise RuntimeError(f"Stage 3 job has invalid sequence/relax counts for {backbone_id}")

    expected: dict[Path, tuple[Path, str]] = {}
    if expected_stage3_mode == "proteinmpnn_only":
        if relax_cycles != 0 or len(tags) != 1:
            raise RuntimeError(f"ProteinMPNN-only job shape is invalid for {backbone_id}")
        for design_index in range(seqs_per_backbone):
            path = output_dir / f"{tags[0]}_dldesign_{design_index}.pdb"
            expected[path.resolve()] = (input_pdbs[0], tags[0])
    else:
        if relax_cycles <= 0 or len(tags) != seqs_per_backbone:
            raise RuntimeError(f"ProteinMPNN-FastRelax job shape is invalid for {backbone_id}")
        for tag, input_pdb in zip(tags, input_pdbs):
            path = output_dir / f"{tag}_dldesign_0_cycle{relax_cycles}.pdb"
            expected[path.resolve()] = (input_pdb, tag)
    return expected


def _validate_stage3_job_row(
    *,
    backbone_id: str,
    backbone_row: Mapping[str, str],
    job_row: Mapping[str, str],
    expected_stage3_mode: str,
    stage3_root: Path,
    stage2_selection_csv: Path,
    stage2_selection_csv_sha256: str,
    aggregate_sources: Mapping[Path, Mapping[str, str]],
) -> dict[Path, tuple[Path, str]]:
    for field in ["global_backbone_id", "design_id", "backbone_id"]:
        if str(job_row.get(field, "")).strip() != backbone_id:
            raise RuntimeError(f"Stage 3 job {field} does not match global backbone ID: {backbone_id}")
    if str(backbone_row.get("global_backbone_id", "")).strip() != backbone_id:
        raise RuntimeError(f"Stage 2.5 row global_backbone_id mismatch: {backbone_id}")
    expected_global_id = (
        f"{str(backbone_row.get('source_run_id', '')).strip()}__"
        f"{str(backbone_row.get('source_local_design_id', '')).strip()}"
    )
    if backbone_id != expected_global_id:
        raise RuntimeError(f"Stage 2.5 global backbone identity formula mismatch: {backbone_id}")
    if str(job_row.get("stage3_mode", "")).strip() != expected_stage3_mode:
        raise RuntimeError(
            f"Stage 3 job mode mismatch for {backbone_id}: expected={expected_stage3_mode}, "
            f"actual={job_row.get('stage3_mode', '')}"
        )
    for field in ["source_local_design_id", "source_batch_label", "backbone_family_id"]:
        if str(job_row.get(field, "")).strip() != str(backbone_row.get(field, "")).strip():
            raise RuntimeError(f"Stage 3 job and Stage 2.5 row disagree for {field}: {backbone_id}")

    _validate_aggregate_source_membership(job_row, aggregate_sources, f"Stage 23 Stage 3 job {backbone_id}")
    _validate_aggregate_source_membership(backbone_row, aggregate_sources, f"Stage 23 Stage 2.5 row {backbone_id}")
    for field in SOURCE_ROUTE_PROVENANCE_FIELDS:
        if str(job_row.get(field, "")).strip() != str(backbone_row.get(field, "")).strip():
            raise RuntimeError(f"Stage 3 job and Stage 2.5 row source provenance disagree for {field}: {backbone_id}")

    recorded_selection = assert_active_route_path(
        _resolve_mixed_path(str(job_row.get("stage2_selection_csv", ""))),
        f"Stage 23 recorded Stage 2.5 selection CSV for {backbone_id}",
    )
    if recorded_selection.resolve() != stage2_selection_csv.resolve():
        raise RuntimeError(f"Stage 3 job references a different Stage 2.5 selection CSV: {backbone_id}")
    if str(job_row.get("stage2_selection_csv_sha256", "")).strip() != stage2_selection_csv_sha256:
        raise RuntimeError(f"Stage 3 job Stage 2.5 selection CSV SHA-256 mismatch: {backbone_id}")

    source_text = str(job_row.get("source_backbone_pdb", "")).strip()
    source_sha256 = str(job_row.get("source_backbone_pdb_sha256", "")).strip().lower()
    job_source = assert_active_route_path(
        _resolve_mixed_path(source_text),
        f"Stage 23 job source PDB for {backbone_id}",
    )
    stage2_source = assert_active_route_path(
        _resolve_mixed_path(str(backbone_row.get("rf_pdb", ""))),
        f"Stage 23 Stage 2.5 source PDB for {backbone_id}",
    )
    if job_source.resolve() != stage2_source.resolve():
        raise RuntimeError(
            f"Stage 3 source backbone path mismatch for {backbone_id}: job={job_source}, stage2={stage2_source}"
        )
    actual_sha256 = _sha256_file(job_source)
    if actual_sha256 != source_sha256 or actual_sha256 != str(backbone_row.get("pdb_sha256", "")).strip():
        raise RuntimeError(f"Stage 3 source backbone SHA-256 mismatch for {backbone_id}")

    for field in ["stage3_job_id", "run_group_id", "protocol_identity_sha256"]:
        if not str(job_row.get(field, "")).strip():
            raise RuntimeError(f"Stage 3 job is missing {field} for {backbone_id}")
    expected_job_id = _safe_token(
        f"{backbone_id}_{expected_stage3_mode}_"
        f"{str(job_row.get('protocol_identity_sha256', '')).strip()[:12]}"
    )
    if str(job_row.get("stage3_job_id", "")).strip() != expected_job_id:
        raise RuntimeError(f"Stage 3 job identity formula mismatch for {backbone_id}")
    return _expected_stage3_outputs(
        job_row=job_row,
        stage3_root=stage3_root,
        expected_stage3_mode=expected_stage3_mode,
        backbone_id=backbone_id,
    )


def _stage3_job_provenance(
    job_row: Mapping[str, str],
    *,
    stage3_jobs_csv: Path,
    stage3_jobs_csv_sha256: str,
    expected_count: int,
    observed_count: int,
) -> dict[str, str | int]:
    return {
        "stage3_job_id": str(job_row.get("stage3_job_id", "")),
        "run_group_id": str(job_row.get("run_group_id", "")),
        "protocol_identity_sha256": str(job_row.get("protocol_identity_sha256", "")),
        "stage3_mode": str(job_row.get("stage3_mode", "")),
        "global_backbone_id": str(job_row.get("global_backbone_id", "")),
        "source_local_design_id": str(job_row.get("source_local_design_id", "")),
        "source_batch_label": str(job_row.get("source_batch_label", "")),
        "backbone_family_id": str(job_row.get("backbone_family_id", "")),
        "source_backbone_pdb": str(job_row.get("source_backbone_pdb", "")),
        "source_backbone_pdb_sha256": str(job_row.get("source_backbone_pdb_sha256", "")),
        "stage3_jobs_csv": str(stage3_jobs_csv),
        "stage3_jobs_csv_sha256": stage3_jobs_csv_sha256,
        "stage2_selection_csv": str(job_row.get("stage2_selection_csv", "")),
        "stage2_selection_csv_sha256": str(job_row.get("stage2_selection_csv_sha256", "")),
        "expected_stage3b_outputs_for_backbone": expected_count,
        "observed_stage3b_outputs_for_backbone": observed_count,
        **{field: str(job_row.get(field, "")) for field in SOURCE_ROUTE_PROVENANCE_FIELDS},
    }


def _residue_label(chain_id: str, residue: Mapping[str, Any]) -> str:
    return f"{chain_id}{residue.get('pdb_residue_number', '')}"


def _residue_atom_coords(residue: Mapping[str, Any]) -> list[tuple[float, float, float]]:
    atoms = residue.get("atoms", {})
    return [coord for coord in atoms.values() if isinstance(coord, tuple) and len(coord) == 3]


def _all_atom_coords(residues: Iterable[Mapping[str, Any]]) -> list[tuple[float, float, float]]:
    coords: list[tuple[float, float, float]] = []
    for residue in residues:
        coords.extend(_residue_atom_coords(residue))
    return coords


def _sq_distance(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return sum((a[idx] - b[idx]) ** 2 for idx in range(3))


def _residue_contact_metrics(
    target_chain_id: str,
    target_residues: Sequence[Mapping[str, Any]],
    peptide_coords: Sequence[tuple[float, float, float]],
    contact_cutoff: float,
) -> tuple[int, float | None, str]:
    if not target_residues or not peptide_coords:
        return 0, None, ""

    cutoff_sq = contact_cutoff * contact_cutoff
    contact_count = 0
    nearest_sq: float | None = None
    nearest_label = ""
    for residue in target_residues:
        residue_coords = _residue_atom_coords(residue)
        if not residue_coords:
            continue
        residue_nearest_sq: float | None = None
        for left in residue_coords:
            for right in peptide_coords:
                current = _sq_distance(left, right)
                if residue_nearest_sq is None or current < residue_nearest_sq:
                    residue_nearest_sq = current
        if residue_nearest_sq is None:
            continue
        if nearest_sq is None or residue_nearest_sq < nearest_sq:
            nearest_sq = residue_nearest_sq
            nearest_label = _residue_label(target_chain_id, residue)
        if residue_nearest_sq <= cutoff_sq:
            contact_count += 1

    return contact_count, math.sqrt(nearest_sq) if nearest_sq is not None else None, nearest_label


def _status_from_distance(
    *,
    contact_count: int,
    min_distance: float | None,
    contact_name: str,
    near_name: str,
    miss_name: str,
    near_distance: float,
    min_contacts: int,
) -> str:
    if min_distance is None:
        return miss_name
    if contact_count >= min_contacts:
        return contact_name
    if min_distance <= near_distance:
        return near_name
    return miss_name


def _target_contact_status(
    *,
    target_contacts: int,
    target_min_distance: float | None,
    contact_cutoff: float,
    near_distance: float,
    min_target_contacts: int,
) -> str:
    if target_contacts >= min_target_contacts:
        return "target_contact_pass"
    if target_contacts > 0:
        return "target_contact_low_count"
    if target_min_distance is not None and target_min_distance <= contact_cutoff:
        return "target_contact_low_count"
    if target_min_distance is not None and target_min_distance <= near_distance:
        return "target_near_only"
    return "detached_from_target_crop"


def _macrocycle_geometry_status(
    peptide_residues: Sequence[Mapping[str, Any]],
    pass_distance: float,
    warn_distance: float,
) -> tuple[str, float | None]:
    if not peptide_residues:
        return "fail_parse_error", None
    first_atoms = peptide_residues[0].get("atoms", {})
    last_atoms = peptide_residues[-1].get("atoms", {})
    n_atom = first_atoms.get("N")
    c_atom = last_atoms.get("C")
    if n_atom is None or c_atom is None:
        return "fail_missing_terminal_atoms", None
    cn_distance = distance(n_atom, c_atom)
    if cn_distance <= pass_distance:
        return "pass_head_to_tail_macrocycle", cn_distance
    if cn_distance <= warn_distance:
        return "warn_cyclic_metadata_missing_but_geometry_close", cn_distance
    return "fail_open_chain_or_no_cyclic_evidence", cn_distance


def _peptide_geometry_summary(peptide_residues: Sequence[Mapping[str, Any]]) -> tuple[tuple[float, float, float], float]:
    ca_coords = [coord for coord in (ca_coord(residue) for residue in peptide_residues) if coord is not None]
    if not ca_coords:
        return (0.0, 0.0, 0.0), 0.0
    center = centroid(ca_coords)
    rog = math.sqrt(sum(distance(coord, center) ** 2 for coord in ca_coords) / len(ca_coords))
    return center, rog


def _load_site_mapping(path: Path) -> tuple[set[str], set[str]]:
    rows = _read_required_csv(path)
    site_numbers: set[str] = set()
    hotspot_numbers: set[str] = set()
    for row in rows:
        residue_number = str(row.get("rfpeptides_residue_number", "")).strip()
        if not residue_number:
            continue
        if str(row.get("is_target_site_residue", "")).strip().lower() == "true":
            site_numbers.add(residue_number)
        if str(row.get("is_selected_hotspot", "")).strip().lower() == "true":
            hotspot_numbers.add(residue_number)
    if not site_numbers:
        raise RuntimeError(f"No target-site residues marked in mapping CSV: {path}")
    if not hotspot_numbers:
        raise RuntimeError(f"No selected hotspots marked in mapping CSV: {path}")
    return site_numbers, hotspot_numbers


def _parse_numeric_residue_number(value: str) -> int | None:
    digits = "".join(ch for ch in str(value) if ch.isdigit() or ch == "-")
    if not digits:
        return None
    try:
        return int(digits)
    except ValueError:
        return None


def _target_number_mapping(
    *,
    source_backbone_pdb: Path,
    target_chain: str,
    relaxed_target_residues: Sequence[Mapping[str, Any]],
) -> tuple[dict[str, str], str, int | str]:
    try:
        source_chains = parse_residues(source_backbone_pdb)
    except Exception:
        return {}, "source_backbone_parse_failed", ""

    source_residues = list(source_chains.get(target_chain, []))
    if not source_residues or not relaxed_target_residues:
        return {}, "target_chain_missing_for_numbering_map", ""
    if len(source_residues) != len(relaxed_target_residues):
        return {}, "target_residue_count_changed_no_numbering_map", ""

    # Stage 0 mapping values are crop sequence positions (1..N), while
    # RFpeptides/Rosetta PDB residue numbers can be offset by peptide length.
    # Map crop positions to the current PDB numbers through target-chain order.
    mapping = {
        str(position): str(relaxed.get("pdb_residue_number", ""))
        for position, relaxed in enumerate(relaxed_target_residues, start=1)
    }
    offsets: set[int] = set()
    for source_number, relaxed_number in mapping.items():
        source_int = _parse_numeric_residue_number(source_number)
        relaxed_int = _parse_numeric_residue_number(relaxed_number)
        if source_int is None or relaxed_int is None:
            continue
        offsets.add(relaxed_int - source_int)

    if not offsets:
        return mapping, "mapped_crop_positions_by_target_order_non_numeric", ""
    if len(offsets) == 1:
        offset = next(iter(offsets))
        if offset == 0:
            return mapping, "mapped_crop_positions_by_target_order_same_pdb_numbers", 0
        return mapping, "mapped_crop_positions_by_target_order_constant_pdb_offset", offset
    return mapping, "mapped_crop_positions_by_target_order_variable_pdb_offset", ""


def _apply_number_mapping(numbers: set[str], mapping: Mapping[str, str]) -> set[str]:
    return {str(mapping.get(number, number)) for number in numbers}


def _forbidden_present(sequence: str, forbidden_aas: str) -> str:
    forbidden = {aa.upper() for aa in forbidden_aas if aa.strip()}
    present = sorted({aa for aa in sequence.upper() if aa in forbidden})
    return ",".join(present)


def _parse_rosetta_energy_table(path: Path, peptide_length: int) -> tuple[str, float | str, float | str, float | str]:
    in_table = False
    total_idx: int | None = None
    pose_total: float | str = ""
    residue_totals: list[float] = []
    try:
        lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    except OSError:
        return "not_readable", "", "", ""

    for line in lines:
        if line.startswith("#BEGIN_POSE_ENERGIES_TABLE"):
            in_table = True
            continue
        if line.startswith("#END_POSE_ENERGIES_TABLE"):
            break
        if not in_table:
            continue
        parts = line.split()
        if not parts:
            continue
        if parts[0] == "label":
            try:
                total_idx = parts.index("total")
            except ValueError:
                total_idx = None
            continue
        if total_idx is None or len(parts) <= total_idx or parts[0] == "weights":
            continue
        try:
            total = float(parts[total_idx])
        except ValueError:
            continue
        if parts[0] == "pose":
            pose_total = total
        else:
            residue_totals.append(total)

    if pose_total == "" and not residue_totals:
        return "not_found", "", "", ""
    peptide_total = sum(residue_totals[:peptide_length]) if residue_totals else ""
    target_total = sum(residue_totals[peptide_length:]) if residue_totals and len(residue_totals) > peptide_length else ""
    return (
        "parsed",
        round(pose_total, 3) if isinstance(pose_total, float) else "",
        round(peptide_total, 3) if isinstance(peptide_total, float) else "",
        round(target_total, 3) if isinstance(target_total, float) else "",
    )


def _validate_complete_stage3_outputs(
    expected_by_backbone: Mapping[str, Mapping[Path, tuple[Path, str]]],
) -> int:
    expected_paths: set[Path] = set()
    output_dirs: set[Path] = set()
    for backbone_id, expected in expected_by_backbone.items():
        if not expected:
            raise RuntimeError(f"Stage 3 job has no expected output paths: {backbone_id}")
        for path in expected:
            resolved = path.resolve()
            if resolved in expected_paths:
                raise RuntimeError(f"Two Stage 3 jobs expect the same output PDB: {resolved}")
            expected_paths.add(resolved)
            output_dirs.add(resolved.parent)

    actual_paths: set[Path] = set()
    for output_dir in output_dirs:
        output_dir = assert_active_route_path(
            output_dir,
            "Stage 23 Stage 3B output directory",
            must_exist=False,
        )
        if output_dir.exists():
            for path in output_dir.glob("*.pdb"):
                actual_paths.add(path.resolve())

    missing = sorted(expected_paths - actual_paths, key=str)
    unexpected = sorted(actual_paths - expected_paths, key=str)
    empty = sorted(
        [path for path in expected_paths & actual_paths if path.stat().st_size == 0],
        key=str,
    )
    if missing or unexpected or empty:
        details = []
        if missing:
            details.append(f"missing={len(missing)} first={missing[:5]}")
        if unexpected:
            details.append(f"unexpected={len(unexpected)} first={unexpected[:5]}")
        if empty:
            details.append(f"empty={len(empty)} first={empty[:5]}")
        raise RuntimeError(
            "Stage 3B output set is incomplete or contaminated; Stage 3C was not started. "
            + "; ".join(details)
        )
    return len(actual_paths)


def _qc_row_for_relaxed_pdb(
    *,
    relaxed_pdb: Path,
    backbone_row: Mapping[str, str],
    site_numbers: set[str],
    hotspot_numbers: set[str],
    contact_cutoff: float,
    site_near_distance: float,
    hotspot_near_distance: float,
    severe_clash_distance: float,
    min_target_contacts: int,
    min_site_contacts: int,
    min_hotspot_contacts: int,
    macrocycle_pass_distance: float,
    macrocycle_warn_distance: float,
    forbidden_aas: str,
) -> dict[str, Any]:
    backbone_id = str(
        backbone_row.get("global_backbone_id", "")
        or backbone_row.get("design_id", "")
    ).strip()
    sequence_design_id = relaxed_pdb.stem
    target_chain = str(backbone_row.get("target_chain", "")).strip() or "A"
    peptide_chain = str(backbone_row.get("peptide_chain", "")).strip() or "B"
    expected_length = _parse_int(backbone_row.get("peptide_length", "0"))
    source_backbone_pdb = _resolve_mixed_path(str(backbone_row.get("rf_pdb", "")))
    base_row: dict[str, Any] = {
        "sequence_design_id": sequence_design_id,
        "backbone_id": backbone_id,
        "site_label": backbone_row.get("site_label", ""),
        "site_id": backbone_row.get("site_id", ""),
        "relaxed_pdb": relaxed_pdb,
        "source_backbone_pdb": source_backbone_pdb,
        "target_chain": target_chain,
        "peptide_chain": peptide_chain,
        "expected_peptide_length": expected_length,
        "forbidden_aas": forbidden_aas,
        "pass_stage3c_qc": "false",
    }

    if not relaxed_pdb.exists() or relaxed_pdb.stat().st_size == 0:
        base_row.update(
            {
                "file_status": "fail_missing_pdb",
                "parse_status": "not_parsed",
                "sequence_status": "not_evaluated",
                "qc_failure_reasons": "missing_relaxed_pdb",
                "qc_notes": "Stage 3B output PDB is missing or empty.",
            }
        )
        return base_row

    base_row["file_status"] = "pass"
    base_row["relaxed_pdb_sha256"] = _sha256_file(relaxed_pdb)
    try:
        chains = parse_residues(relaxed_pdb)
    except Exception as exc:  # pragma: no cover - malformed external PDBs
        base_row.update(
            {
                "parse_status": "fail_parse_error",
                "sequence_status": "not_evaluated",
                "qc_failure_reasons": f"parse_error:{exc.__class__.__name__}",
                "qc_notes": str(exc),
            }
        )
        return base_row

    target_residues = list(chains.get(target_chain, []))
    peptide_residues = list(chains.get(peptide_chain, []))
    if not target_residues or not peptide_residues:
        missing = []
        if not target_residues:
            missing.append(f"target_chain_{target_chain}")
        if not peptide_residues:
            missing.append(f"peptide_chain_{peptide_chain}")
        base_row.update(
            {
                "parse_status": "fail_missing_chain",
                "target_residue_count": len(target_residues),
                "peptide_length": len(peptide_residues),
                "sequence_status": "not_evaluated",
                "qc_failure_reasons": ";".join(missing),
                "qc_notes": "Required chain is missing from relaxed PDB.",
            }
        )
        return base_row

    peptide_sequence = residue_sequence(peptide_residues)
    forbidden_present = _forbidden_present(peptide_sequence, forbidden_aas)
    sequence_failures = []
    if expected_length and len(peptide_residues) != expected_length:
        sequence_failures.append("peptide_length_changed")
    if "X" in peptide_sequence:
        sequence_failures.append("unknown_residue_in_sequence")
    if forbidden_present:
        sequence_failures.append("forbidden_aas_present")
    sequence_status = "pass_sequence" if not sequence_failures else "fail_" + "_and_".join(sequence_failures)

    target_number_map, numbering_status, numbering_offset = _target_number_mapping(
        source_backbone_pdb=source_backbone_pdb,
        target_chain=target_chain,
        relaxed_target_residues=target_residues,
    )
    mapped_site_numbers = _apply_number_mapping(site_numbers, target_number_map)
    mapped_hotspot_numbers = _apply_number_mapping(hotspot_numbers, target_number_map)

    target_by_number = {str(residue.get("pdb_residue_number", "")): residue for residue in target_residues}
    site_residues = [
        target_by_number[number] for number in sorted(mapped_site_numbers, key=_residue_sort_key) if number in target_by_number
    ]
    hotspot_residues = [
        target_by_number[number] for number in sorted(mapped_hotspot_numbers, key=_residue_sort_key) if number in target_by_number
    ]

    peptide_coords = _all_atom_coords(peptide_residues)
    target_contacts, target_min, closest_target = _residue_contact_metrics(target_chain, target_residues, peptide_coords, contact_cutoff)
    site_contacts, site_min, closest_site = _residue_contact_metrics(target_chain, site_residues, peptide_coords, contact_cutoff)
    hotspot_contacts, hotspot_min, closest_hotspot = _residue_contact_metrics(
        target_chain, hotspot_residues, peptide_coords, contact_cutoff
    )

    target_contact_status = _target_contact_status(
        target_contacts=target_contacts,
        target_min_distance=target_min,
        contact_cutoff=contact_cutoff,
        near_distance=max(site_near_distance, hotspot_near_distance),
        min_target_contacts=min_target_contacts,
    )
    site_status = _status_from_distance(
        contact_count=site_contacts,
        min_distance=site_min,
        contact_name="site_contact_pass",
        near_name="site_near_only",
        miss_name="site_missed",
        near_distance=site_near_distance,
        min_contacts=min_site_contacts,
    )
    if target_contact_status in {"target_contact_pass", "target_contact_low_count"} and site_status == "site_missed":
        site_status = "crop_only_contact"
    hotspot_status = _status_from_distance(
        contact_count=hotspot_contacts,
        min_distance=hotspot_min,
        contact_name="hotspot_contact_pass",
        near_name="hotspot_near_only",
        miss_name="hotspot_missed",
        near_distance=hotspot_near_distance,
        min_contacts=min_hotspot_contacts,
    )

    macrocycle_status, cn_distance = _macrocycle_geometry_status(
        peptide_residues,
        pass_distance=macrocycle_pass_distance,
        warn_distance=macrocycle_warn_distance,
    )
    clash_status = "not_evaluated"
    if target_min is not None:
        clash_status = "fail_severe_clash" if target_min < severe_clash_distance else "pass_no_severe_clash"
    center, rog = _peptide_geometry_summary(peptide_residues)
    score_status, pose_total, peptide_total, target_total = _parse_rosetta_energy_table(relaxed_pdb, len(peptide_residues))

    failure_reasons = []
    failure_reasons.extend(sequence_failures)
    if macrocycle_status != "pass_head_to_tail_macrocycle":
        failure_reasons.append(macrocycle_status)
    if target_contact_status not in {"target_contact_pass", "target_contact_low_count"}:
        failure_reasons.append(target_contact_status)
    if site_status != "site_contact_pass":
        failure_reasons.append(site_status)
    if hotspot_status != "hotspot_contact_pass":
        failure_reasons.append(hotspot_status)
    if clash_status.startswith("fail"):
        failure_reasons.append(clash_status)

    pass_qc = not failure_reasons
    notes = []
    if site_status == "crop_only_contact":
        notes.append("Relaxed peptide contacts target crop but misses RFpep_Site_2 residues.")
    if hotspot_status != "hotspot_contact_pass":
        notes.append("Relaxed peptide does not preserve direct hotspot contact.")
    if macrocycle_status != "pass_head_to_tail_macrocycle":
        notes.append("Relaxed peptide does not preserve head-to-tail macrocycle geometry.")
    if score_status != "parsed":
        notes.append("Rosetta pose energy table was not parsed.")

    base_row.update(
        {
            "parse_status": "pass",
            "peptide_sequence": peptide_sequence,
            "peptide_length": len(peptide_residues),
            "forbidden_aas_present": forbidden_present,
            "sequence_status": sequence_status,
            "target_residue_count": len(target_residues),
            "target_residue_numbering_status": numbering_status,
            "target_residue_numbering_offset": numbering_offset,
            "mapped_target_site_pdb_residue_numbers": ",".join(
                sorted(mapped_site_numbers, key=_residue_sort_key)
            ),
            "mapped_hotspot_pdb_residue_numbers": ",".join(
                sorted(mapped_hotspot_numbers, key=_residue_sort_key)
            ),
            "target_site_residue_count": len(site_residues),
            "hotspot_residue_count": len(hotspot_residues),
            "num_target_contacts": target_contacts,
            "num_target_site_contacts": site_contacts,
            "num_hotspot_contacts": hotspot_contacts,
            "peptide_target_min_distance": _round_float(target_min),
            "peptide_site_min_distance": _round_float(site_min),
            "peptide_hotspot_min_distance": _round_float(hotspot_min),
            "closest_target_residue": closest_target,
            "closest_site_residue": closest_site,
            "closest_hotspot_residue": closest_hotspot,
            "target_contact_status": target_contact_status,
            "target_site_recovery_status": site_status,
            "hotspot_recovery_status": hotspot_status,
            "macrocycle_terminal_cn_distance": _round_float(cn_distance),
            "macrocycle_geometry_status": macrocycle_status,
            "clash_status": clash_status,
            "peptide_ca_centroid_x": round(center[0], 3),
            "peptide_ca_centroid_y": round(center[1], 3),
            "peptide_ca_centroid_z": round(center[2], 3),
            "peptide_radius_of_gyration": round(rog, 3),
            "rosetta_score_status": score_status,
            "rosetta_pose_total": pose_total,
            "rosetta_peptide_total_sum": peptide_total,
            "rosetta_target_total_sum": target_total,
            "pass_stage3c_qc": "true" if pass_qc else "false",
            "qc_failure_reasons": ";".join(failure_reasons),
            "qc_notes": "; ".join(notes),
        }
    )
    return base_row


def _write_fastas(rows: Iterable[Mapping[str, Any]], fasta_dir: Path) -> None:
    fasta_dir.mkdir(parents=True, exist_ok=True)
    for row in rows:
        sequence = str(row.get("peptide_sequence", "")).strip()
        if not sequence:
            continue
        header = (
            f"{row.get('sequence_design_id')} "
            f"backbone={row.get('backbone_id')} "
            f"site={row.get('site_label')} "
            f"pass_stage3c_qc={row.get('pass_stage3c_qc')}"
        )
        write_fasta(fasta_dir / f"{_safe_token(str(row.get('sequence_design_id', 'sequence')))}.fa", header, sequence)


def _status_lines(rows: list[Mapping[str, Any]], field: str) -> str:
    counts = Counter(str(row.get(field, "")) for row in rows)
    if not counts:
        return "- none: 0"
    return "\n".join(f"- {key or 'blank'}: {counts[key]}" for key in sorted(counts))


def _stage0_site_inputs(
    *,
    route_manifest: Mapping[str, Any],
    stage0_root: Path,
) -> tuple[str, Path, Path]:
    records = route_manifest.get("stage0_sites")
    if not isinstance(records, list) or len(records) != 1 or not isinstance(records[0], Mapping):
        raise RuntimeError("Stage 23 requires exactly one Stage 0 site record in the aggregate route manifest")
    record = records[0]
    site_label = str(record.get("site_label", "")).strip()
    target_pdb = assert_active_route_path(
        _resolve_mixed_path(str(record.get("target_pdb", ""))),
        "Stage 23 manifest-locked Stage 0 target PDB",
    )
    mapping_csv = assert_active_route_path(
        _resolve_mixed_path(str(record.get("mapping_csv", ""))),
        "Stage 23 manifest-locked Stage 0 mapping CSV",
    )
    _require_within(target_pdb, stage0_root, "Stage 23 Stage 0 target PDB")
    _require_within(mapping_csv, stage0_root, "Stage 23 Stage 0 mapping CSV")
    if _sha256_file(target_pdb) != str(record.get("target_pdb_sha256", "")).strip():
        raise RuntimeError("Stage 23 Stage 0 target PDB SHA-256 does not match the aggregate route manifest")
    if _sha256_file(mapping_csv) != str(record.get("mapping_csv_sha256", "")).strip():
        raise RuntimeError("Stage 23 Stage 0 mapping CSV SHA-256 does not match the aggregate route manifest")
    return site_label, target_pdb, mapping_csv


def _validate_runlist_contract(
    *,
    job_rows: Sequence[Mapping[str, str]],
    stage3_root: Path,
) -> tuple[Path, int]:
    runlist = assert_active_route_path(
        _resolve_mixed_path(_single_value(job_rows, "runlist", "Stage 3 jobs")),
        "Stage 23 Stage 3 runlist",
    )
    _require_within(runlist, stage3_root, "Stage 23 Stage 3 runlist")
    expected_tags: list[str] = []
    for row in job_rows:
        expected_tags.extend(_split_csv(str(row.get("input_tags", ""))))
    if not expected_tags or len(expected_tags) != len(set(expected_tags)):
        raise RuntimeError("Stage 3 jobs contain blank or duplicate input tags")
    actual_tags = [line.strip() for line in runlist.read_text(encoding="utf-8").splitlines() if line.strip()]
    if actual_tags != expected_tags:
        raise RuntimeError(
            "Stage 3 runlist does not exactly match the ordered input_tags in the jobs CSV: "
            f"expected={len(expected_tags)}, actual={len(actual_tags)}"
        )
    return runlist, len(actual_tags)


def _summary_markdown(
    *,
    rows: list[Mapping[str, Any]],
    args: argparse.Namespace,
    output_dir: Path,
    site_numbers: set[str],
    hotspot_numbers: set[str],
    run_group_id: str,
    stage2_selection_csv: Path,
    stage3_jobs_csv: Path,
    job_count: int,
    selected_backbone_count: int,
    expected_output_count: int,
) -> str:
    pass_rows = [row for row in rows if row.get("pass_stage3c_qc") == "true"]
    top_rows = sorted(
        rows,
        key=lambda row: (
            row.get("pass_stage3c_qc") != "true",
            float(row.get("peptide_hotspot_min_distance") or 999.0),
            float(row.get("macrocycle_terminal_cn_distance") or 999.0),
            str(row.get("sequence_design_id", "")),
        ),
    )[: args.top_report]
    columns = [
        "sequence_design_id",
        "peptide_sequence",
        "num_target_site_contacts",
        "num_hotspot_contacts",
        "peptide_site_min_distance",
        "peptide_hotspot_min_distance",
        "macrocycle_terminal_cn_distance",
        "macrocycle_geometry_status",
        "pass_stage3c_qc",
        "qc_failure_reasons",
    ]
    mode_label = "ProteinMPNN-only" if args.stage3_mode == "proteinmpnn_only" else "ProteinMPNN-FastRelax"
    return f"""# FGA RFpeptides Stage 3C {mode_label} Sequence QC

Status: {mode_label} outputs parsed after Stage 3B.

Important rule: Stage 3C does not treat sequence generation alone as success.
A sequence-designed structure must still contact RFpep_Site_2 and the selected
hotspot residues, and it must preserve head-to-tail macrocycle geometry.

Output directory:

```text
{output_dir}
```

Parameters:

```text
stage3_mode: {args.stage3_mode}
stage3_root: {args.stage3_root}
run_group_id: {run_group_id}
stage2_selection_csv: {stage2_selection_csv}
stage3_jobs_csv: {stage3_jobs_csv}
jobs_in_authoritative_table: {job_count}
selected_global_backbones_collected: {selected_backbone_count}
expected_complete_stage3b_outputs: {expected_output_count}
contact_cutoff_A: {args.contact_cutoff}
site_near_distance_A: {args.site_near_distance}
hotspot_near_distance_A: {args.hotspot_near_distance}
macrocycle_pass_distance_A: {args.macrocycle_pass_distance}
macrocycle_warn_distance_A: {args.macrocycle_warn_distance}
forbidden_aas: {args.forbidden_aas}
```

Target-site residue numbers:

```text
{",".join(sorted(site_numbers, key=_residue_sort_key))}
```

These are Stage 0 crop sequence positions. They are mapped to each output PDB
through target-chain residue order before contact calculations and PyMOL
selections; raw RFpeptides/Rosetta residue numbers may include a peptide-length
offset.

Hotspot residue numbers:

```text
{",".join(sorted(hotspot_numbers, key=_residue_sort_key))}
```

## Counts

```text
total_stage3b_outputs: {len(rows)}
pass_stage3c_qc: {len(pass_rows)}
```

Before QC began, Stage 23 required the complete, exact output set derived from
the authoritative Stage 22 jobs table. Missing, empty, or unexpected PDB files
cause a hard failure and no partial Stage 3C table is written.

## Sequence Status

{_status_lines(rows, "sequence_status")}

## Target-Site Recovery Status

{_status_lines(rows, "target_site_recovery_status")}

## Hotspot Recovery Status

{_status_lines(rows, "hotspot_recovery_status")}

## Macrocycle Geometry Status

{_status_lines(rows, "macrocycle_geometry_status")}

## Ranked Stage 3C Rows

{rows_to_markdown(top_rows, columns, "No Stage 3B outputs were parsed.")}

These Stage 3C rows are sequence-designed intermediate structures, not final
peptide candidates.
"""


def _write_pymol_review(
    *,
    rows: list[Mapping[str, Any]],
    output_path: Path,
    target_chain: str,
    peptide_chain: str,
    site_numbers: set[str],
    hotspot_numbers: set[str],
    top_n: int,
    pymol_path_style: str,
) -> None:
    selected = sorted(
        rows,
        key=lambda row: (
            row.get("pass_stage3c_qc") != "true",
            float(row.get("peptide_hotspot_min_distance") or 999.0),
            float(row.get("macrocycle_terminal_cn_distance") or 999.0),
            str(row.get("sequence_design_id", "")),
        ),
    )[:top_n]
    lines = [
        "reinitialize",
        "set retain_order, 1",
        "hide everything, all",
    ]
    for idx, row in enumerate(selected, start=1):
        obj = _safe_token(str(row.get("sequence_design_id", f"seq_{idx}")))
        pdb_path = _format_pymol_path(str(row.get("relaxed_pdb", "")), pymol_path_style)
        peptide_color = "cyan" if row.get("pass_stage3c_qc") == "true" else "magenta"
        row_site_numbers = set(_split_csv(str(row.get("mapped_target_site_pdb_residue_numbers", "")))) or site_numbers
        row_hotspot_numbers = set(_split_csv(str(row.get("mapped_hotspot_pdb_residue_numbers", "")))) or hotspot_numbers
        site_resi = "+".join(sorted(row_site_numbers, key=_residue_sort_key)) or "none"
        hotspot_resi = "+".join(sorted(row_hotspot_numbers, key=_residue_sort_key)) or "none"
        lines.extend(
            [
                f"load \"{pdb_path}\", {obj}",
                f"hide everything, {obj}",
                f"show cartoon, {obj} and chain {target_chain}",
                f"show sticks, {obj} and chain {peptide_chain}",
                f"select {obj}_site, {obj} and chain {target_chain} and resi {site_resi}",
                f"select {obj}_hotspots, {obj} and chain {target_chain} and resi {hotspot_resi}",
                f"show sticks, {obj}_site",
                f"show spheres, {obj}_hotspots",
                f"color gray80, {obj} and chain {target_chain}",
                f"color orange, {obj}_site",
                f"color red, {obj}_hotspots",
                f"color {peptide_color}, {obj} and chain {peptide_chain}",
                f"set sphere_scale, 0.45, {obj}_hotspots",
                f"label {obj}_hotspots and name CA, \"%s%s\" % (chain, resi)",
            ]
        )
    if selected:
        first = _safe_token(str(selected[0].get("sequence_design_id", "seq_1")))
        lines.append(f"zoom {first}_site, 14")
    write_markdown(output_path, "\n".join(lines))


def main() -> int:
    parser = argparse.ArgumentParser(description="Collect and QC Stage 3B ProteinMPNN sequence-design outputs.")
    parser.add_argument("--stage0-root", required=True)
    parser.add_argument("--stage3-root", required=True)
    parser.add_argument("--project-config", required=True)
    parser.add_argument("--output-root", default="", help="Defaults to --stage3-root.")
    parser.add_argument(
        "--stage3-mode",
        choices=["proteinmpnn_fastrelax", "proteinmpnn_only"],
        required=True,
        help="Must exactly match the mode locked in the Stage 22 jobs table.",
    )
    parser.add_argument("--stage2-selection-csv", required=True)
    parser.add_argument("--stage3-jobs-csv", required=True)
    parser.add_argument(
        "--selected-global-backbones",
        default="",
        help=(
            "Optional comma-separated global_backbone_id subset to report. "
            "The complete jobs-table output set is still required before any Stage 3C QC starts."
        ),
    )
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
    parser.add_argument(
        "--pymol-path-style",
        choices=["windows", "wsl", "native"],
        default="windows",
        help="Path style for generated PyMOL review scripts. Use windows for Windows PyMOL.",
    )
    args = parser.parse_args()

    logger = setup_logger("23_collect_proteinmpnn_sequences")
    append_run_header(logger, "23_collect_proteinmpnn_sequences.py")

    if args.contact_cutoff <= 0:
        raise RuntimeError("--contact-cutoff must be > 0")
    if args.site_near_distance < args.contact_cutoff:
        raise RuntimeError("--site-near-distance must be >= --contact-cutoff")
    if args.hotspot_near_distance < args.contact_cutoff:
        raise RuntimeError("--hotspot-near-distance must be >= --contact-cutoff")
    if args.macrocycle_warn_distance < args.macrocycle_pass_distance:
        raise RuntimeError("--macrocycle-warn-distance must be >= --macrocycle-pass-distance")

    stage0_root = assert_active_route_path(
        _resolve_mixed_path(args.stage0_root),
        "Stage 23 Stage 0 root",
    )
    stage3_root = assert_active_route_path(
        _resolve_mixed_path(args.stage3_root),
        "Stage 23 Stage 3 root",
    )
    output_root = assert_active_route_path(
        _resolve_mixed_path(args.output_root) if args.output_root else stage3_root,
        "Stage 23 output root",
        must_exist=False,
    )
    stage2_selection_csv = assert_active_route_path(
        _resolve_mixed_path(args.stage2_selection_csv),
        "Stage 23 Stage 2.5 selection CSV",
    )
    stage3_jobs_csv = assert_active_route_path(
        _resolve_mixed_path(args.stage3_jobs_csv),
        "Stage 23 Stage 3 jobs CSV",
    )
    _require_within(stage2_selection_csv, stage3_root, "Stage 23 Stage 2.5 selection CSV")
    _require_within(stage3_jobs_csv, stage3_root, "Stage 23 Stage 3 jobs CSV")

    route_manifest_path, route_manifest, route_manifest_sha256 = load_route_manifest(stage3_root)
    validate_route_project_config(args.project_config, route_manifest)
    stage3_route_provenance = route_provenance_fields(
        route_manifest_path,
        route_manifest,
        route_manifest_sha256,
    )
    aggregate_sources = _aggregate_source_manifest_lookup(route_manifest)
    manifest_site_label, _, site_mapping_csv = _stage0_site_inputs(
        route_manifest=route_manifest,
        stage0_root=stage0_root,
    )
    site_numbers, hotspot_numbers = _load_site_mapping(site_mapping_csv)

    selection_rows = _read_required_csv(stage2_selection_csv)
    stage2_selection_lookup = _strict_lookup_rows(
        selection_rows,
        "global_backbone_id",
        "Stage 2.5 selection",
    )
    job_rows = _read_required_csv(stage3_jobs_csv)
    stage3_job_lookup = _strict_lookup_rows(
        job_rows,
        "global_backbone_id",
        "Stage 3 job",
    )
    _strict_lookup_rows(job_rows, "stage3_job_id", "Stage 3 job")

    run_group_id = _single_value(job_rows, "run_group_id", "Stage 3 jobs")
    protocol_identity_sha256 = _single_value(
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
    locked_stage3_mode = _single_value(job_rows, "stage3_mode", "Stage 3 jobs")
    if locked_stage3_mode != args.stage3_mode:
        raise RuntimeError(
            f"--stage3-mode={args.stage3_mode} does not match jobs table mode={locked_stage3_mode}"
        )
    _single_value(job_rows, "input_pdb_dir", "Stage 3 jobs")
    _single_value(job_rows, "output_pdb_dir", "Stage 3 jobs")

    stage2_selection_csv_sha256 = _sha256_file(stage2_selection_csv)
    stage3_jobs_csv_sha256 = _sha256_file(stage3_jobs_csv)
    expected_by_backbone: dict[str, dict[Path, tuple[Path, str]]] = {}
    family_ids: set[str] = set()
    for backbone_id, job_row in stage3_job_lookup.items():
        backbone_row = stage2_selection_lookup.get(backbone_id)
        if backbone_row is None:
            raise RuntimeError(
                f"Stage 3 job global backbone is absent from the Stage 2.5 selection CSV: {backbone_id}"
            )
        if str(backbone_row.get("pass_backbone_qc", "")).strip().lower() != "true":
            raise RuntimeError(f"Stage 2.5 global backbone is not pass_backbone_qc=true: {backbone_id}")
        if str(backbone_row.get("stage2_5_selected", "")).strip().lower() != "true":
            raise RuntimeError(f"Stage 2.5 global backbone is not stage2_5_selected=true: {backbone_id}")
        family_id = str(backbone_row.get("backbone_family_id", "")).strip()
        if not family_id:
            raise RuntimeError(f"Stage 2.5 global backbone has no backbone_family_id: {backbone_id}")
        if family_id in family_ids:
            raise RuntimeError(f"Stage 3 jobs contain more than one representative from family {family_id}")
        family_ids.add(family_id)
        if str(backbone_row.get("family_representative_global_backbone_id", "")).strip() != backbone_id:
            raise RuntimeError(f"Stage 2.5 selected row is not its family representative: {backbone_id}")
        if str(backbone_row.get("site_label", "")).strip() != manifest_site_label:
            raise RuntimeError(
                f"Stage 2.5 site label does not match the manifest-locked site for {backbone_id}"
            )
        _validate_cached_route_provenance(
            backbone_row,
            stage3_route_provenance,
            f"Stage 23 Stage 2.5 row {backbone_id}",
        )
        _validate_cached_route_provenance(
            job_row,
            stage3_route_provenance,
            f"Stage 23 Stage 3 job {backbone_id}",
        )
        expected_by_backbone[backbone_id] = _validate_stage3_job_row(
            backbone_id=backbone_id,
            backbone_row=backbone_row,
            job_row=job_row,
            expected_stage3_mode=args.stage3_mode,
            stage3_root=stage3_root,
            stage2_selection_csv=stage2_selection_csv,
            stage2_selection_csv_sha256=stage2_selection_csv_sha256,
            aggregate_sources=aggregate_sources,
        )

    _, runlist_tag_count = _validate_runlist_contract(
        job_rows=job_rows,
        stage3_root=stage3_root,
    )
    expected_output_count = sum(len(paths) for paths in expected_by_backbone.values())
    observed_output_count = _validate_complete_stage3_outputs(expected_by_backbone)
    if observed_output_count != expected_output_count:
        raise RuntimeError(
            "Stage 3B output count changed after completeness validation: "
            f"expected={expected_output_count}, observed={observed_output_count}"
        )

    requested = _split_csv(args.selected_global_backbones)
    if len(requested) != len(set(requested)):
        raise RuntimeError("--selected-global-backbones contains duplicate global backbone IDs")
    if requested:
        missing_requested = sorted(set(requested) - set(stage3_job_lookup))
        if missing_requested:
            raise RuntimeError(
                "Requested global backbone IDs are absent from the authoritative Stage 3 jobs table: "
                f"{missing_requested[:10]}"
            )
        requested_set = set(requested)
        selected_backbones = [
            str(row["global_backbone_id"]).strip()
            for row in job_rows
            if str(row["global_backbone_id"]).strip() in requested_set
        ]
    else:
        selected_backbones = [str(row["global_backbone_id"]).strip() for row in job_rows]
    if not selected_backbones:
        raise RuntimeError("No global backbone IDs were selected from the Stage 3 jobs table")

    if output_root.resolve() != stage3_root.resolve():
        output_manifest_path, output_manifest, output_manifest_sha256 = write_route_manifest(
            output_root,
            route_manifest,
        )
        output_route_provenance = route_provenance_fields(
            output_manifest_path,
            output_manifest,
            output_manifest_sha256,
        )
    else:
        output_route_provenance = stage3_route_provenance

    output_dir = output_root / "05_proteinmpnn_sequences" / "stage3c" / _safe_token(run_group_id)
    fasta_dir = output_dir / "fasta"
    output_prefix = f"FGA_rfpeptides_{_safe_token(run_group_id)}_stage3C"

    all_rows: list[dict[str, Any]] = []
    target_chains: set[str] = set()
    peptide_chains: set[str] = set()
    for backbone_id in selected_backbones:
        backbone_row = stage2_selection_lookup[backbone_id]
        job_row = stage3_job_lookup[backbone_id]
        target_chains.add(str(backbone_row.get("target_chain", "")).strip() or "A")
        peptide_chains.add(str(backbone_row.get("peptide_chain", "")).strip() or "B")
        expected_outputs = expected_by_backbone[backbone_id]
        provenance = _stage3_job_provenance(
            job_row,
            stage3_jobs_csv=stage3_jobs_csv,
            stage3_jobs_csv_sha256=stage3_jobs_csv_sha256,
            expected_count=len(expected_outputs),
            observed_count=len(expected_outputs),
        )
        provenance.update(output_route_provenance)
        for relaxed_pdb, (input_pdb, _) in sorted(
            expected_outputs.items(),
            key=lambda item: str(item[0]),
        ):
            assert_active_route_path(relaxed_pdb, f"Stage 23 Stage 3B PDB for {backbone_id}")
            qc_row = _qc_row_for_relaxed_pdb(
                relaxed_pdb=relaxed_pdb,
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
            qc_row["stage3_input_pdb"] = input_pdb
            qc_row["stage3_input_pdb_sha256"] = _sha256_file(input_pdb)
            qc_row.update(provenance)
            all_rows.append(qc_row)

    if len(target_chains) != 1 or len(peptide_chains) != 1:
        raise RuntimeError(
            "Selected Stage 3 jobs do not share one target/peptide chain contract: "
            f"target={sorted(target_chains)}, peptide={sorted(peptide_chains)}"
        )
    target_chain = next(iter(target_chains))
    peptide_chain = next(iter(peptide_chains))
    pass_rows = [row for row in all_rows if row.get("pass_stage3c_qc") == "true"]
    write_csv(output_dir / f"{output_prefix}_sequences_qc.csv", all_rows, STAGE3C_FIELDS)
    write_csv(output_dir / f"{output_prefix}_sequences_qc_pass.csv", pass_rows, STAGE3C_FIELDS)
    _write_fastas(all_rows, fasta_dir)
    write_markdown(
        output_dir / f"{output_prefix}_sequences_qc.md",
        _summary_markdown(
            rows=all_rows,
            args=args,
            output_dir=output_dir,
            site_numbers=site_numbers,
            hotspot_numbers=hotspot_numbers,
            run_group_id=run_group_id,
            stage2_selection_csv=stage2_selection_csv,
            stage3_jobs_csv=stage3_jobs_csv,
            job_count=len(job_rows),
            selected_backbone_count=len(selected_backbones),
            expected_output_count=expected_output_count,
        ),
    )
    _write_pymol_review(
        rows=all_rows,
        output_path=output_dir
        / f"{_safe_token(manifest_site_label)}_{_safe_token(run_group_id)}_stage3C_sequence_qc_review.pml",
        target_chain=target_chain,
        peptide_chain=peptide_chain,
        site_numbers=site_numbers,
        hotspot_numbers=hotspot_numbers,
        top_n=args.top_pymol,
        pymol_path_style=args.pymol_path_style,
    )

    logger.info("Authoritative Stage 3 jobs: %s", len(job_rows))
    logger.info("Runlist input tags validated: %s", runlist_tag_count)
    logger.info("Complete Stage 3B output set validated: %s", expected_output_count)
    logger.info("Selected global backbones collected: %s", len(selected_backbones))
    logger.info("Parsed Stage 3B PDBs: %s", len(all_rows))
    logger.info("Passed Stage 3C sequence/relax QC: %s", len(pass_rows))
    logger.info("Output directory: %s", output_dir)
    if not pass_rows:
        logger.warning("No Stage 3B outputs passed Stage 3C QC.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
