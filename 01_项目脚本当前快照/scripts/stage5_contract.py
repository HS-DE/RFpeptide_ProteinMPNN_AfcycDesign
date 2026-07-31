from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any, Iterable, Mapping

from common import (
    ROUTE_PROVENANCE_FIELDS,
    SOURCE_ROUTE_PROVENANCE_FIELDS,
    assert_active_route_path,
    load_route_manifest,
    read_csv,
    route_provenance_fields,
    sha256_file,
    validate_route_project_config,
    validate_row_route_provenance,
    validate_source_route_provenance,
)


STAGE4_IDENTITY_FIELDS = [
    "global_backbone_id",
    "backbone_family_id",
    "source_batch_label",
    "source_stage4_design_id",
    "source_stage4_scored_pdb_sha256",
    "stage4_run_id",
    "stage4_protocol_identity_sha256",
    "stage4_run_group_id",
    "stage4_scores_csv",
    "stage4_scores_csv_sha256",
    "stage4_top_candidates_csv",
    "stage4_top_candidates_csv_sha256",
    "stage4_validation_selection_rank",
]

STAGE4_FULL_REQUIRED_FIELDS = {
    "stage4_design_id",
    "stage4_run_id",
    "stage4_protocol_identity_sha256",
    "run_group_id",
    "global_backbone_id",
    "backbone_family_id",
    "source_batch_label",
    "scored_pdb",
    "scored_pdb_sha256",
    "input_pdb",
    "input_pdb_sha256",
    "target_chain",
    "peptide_chain",
    "peptide_sequence",
    "peptide_length",
    "pass_stage4_qc",
    "pyrosetta_score_status",
    "target_site_recovery_status",
    "hotspot_recovery_status",
    "macrocycle_geometry_status",
    "clash_status",
    "detached_or_collapsed_flag",
    "num_target_contacts",
    "num_target_site_contacts",
    "num_hotspot_contacts",
    "stage4_validation_selection_rank",
    "stage4_priority_rank",
    "stage4_priority_class",
    *ROUTE_PROVENANCE_FIELDS,
    *SOURCE_ROUTE_PROVENANCE_FIELDS,
}

STAGE4_TOP_REQUIRED_FIELDS = {
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
    "ddg_proxy_per_peptide_residue",
    "ddg_proxy_per_target_contact",
    "ddg_proxy_per_site_contact",
    "site_contacts",
    "hotspot_contacts",
    "macrocycle_terminal_cn_distance",
    "clash_status",
    "sequence_liability_notes",
    "reason_selected",
}


def resolve_mixed_path(value: str | Path) -> Path:
    text = str(value).strip().replace("\\", "/")
    if os.name == "nt" and text.startswith("/mnt/") and len(text) > 7 and text[6] == "/":
        return Path(f"{text[5].upper()}:{text[6:]}")
    if os.name != "nt" and len(text) >= 3 and text[1] == ":" and text[2] == "/":
        return Path(f"/mnt/{text[0].lower()}{text[2:]}")
    path = Path(text)
    if path.is_absolute():
        return path
    return Path(__file__).resolve().parents[1] / path


def _read_required_csv(path: Path, label: str) -> list[dict[str, str]]:
    rows = read_csv(path)
    if not rows:
        raise RuntimeError(f"Missing or empty {label}: {path}")
    return rows


def _require_fields(rows: Iterable[Mapping[str, Any]], required: set[str], label: str) -> None:
    rows = list(rows)
    observed = set(rows[0]) if rows else set()
    missing = sorted(required - observed)
    if missing:
        raise RuntimeError(f"{label} is missing required fields: {','.join(missing)}")


def _truthy(value: Any) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes", "pass"}


def _require_equal(observed: Any, expected: Any, label: str) -> None:
    observed_text = str(observed).strip()
    expected_text = str(expected).strip()
    if observed_text.lower() in {"true", "false"} and expected_text.lower() in {"true", "false"}:
        observed_text = observed_text.lower()
        expected_text = expected_text.lower()
    if observed_text != expected_text:
        raise RuntimeError(
            f"{label} mismatch: observed={observed_text!r}, "
            f"expected={expected_text!r}"
        )


def _path_key(value: str | Path, label: str) -> str:
    return str(assert_active_route_path(resolve_mixed_path(value), label).resolve())


def _require_under(path: Path, root: Path, label: str) -> None:
    resolved_path = path.resolve()
    resolved_root = root.resolve()
    if resolved_path != resolved_root and resolved_root not in resolved_path.parents:
        raise RuntimeError(f"{label} is outside the declared source root: {resolved_path}")


def _validate_aggregate_route_rows(
    rows: list[dict[str, str]],
    expected_route: Mapping[str, str],
    label: str,
) -> None:
    expected_manifest = _path_key(expected_route["route_manifest_path"], f"{label} expected route manifest")
    source_groups: dict[tuple[str, str, str, str], dict[str, str]] = {}
    for index, row in enumerate(rows, start=1):
        row_label = f"{label} row {index}"
        for field in ROUTE_PROVENANCE_FIELDS:
            if field == "route_manifest_path":
                observed = _path_key(row.get(field, ""), f"{row_label} route manifest")
                if observed != expected_manifest:
                    raise RuntimeError(f"{row_label} route manifest path does not match the source aggregate manifest")
            else:
                _require_equal(row.get(field, ""), expected_route[field], f"{row_label} {field}")
        source_key = tuple(str(row.get(field, "")).strip() for field in SOURCE_ROUTE_PROVENANCE_FIELDS)
        if any(not value for value in source_key):
            raise RuntimeError(f"{row_label} has incomplete source route provenance")
        source_groups.setdefault(source_key, row)

    for source_key, representative in source_groups.items():
        validate_source_route_provenance(
            representative,
            f"{label} source route {source_key[0]}",
        )


def _validate_stage4_hard_gate(row: Mapping[str, Any], label: str) -> None:
    requirements = {
        "pyrosetta_score_status": "success",
        "target_site_recovery_status": "site_contact_pass",
        "hotspot_recovery_status": "hotspot_contact_pass",
        "macrocycle_geometry_status": "pass_head_to_tail_macrocycle",
        "clash_status": "pass_no_severe_clash",
        "detached_or_collapsed_flag": "pass_basic_pose_geometry",
    }
    if not _truthy(row.get("pass_stage4_qc", "")):
        raise RuntimeError(f"{label} is not marked pass_stage4_qc=true")
    for field, expected in requirements.items():
        _require_equal(row.get(field, ""), expected, f"{label} {field}")
    if str(row.get("stage4_failure_reasons", "")).strip():
        raise RuntimeError(f"{label} has non-empty stage4_failure_reasons")
    if str(row.get("target_chain", "")).strip() != "A":
        raise RuntimeError(f"{label} target_chain must be A")
    if str(row.get("peptide_chain", "")).strip() != "B":
        raise RuntimeError(f"{label} peptide_chain must be B")
    if int(str(row.get("num_target_site_contacts", "0"))) < 1:
        raise RuntimeError(f"{label} has no Site_2 contacts")
    if int(str(row.get("num_hotspot_contacts", "0"))) < 1:
        raise RuntimeError(f"{label} has no hotspot contacts")


def _validate_stage4_pdb(row: Mapping[str, Any], label: str) -> Path:
    input_pdb = assert_active_route_path(
        resolve_mixed_path(row.get("input_pdb", "")),
        f"{label} Stage 3D-1 input PDB",
    )
    scored_pdb = assert_active_route_path(
        resolve_mixed_path(row.get("scored_pdb", "")),
        f"{label} Stage 4 scored PDB",
    )
    input_sha256 = sha256_file(input_pdb)
    scored_sha256 = sha256_file(scored_pdb)
    _require_equal(input_sha256, row.get("input_pdb_sha256", ""), f"{label} input PDB SHA-256")
    _require_equal(scored_sha256, row.get("scored_pdb_sha256", ""), f"{label} scored PDB SHA-256")
    _require_equal(input_sha256, scored_sha256, f"{label} score-only structure preservation")
    return scored_pdb


def _top_full_compare_fields() -> tuple[tuple[str, str], ...]:
    return (
        ("stage4_validation_selection_rank", "stage4_validation_selection_rank"),
        ("stage4_priority_rank", "stage4_priority_rank"),
        ("stage4_priority_class", "stage4_priority_class"),
        ("stage4_run_id", "stage4_run_id"),
        ("stage4_protocol_identity_sha256", "stage4_protocol_identity_sha256"),
        ("run_group_id", "run_group_id"),
        ("global_backbone_id", "global_backbone_id"),
        ("backbone_family_id", "backbone_family_id"),
        ("source_batch_label", "source_batch_label"),
        ("peptide_sequence", "peptide_sequence"),
        ("peptide_length", "peptide_length"),
        ("scored_pdb_sha256", "scored_pdb_sha256"),
        ("ddg_proxy_no_repack", "ddg_proxy_no_repack"),
        ("ddg_proxy_per_peptide_residue", "ddg_proxy_per_peptide_residue"),
        ("ddg_proxy_per_target_contact", "ddg_proxy_per_target_contact"),
        ("ddg_proxy_per_site_contact", "ddg_proxy_per_site_contact"),
        ("site_contacts", "num_target_site_contacts"),
        ("hotspot_contacts", "num_hotspot_contacts"),
        ("macrocycle_terminal_cn_distance", "macrocycle_terminal_cn_distance"),
        ("clash_status", "clash_status"),
        ("sequence_liability_notes", "sequence_liability_notes"),
    )


def load_stage4_validation_contract(
    *,
    source_run_root: str | Path,
    stage4_scores_csv: str | Path,
    stage4_top_candidates_csv: str | Path,
    project_config: str | Path,
    candidate_count: int,
) -> dict[str, Any]:
    if candidate_count < 1:
        raise RuntimeError("candidate_count must be positive")

    source_root = assert_active_route_path(
        resolve_mixed_path(source_run_root),
        "Stage 5 source run root",
    )
    manifest_path, manifest, manifest_sha256 = load_route_manifest(source_root)
    validate_route_project_config(project_config, manifest)
    route_provenance = route_provenance_fields(manifest_path, manifest, manifest_sha256)

    scores_csv = assert_active_route_path(
        resolve_mixed_path(stage4_scores_csv),
        "Stage 5 Stage 4 complete score CSV",
    )
    top_csv = assert_active_route_path(
        resolve_mixed_path(stage4_top_candidates_csv),
        "Stage 5 Stage 4 top candidate CSV",
    )
    _require_under(scores_csv, source_root, "Stage 4 complete score CSV")
    _require_under(top_csv, source_root, "Stage 4 top candidate CSV")
    if scores_csv.parent.resolve() != top_csv.parent.resolve():
        raise RuntimeError("Stage 4 complete and top candidate CSVs must come from the same isolated run directory")
    relative_parent = scores_csv.parent.resolve().relative_to(source_root.resolve())
    if not relative_parent.parts or relative_parent.parts[0] != "06_rosetta_scoring":
        raise RuntimeError("Stage 4 tables are not under the source run 06_rosetta_scoring directory")

    score_rows = _read_required_csv(scores_csv, "Stage 4 complete score CSV")
    top_rows = _read_required_csv(top_csv, "Stage 4 top candidate CSV")
    _require_fields(score_rows, STAGE4_FULL_REQUIRED_FIELDS, "Stage 4 complete score CSV")
    _require_fields(top_rows, STAGE4_TOP_REQUIRED_FIELDS, "Stage 4 top candidate CSV")
    if len(top_rows) != candidate_count:
        raise RuntimeError(
            f"Stage 4 top candidate table has {len(top_rows)} rows, expected exactly {candidate_count}"
        )

    _validate_aggregate_route_rows(score_rows, route_provenance, "Stage 4 complete score table")
    stage4_run_ids = {str(row["stage4_run_id"]).strip() for row in score_rows}
    protocol_hashes = {str(row["stage4_protocol_identity_sha256"]).strip() for row in score_rows}
    run_group_ids = {str(row["run_group_id"]).strip() for row in score_rows}
    if len(stage4_run_ids) != 1 or "" in stage4_run_ids:
        raise RuntimeError("Stage 4 complete score table mixes stage4_run_id values")
    if len(protocol_hashes) != 1 or "" in protocol_hashes:
        raise RuntimeError("Stage 4 complete score table mixes protocol identities")
    if len(run_group_ids) != 1 or "" in run_group_ids:
        raise RuntimeError("Stage 4 complete score table mixes Stage 3 run groups")
    stage4_run_id = next(iter(stage4_run_ids))
    protocol_identity = next(iter(protocol_hashes))
    run_group_id = next(iter(run_group_ids))
    if scores_csv.parent.name != stage4_run_id:
        raise RuntimeError("Stage 4 table directory name does not match stage4_run_id")

    full_lookup: dict[str, dict[str, str]] = {}
    for row in score_rows:
        design_id = str(row["stage4_design_id"]).strip()
        if not design_id or design_id in full_lookup:
            raise RuntimeError(f"Stage 4 complete score table has a duplicate/empty design ID: {design_id!r}")
        full_lookup[design_id] = row

    expected_ranks = list(range(1, candidate_count + 1))
    observed_ranks = [int(str(row["stage4_validation_selection_rank"])) for row in top_rows]
    if observed_ranks != expected_ranks:
        raise RuntimeError(
            f"Stage 4 top candidate ranks must be contiguous 1..{candidate_count}; observed={observed_ranks}"
        )
    selected_global_ids: set[str] = set()
    selected_family_ids: set[str] = set()
    selected_rows: list[dict[str, Any]] = []
    scores_sha256 = sha256_file(scores_csv)
    top_sha256 = sha256_file(top_csv)
    for rank, top_row in enumerate(top_rows, start=1):
        design_id = str(top_row["stage4_design_id"]).strip()
        full_row = full_lookup.get(design_id)
        if full_row is None:
            raise RuntimeError(f"Stage 4 top candidate is absent from complete score table: {design_id}")
        label = f"Stage 4 top candidate rank {rank} ({design_id})"
        for top_field, full_field in _top_full_compare_fields():
            _require_equal(top_row.get(top_field, ""), full_row.get(full_field, ""), f"{label} {top_field}")
        _require_equal(
            _path_key(top_row["scored_pdb"], f"{label} top scored PDB"),
            _path_key(full_row["scored_pdb"], f"{label} full scored PDB"),
            f"{label} scored PDB path",
        )
        _validate_stage4_hard_gate(full_row, label)
        scored_pdb = _validate_stage4_pdb(full_row, label)

        global_id = str(full_row["global_backbone_id"]).strip()
        family_id = str(full_row["backbone_family_id"]).strip()
        if not global_id or global_id in selected_global_ids:
            raise RuntimeError(f"Stage 4 top candidates do not have unique global_backbone_id values: {global_id!r}")
        if not family_id or family_id in selected_family_ids:
            raise RuntimeError(f"Stage 4 top candidates do not cover unique backbone families: {family_id!r}")
        selected_global_ids.add(global_id)
        selected_family_ids.add(family_id)

        selected = dict(full_row)
        selected["scored_pdb"] = str(scored_pdb)
        selected["selection_reason"] = str(top_row.get("reason_selected", "")).strip()
        selected["batch"] = str(full_row["source_batch_label"]).strip()
        selected["backbone_id"] = global_id
        selected.update(
            {
                "source_stage4_design_id": design_id,
                "source_stage4_scored_pdb_sha256": str(full_row["scored_pdb_sha256"]).strip(),
                "stage4_run_group_id": run_group_id,
                "stage4_scores_csv": str(scores_csv),
                "stage4_scores_csv_sha256": scores_sha256,
                "stage4_top_candidates_csv": str(top_csv),
                "stage4_top_candidates_csv_sha256": top_sha256,
            }
        )
        selected_rows.append(selected)

    return {
        "source_root": source_root,
        "source_route_manifest_path": manifest_path,
        "source_route_manifest": manifest,
        "source_route_manifest_sha256": manifest_sha256,
        "source_route_provenance": route_provenance,
        "stage4_scores_csv": scores_csv,
        "stage4_scores_csv_sha256": scores_sha256,
        "stage4_top_candidates_csv": top_csv,
        "stage4_top_candidates_csv_sha256": top_sha256,
        "stage4_run_id": stage4_run_id,
        "stage4_protocol_identity_sha256": protocol_identity,
        "stage4_run_group_id": run_group_id,
        "selected_rows": selected_rows,
    }


def stage4_identity_values(row: Mapping[str, Any]) -> dict[str, Any]:
    return {field: row.get(field, "") for field in STAGE4_IDENTITY_FIELDS}


def validate_stage5_identity_link(
    *,
    candidate: Mapping[str, Any],
    job: Mapping[str, Any],
    candidate_id_field: str,
    job_candidate_id_field: str,
    label: str,
) -> None:
    _require_equal(
        job.get(job_candidate_id_field, ""),
        candidate.get(candidate_id_field, ""),
        f"{label} candidate ID",
    )
    fields = [
        "global_backbone_id",
        "backbone_family_id",
        "source_batch_label",
        "source_stage4_design_id",
        "source_stage4_scored_pdb_sha256",
        "stage4_run_id",
        "stage4_protocol_identity_sha256",
        "stage4_run_group_id",
        "stage4_scores_csv_sha256",
        "stage4_top_candidates_csv_sha256",
        "stage4_validation_selection_rank",
        "peptide_sequence_hash",
        "protocol_hash",
    ]
    for field in fields:
        _require_equal(job.get(field, ""), candidate.get(field, ""), f"{label} {field}")


def load_prepared_stage5_contract(
    *,
    stage_output_dir: str | Path,
    candidate_manifest_csv: str | Path,
    jobs_csv: str | Path,
    project_config: str | Path,
    candidate_id_field: str,
    job_id_field: str,
    job_candidate_id_field: str,
) -> dict[str, Any]:
    output_dir = assert_active_route_path(
        resolve_mixed_path(stage_output_dir),
        "Stage 5 prepared output directory",
    )
    route_root = output_dir.parent
    manifest_path, route_manifest, manifest_sha256 = load_route_manifest(route_root)
    validate_route_project_config(project_config, route_manifest)
    route_provenance = route_provenance_fields(manifest_path, route_manifest, manifest_sha256)

    candidate_csv = assert_active_route_path(
        resolve_mixed_path(candidate_manifest_csv),
        "Stage 5 candidate manifest CSV",
    )
    jobs_path = assert_active_route_path(
        resolve_mixed_path(jobs_csv),
        "Stage 5 prediction jobs CSV",
    )
    _require_under(candidate_csv, output_dir, "Stage 5 candidate manifest CSV")
    _require_under(jobs_path, output_dir, "Stage 5 prediction jobs CSV")
    candidates = _read_required_csv(candidate_csv, "Stage 5 candidate manifest CSV")
    jobs = _read_required_csv(jobs_path, "Stage 5 prediction jobs CSV")
    required_common = {
        "global_backbone_id",
        "backbone_family_id",
        "source_batch_label",
        "source_stage4_design_id",
        "source_stage4_scored_pdb_sha256",
        "stage4_run_id",
        "stage4_protocol_identity_sha256",
        "stage4_run_group_id",
        "stage4_scores_csv",
        "stage4_scores_csv_sha256",
        "stage4_top_candidates_csv",
        "stage4_top_candidates_csv_sha256",
        "stage4_validation_selection_rank",
        "peptide_sequence_hash",
        "protocol_hash",
        *ROUTE_PROVENANCE_FIELDS,
        *SOURCE_ROUTE_PROVENANCE_FIELDS,
    }
    _require_fields(candidates, required_common | {candidate_id_field}, "Stage 5 candidate manifest")
    _require_fields(
        jobs,
        required_common
        | {
            job_id_field,
            job_candidate_id_field,
            "seed",
            "job_spec_json",
            "models_per_seed",
        },
        "Stage 5 prediction jobs table",
    )

    candidate_lookup: dict[str, dict[str, str]] = {}
    global_ids: set[str] = set()
    stage4_design_ids: set[str] = set()
    for candidate in candidates:
        candidate_id = str(candidate.get(candidate_id_field, "")).strip()
        if not candidate_id or candidate_id in candidate_lookup:
            raise RuntimeError(f"Stage 5 candidate manifest has duplicate/empty ID: {candidate_id!r}")
        global_id = str(candidate.get("global_backbone_id", "")).strip()
        stage4_design_id = str(candidate.get("source_stage4_design_id", "")).strip()
        if not global_id or global_id in global_ids:
            raise RuntimeError(f"Stage 5 candidate manifest has duplicate/empty global_backbone_id: {global_id!r}")
        if not stage4_design_id or stage4_design_id in stage4_design_ids:
            raise RuntimeError(
                f"Stage 5 candidate manifest has duplicate/empty source_stage4_design_id: {stage4_design_id!r}"
            )
        validate_row_route_provenance(candidate, route_provenance, f"Stage 5 candidate {candidate_id}")
        validate_source_route_provenance(candidate, f"Stage 5 candidate {candidate_id}")
        scored_pdb = candidate.get("staged_design_pdb") or candidate.get("staged_reference_design_pdb")
        if not scored_pdb:
            raise RuntimeError(f"Stage 5 candidate {candidate_id} lacks a staged Stage 4 reference PDB")
        scored_path = assert_active_route_path(
            resolve_mixed_path(scored_pdb),
            f"Stage 5 candidate {candidate_id} staged Stage 4 reference PDB",
        )
        _require_equal(
            sha256_file(scored_path),
            candidate["source_stage4_scored_pdb_sha256"],
            f"Stage 5 candidate {candidate_id} staged reference SHA-256",
        )
        for path_field, sha_field in (
            ("stage4_scores_csv", "stage4_scores_csv_sha256"),
            ("stage4_top_candidates_csv", "stage4_top_candidates_csv_sha256"),
        ):
            source_table = assert_active_route_path(
                resolve_mixed_path(candidate[path_field]),
                f"Stage 5 candidate {candidate_id} {path_field}",
            )
            _require_equal(
                sha256_file(source_table),
                candidate[sha_field],
                f"Stage 5 candidate {candidate_id} {sha_field}",
            )
        candidate_lookup[candidate_id] = candidate
        global_ids.add(global_id)
        stage4_design_ids.add(stage4_design_id)

    job_lookup: dict[str, dict[str, str]] = {}
    seeds_by_candidate: dict[str, set[int]] = {candidate_id: set() for candidate_id in candidate_lookup}
    for job in jobs:
        job_id = str(job.get(job_id_field, "")).strip()
        candidate_id = str(job.get(job_candidate_id_field, "")).strip()
        if not job_id or job_id in job_lookup:
            raise RuntimeError(f"Stage 5 jobs table has duplicate/empty job ID: {job_id!r}")
        candidate = candidate_lookup.get(candidate_id)
        if candidate is None:
            raise RuntimeError(f"Stage 5 job {job_id} references unknown candidate {candidate_id!r}")
        validate_row_route_provenance(job, route_provenance, f"Stage 5 job {job_id}")
        validate_source_route_provenance(job, f"Stage 5 job {job_id}")
        validate_stage5_identity_link(
            candidate=candidate,
            job=job,
            candidate_id_field=candidate_id_field,
            job_candidate_id_field=job_candidate_id_field,
            label=f"Stage 5 job {job_id}",
        )
        seed = int(str(job.get("seed", "")))
        if seed in seeds_by_candidate[candidate_id]:
            raise RuntimeError(f"Stage 5 candidate {candidate_id} has duplicate seed {seed}")
        seeds_by_candidate[candidate_id].add(seed)

        spec_path = assert_active_route_path(
            resolve_mixed_path(job.get("job_spec_json", "")),
            f"Stage 5 job {job_id} spec JSON",
        )
        try:
            spec = json.loads(spec_path.read_text(encoding="utf-8"))
        except (json.JSONDecodeError, UnicodeDecodeError) as exc:
            raise RuntimeError(f"Invalid Stage 5 job spec JSON: {spec_path}") from exc
        if not isinstance(spec, dict):
            raise RuntimeError(f"Stage 5 job spec must be a JSON object: {spec_path}")
        _require_equal(spec.get(job_id_field, ""), job_id, f"Stage 5 job {job_id} spec job ID")
        validate_stage5_identity_link(
            candidate=candidate,
            job=spec,
            candidate_id_field=candidate_id_field,
            job_candidate_id_field=job_candidate_id_field,
            label=f"Stage 5 job {job_id} spec",
        )
        for field in (
            "seed",
            "peptide_sequence",
            "peptide_sequence_hash",
            "protocol_hash",
            "requested_recycles",
            "forward_passes",
            "target_msa_mode",
            "peptide_msa_mode",
            "template_mode",
            "use_initial_guess",
        ):
            _require_equal(spec.get(field, ""), job.get(field, ""), f"Stage 5 job {job_id} spec {field}")
        job_lookup[job_id] = job

    for candidate_id, candidate in candidate_lookup.items():
        expected_seed_count = int(str(candidate.get("seeds_per_candidate", "0")))
        expected_seeds = set(range(expected_seed_count))
        if seeds_by_candidate[candidate_id] != expected_seeds:
            raise RuntimeError(
                f"Stage 5 candidate {candidate_id} seed set mismatch: "
                f"observed={sorted(seeds_by_candidate[candidate_id])}, expected={sorted(expected_seeds)}"
            )
        models_per_seed = int(str(candidate.get("models_per_seed", "0")))
        if models_per_seed < 1:
            raise RuntimeError(f"Stage 5 candidate {candidate_id} has invalid models_per_seed")

    return {
        "output_dir": output_dir,
        "route_root": route_root,
        "route_manifest_path": manifest_path,
        "route_manifest": route_manifest,
        "route_manifest_sha256": manifest_sha256,
        "route_provenance": route_provenance,
        "candidate_manifest_csv": candidate_csv,
        "jobs_csv": jobs_path,
        "candidates": candidates,
        "candidate_lookup": candidate_lookup,
        "jobs": jobs,
        "job_lookup": job_lookup,
    }
