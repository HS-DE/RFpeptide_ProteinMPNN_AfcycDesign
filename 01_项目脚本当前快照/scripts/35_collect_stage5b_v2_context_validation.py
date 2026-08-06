from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import math
import os
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np

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
    sha256_file,
    validate_route_project_config,
    validate_row_route_provenance,
    write_csv,
    write_markdown,
)
from pdb_utils import parse_residues, residue_sequence
from stage5_contract import STAGE4_IDENTITY_FIELDS


SCRIPT_DIR = Path(__file__).resolve().parent
STAGE5B_V1_COLLECTOR = SCRIPT_DIR / "31_collect_stage5b_validation.py"
TOP5_PROTOCOL_VERSION = "stage5B_v2_C1_C3_full_target_template_top5_v1"
ALL_PASS_PROTOCOL_VERSION = "stage5B_v2_C1_C3_full_target_template_allpass_v1"
SUPPORTED_PROTOCOL_VERSIONS = {TOP5_PROTOCOL_VERSION, ALL_PASS_PROTOCOL_VERSION}
PROTOCOL_VERSION = TOP5_PROTOCOL_VERSION
VALIDATION_TEST_TYPE = "target_context_structure_conditioned_recovery"
BACKBONE_ATOMS = ("N", "CA", "C")

MODEL_FIELDS = [
    "stage5B_v2_candidate_context_id",
    "stage5B_v2_candidate_id",
    "stage5B_v2_job_id",
    "context_id",
    "protocol_hash",
    "stage5_campaign_id",
    "stage5_selection_mode",
    "selection_order",
    "batch",
    "backbone_id",
    *STAGE4_IDENTITY_FIELDS,
    "peptide_sequence",
    "seed",
    "model_name",
    "prediction_pdb",
    "prediction_npz",
    "protocol_identity_valid",
    "context_sequence_match",
    "peptide_sequence_match",
    "context_template_coverage",
    "peptide_template_coverage",
    "template_mode",
    "use_initial_guess",
    "reference_design_loaded_by_prediction_runner",
    "fga_global_CA_RMSD_A",
    "stage0_crop_local_CA_RMSD_A",
    "Site_2_RMSD_after_fga_alignment_A",
    "hotspot_RMSD_after_fga_alignment_A",
    "Site_2_RMSD_after_crop_alignment_A",
    "hotspot_RMSD_after_crop_alignment_A",
    "context_after_fga_alignment_CA_RMSD_A",
    "partner_after_fga_alignment_CA_RMSD_A",
    "target_context_recovery_pass",
    "target_context_recovery_failure_reasons",
    "target_aligned_peptide_backbone_RMSD_A",
    "target_aligned_peptide_interface_RMSD_A",
    "reference_design_site2_contact_count",
    "reference_design_hotspot_contact_count",
    "reference_design_hotspot_min_distance_A",
    "num_context_contacts",
    "num_fga_contacts",
    "num_site2_contacts",
    "num_hotspot_contacts",
    "num_partner_contacts",
    "hotspot_min_distance_A",
    "off_site_contact_count",
    "off_site_contact_fraction",
    "same_target_site_flag",
    "off_site_binding_dominant_flag",
    "terminal_C_N_distance_A",
    "macrocycle_geometry_status",
    "severe_clash_count",
    "clash_status",
    "peptide_radius_of_gyration_A",
    "context_mean_pLDDT_100",
    "fga_mean_pLDDT_100",
    "Site_2_mean_pLDDT_100",
    "hotspot_mean_pLDDT_100",
    "peptide_mean_pLDDT_100",
    "interface_pae_context_mean_A",
    "interface_pae_fga_mean_A",
    "interface_pae_site2_mean_A",
    "interface_pae_site2_median_A",
    "interface_pae_hotspot_mean_A",
    "pTM",
    "ipTM",
    "ranking_confidence",
    "pose_recovery_class",
    "model_support_status",
    "prediction_coordinate_hash",
    "prediction_confidence_hash",
    "requested_recycles",
    "forward_passes",
    "loaded_colabdesign_commit",
    "cyclic_topology_encoding",
    "notes",
    *SOURCE_ROUTE_PROVENANCE_FIELDS,
    *ROUTE_PROVENANCE_FIELDS,
]

CONTEXT_SUMMARY_FIELDS = [
    "stage5B_v2_candidate_context_id",
    "stage5B_v2_candidate_id",
    "context_id",
    "stage5_selection_mode",
    "selection_order",
    "batch",
    "backbone_id",
    "global_backbone_id",
    "backbone_family_id",
    "peptide_sequence",
    "models_expected",
    "models_completed",
    "effective_unique_prediction_count",
    "duplicate_prediction_count",
    "target_context_recovery_count",
    "same_target_site_count",
    "strong_pose_recovery_count",
    "moderate_pose_recovery_count",
    "macrocycle_pass_count",
    "no_severe_clash_count",
    "median_fga_global_CA_RMSD_A",
    "best_fga_global_CA_RMSD_A",
    "median_stage0_crop_local_CA_RMSD_A",
    "best_stage0_crop_local_CA_RMSD_A",
    "median_peptide_RMSD_A",
    "best_peptide_RMSD_A",
    "median_hotspot_distance_A",
    "best_hotspot_distance_A",
    "median_peptide_pLDDT_100",
    "median_hotspot_PAE_A",
    "best_model_path",
    "context_support_class",
    "support_reason",
    *SOURCE_ROUTE_PROVENANCE_FIELDS,
    *ROUTE_PROVENANCE_FIELDS,
]

CANDIDATE_SUMMARY_FIELDS = [
    "stage5B_v2_candidate_id",
    "stage5_selection_mode",
    "selection_order",
    "batch",
    "backbone_id",
    "global_backbone_id",
    "backbone_family_id",
    "peptide_sequence",
    "C1_models_completed",
    "C1_target_recovery_count",
    "C1_same_site_count",
    "C1_strong_pose_count",
    "C1_moderate_pose_count",
    "C1_best_peptide_RMSD_A",
    "C1_best_hotspot_distance_A",
    "C1_context_support_class",
    "C3_models_completed",
    "C3_target_recovery_count",
    "C3_same_site_count",
    "C3_strong_pose_count",
    "C3_moderate_pose_count",
    "C3_best_peptide_RMSD_A",
    "C3_best_hotspot_distance_A",
    "C3_context_support_class",
    "cross_context_support_class",
    "cross_context_interpretation",
    *SOURCE_ROUTE_PROVENANCE_FIELDS,
    *ROUTE_PROVENANCE_FIELDS,
]


def _load_stage5b_v1_helpers() -> Any:
    spec = importlib.util.spec_from_file_location("stage5b_v1_collect_helpers", STAGE5B_V1_COLLECTOR)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not import Stage 5B geometry helpers: {STAGE5B_V1_COLLECTOR}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _resolve_mixed_path(value: str | Path) -> Path:
    text = str(value).strip().replace("\\", "/")
    if os.name == "nt" and text.startswith("/mnt/") and len(text) > 7 and text[6] == "/":
        return Path(f"{text[5].upper()}:/{text[7:]}")
    if os.name != "nt" and len(text) >= 3 and text[1] == ":" and text[2] == "/":
        return Path(f"/mnt/{text[0].lower()}{text[2:]}")
    path = Path(text)
    return path if path.is_absolute() else resolve_path(path)


def _truthy(value: Any) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes"}


def _split_ints(value: Any) -> list[int]:
    return [int(item.strip()) for item in str(value).split(",") if item.strip()]


def _split_strings(value: Any) -> list[str]:
    return [item.strip() for item in str(value).split(",") if item.strip()]


def _float(value: Any, default: float = float("nan")) -> float:
    try:
        parsed = float(str(value).strip())
    except (TypeError, ValueError):
        return default
    return parsed if math.isfinite(parsed) else default


def _rounded(value: float | None, digits: int = 3) -> float | str:
    return round(value, digits) if value is not None and math.isfinite(value) else ""


def _median(rows: Sequence[Mapping[str, Any]], field: str) -> float | str:
    values = [_float(row.get(field)) for row in rows]
    finite = [value for value in values if math.isfinite(value)]
    return round(float(np.median(finite)), 3) if finite else ""


def _best(rows: Sequence[Mapping[str, Any]], field: str) -> float | str:
    values = [_float(row.get(field)) for row in rows]
    finite = [value for value in values if math.isfinite(value)]
    return round(min(finite), 3) if finite else ""


def _distribution(rows: Sequence[Mapping[str, Any]], field: str) -> str:
    values = [_float(row.get(field)) for row in rows]
    finite = [value for value in values if math.isfinite(value)]
    if not finite:
        return "not_available"
    return f"{min(finite):.3f} / {float(np.median(finite)):.3f} / {max(finite):.3f}"


def _sha1_arrays(*arrays: np.ndarray) -> str:
    digest = hashlib.sha1()
    for array in arrays:
        contiguous = np.ascontiguousarray(array)
        digest.update(str(contiguous.shape).encode("ascii"))
        digest.update(str(contiguous.dtype).encode("ascii"))
        digest.update(contiguous.tobytes())
    return digest.hexdigest()


def _validate_manifest_and_jobs(
    *,
    stage5_root: Path,
    manifest_csv: Path,
    jobs_csv: Path,
    project_config: str,
) -> tuple[list[dict[str, str]], list[dict[str, str]], dict[str, str], dict[str, Any]]:
    route_path, route_manifest, route_sha256 = load_route_manifest(stage5_root.parent)
    validate_route_project_config(project_config, route_manifest)
    protocol = route_manifest.get("stage5B_v2_protocol", {})
    protocol_version = str(protocol.get("protocol_version", ""))
    selection_mode = str(
        route_manifest.get("stage5_selection_mode", protocol.get("candidate_selection_mode", "top_validation"))
    )
    expected_protocol_version = (
        TOP5_PROTOCOL_VERSION if selection_mode == "top_validation" else ALL_PASS_PROTOCOL_VERSION
    )
    if (
        protocol_version not in SUPPORTED_PROTOCOL_VERSIONS
        or selection_mode not in {"top_validation", "all_stage4_pass"}
        or protocol_version != expected_protocol_version
    ):
        raise RuntimeError("Stage 5B-v2 route manifest protocol version mismatch")
    if protocol.get("contexts") != ["C1_crop86_full_template", "C3_native_GHI301_full_template"]:
        raise RuntimeError("Stage 5B-v2 route manifest does not contain ordered C1/C3 contexts")
    if any(protocol.get(field) for field in (
        "template_sequence_masked",
        "template_sidechains_masked",
        "template_interchain_features_masked",
        "use_initial_guess",
        "reference_design_loaded_by_prediction_runner",
        "use_mlm",
        "use_dropout",
    )):
        raise RuntimeError("Stage 5B-v2 route manifest enables a forbidden protocol option")
    provenance = route_provenance_fields(route_path, route_manifest, route_sha256)
    manifests = read_csv(manifest_csv)
    jobs = read_csv(jobs_csv)
    for row in manifests:
        row.setdefault("stage5_selection_mode", selection_mode)
    for row in jobs:
        row.setdefault("stage5_selection_mode", selection_mode)
    candidate_count = int(route_manifest.get("stage5_candidate_count", 0))
    expected_manifest_rows = candidate_count * 2
    if candidate_count < 1 or len(manifests) != expected_manifest_rows:
        raise RuntimeError(
            f"Stage 5B-v2 requires {expected_manifest_rows} candidate-context rows, observed {len(manifests)}"
        )
    if not jobs:
        raise RuntimeError("Stage 5B-v2 job table is empty")
    expected_jobs = int(route_manifest.get("stage5_seed_job_count", 0))
    if expected_jobs < 1 or len(jobs) != expected_jobs:
        raise RuntimeError(f"Stage 5B-v2 requires {expected_jobs} jobs, observed {len(jobs)}")
    for row in manifests:
        validate_row_route_provenance(row, provenance, f"Stage 35 manifest {row.get('stage5B_v2_candidate_context_id')}")
    for row in jobs:
        validate_row_route_provenance(row, provenance, f"Stage 35 job {row.get('stage5B_v2_job_id')}")
    pair_ids = [row["stage5B_v2_candidate_context_id"] for row in manifests]
    if len(set(pair_ids)) != len(pair_ids):
        raise RuntimeError("Stage 5B-v2 candidate-context IDs are not unique")
    grouped_contexts: dict[str, set[str]] = defaultdict(set)
    for row in manifests:
        grouped_contexts[row["stage5B_v2_candidate_id"]].add(row["context_id"])
    expected_contexts = {"C1_crop86_full_template", "C3_native_GHI301_full_template"}
    if len(grouped_contexts) != candidate_count or any(value != expected_contexts for value in grouped_contexts.values()):
        raise RuntimeError("Each selected candidate must have exactly C1 and C3 rows")
    manifest_lookup = {row["stage5B_v2_candidate_context_id"]: row for row in manifests}
    job_ids: set[str] = set()
    for job in jobs:
        job_id = job["stage5B_v2_job_id"]
        if job_id in job_ids:
            raise RuntimeError(f"Duplicate Stage 5B-v2 job ID: {job_id}")
        job_ids.add(job_id)
        candidate = manifest_lookup.get(job["stage5B_v2_candidate_context_id"])
        if candidate is None:
            raise RuntimeError(f"Job {job_id} has no candidate-context manifest row")
        for field in (
            "stage5B_v2_candidate_id",
            "context_id",
            "peptide_sequence_hash",
            "protocol_hash",
            "stage5_campaign_id",
            "stage5_selection_mode",
        ):
            if str(job.get(field)) != str(candidate.get(field)):
                raise RuntimeError(f"Job {job_id} differs from manifest on {field}")
        spec_path = assert_active_route_path(job["job_spec_json"], f"Stage 35 job spec {job_id}")
        if sha256_file(spec_path) != str(job["job_spec_sha256"]):
            raise RuntimeError(f"Job spec SHA-256 mismatch: {job_id}")
        spec = json.loads(spec_path.read_text(encoding="utf-8"))
        for field in (
            "stage5B_v2_job_id",
            "stage5B_v2_candidate_context_id",
            "stage5B_v2_candidate_id",
            "context_id",
            "peptide_sequence_hash",
            "protocol_hash",
            "stage5_campaign_id",
            "stage5_selection_mode",
        ):
            if str(spec.get(field)) != str(job.get(field)):
                raise RuntimeError(f"Job spec {job_id} differs from job table on {field}")
        if spec.get("protocol_version") != protocol_version or spec.get("validation_test_type") != VALIDATION_TEST_TYPE:
            raise RuntimeError(f"Job spec {job_id} has the wrong protocol identity")
        if spec.get("template_mode") != "target_context_full" or any(bool(spec.get(field)) for field in (
            "template_sequence_masked",
            "template_sidechains_masked",
            "template_interchain_features_masked",
            "use_initial_guess",
            "reference_design_loaded_by_prediction_runner",
            "use_mlm",
            "use_dropout",
        )):
            raise RuntimeError(f"Job spec {job_id} violates Stage 5B-v2 protocol invariants")
        if int(spec.get("peptide_template_expected_coverage", -1)) != 0:
            raise RuntimeError(f"Job spec {job_id} has nonzero peptide template coverage")
        context_pdb = assert_active_route_path(spec["context_pdb"], f"Stage 35 context PDB {job_id}")
        mapping_csv_path = assert_active_route_path(spec["context_mapping_csv"], f"Stage 35 mapping {job_id}")
        reference_pdb = assert_active_route_path(
            spec["reference_design_pdb_for_posthoc"], f"Stage 35 posthoc reference {job_id}"
        )
        if sha256_file(context_pdb) != spec["context_pdb_sha256"]:
            raise RuntimeError(f"Context PDB SHA-256 mismatch: {job_id}")
        if sha256_file(mapping_csv_path) != spec["context_mapping_csv_sha256"]:
            raise RuntimeError(f"Context mapping SHA-256 mismatch: {job_id}")
        if sha256_file(reference_pdb) != spec["reference_design_pdb_sha256"]:
            raise RuntimeError(f"Posthoc reference SHA-256 mismatch: {job_id}")
    return manifests, jobs, provenance, {
        "selection_mode": selection_mode,
        "protocol_version": protocol_version,
        "candidate_count": candidate_count,
    }


def _metadata_valid(metadata: Mapping[str, Any], job: Mapping[str, str], candidate: Mapping[str, str]) -> bool:
    expected_protocol_version = (
        TOP5_PROTOCOL_VERSION
        if candidate["stage5_selection_mode"] == "top_validation"
        else ALL_PASS_PROTOCOL_VERSION
    )
    checks = [
        metadata.get("stage5B_v2_job_id") == job["stage5B_v2_job_id"],
        metadata.get("stage5B_v2_candidate_context_id") == candidate["stage5B_v2_candidate_context_id"],
        metadata.get("stage5B_v2_candidate_id") == candidate["stage5B_v2_candidate_id"],
        metadata.get("context_id") == candidate["context_id"],
        metadata.get("protocol_hash") == candidate["protocol_hash"],
        metadata.get("protocol_version") == expected_protocol_version,
        metadata.get("stage5_selection_mode") == candidate["stage5_selection_mode"],
        metadata.get("validation_test_type") == VALIDATION_TEST_TYPE,
        metadata.get("template_mode") == "target_context_full",
        metadata.get("template_sequence_masked") is False,
        metadata.get("template_sidechains_masked") is False,
        metadata.get("template_interchain_features_masked") is False,
        metadata.get("template_input_verified") is True,
        int(metadata.get("context_template_residues_covered", 0)) == int(candidate["context_total_length"]),
        int(metadata.get("peptide_template_residues_covered", -1)) == 0,
        metadata.get("peptide_template_loaded_by_prediction_runner") is False,
        metadata.get("reference_design_loaded_by_prediction_runner") is False,
        metadata.get("initial_guess_loaded_by_prediction_runner") is False,
        metadata.get("cyclic_offset_applied") is True,
        int(metadata.get("cyclic_offset_chain_index", -1)) == int(candidate["cyclic_chain_index"]),
        metadata.get("loaded_colabdesign_commit") == "5ab4efaba2321a6c3c314b82d2fff8e0241f5c2d",
    ]
    return all(checks)


def _model_row(
    *,
    helper: Any,
    candidate: Mapping[str, str],
    job: Mapping[str, str],
    metadata: Mapping[str, Any],
    metric: Mapping[str, str],
    args: argparse.Namespace,
) -> dict[str, Any]:
    context_pdb = _resolve_mixed_path(candidate["context_pdb"])
    reference_pdb = _resolve_mixed_path(candidate["staged_reference_design_pdb"])
    prediction_pdb = _resolve_mixed_path(metric["prediction_pdb"])
    prediction_npz = _resolve_mixed_path(metric["prediction_npz"])
    if not prediction_pdb.is_file() or not prediction_npz.is_file():
        raise RuntimeError(f"Missing Stage 5B-v2 prediction files for {job['stage5B_v2_job_id']}")

    context_chain_ids = _split_strings(candidate["context_chain_ids"])
    context_lengths = _split_ints(candidate["context_chain_lengths"])
    peptide_chain = str(candidate["peptide_chain_id"])
    reference_context = parse_residues(context_pdb)
    reference_design = parse_residues(reference_pdb)
    predicted = parse_residues(prediction_pdb)
    expected_prediction_chains = context_chain_ids + [peptide_chain]
    protocol_valid = _metadata_valid(metadata, job, candidate)
    if list(predicted) != expected_prediction_chains:
        raise RuntimeError(
            f"Prediction chain order {list(predicted)} != expected {expected_prediction_chains}: {prediction_pdb}"
        )
    context_sequence_match = all(
        residue_sequence(predicted[chain]) == residue_sequence(reference_context[chain])
        for chain in context_chain_ids
    )
    peptide_sequence_match = residue_sequence(predicted[peptide_chain]) == candidate["peptide_sequence"]
    if not context_sequence_match or not peptide_sequence_match:
        raise RuntimeError(f"Prediction sequence mismatch: {prediction_pdb}")

    reference_fga = reference_context["A"]
    predicted_fga = predicted["A"]
    reference_target_chain = str(candidate["reference_target_chain"])
    reference_peptide_chain = str(candidate["reference_peptide_chain"])
    if set(reference_design) != {reference_target_chain, reference_peptide_chain}:
        raise RuntimeError(
            f"Stage 4 reference chains {list(reference_design)} do not match "
            f"target={reference_target_chain}, peptide={reference_peptide_chain}"
        )
    reference_design_target = reference_design[reference_target_chain]
    reference_design_peptide = reference_design[reference_peptide_chain]
    predicted_peptide = predicted[peptide_chain]
    if residue_sequence(reference_design_target) != "".join(
        residue_sequence(reference_fga)[index - 1]
        for index in _split_ints(candidate["stage0_crop_context_indices_1based"])
    ):
        raise RuntimeError("Stage 4 reference target does not match the mapped context crop")
    if residue_sequence(reference_design_peptide) != candidate["peptide_sequence"]:
        raise RuntimeError("Stage 4 reference peptide does not match candidate sequence")

    fga_mobile = helper._ca_coords(predicted_fga)
    fga_reference = helper._ca_coords(reference_fga)
    fga_rotation, fga_translation, fga_rmsd = helper._kabsch(fga_mobile, fga_reference)
    crop_indices_0 = np.asarray([value - 1 for value in _split_ints(candidate["stage0_crop_context_indices_1based"])], dtype=int)
    predicted_crop = fga_mobile[crop_indices_0]
    design_target_ca = helper._ca_coords(reference_design_target)
    crop_rotation, crop_translation, crop_rmsd = helper._kabsch(predicted_crop, design_target_ca)
    site = _split_ints(candidate["site2_fga_indices_1based"])
    hotspots = _split_ints(candidate["hotspot_fga_indices_1based"])
    site_after_fga = helper._subset_ca_rmsd(fga_mobile, fga_reference, site, fga_rotation, fga_translation)
    hotspot_after_fga = helper._subset_ca_rmsd(fga_mobile, fga_reference, hotspots, fga_rotation, fga_translation)
    crop_position = {context_index: position + 1 for position, context_index in enumerate(crop_indices_0 + 1)}
    site_crop = [crop_position[value] for value in site]
    hotspot_crop = [crop_position[value] for value in hotspots]
    site_after_crop = helper._subset_ca_rmsd(predicted_crop, design_target_ca, site_crop, crop_rotation, crop_translation)
    hotspot_after_crop = helper._subset_ca_rmsd(predicted_crop, design_target_ca, hotspot_crop, crop_rotation, crop_translation)

    all_reference_target = [residue for chain in context_chain_ids for residue in reference_context[chain]]
    all_predicted_target = [residue for chain in context_chain_ids for residue in predicted[chain]]
    all_reference_ca = helper._ca_coords(all_reference_target)
    all_predicted_ca = helper._ca_coords(all_predicted_target)
    context_after_fga = helper._rmsd(all_predicted_ca, all_reference_ca, fga_rotation, fga_translation)
    partner_after_fga: float | None = None
    if len(context_chain_ids) > 1:
        partner_reference = helper._ca_coords(
            [residue for chain in context_chain_ids[1:] for residue in reference_context[chain]]
        )
        partner_predicted = helper._ca_coords(
            [residue for chain in context_chain_ids[1:] for residue in predicted[chain]]
        )
        partner_after_fga = helper._rmsd(partner_predicted, partner_reference, fga_rotation, fga_translation)

    peptide_mobile, peptide_reference = helper._paired_atom_coords(
        predicted_peptide, reference_design_peptide, BACKBONE_ATOMS
    )
    peptide_rmsd = helper._rmsd(peptide_mobile, peptide_reference, crop_rotation, crop_translation)
    design_contacts, design_peptide_contacts, _, _ = helper._contact_sets(
        reference_design_target, reference_design_peptide, args.contact_cutoff
    )
    reference_site_contacts = design_contacts & set(site_crop)
    reference_hotspot_contacts = design_contacts & set(hotspot_crop)
    reference_hotspot_min = helper._minimum_to_target_indices(
        reference_design_target, reference_design_peptide, hotspot_crop
    )
    if design_peptide_contacts:
        interface_indices = sorted(value - 1 for value in design_peptide_contacts)
        interface_mobile, interface_reference = helper._paired_atom_coords(
            predicted_peptide,
            reference_design_peptide,
            BACKBONE_ATOMS,
            interface_indices,
        )
        interface_rmsd = helper._rmsd(interface_mobile, interface_reference, crop_rotation, crop_translation)
    else:
        interface_rmsd = float("nan")

    context_contacts, _, _, _ = helper._contact_sets(all_predicted_target, predicted_peptide, args.contact_cutoff)
    fga_contacts, _, _, _ = helper._contact_sets(predicted_fga, predicted_peptide, args.contact_cutoff)
    site_contacts = fga_contacts & set(site)
    hotspot_contacts = fga_contacts & set(hotspots)
    partner_contacts = {value for value in context_contacts if value > context_lengths[0]}
    offsite = context_contacts - set(site)
    hotspot_min = helper._minimum_to_target_indices(predicted_fga, predicted_peptide, hotspots)
    offsite_fraction = len(offsite) / max(1, len(context_contacts))
    same_site = bool(site_contacts and hotspot_contacts and hotspot_min <= args.hotspot_contact_distance)
    offsite_dominant = bool(not same_site or (len(site_contacts) / max(1, len(context_contacts))) < args.min_site_contact_fraction)
    terminal_distance = helper._terminal_cn_distance(predicted_peptide)
    macrocycle_status = helper._macrocycle_status(
        terminal_distance, args.macrocycle_pass_distance, args.macrocycle_warn_distance
    )
    severe_clashes = helper._severe_clash_count(all_predicted_target, predicted_peptide, args.severe_clash_distance)
    clash_status = "pass_no_severe_clash" if severe_clashes == 0 else "fail_severe_clash"
    peptide_rog = helper._radius_of_gyration(predicted_peptide)

    target_failures: list[str] = []
    if fga_rmsd > (args.max_c1_fga_rmsd if candidate["context_id"].startswith("C1_") else args.max_c3_fga_rmsd):
        target_failures.append("high_fga_global_CA_RMSD")
    if crop_rmsd > args.max_crop_local_rmsd:
        target_failures.append("high_stage0_crop_local_CA_RMSD")
    if site_after_crop > args.max_site2_local_rmsd:
        target_failures.append("high_Site_2_local_CA_RMSD")
    if hotspot_after_crop > args.max_hotspot_local_rmsd:
        target_failures.append("high_hotspot_local_CA_RMSD")
    if candidate["context_id"].startswith("C3_") and (
        partner_after_fga is None or partner_after_fga > args.max_partner_after_fga_rmsd
    ):
        target_failures.append("native_partner_assembly_not_recovered")
    target_pass = not target_failures
    if peptide_rmsd <= args.strong_pose_rmsd:
        pose_class = "strong_pose_recovery"
    elif peptide_rmsd <= args.moderate_pose_rmsd:
        pose_class = "moderate_pose_recovery"
    else:
        pose_class = "poor_pose_recovery"
    hard_geometry = macrocycle_status == "pass_head_to_tail_geometry" and severe_clashes == 0
    metric_identity = all(
        str(metric.get(field, "")) == str(expected)
        for field, expected in {
            "stage5B_v2_job_id": job["stage5B_v2_job_id"],
            "stage5B_v2_candidate_context_id": candidate["stage5B_v2_candidate_context_id"],
            "stage5B_v2_candidate_id": candidate["stage5B_v2_candidate_id"],
            "context_id": candidate["context_id"],
            "protocol_hash": candidate["protocol_hash"],
            "stage5_selection_mode": candidate["stage5_selection_mode"],
        }.items()
    )
    protocol_valid = protocol_valid and metric_identity
    if not protocol_valid:
        model_support = "stage5B_v2_protocol_failure"
    elif not target_pass:
        model_support = "stage5B_v2_target_context_not_recovered"
    elif hard_geometry and same_site and pose_class == "strong_pose_recovery":
        model_support = "stage5B_v2_strong_support"
    elif hard_geometry and (same_site or pose_class == "moderate_pose_recovery"):
        model_support = "stage5B_v2_partial_support"
    else:
        model_support = "stage5B_v2_peptide_not_recovered"

    prediction_coords = helper._ca_coords(all_predicted_target + predicted_peptide)
    coordinate_hash = _sha1_arrays(np.round(prediction_coords, 3))
    confidence_hash = helper._confidence_hash(prediction_npz)
    return {
        "stage5B_v2_candidate_context_id": candidate["stage5B_v2_candidate_context_id"],
        "stage5B_v2_candidate_id": candidate["stage5B_v2_candidate_id"],
        "stage5B_v2_job_id": job["stage5B_v2_job_id"],
        "context_id": candidate["context_id"],
        "protocol_hash": candidate["protocol_hash"],
        "stage5_campaign_id": candidate["stage5_campaign_id"],
        "stage5_selection_mode": candidate["stage5_selection_mode"],
        "selection_order": candidate["selection_order"],
        "batch": candidate["batch"],
        "backbone_id": candidate["backbone_id"],
        **{field: candidate.get(field, "") for field in STAGE4_IDENTITY_FIELDS},
        "peptide_sequence": candidate["peptide_sequence"],
        "seed": job["seed"],
        "model_name": metric["model_name"],
        "prediction_pdb": prediction_pdb,
        "prediction_npz": prediction_npz,
        "protocol_identity_valid": str(protocol_valid).lower(),
        "context_sequence_match": str(context_sequence_match).lower(),
        "peptide_sequence_match": str(peptide_sequence_match).lower(),
        "context_template_coverage": metadata.get("context_template_residues_covered", ""),
        "peptide_template_coverage": metadata.get("peptide_template_residues_covered", ""),
        "template_mode": metadata.get("template_mode", ""),
        "use_initial_guess": str(metadata.get("use_initial_guess", "")).lower(),
        "reference_design_loaded_by_prediction_runner": str(
            metadata.get("reference_design_loaded_by_prediction_runner", "")
        ).lower(),
        "fga_global_CA_RMSD_A": _rounded(fga_rmsd),
        "stage0_crop_local_CA_RMSD_A": _rounded(crop_rmsd),
        "Site_2_RMSD_after_fga_alignment_A": _rounded(site_after_fga),
        "hotspot_RMSD_after_fga_alignment_A": _rounded(hotspot_after_fga),
        "Site_2_RMSD_after_crop_alignment_A": _rounded(site_after_crop),
        "hotspot_RMSD_after_crop_alignment_A": _rounded(hotspot_after_crop),
        "context_after_fga_alignment_CA_RMSD_A": _rounded(context_after_fga),
        "partner_after_fga_alignment_CA_RMSD_A": _rounded(partner_after_fga),
        "target_context_recovery_pass": str(target_pass).lower(),
        "target_context_recovery_failure_reasons": ";".join(target_failures),
        "target_aligned_peptide_backbone_RMSD_A": _rounded(peptide_rmsd),
        "target_aligned_peptide_interface_RMSD_A": _rounded(interface_rmsd),
        "reference_design_site2_contact_count": len(reference_site_contacts),
        "reference_design_hotspot_contact_count": len(reference_hotspot_contacts),
        "reference_design_hotspot_min_distance_A": _rounded(reference_hotspot_min),
        "num_context_contacts": len(context_contacts),
        "num_fga_contacts": len(fga_contacts),
        "num_site2_contacts": len(site_contacts),
        "num_hotspot_contacts": len(hotspot_contacts),
        "num_partner_contacts": len(partner_contacts),
        "hotspot_min_distance_A": _rounded(hotspot_min),
        "off_site_contact_count": len(offsite),
        "off_site_contact_fraction": _rounded(offsite_fraction, 4),
        "same_target_site_flag": str(same_site).lower(),
        "off_site_binding_dominant_flag": str(offsite_dominant).lower(),
        "terminal_C_N_distance_A": _rounded(terminal_distance),
        "macrocycle_geometry_status": macrocycle_status,
        "severe_clash_count": severe_clashes,
        "clash_status": clash_status,
        "peptide_radius_of_gyration_A": _rounded(peptide_rog),
        "context_mean_pLDDT_100": metric.get("plddt_context_mean_100", ""),
        "fga_mean_pLDDT_100": metric.get("plddt_fga_mean_100", ""),
        "Site_2_mean_pLDDT_100": metric.get("plddt_site2_mean_100", ""),
        "hotspot_mean_pLDDT_100": metric.get("plddt_hotspot_mean_100", ""),
        "peptide_mean_pLDDT_100": metric.get("plddt_peptide_mean_100", ""),
        "interface_pae_context_mean_A": metric.get("interface_pae_context_mean_A", ""),
        "interface_pae_fga_mean_A": metric.get("interface_pae_fga_mean_A", ""),
        "interface_pae_site2_mean_A": metric.get("interface_pae_site2_mean_A", ""),
        "interface_pae_site2_median_A": metric.get("interface_pae_site2_median_A", ""),
        "interface_pae_hotspot_mean_A": metric.get("interface_pae_hotspot_mean_A", ""),
        "pTM": metric.get("ptm", ""),
        "ipTM": metric.get("iptm", ""),
        "ranking_confidence": metric.get("ranking_confidence", ""),
        "pose_recovery_class": pose_class,
        "model_support_status": model_support,
        "prediction_coordinate_hash": coordinate_hash,
        "prediction_confidence_hash": confidence_hash,
        "requested_recycles": metadata.get("requested_recycles", ""),
        "forward_passes": metadata.get("forward_passes", ""),
        "loaded_colabdesign_commit": metadata.get("loaded_colabdesign_commit", ""),
        "cyclic_topology_encoding": metadata.get("cyclic_topology_encoding", ""),
        "notes": (
            "C1 is a crop-constrained protocol sanity test; C3 is the native-context test. "
            "Stage 4 coordinates were used only after prediction for alignment and recovery metrics."
        ),
        **{field: candidate.get(field, "") for field in SOURCE_ROUTE_PROVENANCE_FIELDS},
        **{field: candidate.get(field, "") for field in ROUTE_PROVENANCE_FIELDS},
    }


def _best_model(rows: Sequence[Mapping[str, Any]]) -> Mapping[str, Any] | None:
    if not rows:
        return None
    return min(
        rows,
        key=lambda row: (
            row.get("model_support_status") != "stage5B_v2_strong_support",
            row.get("model_support_status") != "stage5B_v2_partial_support",
            _float(row.get("target_aligned_peptide_backbone_RMSD_A"), 999.0),
            _float(row.get("hotspot_min_distance_A"), 999.0),
            -_float(row.get("peptide_mean_pLDDT_100"), -1.0),
        ),
    )


def _context_summaries(
    manifests: Sequence[Mapping[str, str]], models: Sequence[Mapping[str, Any]], expected_models: int
) -> list[dict[str, Any]]:
    grouped: dict[str, list[Mapping[str, Any]]] = defaultdict(list)
    for row in models:
        grouped[str(row["stage5B_v2_candidate_context_id"])].append(row)
    summaries = []
    for manifest in manifests:
        rows = grouped.get(manifest["stage5B_v2_candidate_context_id"], [])
        best_model = _best_model(rows)
        unique = len({(row["prediction_coordinate_hash"], row["prediction_confidence_hash"]) for row in rows})
        strong = sum(row["model_support_status"] == "stage5B_v2_strong_support" for row in rows)
        partial = sum(row["model_support_status"] == "stage5B_v2_partial_support" for row in rows)
        target_count = sum(_truthy(row["target_context_recovery_pass"]) for row in rows)
        if len(rows) != expected_models:
            support_class = "stage5B_v2_incomplete"
            reason = f"completed {len(rows)}/{expected_models} models"
        elif strong:
            support_class = "stage5B_v2_context_strong_support"
            reason = f"{strong}/{len(rows)} models met target, same-site, geometry, and <=3 A pose criteria"
        elif partial:
            support_class = "stage5B_v2_context_partial_support"
            reason = f"{partial}/{len(rows)} models gave partial context-conditioned recovery"
        elif not target_count:
            support_class = "stage5B_v2_context_not_evaluable"
            reason = "target context did not meet internal recovery thresholds"
        else:
            support_class = "stage5B_v2_context_not_recovered"
            reason = "target context recovered in at least one model, but peptide pose/site did not"
        summaries.append(
            {
                "stage5B_v2_candidate_context_id": manifest["stage5B_v2_candidate_context_id"],
                "stage5B_v2_candidate_id": manifest["stage5B_v2_candidate_id"],
                "context_id": manifest["context_id"],
                "stage5_selection_mode": manifest["stage5_selection_mode"],
                "selection_order": manifest["selection_order"],
                "batch": manifest["batch"],
                "backbone_id": manifest["backbone_id"],
                "global_backbone_id": manifest["global_backbone_id"],
                "backbone_family_id": manifest["backbone_family_id"],
                "peptide_sequence": manifest["peptide_sequence"],
                "models_expected": expected_models,
                "models_completed": len(rows),
                "effective_unique_prediction_count": unique,
                "duplicate_prediction_count": len(rows) - unique,
                "target_context_recovery_count": target_count,
                "same_target_site_count": sum(_truthy(row["same_target_site_flag"]) for row in rows),
                "strong_pose_recovery_count": sum(row["pose_recovery_class"] == "strong_pose_recovery" for row in rows),
                "moderate_pose_recovery_count": sum(row["pose_recovery_class"] == "moderate_pose_recovery" for row in rows),
                "macrocycle_pass_count": sum(row["macrocycle_geometry_status"] == "pass_head_to_tail_geometry" for row in rows),
                "no_severe_clash_count": sum(row["clash_status"] == "pass_no_severe_clash" for row in rows),
                "median_fga_global_CA_RMSD_A": _median(rows, "fga_global_CA_RMSD_A"),
                "best_fga_global_CA_RMSD_A": _best(rows, "fga_global_CA_RMSD_A"),
                "median_stage0_crop_local_CA_RMSD_A": _median(rows, "stage0_crop_local_CA_RMSD_A"),
                "best_stage0_crop_local_CA_RMSD_A": _best(rows, "stage0_crop_local_CA_RMSD_A"),
                "median_peptide_RMSD_A": _median(rows, "target_aligned_peptide_backbone_RMSD_A"),
                "best_peptide_RMSD_A": _best(rows, "target_aligned_peptide_backbone_RMSD_A"),
                "median_hotspot_distance_A": _median(rows, "hotspot_min_distance_A"),
                "best_hotspot_distance_A": _best(rows, "hotspot_min_distance_A"),
                "median_peptide_pLDDT_100": _median(rows, "peptide_mean_pLDDT_100"),
                "median_hotspot_PAE_A": _median(rows, "interface_pae_hotspot_mean_A"),
                "best_model_path": best_model["prediction_pdb"] if best_model else "",
                "context_support_class": support_class,
                "support_reason": reason,
                **{field: manifest.get(field, "") for field in SOURCE_ROUTE_PROVENANCE_FIELDS},
                **{field: manifest.get(field, "") for field in ROUTE_PROVENANCE_FIELDS},
            }
        )
    return summaries


def _candidate_summaries(
    manifests: Sequence[Mapping[str, str]], context_rows: Sequence[Mapping[str, Any]]
) -> list[dict[str, Any]]:
    manifest_by_candidate: dict[str, Mapping[str, str]] = {}
    grouped: dict[str, dict[str, Mapping[str, Any]]] = defaultdict(dict)
    for manifest in manifests:
        manifest_by_candidate.setdefault(manifest["stage5B_v2_candidate_id"], manifest)
    for row in context_rows:
        grouped[str(row["stage5B_v2_candidate_id"])][str(row["context_id"])] = row
    output = []
    c1_id, c3_id = "C1_crop86_full_template", "C3_native_GHI301_full_template"
    for candidate_id, manifest in sorted(
        manifest_by_candidate.items(), key=lambda item: int(item[1]["selection_order"])
    ):
        c1 = grouped[candidate_id][c1_id]
        c3 = grouped[candidate_id][c3_id]
        if c3["context_support_class"] == "stage5B_v2_context_strong_support":
            support = "stage5B_v2_native_context_strong_support"
            interpretation = "At least one C3 model recovered the target context, Site_2/hotspot, macrocycle, and strong peptide pose."
        elif c3["context_support_class"] == "stage5B_v2_context_partial_support":
            support = "stage5B_v2_native_context_partial_support"
            interpretation = "C3 provides partial native-context support; further seeds or an orthogonal method are required."
        elif c3["target_context_recovery_count"] == 0:
            support = "stage5B_v2_native_context_not_evaluable"
            interpretation = "C3 target context was not recovered under the internal thresholds, so peptide rejection is not justified."
        elif c1["context_support_class"] in {
            "stage5B_v2_context_strong_support",
            "stage5B_v2_context_partial_support",
        }:
            support = "stage5B_v2_crop_only_support"
            interpretation = "Support was limited to the spatially restricted C1 crop and did not transfer to C3."
        else:
            support = "stage5B_v2_not_recovered"
            interpretation = "The target was evaluable, but no qualifying peptide pose/site recovery was observed."
        output.append(
            {
                "stage5B_v2_candidate_id": candidate_id,
                "stage5_selection_mode": manifest["stage5_selection_mode"],
                "selection_order": manifest["selection_order"],
                "batch": manifest["batch"],
                "backbone_id": manifest["backbone_id"],
                "global_backbone_id": manifest["global_backbone_id"],
                "backbone_family_id": manifest["backbone_family_id"],
                "peptide_sequence": manifest["peptide_sequence"],
                "C1_models_completed": c1["models_completed"],
                "C1_target_recovery_count": c1["target_context_recovery_count"],
                "C1_same_site_count": c1["same_target_site_count"],
                "C1_strong_pose_count": c1["strong_pose_recovery_count"],
                "C1_moderate_pose_count": c1["moderate_pose_recovery_count"],
                "C1_best_peptide_RMSD_A": c1["best_peptide_RMSD_A"],
                "C1_best_hotspot_distance_A": c1["best_hotspot_distance_A"],
                "C1_context_support_class": c1["context_support_class"],
                "C3_models_completed": c3["models_completed"],
                "C3_target_recovery_count": c3["target_context_recovery_count"],
                "C3_same_site_count": c3["same_target_site_count"],
                "C3_strong_pose_count": c3["strong_pose_recovery_count"],
                "C3_moderate_pose_count": c3["moderate_pose_recovery_count"],
                "C3_best_peptide_RMSD_A": c3["best_peptide_RMSD_A"],
                "C3_best_hotspot_distance_A": c3["best_hotspot_distance_A"],
                "C3_context_support_class": c3["context_support_class"],
                "cross_context_support_class": support,
                "cross_context_interpretation": interpretation,
                **{field: manifest.get(field, "") for field in SOURCE_ROUTE_PROVENANCE_FIELDS},
                **{field: manifest.get(field, "") for field in ROUTE_PROVENANCE_FIELDS},
            }
        )
    return output


def _report(
    *,
    models: Sequence[Mapping[str, Any]],
    context_rows: Sequence[Mapping[str, Any]],
    candidates: Sequence[Mapping[str, Any]],
    selection_mode: str,
    args: argparse.Namespace,
) -> str:
    candidate_columns = [
        "selection_order",
        "peptide_sequence",
        "C1_target_recovery_count",
        "C1_same_site_count",
        "C1_best_peptide_RMSD_A",
        "C3_target_recovery_count",
        "C3_same_site_count",
        "C3_best_peptide_RMSD_A",
        "cross_context_support_class",
        "cross_context_interpretation",
    ]
    context_columns = [
        "selection_order",
        "context_id",
        "target_context_recovery_count",
        "same_target_site_count",
        "strong_pose_recovery_count",
        "moderate_pose_recovery_count",
        "best_fga_global_CA_RMSD_A",
        "best_stage0_crop_local_CA_RMSD_A",
        "best_peptide_RMSD_A",
        "best_hotspot_distance_A",
        "context_support_class",
    ]
    support_counts = Counter(str(row["cross_context_support_class"]) for row in candidates)
    support_order = {
        "stage5B_v2_native_context_strong_support": 0,
        "stage5B_v2_native_context_partial_support": 1,
        "stage5B_v2_crop_only_support": 2,
        "stage5B_v2_not_recovered": 3,
        "stage5B_v2_native_context_not_evaluable": 4,
    }
    displayed_candidates = sorted(
        candidates,
        key=lambda row: (
            support_order.get(str(row["cross_context_support_class"]), 99),
            int(row["selection_order"]),
        ),
    )[: args.top_report]
    displayed_ids = {str(row["stage5B_v2_candidate_id"]) for row in displayed_candidates}
    displayed_contexts = [
        row for row in context_rows if str(row["stage5B_v2_candidate_id"]) in displayed_ids
    ]
    support_count_lines = "\n".join(
        f"{key}: {support_counts[key]}" for key in sorted(support_counts, key=lambda key: support_order.get(key, 99))
    ) or "none: 0"
    return f"""# Stage 5B-v2 C1/C3 Validation Report

This report evaluates corrected-route Stage 4 candidates under two full
target-template contexts. C1 is a restricted 86-aa protocol sanity test. C3 is
the native G/H/I context and carries more biological weight for Site_2 recovery.

Prediction inputs retained target sequence, backbone, sidechains, and
target-target interchain template information. Peptide template coverage was
zero, no peptide design coordinates or initial guess were loaded, and only the
peptide chain received the cyclic positional offset. Stage 4 complex coordinates
were used only after prediction for alignment and pose/contact recovery metrics.

## Completion

```text
models_completed: {len(models)}
candidate_context_pairs: {len(context_rows)}
candidates: {len(candidates)}
selection_mode: {selection_mode}
fga_global_CA_RMSD_A min/median/max: {_distribution(models, 'fga_global_CA_RMSD_A')}
stage0_crop_local_CA_RMSD_A min/median/max: {_distribution(models, 'stage0_crop_local_CA_RMSD_A')}
peptide_backbone_RMSD_A min/median/max: {_distribution(models, 'target_aligned_peptide_backbone_RMSD_A')}
hotspot_min_distance_A min/median/max: {_distribution(models, 'hotspot_min_distance_A')}
```

Internal diagnostic thresholds are not experimental standards:

```text
C1_FGA_RMSD_A: <= {args.max_c1_fga_rmsd:g}
C3_FGA_RMSD_A: <= {args.max_c3_fga_rmsd:g}
crop_local_RMSD_A: <= {args.max_crop_local_rmsd:g}
Site_2_local_RMSD_A: <= {args.max_site2_local_rmsd:g}
hotspot_local_RMSD_A: <= {args.max_hotspot_local_rmsd:g}
C3_partner_after_FGA_alignment_RMSD_A: <= {args.max_partner_after_fga_rmsd:g}
same_site_hotspot_distance_A: <= {args.hotspot_contact_distance:g}
strong_peptide_pose_RMSD_A: <= {args.strong_pose_rmsd:g}
moderate_peptide_pose_RMSD_A: <= {args.moderate_pose_rmsd:g}
```

## Cross-context Candidate Interpretation

```text
{support_count_lines}
```

The table below is limited to the first {args.top_report} rows after sorting by
support class and Stage 4 selection order. The CSV contains every candidate.

{rows_to_markdown(displayed_candidates, candidate_columns, "No candidate results were collected.")}

## Candidate-context Results

{rows_to_markdown(displayed_contexts, context_columns, "No context results were collected.")}

These are computational recovery tests, not final peptide candidates and not
experimental affinity measurements. A C1-only recovery can be caused by the
restricted crop. A C3 failure is not a peptide failure when the C3 target
context itself does not satisfy the recovery checks.
"""


def main() -> int:
    parser = argparse.ArgumentParser(description="Collect Stage 5B-v2 C1/C3 target-context-conditioned predictions.")
    parser.add_argument("--stage5b-v2-root", required=True)
    parser.add_argument("--candidate-context-manifest-csv", default="")
    parser.add_argument("--jobs-csv", default="")
    parser.add_argument("--project-config", required=True)
    parser.add_argument("--validate-inputs-only", action="store_true")
    parser.add_argument("--contact-cutoff", type=float, default=5.0)
    parser.add_argument("--hotspot-contact-distance", type=float, default=5.0)
    parser.add_argument("--min-site-contact-fraction", type=float, default=0.20)
    parser.add_argument("--severe-clash-distance", type=float, default=1.2)
    parser.add_argument("--macrocycle-pass-distance", type=float, default=2.0)
    parser.add_argument("--macrocycle-warn-distance", type=float, default=3.0)
    parser.add_argument("--max-c1-fga-rmsd", type=float, default=3.0)
    parser.add_argument("--max-c3-fga-rmsd", type=float, default=4.0)
    parser.add_argument("--max-crop-local-rmsd", type=float, default=2.0)
    parser.add_argument("--max-site2-local-rmsd", type=float, default=2.0)
    parser.add_argument("--max-hotspot-local-rmsd", type=float, default=2.0)
    parser.add_argument("--max-partner-after-fga-rmsd", type=float, default=4.0)
    parser.add_argument("--strong-pose-rmsd", type=float, default=3.0)
    parser.add_argument("--moderate-pose-rmsd", type=float, default=5.0)
    parser.add_argument("--top-report", type=int, default=100)
    args = parser.parse_args()

    logger = setup_logger("35_collect_stage5b_v2_context_validation")
    append_run_header(logger, "35_collect_stage5b_v2_context_validation.py")
    if not (0.0 <= args.min_site_contact_fraction <= 1.0):
        raise RuntimeError("--min-site-contact-fraction must be between 0 and 1")
    if args.moderate_pose_rmsd < args.strong_pose_rmsd:
        raise RuntimeError("--moderate-pose-rmsd must be >= --strong-pose-rmsd")
    if args.top_report < 1:
        raise RuntimeError("--top-report must be positive")
    stage5_root = assert_active_route_path(
        _resolve_mixed_path(args.stage5b_v2_root), "Stage 35 Stage 5B-v2 root"
    )
    manifest_csv = assert_active_route_path(
        _resolve_mixed_path(args.candidate_context_manifest_csv)
        if args.candidate_context_manifest_csv
        else stage5_root / "FGA_rfpeptides_stage5B_v2_candidate_context_manifest.csv",
        "Stage 35 candidate-context manifest",
    )
    jobs_csv = assert_active_route_path(
        _resolve_mixed_path(args.jobs_csv)
        if args.jobs_csv
        else stage5_root / "FGA_rfpeptides_stage5B_v2_prediction_jobs.csv",
        "Stage 35 jobs CSV",
    )
    manifests, jobs, _, campaign = _validate_manifest_and_jobs(
        stage5_root=stage5_root,
        manifest_csv=manifest_csv,
        jobs_csv=jobs_csv,
        project_config=args.project_config,
    )
    if args.validate_inputs_only:
        logger.info("Stage 5B-v2 candidate-context rows validated: %s", len(manifests))
        logger.info("Stage 5B-v2 seed jobs validated: %s", len(jobs))
        logger.info("No prediction output was required or parsed.")
        return 0

    helper = _load_stage5b_v1_helpers()
    manifest_lookup = {row["stage5B_v2_candidate_context_id"]: row for row in manifests}
    model_rows: list[dict[str, Any]] = []
    expected_models_per_pair: dict[str, int] = defaultdict(int)
    for job in jobs:
        candidate = manifest_lookup[job["stage5B_v2_candidate_context_id"]]
        output_dir = assert_active_route_path(job["prediction_output_dir"], f"Stage 35 output {job['stage5B_v2_job_id']}")
        metadata_path = output_dir / "run_metadata.json"
        metrics_path = output_dir / "model_metrics.csv"
        if not metadata_path.is_file() or not metrics_path.is_file():
            raise RuntimeError(f"Incomplete Stage 5B-v2 output: {output_dir}")
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        with metrics_path.open("r", encoding="utf-8", newline="") as handle:
            metrics = list(csv.DictReader(handle))
        if len(metrics) != int(job["models_per_seed"]):
            raise RuntimeError(f"Job {job['stage5B_v2_job_id']} has {len(metrics)} model rows")
        expected_models_per_pair[candidate["stage5B_v2_candidate_context_id"]] += int(job["models_per_seed"])
        for metric in metrics:
            model_rows.append(
                _model_row(
                    helper=helper,
                    candidate=candidate,
                    job=job,
                    metadata=metadata,
                    metric=metric,
                    args=args,
                )
            )

    expected_set = set(manifest_lookup)
    if {row["stage5B_v2_candidate_context_id"] for row in model_rows} != expected_set:
        raise RuntimeError("Collected model rows do not cover every candidate-context pair")
    if len(set(expected_models_per_pair.values())) != 1:
        raise RuntimeError("Candidate-context pairs have inconsistent expected model counts")
    expected_models = next(iter(expected_models_per_pair.values()))
    context_rows = _context_summaries(manifests, model_rows, expected_models)
    candidate_rows = _candidate_summaries(manifests, context_rows)
    write_csv(stage5_root / "FGA_rfpeptides_stage5B_v2_model_results.csv", model_rows, MODEL_FIELDS)
    write_csv(
        stage5_root / "FGA_rfpeptides_stage5B_v2_candidate_context_summary.csv",
        context_rows,
        CONTEXT_SUMMARY_FIELDS,
    )
    write_csv(
        stage5_root / "FGA_rfpeptides_stage5B_v2_candidate_summary.csv",
        candidate_rows,
        CANDIDATE_SUMMARY_FIELDS,
    )
    write_markdown(
        stage5_root / "FGA_rfpeptides_stage5B_v2_validation_report.md",
        _report(
            models=model_rows,
            context_rows=context_rows,
            candidates=candidate_rows,
            selection_mode=str(campaign["selection_mode"]),
            args=args,
        ),
    )
    logger.info("Stage 5B-v2 model predictions parsed: %s", len(model_rows))
    logger.info("Stage 5B-v2 candidate-context summaries: %s", len(context_rows))
    logger.info("Stage 5B-v2 cross-context candidate summaries: %s", len(candidate_rows))
    logger.info("Output directory: %s", stage5_root)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
