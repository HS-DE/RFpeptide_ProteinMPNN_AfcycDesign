from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np

from common import (
    ROUTE_PROVENANCE_FIELDS,
    SOURCE_ROUTE_PROVENANCE_FIELDS,
    add_route_provenance,
    assert_active_route_path,
    append_run_header,
    load_route_manifest,
    read_csv,
    resolve_path,
    route_provenance_fields,
    rows_to_markdown,
    setup_logger,
    validate_route_project_config,
    validate_row_route_provenance,
    validate_source_route_provenance,
    write_csv,
    write_markdown,
)
from pdb_utils import parse_residues, residue_sequence


BACKBONE_ATOMS = ("N", "CA", "C")

MODEL_FIELDS = [
    "stage5B_candidate_id",
    "stage5B_job_id",
    "protocol_hash",
    "batch",
    "backbone_id",
    "peptide_sequence",
    "seed",
    "model_name",
    "prediction_pdb",
    "prediction_npz",
    "target_template_path",
    "target_template_sha1",
    "target_template_coverage",
    "target_template_coverage_fraction",
    "peptide_template_coverage",
    "template_mode",
    "template_input_verified",
    "template_effect_observed",
    "use_initial_guess",
    "target_sequence_match",
    "peptide_sequence_match",
    "reference_target_to_template_CA_RMSD_A",
    "target_global_CA_RMSD_A",
    "Site_2_local_CA_RMSD_A",
    "hotspot_local_CA_RMSD_A",
    "target_template_recovery_pass",
    "target_template_recovery_failure_reasons",
    "target_aligned_peptide_backbone_RMSD_A",
    "target_aligned_peptide_interface_RMSD_A",
    "design_target_contact_count",
    "reference_design_site2_contact_count",
    "reference_design_hotspot_contact_count",
    "reference_design_hotspot_min_distance_A",
    "reference_design_same_target_site_flag",
    "target_contact_recovery_count",
    "num_target_contacts",
    "num_site2_contacts",
    "num_hotspot_contacts",
    "hotspot_min_distance_A",
    "off_site_contact_count",
    "off_site_contact_fraction",
    "same_target_site_flag",
    "terminal_C_N_distance_A",
    "macrocycle_geometry_status",
    "severe_clash_count",
    "clash_status",
    "peptide_radius_of_gyration_A",
    "target_mean_pLDDT_100",
    "Site_2_mean_pLDDT_100",
    "hotspot_mean_pLDDT_100",
    "peptide_mean_pLDDT_100",
    "interface_pae_global_mean_A",
    "interface_pae_site2_mean_A",
    "interface_pae_site2_median_A",
    "interface_pae_hotspot_mean_A",
    "pTM",
    "ipTM",
    "ranking_confidence",
    "ranking_confidence_source",
    "confidence_support_status",
    "pose_recovery_class",
    "model_support_status",
    "prediction_coordinate_hash",
    "prediction_confidence_hash",
    "prediction_identity_hash",
    "validation_test_type",
    "target_msa_mode",
    "peptide_msa_mode",
    "use_mlm",
    "use_dropout",
    "requested_recycles",
    "forward_passes",
    "loaded_colabdesign_commit",
    "cyclic_topology_encoding",
    "protocol_identity_valid",
    "notes",
] + ROUTE_PROVENANCE_FIELDS + SOURCE_ROUTE_PROVENANCE_FIELDS

CANDIDATE_FIELDS = [
    "stage5B_candidate_id",
    "batch",
    "backbone_id",
    "peptide_sequence",
    "models_expected",
    "models_completed",
    "protocol_valid_model_count",
    "effective_unique_prediction_count",
    "duplicate_prediction_count",
    "seed_diversity_status",
    "reference_design_site2_contact_count",
    "reference_design_hotspot_contact_count",
    "reference_design_hotspot_min_distance_A",
    "reference_design_same_target_site_flag",
    "target_template_recovery_count",
    "Site_2_recovery_count",
    "hotspot_recovery_count",
    "same_target_site_count",
    "strong_pose_recovery_count",
    "moderate_pose_recovery_count",
    "macrocycle_pass_count",
    "confidence_pass_count",
    "median_peptide_RMSD_A",
    "best_peptide_RMSD_A",
    "median_Site_2_PAE_A",
    "best_Site_2_PAE_A",
    "median_hotspot_distance_A",
    "best_hotspot_distance_A",
    "median_peptide_pLDDT_100",
    "best_model_path",
    "best_model_name",
    "best_seed",
    "stage5B_support_class",
    "support_reason",
] + ROUTE_PROVENANCE_FIELDS + SOURCE_ROUTE_PROVENANCE_FIELDS


def _resolve_mixed_path(value: str | Path) -> Path:
    text = str(value).strip().replace("\\", "/")
    if os.name == "nt" and text.startswith("/mnt/") and len(text) > 7 and text[6] == "/":
        return Path(f"{text[5].upper()}:/{text[7:]}")
    if os.name != "nt" and len(text) >= 3 and text[1] == ":" and text[2] == "/":
        return Path(f"/mnt/{text[0].lower()}{text[2:]}")
    path = Path(text)
    return path if path.is_absolute() else resolve_path(path)


def _float(value: Any, default: float = float("nan")) -> float:
    try:
        parsed = float(str(value).strip())
    except (TypeError, ValueError):
        return default
    return parsed if math.isfinite(parsed) else default


def _rounded(value: float | None, digits: int = 3) -> float | str:
    return round(value, digits) if value is not None and math.isfinite(value) else ""


def _median(rows: Sequence[Mapping[str, Any]], field: str) -> float | str:
    values = [_float(row.get(field, "")) for row in rows]
    finite = [value for value in values if math.isfinite(value)]
    return round(float(np.median(finite)), 3) if finite else ""


def _best_numeric(rows: Sequence[Mapping[str, Any]], field: str, higher: bool = False) -> float | str:
    values = [_float(row.get(field, "")) for row in rows]
    finite = [value for value in values if math.isfinite(value)]
    if not finite:
        return ""
    return round(max(finite) if higher else min(finite), 3)


def _distribution_text(rows: Sequence[Mapping[str, Any]], field: str) -> str:
    values = [_float(row.get(field, "")) for row in rows]
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


def _ca_coords(residues: Sequence[Mapping[str, Any]]) -> np.ndarray:
    points = []
    for residue in residues:
        atoms = residue.get("atoms", {})
        if "CA" not in atoms:
            raise RuntimeError("A residue lacks a CA atom")
        points.append(atoms["CA"])
    return np.asarray(points, dtype=float)


def _paired_atom_coords(
    mobile: Sequence[Mapping[str, Any]],
    reference: Sequence[Mapping[str, Any]],
    atom_names: Sequence[str],
    residue_indices_0based: Iterable[int] | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    if len(mobile) != len(reference):
        raise RuntimeError(f"Chain length mismatch: {len(mobile)} != {len(reference)}")
    indices = range(len(mobile)) if residue_indices_0based is None else residue_indices_0based
    mobile_points = []
    reference_points = []
    for index in indices:
        mobile_atoms = mobile[index].get("atoms", {})
        reference_atoms = reference[index].get("atoms", {})
        for atom_name in atom_names:
            if atom_name in mobile_atoms and atom_name in reference_atoms:
                mobile_points.append(mobile_atoms[atom_name])
                reference_points.append(reference_atoms[atom_name])
    if len(mobile_points) < 3:
        raise RuntimeError("Too few paired atoms for RMSD")
    return np.asarray(mobile_points, dtype=float), np.asarray(reference_points, dtype=float)


def _kabsch(mobile: np.ndarray, reference: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    if mobile.shape != reference.shape or mobile.shape[0] < 3:
        raise RuntimeError(f"Cannot align coordinate arrays with shapes {mobile.shape} and {reference.shape}")
    mobile_center = mobile.mean(axis=0)
    reference_center = reference.mean(axis=0)
    covariance = (mobile - mobile_center).T @ (reference - reference_center)
    u_matrix, _, vt_matrix = np.linalg.svd(covariance)
    rotation = vt_matrix.T @ u_matrix.T
    if np.linalg.det(rotation) < 0:
        vt_matrix[-1, :] *= -1
        rotation = vt_matrix.T @ u_matrix.T
    translation = reference_center - mobile_center @ rotation.T
    aligned = mobile @ rotation.T + translation
    rmsd = float(np.sqrt(np.mean(np.sum((aligned - reference) ** 2, axis=1))))
    return rotation, translation, rmsd


def _transform(coords: np.ndarray, rotation: np.ndarray, translation: np.ndarray) -> np.ndarray:
    return coords @ rotation.T + translation


def _rmsd(mobile: np.ndarray, reference: np.ndarray, rotation: np.ndarray, translation: np.ndarray) -> float:
    aligned = _transform(mobile, rotation, translation)
    return float(np.sqrt(np.mean(np.sum((aligned - reference) ** 2, axis=1))))


def _subset_ca_rmsd(
    mobile_ca: np.ndarray,
    reference_ca: np.ndarray,
    indices_1based: Sequence[int],
    rotation: np.ndarray,
    translation: np.ndarray,
) -> float:
    indices = np.asarray([int(index) - 1 for index in indices_1based], dtype=int)
    delta = _transform(mobile_ca[indices], rotation, translation) - reference_ca[indices]
    return float(np.sqrt(np.mean(np.sum(delta * delta, axis=1))))


def _residue_atoms(residue: Mapping[str, Any]) -> np.ndarray:
    return np.asarray(list(residue.get("atoms", {}).values()), dtype=float)


def _contact_sets(
    target: Sequence[Mapping[str, Any]],
    peptide: Sequence[Mapping[str, Any]],
    cutoff: float,
) -> tuple[set[int], set[int], set[tuple[int, int]], float]:
    target_contacts: set[int] = set()
    peptide_contacts: set[int] = set()
    pairs: set[tuple[int, int]] = set()
    minimum = float("inf")
    cutoff_sq = cutoff * cutoff
    for target_index, target_residue in enumerate(target, start=1):
        target_atoms = _residue_atoms(target_residue)
        if not target_atoms.size:
            continue
        for peptide_index, peptide_residue in enumerate(peptide, start=1):
            peptide_atoms = _residue_atoms(peptide_residue)
            if not peptide_atoms.size:
                continue
            squared = np.sum((target_atoms[:, None, :] - peptide_atoms[None, :, :]) ** 2, axis=2)
            pair_min_sq = float(squared.min())
            minimum = min(minimum, math.sqrt(pair_min_sq))
            if pair_min_sq <= cutoff_sq:
                target_contacts.add(target_index)
                peptide_contacts.add(peptide_index)
                pairs.add((target_index, peptide_index))
    return target_contacts, peptide_contacts, pairs, minimum


def _minimum_to_target_indices(
    target: Sequence[Mapping[str, Any]],
    peptide: Sequence[Mapping[str, Any]],
    indices_1based: Sequence[int],
) -> float:
    peptide_atoms = np.concatenate([_residue_atoms(residue) for residue in peptide if _residue_atoms(residue).size], axis=0)
    values = []
    for index in indices_1based:
        atoms = _residue_atoms(target[int(index) - 1])
        if atoms.size:
            values.append(float(np.sqrt(np.sum((atoms[:, None, :] - peptide_atoms[None, :, :]) ** 2, axis=2)).min()))
    return min(values, default=float("inf"))


def _severe_clash_count(
    target: Sequence[Mapping[str, Any]],
    peptide: Sequence[Mapping[str, Any]],
    cutoff: float,
) -> int:
    target_atoms = np.concatenate([_residue_atoms(residue) for residue in target if _residue_atoms(residue).size], axis=0)
    peptide_atoms = np.concatenate([_residue_atoms(residue) for residue in peptide if _residue_atoms(residue).size], axis=0)
    squared = np.sum((target_atoms[:, None, :] - peptide_atoms[None, :, :]) ** 2, axis=2)
    return int(np.count_nonzero(squared < cutoff * cutoff))


def _terminal_cn_distance(peptide: Sequence[Mapping[str, Any]]) -> float:
    if not peptide:
        return float("nan")
    first = peptide[0].get("atoms", {})
    last = peptide[-1].get("atoms", {})
    if "N" not in first or "C" not in last:
        return float("nan")
    return float(np.linalg.norm(np.asarray(last["C"]) - np.asarray(first["N"])))


def _radius_of_gyration(peptide: Sequence[Mapping[str, Any]]) -> float:
    coords = _ca_coords(peptide)
    center = coords.mean(axis=0)
    return float(np.sqrt(np.mean(np.sum((coords - center) ** 2, axis=1))))


def _macrocycle_status(distance: float, pass_distance: float, warn_distance: float) -> str:
    if not math.isfinite(distance):
        return "not_evaluable"
    if distance <= pass_distance:
        return "pass_head_to_tail_geometry"
    if distance <= warn_distance:
        return "warn_head_to_tail_geometry"
    return "fail_open_chain_geometry"


def _mapping_indices(path: Path, expected_length: int) -> tuple[list[int], list[int]]:
    rows = read_csv(path)
    if len(rows) != expected_length:
        raise RuntimeError(f"Stage 0 mapping has {len(rows)} rows, expected {expected_length}")
    site = [int(row["rfpeptides_residue_number"]) for row in rows if str(row.get("is_target_site_residue", "")).lower() == "true"]
    hotspots = [int(row["rfpeptides_residue_number"]) for row in rows if str(row.get("is_selected_hotspot", "")).lower() == "true"]
    if not site or not hotspots:
        raise RuntimeError("Stage 0 mapping lacks Site_2 or hotspot indices")
    return site, hotspots


def _protocol_identity_valid(
    job: Mapping[str, str],
    manifest: Mapping[str, str],
    metadata: Mapping[str, Any],
    metric: Mapping[str, str],
) -> bool:
    checks = [
        metadata.get("stage5B_job_id") == job.get("stage5B_job_id"),
        metadata.get("stage5B_candidate_id") == manifest.get("stage5B_candidate_id"),
        metadata.get("protocol_hash") == manifest.get("protocol_hash"),
        metadata.get("target_template_sha1") == manifest.get("target_template_sha1"),
        metadata.get("peptide_sequence") == manifest.get("peptide_sequence"),
        metadata.get("template_mode") == "target_only",
        metadata.get("use_initial_guess") is False,
        metadata.get("reference_design_loaded_by_prediction_runner") is False,
        metadata.get("template_input_verified") is True,
        int(metadata.get("target_template_residues_covered", 0)) == 86,
        int(metadata.get("peptide_template_residues_covered", -1)) == 0,
        metric.get("stage5B_job_id") == job.get("stage5B_job_id"),
        metric.get("target_template_sha1") == manifest.get("target_template_sha1"),
    ]
    return all(checks)


def _confidence_hash(npz_path: Path) -> str:
    with np.load(npz_path, allow_pickle=False) as data:
        plddt = np.round(np.asarray(data["plddt_fraction"], dtype=np.float32), 4)
        pae = np.round(np.asarray(data["pae"], dtype=np.float32), 2)
    return _sha1_arrays(plddt, pae)


def _model_row(
    manifest: Mapping[str, str],
    job: Mapping[str, str],
    metric: Mapping[str, str],
    metadata: Mapping[str, Any],
    prediction_pdb: Path,
    prediction_npz: Path,
    reference_chains: Mapping[str, Sequence[Mapping[str, Any]]],
    template_target: Sequence[Mapping[str, Any]],
    site_indices: Sequence[int],
    hotspot_indices: Sequence[int],
    args: argparse.Namespace,
) -> dict[str, Any]:
    predicted_chains = parse_residues(prediction_pdb)
    if "A" not in predicted_chains or "B" not in predicted_chains:
        raise RuntimeError(f"Predicted complex lacks chain A or B: {prediction_pdb}")
    predicted_target = predicted_chains["A"]
    predicted_peptide = predicted_chains["B"]
    reference_target = list(reference_chains["A"])
    reference_peptide = list(reference_chains["B"])
    if len(predicted_target) != len(template_target) or len(predicted_peptide) != len(reference_peptide):
        raise RuntimeError("Predicted chain lengths do not match target template/reference peptide")
    if len(reference_target) != len(template_target):
        raise RuntimeError("Stage 4 reference target length does not match the Stage 0 target template")

    target_sequence_match = residue_sequence(predicted_target) == residue_sequence(template_target)
    peptide_sequence_match = residue_sequence(predicted_peptide) == str(manifest["peptide_sequence"])
    predicted_target_ca = _ca_coords(predicted_target)
    template_target_ca = _ca_coords(template_target)
    rotation, translation, target_rmsd = _kabsch(predicted_target_ca, template_target_ca)
    reference_rotation, reference_translation, reference_target_rmsd = _kabsch(
        _ca_coords(reference_target),
        template_target_ca,
    )
    site_rmsd = _subset_ca_rmsd(predicted_target_ca, template_target_ca, site_indices, rotation, translation)
    hotspot_rmsd = _subset_ca_rmsd(predicted_target_ca, template_target_ca, hotspot_indices, rotation, translation)

    peptide_mobile, peptide_reference = _paired_atom_coords(predicted_peptide, reference_peptide, BACKBONE_ATOMS)
    peptide_reference_in_template_frame = _transform(peptide_reference, reference_rotation, reference_translation)
    peptide_rmsd = _rmsd(peptide_mobile, peptide_reference_in_template_frame, rotation, translation)
    design_target_contacts, design_peptide_contacts, _, _ = _contact_sets(reference_target, reference_peptide, args.contact_cutoff)
    reference_site_contacts = design_target_contacts & set(site_indices)
    reference_hotspot_contacts = design_target_contacts & set(hotspot_indices)
    reference_hotspot_min = _minimum_to_target_indices(reference_target, reference_peptide, hotspot_indices)
    reference_same_site = bool(
        reference_site_contacts
        and reference_hotspot_contacts
        and reference_hotspot_min <= args.hotspot_contact_distance
    )
    interface_indices = sorted(index - 1 for index in design_peptide_contacts)
    interface_mobile, interface_reference = _paired_atom_coords(
        predicted_peptide,
        reference_peptide,
        BACKBONE_ATOMS,
        interface_indices,
    )
    interface_reference_in_template_frame = _transform(interface_reference, reference_rotation, reference_translation)
    interface_rmsd = _rmsd(interface_mobile, interface_reference_in_template_frame, rotation, translation)

    predicted_contacts, _, _, _ = _contact_sets(predicted_target, predicted_peptide, args.contact_cutoff)
    site_contacts = predicted_contacts & set(site_indices)
    hotspot_contacts = predicted_contacts & set(hotspot_indices)
    offsite_contacts = predicted_contacts - set(site_indices)
    hotspot_min = _minimum_to_target_indices(predicted_target, predicted_peptide, hotspot_indices)
    same_site = bool(site_contacts and hotspot_contacts and hotspot_min <= args.hotspot_contact_distance)
    offsite_fraction = len(offsite_contacts) / max(1, len(predicted_contacts))

    terminal_cn = _terminal_cn_distance(predicted_peptide)
    macrocycle_status = _macrocycle_status(terminal_cn, args.macrocycle_pass_distance, args.macrocycle_warn_distance)
    severe_clashes = _severe_clash_count(predicted_target, predicted_peptide, args.severe_clash_distance)
    clash_status = "pass_no_severe_clash" if severe_clashes == 0 else "fail_severe_interchain_clash"
    peptide_rog = _radius_of_gyration(predicted_peptide)

    protocol_valid = _protocol_identity_valid(job, manifest, metadata, metric)
    target_failures = []
    if not protocol_valid:
        target_failures.append("protocol_identity_invalid")
    if not target_sequence_match:
        target_failures.append("target_sequence_mismatch")
    if target_rmsd > args.max_target_global_rmsd:
        target_failures.append("high_target_global_CA_RMSD")
    if site_rmsd > args.max_site2_local_rmsd:
        target_failures.append("high_Site_2_local_CA_RMSD")
    if hotspot_rmsd > args.max_hotspot_local_rmsd:
        target_failures.append("high_hotspot_local_CA_RMSD")
    target_recovery = not target_failures

    if peptide_rmsd <= args.strong_pose_rmsd:
        pose_class = "strong"
    elif peptide_rmsd <= args.moderate_pose_rmsd:
        pose_class = "moderate"
    else:
        pose_class = "poor"

    peptide_plddt = _float(metric.get("plddt_peptide_mean_100", ""))
    site_pae = _float(metric.get("interface_pae_site2_mean_A", ""))
    hotspot_pae = _float(metric.get("interface_pae_hotspot_mean_A", ""))
    ptm_value = _float(metric.get("ptm", ""))
    iptm_value = _float(metric.get("iptm", ""))
    ranking_value = _float(metric.get("ranking_confidence", ""))
    ranking_source = str(metric.get("ranking_confidence_source", "")).strip()
    if not math.isfinite(ranking_value) and math.isfinite(ptm_value) and math.isfinite(iptm_value):
        ranking_value = (0.8 * iptm_value) + (0.2 * ptm_value)
        ranking_source = "derived_0.8_iPTM_plus_0.2_pTM"
    if peptide_plddt >= 70.0 and site_pae <= 15.0 and hotspot_pae <= 15.0:
        confidence_status = "supportive"
    elif peptide_plddt >= 60.0 and min(site_pae, hotspot_pae) <= 20.0:
        confidence_status = "mixed"
    else:
        confidence_status = "weak"

    if not protocol_valid:
        model_support = "protocol_failure"
    elif not target_recovery:
        model_support = "target_conditioning_not_recovered"
    elif macrocycle_status != "pass_head_to_tail_geometry" or severe_clashes:
        model_support = "pose_not_recovered"
    elif same_site and pose_class == "strong":
        model_support = "strong_pose_support"
    elif same_site and pose_class == "moderate":
        model_support = "moderate_pose_support"
    else:
        model_support = "pose_not_recovered"

    aligned_target_ca = _transform(predicted_target_ca, rotation, translation)
    aligned_peptide_ca = _transform(_ca_coords(predicted_peptide), rotation, translation)
    coordinate_hash = _sha1_arrays(np.round(np.concatenate([aligned_target_ca, aligned_peptide_ca]), 2).astype(np.float32))
    confidence_hash = _confidence_hash(prediction_npz)
    identity_hash = hashlib.sha1(f"{coordinate_hash}:{confidence_hash}".encode("ascii")).hexdigest()

    return {
        "stage5B_candidate_id": manifest["stage5B_candidate_id"],
        "stage5B_job_id": job["stage5B_job_id"],
        "protocol_hash": manifest["protocol_hash"],
        "batch": manifest["batch"],
        "backbone_id": manifest["backbone_id"],
        "peptide_sequence": manifest["peptide_sequence"],
        "seed": job["seed"],
        "model_name": metric.get("model_name", prediction_pdb.stem),
        "prediction_pdb": prediction_pdb,
        "prediction_npz": prediction_npz,
        "target_template_path": metadata.get("target_template_path", manifest.get("target_template_pdb", "")),
        "target_template_sha1": metadata.get("target_template_sha1", ""),
        "target_template_coverage": f"{metadata.get('target_template_residues_covered', 0)}/86",
        "target_template_coverage_fraction": metadata.get("target_template_coverage_fraction", ""),
        "peptide_template_coverage": metadata.get("peptide_template_residues_covered", ""),
        "template_mode": metadata.get("template_mode", ""),
        "template_input_verified": str(bool(metadata.get("template_input_verified"))).lower(),
        "template_effect_observed": str(target_recovery).lower(),
        "use_initial_guess": str(bool(metadata.get("use_initial_guess"))).lower(),
        "target_sequence_match": str(target_sequence_match).lower(),
        "peptide_sequence_match": str(peptide_sequence_match).lower(),
        "reference_target_to_template_CA_RMSD_A": _rounded(reference_target_rmsd),
        "target_global_CA_RMSD_A": _rounded(target_rmsd),
        "Site_2_local_CA_RMSD_A": _rounded(site_rmsd),
        "hotspot_local_CA_RMSD_A": _rounded(hotspot_rmsd),
        "target_template_recovery_pass": str(target_recovery).lower(),
        "target_template_recovery_failure_reasons": ";".join(target_failures),
        "target_aligned_peptide_backbone_RMSD_A": _rounded(peptide_rmsd),
        "target_aligned_peptide_interface_RMSD_A": _rounded(interface_rmsd),
        "design_target_contact_count": len(design_target_contacts),
        "reference_design_site2_contact_count": len(reference_site_contacts),
        "reference_design_hotspot_contact_count": len(reference_hotspot_contacts),
        "reference_design_hotspot_min_distance_A": _rounded(reference_hotspot_min),
        "reference_design_same_target_site_flag": str(reference_same_site).lower(),
        "target_contact_recovery_count": len(design_target_contacts & predicted_contacts),
        "num_target_contacts": len(predicted_contacts),
        "num_site2_contacts": len(site_contacts),
        "num_hotspot_contacts": len(hotspot_contacts),
        "hotspot_min_distance_A": _rounded(hotspot_min),
        "off_site_contact_count": len(offsite_contacts),
        "off_site_contact_fraction": _rounded(offsite_fraction, 4),
        "same_target_site_flag": str(same_site).lower(),
        "terminal_C_N_distance_A": _rounded(terminal_cn),
        "macrocycle_geometry_status": macrocycle_status,
        "severe_clash_count": severe_clashes,
        "clash_status": clash_status,
        "peptide_radius_of_gyration_A": _rounded(peptide_rog),
        "target_mean_pLDDT_100": metric.get("plddt_target_mean_100", ""),
        "Site_2_mean_pLDDT_100": metric.get("plddt_site2_mean_100", ""),
        "hotspot_mean_pLDDT_100": metric.get("plddt_hotspot_mean_100", ""),
        "peptide_mean_pLDDT_100": metric.get("plddt_peptide_mean_100", ""),
        "interface_pae_global_mean_A": metric.get("interface_pae_global_mean_A", ""),
        "interface_pae_site2_mean_A": metric.get("interface_pae_site2_mean_A", ""),
        "interface_pae_site2_median_A": metric.get("interface_pae_site2_median_A", ""),
        "interface_pae_hotspot_mean_A": metric.get("interface_pae_hotspot_mean_A", ""),
        "pTM": _rounded(ptm_value, 6),
        "ipTM": _rounded(iptm_value, 6),
        "ranking_confidence": _rounded(ranking_value, 6),
        "ranking_confidence_source": ranking_source,
        "confidence_support_status": confidence_status,
        "pose_recovery_class": pose_class,
        "model_support_status": model_support,
        "prediction_coordinate_hash": coordinate_hash,
        "prediction_confidence_hash": confidence_hash,
        "prediction_identity_hash": identity_hash,
        "validation_test_type": metadata.get("validation_test_type", ""),
        "target_msa_mode": metadata.get("target_msa_mode", ""),
        "peptide_msa_mode": metadata.get("peptide_msa_mode", ""),
        "use_mlm": str(bool(metadata.get("use_mlm"))).lower(),
        "use_dropout": str(bool(metadata.get("use_dropout"))).lower(),
        "requested_recycles": metadata.get("requested_recycles", ""),
        "forward_passes": metadata.get("forward_passes", ""),
        "loaded_colabdesign_commit": metadata.get("loaded_colabdesign_commit", ""),
        "cyclic_topology_encoding": metadata.get("cyclic_topology_encoding", ""),
        "protocol_identity_valid": str(protocol_valid).lower(),
        "notes": (
            "Stage 4 reference peptide coordinates were transformed through a reference-target-to-Stage-0 alignment before pose RMSD. "
            "Off-site contact metrics are descriptive and are not a hard gate in Stage 5B v1."
        ),
        **{field: manifest.get(field, "") for field in SOURCE_ROUTE_PROVENANCE_FIELDS},
    }


def _best_model(rows: Sequence[Mapping[str, Any]]) -> Mapping[str, Any] | None:
    if not rows:
        return None
    order = {"strong_pose_support": 0, "moderate_pose_support": 1, "pose_not_recovered": 2, "target_conditioning_not_recovered": 3, "protocol_failure": 4}
    return min(
        rows,
        key=lambda row: (
            order.get(str(row.get("model_support_status", "")), 9),
            _float(row.get("target_aligned_peptide_backbone_RMSD_A", ""), 999.0),
            _float(row.get("hotspot_min_distance_A", ""), 999.0),
            _float(row.get("interface_pae_site2_mean_A", ""), 999.0),
            -_float(row.get("peptide_mean_pLDDT_100", ""), -1.0),
        ),
    )


def _seed_diversity(rows: Sequence[Mapping[str, Any]], expected_seeds: int) -> str:
    by_model: dict[str, list[Mapping[str, Any]]] = defaultdict(list)
    for row in rows:
        by_model[str(row.get("model_name", ""))].append(row)
    if len(by_model) != 5 or any(len({str(row.get("seed", "")) for row in group}) != expected_seeds for group in by_model.values()):
        return "incomplete_seed_model_matrix"
    unique_per_model = [len({str(row.get("prediction_identity_hash", "")) for row in group}) for group in by_model.values()]
    if all(count == expected_seeds for count in unique_per_model):
        return "seed_diverse"
    if all(count == 1 for count in unique_per_model):
        return "fully_redundant_across_seeds"
    return "partially_redundant_across_seeds"


def _candidate_summaries(
    manifest_rows: Sequence[Mapping[str, str]],
    model_rows: Sequence[Mapping[str, Any]],
    expected_models: int,
    expected_seeds: int,
) -> list[dict[str, Any]]:
    grouped: dict[str, list[Mapping[str, Any]]] = defaultdict(list)
    for row in model_rows:
        grouped[str(row["stage5B_candidate_id"])].append(row)
    summaries = []
    for manifest in manifest_rows:
        candidate_id = str(manifest["stage5B_candidate_id"])
        rows = grouped.get(candidate_id, [])
        identity_hashes = {str(row.get("prediction_identity_hash", "")) for row in rows if row.get("prediction_identity_hash")}
        protocol_valid = [row for row in rows if row.get("protocol_identity_valid") == "true"]
        strong = [row for row in rows if row.get("model_support_status") == "strong_pose_support"]
        moderate = [row for row in rows if row.get("model_support_status") == "moderate_pose_support"]
        strong_unique = {str(row["prediction_identity_hash"]) for row in strong}
        strong_seeds = {str(row["seed"]) for row in strong}
        strong_models = {str(row["model_name"]) for row in strong}
        completed = len(rows)
        reference_row = rows[0] if rows else {}
        reference_same_site = reference_row.get("reference_design_same_target_site_flag") == "true"
        if completed != expected_models or len(protocol_valid) != completed:
            support_class = "stage5B_protocol_failure"
            reason = f"Completed {completed}/{expected_models} models or one or more protocol identities were invalid."
        elif not reference_same_site:
            support_class = "stage5B_protocol_failure"
            reason = (
                "The staged Stage 4 reference does not contact the corrected Site_2/hotspot mapping, "
                "so target-site pose recovery cannot be interpreted for this candidate."
            )
        elif len(strong_unique) >= 3 and len(strong_seeds) >= 2 and len(strong_models) >= 2:
            support_class = "stage5B_strong_support"
            reason = "At least three unique strong same-site pose recoveries span at least two seeds and two model parameter sets."
        elif strong or moderate:
            support_class = "stage5B_partial_support"
            reason = "At least one protocol-valid strong or moderate same-site pose recovery was observed, but strong-support replication was not met."
        else:
            support_class = "stage5B_not_recovered"
            reason = "No protocol-valid strong or moderate same-site pose recovery was observed."
        best = _best_model(rows)
        summaries.append(
            {
                "stage5B_candidate_id": candidate_id,
                "batch": manifest.get("batch", ""),
                "backbone_id": manifest.get("backbone_id", ""),
                "peptide_sequence": manifest.get("peptide_sequence", ""),
                "models_expected": expected_models,
                "models_completed": completed,
                "protocol_valid_model_count": len(protocol_valid),
                "effective_unique_prediction_count": len(identity_hashes),
                "duplicate_prediction_count": completed - len(identity_hashes),
                "seed_diversity_status": _seed_diversity(rows, expected_seeds),
                "reference_design_site2_contact_count": reference_row.get("reference_design_site2_contact_count", ""),
                "reference_design_hotspot_contact_count": reference_row.get("reference_design_hotspot_contact_count", ""),
                "reference_design_hotspot_min_distance_A": reference_row.get("reference_design_hotspot_min_distance_A", ""),
                "reference_design_same_target_site_flag": reference_row.get("reference_design_same_target_site_flag", ""),
                "target_template_recovery_count": sum(row.get("target_template_recovery_pass") == "true" for row in rows),
                "Site_2_recovery_count": sum(int(row.get("num_site2_contacts", 0)) >= 1 for row in rows),
                "hotspot_recovery_count": sum(int(row.get("num_hotspot_contacts", 0)) >= 1 for row in rows),
                "same_target_site_count": sum(row.get("same_target_site_flag") == "true" for row in rows),
                "strong_pose_recovery_count": len(strong),
                "moderate_pose_recovery_count": len(moderate),
                "macrocycle_pass_count": sum(row.get("macrocycle_geometry_status") == "pass_head_to_tail_geometry" for row in rows),
                "confidence_pass_count": sum(row.get("confidence_support_status") == "supportive" for row in rows),
                "median_peptide_RMSD_A": _median(rows, "target_aligned_peptide_backbone_RMSD_A"),
                "best_peptide_RMSD_A": _best_numeric(rows, "target_aligned_peptide_backbone_RMSD_A"),
                "median_Site_2_PAE_A": _median(rows, "interface_pae_site2_mean_A"),
                "best_Site_2_PAE_A": _best_numeric(rows, "interface_pae_site2_mean_A"),
                "median_hotspot_distance_A": _median(rows, "hotspot_min_distance_A"),
                "best_hotspot_distance_A": _best_numeric(rows, "hotspot_min_distance_A"),
                "median_peptide_pLDDT_100": _median(rows, "peptide_mean_pLDDT_100"),
                "best_model_path": best.get("prediction_pdb", "") if best else "",
                "best_model_name": best.get("model_name", "") if best else "",
                "best_seed": best.get("seed", "") if best else "",
                "stage5B_support_class": support_class,
                "support_reason": reason,
                **{field: manifest.get(field, "") for field in SOURCE_ROUTE_PROVENANCE_FIELDS},
            }
        )
    return summaries


def _report(
    models: Sequence[Mapping[str, Any]],
    candidates: Sequence[Mapping[str, Any]],
    jobs: Sequence[Mapping[str, Any]],
    args: argparse.Namespace,
    stage5b_dir: Path,
) -> str:
    candidate_columns = [
        "stage5B_candidate_id", "backbone_id", "models_completed", "effective_unique_prediction_count",
        "reference_design_site2_contact_count", "reference_design_hotspot_contact_count",
        "reference_design_hotspot_min_distance_A", "reference_design_same_target_site_flag",
        "target_template_recovery_count", "Site_2_recovery_count", "hotspot_recovery_count",
        "same_target_site_count", "strong_pose_recovery_count", "moderate_pose_recovery_count",
        "macrocycle_pass_count", "seed_diversity_status", "median_peptide_RMSD_A",
        "best_hotspot_distance_A", "median_peptide_pLDDT_100", "stage5B_support_class",
    ]
    supported = [
        row
        for row in candidates
        if row.get("stage5B_support_class") in {"stage5B_strong_support", "stage5B_partial_support"}
    ]
    valid_references = [row for row in candidates if row.get("reference_design_same_target_site_flag") == "true"]
    models_by_job: dict[str, set[str]] = defaultdict(set)
    for row in models:
        models_by_job[str(row.get("stage5B_job_id", ""))].add(str(row.get("model_name", "")))
    completed_jobs = sum(len(model_names) == 5 for model_names in models_by_job.values())
    expected_jobs = len(jobs)
    expected_models = expected_jobs * 5
    template_verified = sum(row.get("template_input_verified") == "true" for row in models)
    peptide_zero = sum(str(row.get("peptide_template_coverage", "")) == "0" for row in models)
    return f"""# FGA RFpeptides Stage 5B Validation Report

Stage 5B is a target-structure-conditioned peptide pose-recovery test. It is
not sequence-only independent recovery, experimental binding validation, or a
final candidate decision. The Stage 0 86-aa target backbone is provided as a
target-only template; peptide template coverage and initial-guess usage remain
zero.

## Completion

```text
seed_jobs_completed: {completed_jobs}/{expected_jobs}
models_parsed: {len(models)}/{expected_models}
target_template_input_verified: {template_verified}/{len(models)}
peptide_template_coverage_zero: {peptide_zero}/{len(models)}
candidates_with_strong_or_partial_support: {len(supported)}/{len(candidates)}
reference_designs_valid_for_corrected_Site_2: {len(valid_references)}/{len(candidates)}
output_directory: {stage5b_dir}
```

## Target/template Recovery Distributions

Values are reported as `minimum / median / maximum` in angstrom.

```text
target_global_CA_RMSD_A: {_distribution_text(models, "target_global_CA_RMSD_A")}
Site_2_local_CA_RMSD_A: {_distribution_text(models, "Site_2_local_CA_RMSD_A")}
hotspot_local_CA_RMSD_A: {_distribution_text(models, "hotspot_local_CA_RMSD_A")}
```

## Internal Reporting Rules

```text
target_global_CA_RMSD_pass_A: <= {args.max_target_global_rmsd:g}
Site_2_local_CA_RMSD_pass_A: <= {args.max_site2_local_rmsd:g}
hotspot_local_CA_RMSD_pass_A: <= {args.max_hotspot_local_rmsd:g}
contact_cutoff_A: <= {args.contact_cutoff:g}
hotspot_same_site_distance_A: <= {args.hotspot_contact_distance:g}
macrocycle_pass_terminal_C_N_A: <= {args.macrocycle_pass_distance:g}
severe_clash_heavy_atom_pair_A: < {args.severe_clash_distance:g}
strong_pose_RMSD_A: <= {args.strong_pose_rmsd:g}
moderate_pose_RMSD_A: > {args.strong_pose_rmsd:g} and <= {args.moderate_pose_rmsd:g}
```

These thresholds are internal computational reporting rules, not experimental
standards. Confidence is reported separately and does not veto a structure on
one PAE field alone. Off-site contact counts are descriptive rather than a hard
gate because the Stage 4 reference poses themselves can contact neighboring
target-crop residues outside the mapped Site_2 core.

The collector independently checks whether every staged Stage 4 reference pose
contacts Site_2 and a selected hotspot after mapping Stage 0 crop positions by
target-chain order. A reference that fails this check is an upstream input
premise failure and receives `stage5B_protocol_failure`; its predictions cannot
be interpreted as a valid recovery test for the intended site.

## Candidate Summary

{rows_to_markdown(candidates, candidate_columns, "No Stage 5B predictions have been parsed.")}
"""


def main() -> int:
    parser = argparse.ArgumentParser(description="Collect Stage 5B target-conditioned AfCycDesign predictions.")
    parser.add_argument("--stage5b-root", required=True)
    parser.add_argument("--stage0-root", required=True)
    parser.add_argument("--project-config", required=True)
    parser.add_argument("--contact-cutoff", type=float, default=5.0)
    parser.add_argument("--hotspot-contact-distance", type=float, default=5.0)
    parser.add_argument("--severe-clash-distance", type=float, default=1.2)
    parser.add_argument("--macrocycle-pass-distance", type=float, default=2.0)
    parser.add_argument("--macrocycle-warn-distance", type=float, default=3.0)
    parser.add_argument("--max-target-global-rmsd", type=float, default=3.0)
    parser.add_argument("--max-site2-local-rmsd", type=float, default=2.0)
    parser.add_argument("--max-hotspot-local-rmsd", type=float, default=2.0)
    parser.add_argument("--strong-pose-rmsd", type=float, default=3.0)
    parser.add_argument("--moderate-pose-rmsd", type=float, default=5.0)
    args = parser.parse_args()
    if args.macrocycle_warn_distance < args.macrocycle_pass_distance:
        raise RuntimeError("--macrocycle-warn-distance must be >= --macrocycle-pass-distance")
    if args.moderate_pose_rmsd < args.strong_pose_rmsd:
        raise RuntimeError("--moderate-pose-rmsd must be >= --strong-pose-rmsd")

    logger = setup_logger("31_collect_stage5b_validation")
    append_run_header(logger, "31_collect_stage5b_validation.py")
    stage5b_dir = _resolve_mixed_path(args.stage5b_root)
    stage0_root = _resolve_mixed_path(args.stage0_root)
    assert_active_route_path(stage5b_dir, "Stage 31 Stage 5B root")
    assert_active_route_path(stage0_root, "Stage 31 Stage 0 root")
    route_manifest_path, route_manifest, route_manifest_sha256 = load_route_manifest(stage5b_dir.parent)
    validate_route_project_config(args.project_config, route_manifest)
    route_provenance = route_provenance_fields(route_manifest_path, route_manifest, route_manifest_sha256)
    stage0_dir = stage0_root / "00_target_inputs"
    manifest_csv = stage5b_dir / "FGA_rfpeptides_stage5B_candidate_manifest.csv"
    jobs_csv = stage5b_dir / "FGA_rfpeptides_stage5B_prediction_jobs.csv"
    assert_active_route_path(manifest_csv, "Stage 31 Stage 5B candidate manifest CSV")
    assert_active_route_path(jobs_csv, "Stage 31 Stage 5B prediction jobs CSV")
    manifest_rows = read_csv(manifest_csv)
    job_rows = read_csv(jobs_csv)
    if not manifest_rows or not job_rows:
        raise RuntimeError(f"Missing Stage 5B manifest or jobs CSV: {stage5b_dir}")
    for row in manifest_rows:
        validate_row_route_provenance(row, route_provenance, "Stage 31 Stage 5B candidate row")
        validate_source_route_provenance(row, "Stage 31 Stage 5B candidate row")
    for row in job_rows:
        validate_row_route_provenance(row, route_provenance, "Stage 31 Stage 5B job row")
        validate_source_route_provenance(row, "Stage 31 Stage 5B job row")
    manifest_lookup = {row["stage5B_candidate_id"]: row for row in manifest_rows}
    target_template = stage0_dir / "RFpep_Site_2_target.pdb"
    assert_active_route_path(target_template, "Stage 31 Stage 0 target template")
    template_chains = parse_residues(target_template)
    if set(template_chains) != {"A"} or len(template_chains["A"]) != 86:
        raise RuntimeError("Stage 0 target template must contain only chain A with 86 residues")
    template_target = template_chains["A"]
    mapping_csv = stage0_dir / "RFpep_Site_2_crop_renumbering_mapping.csv"
    assert_active_route_path(mapping_csv, "Stage 31 Stage 0 mapping CSV")
    site_indices, hotspot_indices = _mapping_indices(
        mapping_csv,
        86,
    )

    model_rows: list[dict[str, Any]] = []
    for job in job_rows:
        manifest = manifest_lookup.get(job["stage5B_candidate_id"])
        if manifest is None:
            logger.error("No manifest row for job %s", job.get("stage5B_job_id", ""))
            continue
        prediction_dir = _resolve_mixed_path(job["prediction_output_dir"])
        assert_active_route_path(prediction_dir, "Stage 31 prediction directory", must_exist=False)
        metadata_path = prediction_dir / "run_metadata.json"
        metrics_path = prediction_dir / "model_metrics.csv"
        if not metadata_path.is_file() or not metrics_path.is_file():
            continue
        assert_active_route_path(metadata_path, "Stage 31 run metadata JSON")
        assert_active_route_path(metrics_path, "Stage 31 model metrics CSV")
        try:
            metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
            with metrics_path.open("r", encoding="utf-8", newline="") as handle:
                metrics = list(csv.DictReader(handle))
            reference_pdb = _resolve_mixed_path(manifest["staged_reference_design_pdb"])
            assert_active_route_path(reference_pdb, "Stage 31 staged reference PDB")
            reference_chains = parse_residues(reference_pdb)
            if "A" not in reference_chains or "B" not in reference_chains:
                raise RuntimeError("Stage 4 reference PDB lacks chain A or B")
            for metric in metrics:
                prediction_pdb = _resolve_mixed_path(metric.get("prediction_pdb", ""))
                prediction_npz = _resolve_mixed_path(metric.get("prediction_npz", ""))
                if not prediction_pdb.is_file() or not prediction_npz.is_file():
                    logger.error("Missing Stage 5B PDB/NPZ pair for %s", metric.get("model_name", ""))
                    continue
                assert_active_route_path(prediction_pdb, "Stage 31 prediction PDB")
                assert_active_route_path(prediction_npz, "Stage 31 prediction NPZ")
                model_rows.append(
                    _model_row(
                        manifest,
                        job,
                        metric,
                        metadata,
                        prediction_pdb,
                        prediction_npz,
                        reference_chains,
                        template_target,
                        site_indices,
                        hotspot_indices,
                        args,
                    )
                )
        except Exception as exc:
            logger.error("Failed to parse job %s: %s: %s", job.get("stage5B_job_id", ""), exc.__class__.__name__, exc)

    expected_models_per_candidate = 25
    add_route_provenance(model_rows, route_provenance)
    candidate_rows = _candidate_summaries(manifest_rows, model_rows, expected_models_per_candidate, 5)
    add_route_provenance(candidate_rows, route_provenance)
    write_csv(stage5b_dir / "FGA_rfpeptides_stage5B_model_results.csv", model_rows, MODEL_FIELDS)
    write_csv(stage5b_dir / "FGA_rfpeptides_stage5B_candidate_summary.csv", candidate_rows, CANDIDATE_FIELDS)
    write_markdown(
        stage5b_dir / "FGA_rfpeptides_stage5B_validation_report.md",
        _report(model_rows, candidate_rows, job_rows, args, stage5b_dir),
    )
    logger.info("Stage 5B model predictions parsed: %s/125", len(model_rows))
    logger.info("Stage 5B candidate summaries: %s", len(candidate_rows))
    logger.info("Output directory: %s", stage5b_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
