from __future__ import annotations

import argparse
import math
import os
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np

from common import append_run_header, read_csv, resolve_path, rows_to_markdown, setup_logger, write_csv, write_markdown
from pdb_utils import parse_residues, residue_sequence


BACKBONE_ATOMS = ("N", "CA", "C")

MODEL_FIELDS = [
    "stage5_candidate_id",
    "stage5_job_id",
    "protocol_hash",
    "peptide_sequence_hash",
    "batch",
    "backbone_id",
    "peptide_sequence",
    "seed",
    "model_name",
    "prediction_pdb",
    "prediction_npz",
    "target_chain",
    "peptide_chain",
    "target_sequence_match",
    "peptide_sequence_match",
    "target_alignment_rmsd",
    "target_aligned_peptide_backbone_rmsd",
    "design_site_contact_count",
    "predicted_site_contact_count",
    "site_contact_recovery_count",
    "site_contact_recovery_fraction",
    "site_contact_recovery_status",
    "design_hotspot_contact_count",
    "predicted_hotspot_contact_count",
    "hotspot_contact_recovery_count",
    "hotspot_contact_recovery_fraction",
    "hotspot_min_distance",
    "off_site_contact_count",
    "predicted_site_contact_residues",
    "predicted_hotspot_contact_residues",
    "predicted_off_site_contact_residues",
    "same_target_site_flag",
    "macrocycle_terminal_cn_distance",
    "cyclic_topology_status",
    "cyclic_topology_note",
    "peptide_target_min_distance",
    "severe_clash",
    "clash_status",
    "plddt_mean_fraction",
    "plddt_peptide_mean_fraction",
    "plddt_mean_100",
    "plddt_peptide_mean_100",
    "interface_pae_global_mean_A",
    "interface_pae_site2_mean_A",
    "interface_pae_site2_median_A",
    "interface_pae_hotspot_mean_A",
    "iptm",
    "ptm",
    "geometry_recovery_flag",
    "confidence_pass",
    "recovery_success",
    "recovery_failure_reasons",
]

CANDIDATE_FIELDS = [
    "stage5_candidate_id",
    "batch",
    "backbone_id",
    "peptide_sequence",
    "seeds_expected",
    "seeds_completed",
    "models_completed",
    "seeds_with_geometry_recovery",
    "recovery_success_count",
    "best_model_name",
    "best_seed",
    "best_prediction_pdb",
    "best_target_aligned_peptide_backbone_rmsd",
    "best_predicted_site_contact_count",
    "best_predicted_hotspot_contact_count",
    "best_hotspot_min_distance",
    "best_macrocycle_terminal_cn_distance",
    "best_plddt_peptide_mean_100",
    "best_interface_pae_global_mean_A",
    "best_interface_pae_site2_mean_A",
    "best_interface_pae_hotspot_mean_A",
    "best_iptm",
    "candidate_validation_status",
    "notes",
]


def _split_csv(value: str) -> list[str]:
    return [item.strip() for item in str(value or "").split(",") if item.strip()]


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
        out = float(str(value).strip())
    except (TypeError, ValueError):
        return default
    return out if math.isfinite(out) else default


def _round(value: float | None, digits: int = 3) -> float | str:
    if value is None or not math.isfinite(value):
        return ""
    return round(value, digits)


def _coords(residues: Sequence[Mapping[str, Any]], atom_names: Sequence[str]) -> np.ndarray:
    points = []
    for residue in residues:
        atoms = residue.get("atoms", {})
        for atom_name in atom_names:
            if atom_name in atoms:
                points.append(atoms[atom_name])
    return np.asarray(points, dtype=float)


def _paired_backbone_coords(
    mobile: Sequence[Mapping[str, Any]],
    reference: Sequence[Mapping[str, Any]],
) -> tuple[np.ndarray, np.ndarray]:
    if len(mobile) != len(reference):
        raise RuntimeError(f"Chain length mismatch for alignment: {len(mobile)} != {len(reference)}")
    mobile_points = []
    reference_points = []
    for mobile_residue, reference_residue in zip(mobile, reference):
        mobile_atoms = mobile_residue.get("atoms", {})
        reference_atoms = reference_residue.get("atoms", {})
        for atom_name in BACKBONE_ATOMS:
            if atom_name in mobile_atoms and atom_name in reference_atoms:
                mobile_points.append(mobile_atoms[atom_name])
                reference_points.append(reference_atoms[atom_name])
    if len(mobile_points) < 3:
        raise RuntimeError("Too few paired target backbone atoms for alignment")
    return np.asarray(mobile_points, dtype=float), np.asarray(reference_points, dtype=float)


def _kabsch_transform(mobile: np.ndarray, reference: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
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


def _rmsd(left: np.ndarray, right: np.ndarray) -> float:
    if left.shape != right.shape or not left.size:
        return float("nan")
    return float(np.sqrt(np.mean(np.sum((left - right) ** 2, axis=1))))


def _residue_min_distance(left: Mapping[str, Any], right_atoms: np.ndarray) -> float:
    atoms = np.asarray(list(left.get("atoms", {}).values()), dtype=float)
    if not atoms.size or not right_atoms.size:
        return float("inf")
    delta = atoms[:, None, :] - right_atoms[None, :, :]
    return float(np.sqrt(np.sum(delta * delta, axis=2)).min())


def _contact_indices(
    target_residues: Sequence[Mapping[str, Any]],
    peptide_residues: Sequence[Mapping[str, Any]],
    cutoff: float,
) -> tuple[set[int], float]:
    peptide_atoms = _coords(peptide_residues, tuple({atom for residue in peptide_residues for atom in residue.get("atoms", {})}))
    contacts: set[int] = set()
    minimum = float("inf")
    for index, residue in enumerate(target_residues, start=1):
        residue_distance = _residue_min_distance(residue, peptide_atoms)
        minimum = min(minimum, residue_distance)
        if residue_distance <= cutoff:
            contacts.add(index)
    return contacts, minimum


def _minimum_for_indices(
    target_residues: Sequence[Mapping[str, Any]],
    peptide_residues: Sequence[Mapping[str, Any]],
    indices: set[int],
) -> float:
    peptide_atoms = _coords(peptide_residues, tuple({atom for residue in peptide_residues for atom in residue.get("atoms", {})}))
    values = [
        _residue_min_distance(target_residues[index - 1], peptide_atoms)
        for index in sorted(indices)
        if 1 <= index <= len(target_residues)
    ]
    return min(values, default=float("inf"))


def _terminal_cn_distance(peptide_residues: Sequence[Mapping[str, Any]]) -> float:
    if not peptide_residues:
        return float("nan")
    first_atoms = peptide_residues[0].get("atoms", {})
    last_atoms = peptide_residues[-1].get("atoms", {})
    if "N" not in first_atoms or "C" not in last_atoms:
        return float("nan")
    return float(np.linalg.norm(np.asarray(last_atoms["C"]) - np.asarray(first_atoms["N"])))


def _load_site_indices(mapping_csv: Path) -> tuple[set[int], set[int]]:
    rows = read_csv(mapping_csv)
    if not rows:
        raise RuntimeError(f"Missing or empty Stage 0 mapping CSV: {mapping_csv}")
    site = {
        int(row["rfpeptides_residue_number"])
        for row in rows
        if str(row.get("is_target_site_residue", "")).strip().lower() == "true"
    }
    hotspots = {
        int(row["rfpeptides_residue_number"])
        for row in rows
        if str(row.get("is_selected_hotspot", "")).strip().lower() == "true"
    }
    if not site or not hotspots:
        raise RuntimeError(f"Stage 0 mapping lacks site or hotspot residues: {mapping_csv}")
    return site, hotspots


def _metric_lookup(seed_dir: Path) -> dict[str, dict[str, str]]:
    rows = read_csv(seed_dir / "model_metrics.csv")
    lookup: dict[str, dict[str, str]] = {}
    for row in rows:
        path = _resolve_mixed_path(row.get("prediction_pdb", ""))
        lookup[path.name] = row
    return lookup


def _cyclic_status(distance: float, pass_distance: float, warn_distance: float) -> str:
    if not math.isfinite(distance):
        return "not_evaluable"
    if distance <= pass_distance:
        return "pass_head_to_tail_geometry"
    if distance <= warn_distance:
        return "warn_head_to_tail_geometry"
    return "fail_open_chain_geometry"


def _contact_recovery_status(design_contacts: set[int], predicted_contacts: set[int]) -> tuple[int, float, str]:
    overlap = len(design_contacts & predicted_contacts)
    if not design_contacts:
        return overlap, float("nan"), "design_has_no_contacts"
    fraction = overlap / len(design_contacts)
    if fraction >= 1.0:
        status = "full_contact_recovery"
    elif overlap:
        status = "partial_contact_recovery"
    else:
        status = "no_exact_contact_recovery"
    return overlap, fraction, status


def _labels(indices: Iterable[int]) -> str:
    return ",".join(f"A{index}" for index in sorted(indices))


def _model_row(
    *,
    candidate: Mapping[str, str],
    prediction_pdb: Path,
    metric: Mapping[str, str],
    reference_chains: Mapping[str, Sequence[Mapping[str, Any]]],
    site_indices: set[int],
    hotspot_indices: set[int],
    args: argparse.Namespace,
) -> dict[str, Any]:
    target_chain = str(candidate.get("target_chain", "A")) or "A"
    peptide_chain = str(candidate.get("peptide_chain", "B")) or "B"
    prediction_chains = parse_residues(prediction_pdb)
    if target_chain not in prediction_chains or peptide_chain not in prediction_chains:
        raise RuntimeError(f"Expected chains {target_chain}/{peptide_chain} in {prediction_pdb}")

    reference_target = list(reference_chains[target_chain])
    reference_peptide = list(reference_chains[peptide_chain])
    predicted_target = prediction_chains[target_chain]
    predicted_peptide = prediction_chains[peptide_chain]
    target_match = residue_sequence(predicted_target) == residue_sequence(reference_target)
    peptide_match = residue_sequence(predicted_peptide) == str(candidate["peptide_sequence"])

    mobile_target, reference_target_coords = _paired_backbone_coords(predicted_target, reference_target)
    rotation, translation, target_alignment_rmsd = _kabsch_transform(mobile_target, reference_target_coords)
    mobile_peptide, reference_peptide_coords = _paired_backbone_coords(predicted_peptide, reference_peptide)
    peptide_rmsd = _rmsd(_transform(mobile_peptide, rotation, translation), reference_peptide_coords)

    design_contacts, _ = _contact_indices(reference_target, reference_peptide, args.contact_cutoff)
    predicted_contacts, peptide_target_min = _contact_indices(predicted_target, predicted_peptide, args.contact_cutoff)
    design_site = design_contacts & site_indices
    predicted_site = predicted_contacts & site_indices
    design_hotspots = design_contacts & hotspot_indices
    predicted_hotspots = predicted_contacts & hotspot_indices
    predicted_offsite = predicted_contacts - site_indices
    site_overlap, site_fraction, site_status = _contact_recovery_status(design_site, predicted_site)
    hotspot_overlap, hotspot_fraction, _ = _contact_recovery_status(design_hotspots, predicted_hotspots)
    hotspot_min = _minimum_for_indices(predicted_target, predicted_peptide, hotspot_indices)
    terminal_cn = _terminal_cn_distance(predicted_peptide)
    cyclic_status = _cyclic_status(terminal_cn, args.macrocycle_pass_distance, args.macrocycle_warn_distance)
    severe_clash = peptide_target_min < args.severe_clash_distance
    same_site = bool(predicted_site and predicted_hotspots)

    plddt_mean_fraction = _float(metric.get("plddt_mean_fraction", metric.get("plddt_mean", "")))
    plddt_peptide_mean_fraction = _float(
        metric.get("plddt_peptide_mean_fraction", metric.get("plddt_peptide_mean", ""))
    )
    plddt_mean_100 = _float(metric.get("plddt_mean_100", ""))
    if not math.isfinite(plddt_mean_100) and math.isfinite(plddt_mean_fraction):
        plddt_mean_100 = plddt_mean_fraction * 100.0
    plddt_peptide_mean_100 = _float(metric.get("plddt_peptide_mean_100", ""))
    if not math.isfinite(plddt_peptide_mean_100) and math.isfinite(plddt_peptide_mean_fraction):
        plddt_peptide_mean_100 = plddt_peptide_mean_fraction * 100.0

    interface_pae_global_mean = _float(
        metric.get("interface_pae_global_mean_A", metric.get("interface_pae_normalized", ""))
    )
    interface_pae_site2_mean = _float(
        metric.get("interface_pae_site2_mean_A", metric.get("interface_pae_normalized", ""))
    )
    interface_pae_site2_median = _float(metric.get("interface_pae_site2_median_A", ""))
    interface_pae_hotspot_mean = _float(metric.get("interface_pae_hotspot_mean_A", ""))
    iptm = _float(metric.get("iptm", ""))
    ptm = _float(metric.get("ptm", ""))
    geometry_recovery = bool(
        target_match
        and peptide_match
        and same_site
        and cyclic_status == "pass_head_to_tail_geometry"
        and not severe_clash
        and peptide_rmsd <= args.max_peptide_backbone_rmsd
    )
    confidence_pass = bool(
        math.isfinite(plddt_peptide_mean_fraction)
        and plddt_peptide_mean_fraction >= args.min_peptide_plddt
        and math.isfinite(interface_pae_site2_mean)
        and interface_pae_site2_mean <= args.max_interface_pae
    )

    failure_reasons = []
    if not target_match:
        failure_reasons.append("target_sequence_mismatch")
    if not peptide_match:
        failure_reasons.append("peptide_sequence_mismatch")
    if not predicted_site:
        failure_reasons.append("no_Site_2_contact")
    if not predicted_hotspots:
        failure_reasons.append("no_hotspot_contact")
    if cyclic_status != "pass_head_to_tail_geometry":
        failure_reasons.append(cyclic_status)
    if severe_clash:
        failure_reasons.append("severe_interchain_clash")
    if peptide_rmsd > args.max_peptide_backbone_rmsd:
        failure_reasons.append("peptide_backbone_not_recovered")
    if not math.isfinite(plddt_peptide_mean_fraction) or plddt_peptide_mean_fraction < args.min_peptide_plddt:
        failure_reasons.append("low_peptide_plddt")
    if not math.isfinite(interface_pae_site2_mean) or interface_pae_site2_mean > args.max_interface_pae:
        failure_reasons.append("high_Site_2_interface_pae")

    return {
        "stage5_candidate_id": candidate["stage5_candidate_id"],
        "stage5_job_id": metric.get("stage5_job_id", ""),
        "protocol_hash": metric.get("protocol_hash", candidate.get("protocol_hash", "")),
        "peptide_sequence_hash": metric.get(
            "peptide_sequence_hash", candidate.get("peptide_sequence_hash", "")
        ),
        "batch": candidate.get("batch", ""),
        "backbone_id": candidate.get("backbone_id", ""),
        "peptide_sequence": candidate.get("peptide_sequence", ""),
        "seed": metric.get("seed", ""),
        "model_name": metric.get("model_name", prediction_pdb.stem),
        "prediction_pdb": prediction_pdb,
        "prediction_npz": _resolve_mixed_path(metric.get("prediction_npz", "")) if metric.get("prediction_npz") else "",
        "target_chain": target_chain,
        "peptide_chain": peptide_chain,
        "target_sequence_match": str(target_match).lower(),
        "peptide_sequence_match": str(peptide_match).lower(),
        "target_alignment_rmsd": _round(target_alignment_rmsd),
        "target_aligned_peptide_backbone_rmsd": _round(peptide_rmsd),
        "design_site_contact_count": len(design_site),
        "predicted_site_contact_count": len(predicted_site),
        "site_contact_recovery_count": site_overlap,
        "site_contact_recovery_fraction": _round(site_fraction),
        "site_contact_recovery_status": site_status,
        "design_hotspot_contact_count": len(design_hotspots),
        "predicted_hotspot_contact_count": len(predicted_hotspots),
        "hotspot_contact_recovery_count": hotspot_overlap,
        "hotspot_contact_recovery_fraction": _round(hotspot_fraction),
        "hotspot_min_distance": _round(hotspot_min),
        "off_site_contact_count": len(predicted_offsite),
        "predicted_site_contact_residues": _labels(predicted_site),
        "predicted_hotspot_contact_residues": _labels(predicted_hotspots),
        "predicted_off_site_contact_residues": _labels(predicted_offsite),
        "same_target_site_flag": str(same_site).lower(),
        "macrocycle_terminal_cn_distance": _round(terminal_cn),
        "cyclic_topology_status": cyclic_status,
        "cyclic_topology_note": "relative-position cyclic offset; no explicit terminal covalent bond in PDB",
        "peptide_target_min_distance": _round(peptide_target_min),
        "severe_clash": str(severe_clash).lower(),
        "clash_status": "fail_severe_interchain_clash" if severe_clash else "pass_no_severe_interchain_clash",
        "plddt_mean_fraction": _round(plddt_mean_fraction, 5),
        "plddt_peptide_mean_fraction": _round(plddt_peptide_mean_fraction, 5),
        "plddt_mean_100": _round(plddt_mean_100, 2),
        "plddt_peptide_mean_100": _round(plddt_peptide_mean_100, 2),
        "interface_pae_global_mean_A": _round(interface_pae_global_mean, 3),
        "interface_pae_site2_mean_A": _round(interface_pae_site2_mean, 3),
        "interface_pae_site2_median_A": _round(interface_pae_site2_median, 3),
        "interface_pae_hotspot_mean_A": _round(interface_pae_hotspot_mean, 3),
        "iptm": _round(iptm, 5),
        "ptm": _round(ptm, 5),
        "geometry_recovery_flag": str(geometry_recovery).lower(),
        "confidence_pass": str(confidence_pass).lower(),
        "recovery_success": str(geometry_recovery and confidence_pass).lower(),
        "recovery_failure_reasons": ";".join(failure_reasons),
    }


def _best_row(rows: Sequence[Mapping[str, Any]]) -> Mapping[str, Any]:
    return min(
        rows,
        key=lambda row: (
            row.get("recovery_success") != "true",
            row.get("geometry_recovery_flag") != "true",
            row.get("same_target_site_flag") != "true",
            -_float(row.get("predicted_hotspot_contact_count", ""), 0.0),
            _float(row.get("hotspot_min_distance", ""), 999.0),
            -_float(row.get("predicted_site_contact_count", ""), 0.0),
            _float(row.get("target_aligned_peptide_backbone_rmsd", ""), 999.0),
            _float(row.get("interface_pae_site2_mean_A", ""), 999.0),
            -_float(row.get("plddt_peptide_mean_100", ""), -1.0),
        ),
    )


def _candidate_rows(
    candidates: Sequence[Mapping[str, str]],
    model_rows: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    grouped: dict[str, list[Mapping[str, Any]]] = defaultdict(list)
    for row in model_rows:
        grouped[str(row["stage5_candidate_id"])].append(row)
    summaries = []
    for candidate in candidates:
        candidate_id = candidate["stage5_candidate_id"]
        rows = grouped.get(candidate_id, [])
        seeds_completed = sorted({str(row.get("seed", "")) for row in rows})
        geometry_seeds = {
            str(row.get("seed", "")) for row in rows if row.get("geometry_recovery_flag") == "true"
        }
        success_seeds = {str(row.get("seed", "")) for row in rows if row.get("recovery_success") == "true"}
        if rows:
            best = _best_row(rows)
            expected_seeds = int(candidate.get("seeds_per_candidate", "5"))
            if len(success_seeds) >= 3:
                status = "validation_pass"
            elif len(seeds_completed) < expected_seeds and not geometry_seeds:
                status = "smoke_no_recovery_hold"
            elif len(seeds_completed) < expected_seeds and not success_seeds:
                status = "smoke_geometry_only_hold"
            else:
                status = "partial_pending_more_seeds"
            if len(seeds_completed) >= expected_seeds and len(success_seeds) < 3:
                status = "validation_not_recovered"
            notes = (
                "Independent recovery; success requires Site_2 and hotspot contact, peptide backbone recovery, "
                "head-to-tail geometry, no severe clash, peptide pLDDT and interface PAE thresholds."
            )
        else:
            best = {}
            status = "not_run"
            notes = "No prediction models found."
        summaries.append(
            {
                "stage5_candidate_id": candidate_id,
                "batch": candidate.get("batch", ""),
                "backbone_id": candidate.get("backbone_id", ""),
                "peptide_sequence": candidate.get("peptide_sequence", ""),
                "seeds_expected": candidate.get("seeds_per_candidate", "5"),
                "seeds_completed": len(seeds_completed),
                "models_completed": len(rows),
                "seeds_with_geometry_recovery": len(geometry_seeds),
                "recovery_success_count": len(success_seeds),
                "best_model_name": best.get("model_name", ""),
                "best_seed": best.get("seed", ""),
                "best_prediction_pdb": best.get("prediction_pdb", ""),
                "best_target_aligned_peptide_backbone_rmsd": best.get("target_aligned_peptide_backbone_rmsd", ""),
                "best_predicted_site_contact_count": best.get("predicted_site_contact_count", ""),
                "best_predicted_hotspot_contact_count": best.get("predicted_hotspot_contact_count", ""),
                "best_hotspot_min_distance": best.get("hotspot_min_distance", ""),
                "best_macrocycle_terminal_cn_distance": best.get("macrocycle_terminal_cn_distance", ""),
                "best_plddt_peptide_mean_100": best.get("plddt_peptide_mean_100", ""),
                "best_interface_pae_global_mean_A": best.get("interface_pae_global_mean_A", ""),
                "best_interface_pae_site2_mean_A": best.get("interface_pae_site2_mean_A", ""),
                "best_interface_pae_hotspot_mean_A": best.get("interface_pae_hotspot_mean_A", ""),
                "best_iptm": best.get("iptm", ""),
                "candidate_validation_status": status,
                "notes": notes,
            }
        )
    return summaries


def _summary_markdown(
    model_rows: Sequence[Mapping[str, Any]],
    candidate_rows: Sequence[Mapping[str, Any]],
    args: argparse.Namespace,
) -> str:
    completed = [row for row in candidate_rows if row.get("models_completed", 0)]
    model_passes = [row for row in model_rows if row.get("recovery_success") == "true"]
    geometry_passes = [row for row in model_rows if row.get("geometry_recovery_flag") == "true"]
    model_columns = [
        "stage5_candidate_id",
        "seed",
        "model_name",
        "target_aligned_peptide_backbone_rmsd",
        "predicted_site_contact_count",
        "predicted_hotspot_contact_count",
        "hotspot_min_distance",
        "off_site_contact_count",
        "same_target_site_flag",
        "macrocycle_terminal_cn_distance",
        "cyclic_topology_status",
        "clash_status",
        "plddt_peptide_mean_100",
        "interface_pae_global_mean_A",
        "interface_pae_site2_mean_A",
        "interface_pae_hotspot_mean_A",
        "iptm",
        "geometry_recovery_flag",
        "recovery_success",
        "recovery_failure_reasons",
    ]
    candidate_columns = [
        "stage5_candidate_id",
        "seeds_completed",
        "models_completed",
        "seeds_with_geometry_recovery",
        "recovery_success_count",
        "best_model_name",
        "best_target_aligned_peptide_backbone_rmsd",
        "best_predicted_site_contact_count",
        "best_predicted_hotspot_contact_count",
        "best_plddt_peptide_mean_100",
        "best_interface_pae_global_mean_A",
        "best_interface_pae_site2_mean_A",
        "best_interface_pae_hotspot_mean_A",
        "candidate_validation_status",
    ]
    return f"""# FGA RFpeptides Stage 5 AfCycDesign Validation

This report parses sequence/target-based independent-recovery predictions.
The Stage 4 design pose is used only after prediction for target alignment and
recovery measurements; it was not supplied as a template or initial guess.

AfCycDesign encodes peptide cyclicity through a relative-position cyclic
offset. The output PDB does not contain an explicit terminal C-N covalent bond,
so terminal C-N geometry is checked directly.

## Thresholds

```text
contact_cutoff_A: {args.contact_cutoff:g}
severe_clash_distance_A: {args.severe_clash_distance:g}
macrocycle_pass_distance_A: {args.macrocycle_pass_distance:g}
macrocycle_warn_distance_A: {args.macrocycle_warn_distance:g}
max_target_aligned_peptide_backbone_rmsd_A: {args.max_peptide_backbone_rmsd:g}
min_peptide_plddt_fraction: {args.min_peptide_plddt:g}
max_Site_2_interface_pae_mean_A: {args.max_interface_pae:g}
```

`geometry_recovery_flag` requires sequence identity, Site_2 and hotspot
contact, target-aligned peptide backbone RMSD, head-to-tail terminal geometry,
and no severe inter-chain clash. `recovery_success` additionally requires the
peptide pLDDT and Site_2-local interface PAE thresholds. pLDDT is reported on
the 0-100 scale in this report; PAE fields are in angstrom. The global PAE is
the mean across both target-to-peptide directions for all target residues and
is not treated as a Site_2-specific confidence metric.

## Current Counts

```text
candidates_with_any_output: {len(completed)}
models_parsed: {len(model_rows)}
models_passing_geometry_recovery: {len(geometry_passes)}
models_passing_strict_recovery: {len(model_passes)}
```

## Candidate Summary

{rows_to_markdown(candidate_rows, candidate_columns, "No candidate outputs were found.")}

## Model Results

{rows_to_markdown(model_rows, model_columns, "No prediction models were found.")}
"""


def main() -> int:
    parser = argparse.ArgumentParser(description="Parse Stage 5 AfCycDesign independent-recovery predictions.")
    parser.add_argument(
        "--stage5-root",
        default="results/rfpeptides_article_route_clean_20260623_stage5A_v2_batch01_batch02/07_structure_validation",
    )
    parser.add_argument("--stage0-root", default="results/rfpeptides_article_route_clean_20260615_fpocket")
    parser.add_argument("--selected-candidates", default="")
    parser.add_argument("--contact-cutoff", type=float, default=5.0)
    parser.add_argument("--severe-clash-distance", type=float, default=1.2)
    parser.add_argument("--macrocycle-pass-distance", type=float, default=2.0)
    parser.add_argument("--macrocycle-warn-distance", type=float, default=3.0)
    parser.add_argument("--max-peptide-backbone-rmsd", type=float, default=4.0)
    parser.add_argument("--min-peptide-plddt", type=float, default=0.60)
    parser.add_argument("--max-interface-pae", type=float, default=15.0)
    args = parser.parse_args()

    logger = setup_logger("27_collect_afcycdesign_validation")
    append_run_header(logger, "27_collect_afcycdesign_validation.py")
    if args.contact_cutoff <= 0 or args.severe_clash_distance <= 0:
        raise RuntimeError("Distance cutoffs must be > 0")
    if args.macrocycle_warn_distance < args.macrocycle_pass_distance:
        raise RuntimeError("--macrocycle-warn-distance must be >= --macrocycle-pass-distance")
    if not 0.0 <= args.min_peptide_plddt <= 1.0:
        raise RuntimeError("--min-peptide-plddt must be between 0 and 1")

    stage5_root = _resolve_mixed_path(args.stage5_root)
    stage0_root = _resolve_mixed_path(args.stage0_root)
    candidates = read_csv(stage5_root / "FGA_rfpeptides_stage5_candidate_manifest.csv")
    if not candidates:
        raise RuntimeError(f"Missing Stage 5 candidate manifest: {stage5_root}")
    selected = set(_split_csv(args.selected_candidates))
    if selected:
        candidates = [row for row in candidates if row.get("stage5_candidate_id", "") in selected]
    site_indices, hotspot_indices = _load_site_indices(
        stage0_root / "00_target_inputs" / "RFpep_Site_2_crop_renumbering_mapping.csv"
    )

    model_rows: list[dict[str, Any]] = []
    for candidate in candidates:
        candidate_id = candidate["stage5_candidate_id"]
        reference_pdb = _resolve_mixed_path(candidate["staged_design_pdb"])
        reference_chains = parse_residues(reference_pdb)
        target_chain = candidate.get("target_chain", "A") or "A"
        peptide_chain = candidate.get("peptide_chain", "B") or "B"
        if target_chain not in reference_chains or peptide_chain not in reference_chains:
            raise RuntimeError(f"Reference design lacks expected chains: {reference_pdb}")
        candidate_dir = stage5_root / "predictions" / candidate_id
        for seed_dir in sorted(candidate_dir.rglob("seed_*")) if candidate_dir.exists() else []:
            metrics = _metric_lookup(seed_dir)
            for prediction_pdb in sorted(seed_dir.glob("*.pdb")):
                metric = metrics.get(prediction_pdb.name)
                if metric is None:
                    logger.warning("Skipping PDB without model_metrics.csv row: %s", prediction_pdb)
                    continue
                try:
                    model_rows.append(
                        _model_row(
                            candidate=candidate,
                            prediction_pdb=prediction_pdb,
                            metric=metric,
                            reference_chains=reference_chains,
                            site_indices=site_indices,
                            hotspot_indices=hotspot_indices,
                            args=args,
                        )
                    )
                except Exception as exc:
                    logger.error("Failed to parse %s: %s: %s", prediction_pdb, exc.__class__.__name__, exc)

    summaries = _candidate_rows(candidates, model_rows)
    write_csv(stage5_root / "FGA_rfpeptides_stage5_model_validation.csv", model_rows, MODEL_FIELDS)
    write_csv(stage5_root / "FGA_rfpeptides_stage5_candidate_validation_summary.csv", summaries, CANDIDATE_FIELDS)
    write_markdown(
        stage5_root / "FGA_rfpeptides_stage5_validation.md",
        _summary_markdown(model_rows, summaries, args),
    )
    logger.info("Stage 5 models parsed: %s", len(model_rows))
    logger.info("Geometry recovery passes: %s", sum(row["geometry_recovery_flag"] == "true" for row in model_rows))
    logger.info("Strict recovery passes: %s", sum(row["recovery_success"] == "true" for row in model_rows))
    logger.info("Output directory: %s", stage5_root)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
