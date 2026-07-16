from __future__ import annotations

from typing import Any, Dict, Iterable, List, Mapping


FINAL_COLUMNS = [
    "peptide_id",
    "target",
    "target_uniprot",
    "target_patch",
    "uniprot_region",
    "hotspot_residues",
    "core_sequence",
    "core_length",
    "final_synthesis_format",
    "cyclization",
    "mean_iptm",
    "best_iptm",
    "mean_interface_pae",
    "best_interface_pae",
    "peptide_plddt",
    "interface_contacts",
    "pose_consistency_rmsd",
    "cys_cys_geometry",
    "net_charge",
    "hydrophobicity_flag",
    "sequence_filter_pass",
    "negative_screen_flag",
    "final_score",
    "priority",
    "notes",
]


def to_float(value: Any, default: float = 0.0) -> float:
    try:
        if value in ("", None):
            return default
        return float(value)
    except (TypeError, ValueError):
        return default


def normalized_complex_score(row: Mapping[str, Any]) -> float:
    iptm = max(0.0, min(1.0, to_float(row.get("best_iptm"))))
    pae = to_float(row.get("best_interface_pae"), 30.0)
    pae_score = max(0.0, min(1.0, (30.0 - pae) / 30.0))
    plddt = max(0.0, min(1.0, to_float(row.get("mean_peptide_plddt")) / 100.0))
    contacts = max(0.0, min(1.0, to_float(row.get("interface_contacts")) / 20.0))
    return 0.35 * iptm + 0.25 * pae_score + 0.20 * plddt + 0.20 * contacts


def synthesis_score(core_sequence: str) -> float:
    length = len(core_sequence)
    if length in (12, 14, 16):
        length_score = 1.0
    elif 10 <= length <= 18:
        length_score = 0.8
    else:
        length_score = 0.0
    oxidation_penalty = 0.1 * core_sequence.count("W") + 0.1 * core_sequence.count("M")
    return max(0.0, length_score - oxidation_penalty)


def score_candidate(row: Mapping[str, Any], negative_pass: bool = True) -> float:
    patch_consistency_score = 1.0 if str(row.get("patch_consistency_flag", "")).lower() == "pass" else 0.0
    cys_geometry_score = 1.0 if str(row.get("cys_cys_geometry", "")).lower() == "pass" else 0.0
    neg_score = 1.0 if negative_pass else 0.0
    core = str(row.get("core_sequence", ""))
    final = (
        0.35 * normalized_complex_score(row)
        + 0.20 * patch_consistency_score
        + 0.15 * cys_geometry_score
        + 0.15 * synthesis_score(core)
        + 0.15 * neg_score
    )
    return round(final, 4)


def rank_rows(rows: Iterable[Dict[str, Any]]) -> List[Dict[str, Any]]:
    usable = []
    for row in rows:
        neg_pass = str(row.get("negative_screen_flag", "")).lower() == "pass"
        row = dict(row)
        row["final_score"] = score_candidate(row, neg_pass)
        row["priority"] = "exclude" if row.get("cys_cys_geometry") == "fail" else row.get("priority", "ranked")
        usable.append(row)
    return sorted(usable, key=lambda r: (r.get("priority") == "exclude", -to_float(r.get("final_score"))))
