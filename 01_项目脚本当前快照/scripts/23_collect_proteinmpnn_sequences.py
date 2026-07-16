from __future__ import annotations

import argparse
import math
import os
from collections import Counter
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

from common import append_run_header, read_csv, resolve_path, rows_to_markdown, setup_logger, write_csv, write_fasta, write_markdown
from pdb_utils import ca_coord, centroid, distance, parse_residues, residue_sequence


STAGE3C_FIELDS = [
    "sequence_design_id",
    "backbone_id",
    "site_label",
    "site_id",
    "relaxed_pdb",
    "source_backbone_pdb",
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
]


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


def _lookup_rows(rows: Iterable[Mapping[str, str]], keys: Sequence[str]) -> dict[str, dict[str, str]]:
    lookup: dict[str, dict[str, str]] = {}
    for row in rows:
        item = dict(row)
        for key in keys:
            value = item.get(key, "")
            if value:
                lookup[value] = item
    return lookup


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


def _relaxed_pdbs_for_backbone(output_pdb_dir: Path, backbone_id: str) -> list[Path]:
    patterns = [
        f"{backbone_id}_mpnnonly_*_dldesign_*.pdb",
        f"{backbone_id}_mpnnfr_*_dldesign_*_cycle*.pdb",
        f"{backbone_id}*_dldesign_*.pdb",
        f"{backbone_id}*_dldesign_*_cycle*.pdb",
        f"{backbone_id}*.pdb",
    ]
    seen: set[Path] = set()
    files: list[Path] = []
    for pattern in patterns:
        for path in sorted(output_pdb_dir.glob(pattern)):
            if path not in seen:
                seen.add(path)
                files.append(path)
    return files


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
    backbone_id = str(backbone_row.get("design_id", "")).strip()
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


def _summary_markdown(
    *,
    rows: list[Mapping[str, Any]],
    args: argparse.Namespace,
    output_dir: Path,
    site_numbers: set[str],
    hotspot_numbers: set[str],
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
selected_backbones: {args.selected_backbones}
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
    parser.add_argument("--stage0-root", default="results/rfpeptides_article_route_clean_20260615_fpocket")
    parser.add_argument("--stage3-root", default="results/rfpeptides_article_route_clean_20260615_fpocket_stage1_N1000_no_traj")
    parser.add_argument("--output-root", default="", help="Defaults to --stage3-root.")
    parser.add_argument(
        "--stage3-mode",
        choices=["proteinmpnn_fastrelax", "proteinmpnn_only"],
        default="proteinmpnn_fastrelax",
        help="Select mode-specific Stage 3B job/output defaults.",
    )
    parser.add_argument("--selected-backbones", default="RFpep_Site_2_0007")
    parser.add_argument("--stage2-pass-csv", default="")
    parser.add_argument("--stage3-jobs-csv", default="")
    parser.add_argument("--input-pdb-dir", default="")
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

    stage0_root = _resolve_mixed_path(args.stage0_root)
    stage3_root = _resolve_mixed_path(args.stage3_root)
    output_root = _resolve_mixed_path(args.output_root) if args.output_root else stage3_root
    output_dir = output_root / "05_proteinmpnn_sequences"
    output_name_suffix = "" if args.stage3_mode == "proteinmpnn_fastrelax" else f"_{args.stage3_mode}"
    fasta_dir = output_dir / ("fasta" if args.stage3_mode == "proteinmpnn_fastrelax" else f"fasta_{args.stage3_mode}")
    default_input_subdir = "proteinmpnn_only_pdbs" if args.stage3_mode == "proteinmpnn_only" else "fastrelax_pdbs"
    input_pdb_dir = _resolve_mixed_path(args.input_pdb_dir) if args.input_pdb_dir else output_dir / default_input_subdir
    stage2_pass_csv = (
        _resolve_mixed_path(args.stage2_pass_csv)
        if args.stage2_pass_csv
        else stage3_root / "03_backbone_qc" / "FGA_rfpeptides_backbones_qc_pass.csv"
    )
    stage3_jobs_csv = (
        _resolve_mixed_path(args.stage3_jobs_csv)
        if args.stage3_jobs_csv
        else stage3_root / "04_proteinmpnn_inputs" / f"FGA_rfpeptides_stage3B_{args.stage3_mode}_jobs.csv"
    )
    old_fastrelax_jobs_csv = stage3_root / "04_proteinmpnn_inputs" / "FGA_rfpeptides_stage3_proteinmpnn_jobs.csv"
    if not args.stage3_jobs_csv and not stage3_jobs_csv.exists() and args.stage3_mode == "proteinmpnn_fastrelax":
        stage3_jobs_csv = old_fastrelax_jobs_csv

    stage2_pass_lookup = _lookup_rows(_read_required_csv(stage2_pass_csv), ["design_id"])
    stage3_job_rows = read_csv(stage3_jobs_csv)
    stage3_job_lookup = _lookup_rows(stage3_job_rows, ["design_id"]) if stage3_job_rows else {}

    selected_backbones = _split_csv(args.selected_backbones)
    if not selected_backbones:
        raise RuntimeError("--selected-backbones must not be empty")

    all_rows: list[dict[str, Any]] = []
    last_site_numbers: set[str] = set()
    last_hotspot_numbers: set[str] = set()
    last_target_chain = "A"
    last_peptide_chain = "B"
    for backbone_id in selected_backbones:
        backbone_row = stage2_pass_lookup.get(backbone_id)
        if backbone_row is None:
            raise RuntimeError(f"Selected backbone not found in Stage 2 pass CSV: {backbone_id}")
        if str(backbone_row.get("pass_backbone_qc", "")).strip().lower() != "true":
            raise RuntimeError(f"Selected backbone is not marked pass_backbone_qc=true: {backbone_id}")

        site_label = str(backbone_row.get("site_label", ""))
        site_mapping_csv = stage0_root / "00_target_inputs" / f"{_safe_token(site_label)}_crop_renumbering_mapping.csv"
        site_numbers, hotspot_numbers = _load_site_mapping(site_mapping_csv)
        last_site_numbers = site_numbers
        last_hotspot_numbers = hotspot_numbers
        last_target_chain = str(backbone_row.get("target_chain", "")).strip() or "A"
        last_peptide_chain = str(backbone_row.get("peptide_chain", "")).strip() or "B"

        job_row = stage3_job_lookup.get(backbone_id, {})
        output_pdb_dir = _resolve_mixed_path(str(job_row.get("output_pdb_dir", ""))) if job_row.get("output_pdb_dir") else input_pdb_dir
        relaxed_pdbs = _relaxed_pdbs_for_backbone(output_pdb_dir, backbone_id)
        if not relaxed_pdbs:
            missing_row = dict(backbone_row)
            missing_row["design_id"] = backbone_id
            all_rows.append(
                {
                    "sequence_design_id": f"{backbone_id}_missing_stage3b_output",
                    "backbone_id": backbone_id,
                    "site_label": backbone_row.get("site_label", ""),
                    "site_id": backbone_row.get("site_id", ""),
                    "source_backbone_pdb": backbone_row.get("rf_pdb", ""),
                    "file_status": "fail_missing_stage3b_outputs",
                    "parse_status": "not_parsed",
                    "target_chain": last_target_chain,
                    "peptide_chain": last_peptide_chain,
                    "expected_peptide_length": backbone_row.get("peptide_length", ""),
                    "forbidden_aas": args.forbidden_aas,
                    "sequence_status": "not_evaluated",
                    "pass_stage3c_qc": "false",
                    "qc_failure_reasons": "missing_stage3b_outputs",
                    "qc_notes": f"No Stage 3B PDBs found in {output_pdb_dir}",
                }
            )
            continue

        for relaxed_pdb in relaxed_pdbs:
            all_rows.append(
                _qc_row_for_relaxed_pdb(
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
            )

    pass_rows = [row for row in all_rows if row.get("pass_stage3c_qc") == "true"]
    write_csv(output_dir / f"FGA_rfpeptides_stage3{output_name_suffix}_sequences_qc.csv", all_rows, STAGE3C_FIELDS)
    write_csv(output_dir / f"FGA_rfpeptides_stage3{output_name_suffix}_sequences_qc_pass.csv", pass_rows, STAGE3C_FIELDS)
    _write_fastas(all_rows, fasta_dir)
    write_markdown(
        output_dir / f"FGA_rfpeptides_stage3{output_name_suffix}_sequences_qc.md",
        _summary_markdown(
            rows=all_rows,
            args=args,
            output_dir=output_dir,
            site_numbers=last_site_numbers,
            hotspot_numbers=last_hotspot_numbers,
        ),
    )
    _write_pymol_review(
        rows=all_rows,
        output_path=output_dir / f"RFpep_Site_2_stage3C{output_name_suffix}_sequence_qc_review.pml",
        target_chain=last_target_chain,
        peptide_chain=last_peptide_chain,
        site_numbers=last_site_numbers,
        hotspot_numbers=last_hotspot_numbers,
        top_n=args.top_pymol,
        pymol_path_style=args.pymol_path_style,
    )

    logger.info("Parsed Stage 3B relaxed PDBs: %s", len(all_rows))
    logger.info("Passed Stage 3C sequence/relax QC: %s", len(pass_rows))
    logger.info("Output directory: %s", output_dir)
    if not pass_rows:
        logger.warning("No Stage 3B outputs passed Stage 3C QC.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
