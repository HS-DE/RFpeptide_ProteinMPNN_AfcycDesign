from __future__ import annotations

import argparse
import math
import os
from collections import Counter
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

from common import append_run_header, read_csv, resolve_path, rows_to_markdown, setup_logger, write_csv, write_markdown
from pdb_utils import ca_coord, centroid, distance, parse_residues


QC_FIELDS = [
    "design_id",
    "site_label",
    "site_id",
    "rf_pdb",
    "trb",
    "file_status",
    "parse_status",
    "target_chain",
    "peptide_chain",
    "peptide_length",
    "requested_length_min",
    "requested_length_max",
    "target_residue_count",
    "target_position_mapping_status",
    "target_pdb_residue_number_start",
    "target_pdb_residue_number_end",
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
    "pass_backbone_qc",
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
        drive = text[5].upper()
        return Path(f"{drive}:/{text[7:]}")
    if os.name != "nt" and len(text) >= 3 and text[1] == ":" and text[2] == "/":
        drive = text[0].lower()
        return Path(f"/mnt/{drive}{text[2:]}")
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


def _read_required_csv(path: Path) -> list[dict[str, str]]:
    rows = read_csv(path)
    if not rows:
        raise RuntimeError(f"Missing or empty CSV: {path}")
    return rows


def _stage0_lookup(rows: Iterable[Mapping[str, str]]) -> dict[str, dict[str, str]]:
    lookup: dict[str, dict[str, str]] = {}
    for row in rows:
        item = dict(row)
        for key in [item.get("site_label", ""), item.get("site_id", "")]:
            if key:
                lookup[key] = item
    return lookup


def _stage1_job_lookup(rows: Iterable[Mapping[str, str]]) -> dict[str, dict[str, str]]:
    lookup: dict[str, dict[str, str]] = {}
    for row in rows:
        item = dict(row)
        for key in [item.get("site_label", ""), item.get("site_id", ""), item.get("rfpeptides_job_id", "")]:
            if key:
                lookup[key] = item
    return lookup


def _parse_rf_range(value: str) -> tuple[str, int, int]:
    text = str(value).strip()
    if "-" not in text:
        raise RuntimeError(f"Invalid rfpeptides_residue_range: {value}")
    left, right = text.split("-", 1)
    if not left or not right or left[0] != right[0]:
        raise RuntimeError(f"Invalid rfpeptides_residue_range: {value}")
    try:
        start = int(left[1:])
        end = int(right[1:])
    except ValueError as exc:
        raise RuntimeError(f"Invalid rfpeptides_residue_range: {value}") from exc
    return left[0], start, end


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


def _residue_label(chain_id: str, residue: Mapping[str, Any]) -> str:
    return f"{chain_id}{residue.get('pdb_residue_number', '')}"


def _residues_at_crop_positions(
    target_residues: Sequence[Mapping[str, Any]],
    crop_positions: Iterable[str],
) -> list[Mapping[str, Any]]:
    selected = []
    for value in sorted(crop_positions, key=_residue_sort_key):
        position = _parse_int(value)
        if position < 1 or position > len(target_residues):
            raise RuntimeError(
                f"Crop position {value} is outside target-chain length {len(target_residues)}"
            )
        selected.append(target_residues[position - 1])
    return selected


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


def _macrocycle_geometry_status(peptide_residues: Sequence[Mapping[str, Any]]) -> tuple[str, float | None]:
    if not peptide_residues:
        return "fail_parse_error", None
    first_atoms = peptide_residues[0].get("atoms", {})
    last_atoms = peptide_residues[-1].get("atoms", {})
    n_atom = first_atoms.get("N")
    c_atom = last_atoms.get("C")
    if n_atom is None or c_atom is None:
        return "warn_chain_or_residue_numbering_unclear", None
    cn_distance = distance(n_atom, c_atom)
    if cn_distance <= 2.0:
        return "pass_head_to_tail_macrocycle", cn_distance
    if cn_distance <= 3.0:
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


def _find_design_files(output_prefix: Path, num_designs: int) -> list[tuple[int, Path, Path]]:
    files: list[tuple[int, Path, Path]] = []
    for design_index in range(num_designs):
        pdb = Path(f"{output_prefix}_{design_index}.pdb")
        trb = Path(f"{output_prefix}_{design_index}.trb")
        files.append((design_index, pdb, trb))
    return files


def _chain_candidates(
    chains: Mapping[str, Sequence[Mapping[str, Any]]],
    target_chain: str,
    length_min: int,
    length_max: int,
) -> list[tuple[str, Sequence[Mapping[str, Any]]]]:
    candidates = []
    for chain_id, residues in chains.items():
        if chain_id == target_chain:
            continue
        if length_min <= len(residues) <= length_max:
            candidates.append((chain_id, residues))
    return candidates


def _choose_peptide_chain(
    chains: Mapping[str, Sequence[Mapping[str, Any]]],
    target_chain: str,
    length_min: int,
    length_max: int,
) -> tuple[str, Sequence[Mapping[str, Any]], str]:
    candidates = _chain_candidates(chains, target_chain, length_min, length_max)
    if not candidates:
        return "", [], "fail_no_designed_chain_in_length_range"
    if len(candidates) == 1:
        return candidates[0][0], candidates[0][1], "pass"
    target_coords = _all_atom_coords(chains.get(target_chain, []))
    ranked = []
    for chain_id, residues in candidates:
        peptide_coords = _all_atom_coords(residues)
        min_dist = min((distance(a, b) for a in target_coords for b in peptide_coords), default=999.0)
        ranked.append((min_dist, chain_id, residues))
    ranked.sort(key=lambda item: (item[0], item[1]))
    return ranked[0][1], ranked[0][2], "warn_multiple_candidate_peptide_chains"


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


def _qc_row_for_design(
    *,
    design_index: int,
    pdb_path: Path,
    trb_path: Path,
    site_label: str,
    site_id: str,
    target_chain: str,
    length_min: int,
    length_max: int,
    site_residue_numbers: set[str],
    hotspot_residue_numbers: set[str],
    contact_cutoff: float,
    site_near_distance: float,
    hotspot_near_distance: float,
    severe_clash_distance: float,
    min_target_contacts: int,
    min_site_contacts: int,
    min_hotspot_contacts: int,
    allow_near_pass: bool,
) -> dict[str, Any]:
    design_id = f"{site_label}_{design_index:04d}"
    file_missing = []
    if not pdb_path.exists() or pdb_path.stat().st_size == 0:
        file_missing.append("pdb")
    if not trb_path.exists() or trb_path.stat().st_size == 0:
        file_missing.append("trb")
    file_status = "pass" if not file_missing else "fail_missing_" + "_".join(file_missing)
    base_row: dict[str, Any] = {
        "design_id": design_id,
        "site_label": site_label,
        "site_id": site_id,
        "rf_pdb": pdb_path,
        "trb": trb_path,
        "file_status": file_status,
        "requested_length_min": length_min,
        "requested_length_max": length_max,
        "target_chain": target_chain,
        "pass_backbone_qc": "false",
    }
    if file_status != "pass":
        base_row.update(
            {
                "parse_status": "not_parsed",
                "qc_failure_reasons": file_status,
                "qc_notes": "Output PDB/TRB pair is incomplete.",
            }
        )
        return base_row

    try:
        chains = parse_residues(pdb_path)
    except Exception as exc:  # pragma: no cover - malformed external PDBs
        base_row.update(
            {
                "parse_status": "fail_parse_error",
                "qc_failure_reasons": f"parse_error:{exc.__class__.__name__}",
                "qc_notes": str(exc),
            }
        )
        return base_row

    target_residues = list(chains.get(target_chain, []))
    if not target_residues:
        base_row.update(
            {
                "parse_status": "fail_missing_target_chain",
                "qc_failure_reasons": "missing_target_chain",
                "qc_notes": f"Target chain {target_chain} was not found in RFpeptides output.",
            }
        )
        return base_row

    peptide_chain, peptide_residues, peptide_parse_status = _choose_peptide_chain(chains, target_chain, length_min, length_max)
    if not peptide_residues:
        base_row.update(
            {
                "parse_status": peptide_parse_status,
                "target_residue_count": len(target_residues),
                "peptide_chain": peptide_chain,
                "peptide_length": 0,
                "qc_failure_reasons": "missing_designed_chain",
                "qc_notes": f"No non-target chain with length {length_min}-{length_max} was found.",
            }
        )
        return base_row

    try:
        site_residues = _residues_at_crop_positions(target_residues, site_residue_numbers)
        hotspot_residues = _residues_at_crop_positions(target_residues, hotspot_residue_numbers)
    except RuntimeError as exc:
        base_row.update(
            {
                "parse_status": "fail_target_position_mapping",
                "target_residue_count": len(target_residues),
                "peptide_chain": peptide_chain,
                "peptide_length": len(peptide_residues),
                "target_position_mapping_status": "failed",
                "qc_failure_reasons": "target_crop_position_mapping_failed",
                "qc_notes": str(exc),
            }
        )
        return base_row

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
        near_name="site_near_pass",
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
        near_name="hotspot_near_pass",
        miss_name="hotspot_missed",
        near_distance=hotspot_near_distance,
        min_contacts=min_hotspot_contacts,
    )

    macrocycle_status, cn_distance = _macrocycle_geometry_status(peptide_residues)
    clash_status = "not_evaluated"
    if target_min is not None:
        clash_status = "fail_severe_clash" if target_min < severe_clash_distance else "pass_no_severe_clash"

    center, rog = _peptide_geometry_summary(peptide_residues)
    length_status = "pass_length" if length_min <= len(peptide_residues) <= length_max else "fail_length_out_of_range"
    failure_reasons = []
    if peptide_parse_status.startswith("fail"):
        failure_reasons.append(peptide_parse_status)
    if length_status.startswith("fail"):
        failure_reasons.append(length_status)
    if macrocycle_status.startswith("fail"):
        failure_reasons.append(macrocycle_status)
    target_ok = target_contact_status in {"target_contact_pass", "target_contact_low_count"}
    site_ok = site_status == "site_contact_pass" or (allow_near_pass and site_status == "site_near_pass")
    hotspot_ok = hotspot_status == "hotspot_contact_pass" or (
        allow_near_pass and hotspot_status == "hotspot_near_pass"
    )
    if not target_ok:
        failure_reasons.append(target_contact_status)
    if not site_ok:
        failure_reasons.append(site_status)
    if not hotspot_ok:
        failure_reasons.append(hotspot_status)
    if clash_status.startswith("fail"):
        failure_reasons.append(clash_status)

    pass_qc = not failure_reasons
    notes = []
    if peptide_parse_status.startswith("warn"):
        notes.append(peptide_parse_status)
    if macrocycle_status.startswith("warn"):
        notes.append(macrocycle_status)
    if target_contact_status == "target_contact_low_count":
        notes.append("Target crop contact count is below the strict count threshold, but direct target-site/hotspot recovery is evaluated separately.")
    if site_status == "crop_only_contact":
        notes.append("Peptide contacts the crop but misses RFpep_Site_2 residues.")
    if hotspot_status == "hotspot_missed":
        notes.append("Peptide does not recover proximity to selected hotspots.")

    base_row.update(
        {
            "parse_status": peptide_parse_status,
            "peptide_chain": peptide_chain,
            "peptide_length": len(peptide_residues),
            "target_residue_count": len(target_residues),
            "target_position_mapping_status": "mapped_crop_positions_by_target_chain_order",
            "target_pdb_residue_number_start": target_residues[0].get("pdb_residue_number", ""),
            "target_pdb_residue_number_end": target_residues[-1].get("pdb_residue_number", ""),
            "mapped_target_site_pdb_residue_numbers": ",".join(
                str(residue.get("pdb_residue_number", "")) for residue in site_residues
            ),
            "mapped_hotspot_pdb_residue_numbers": ",".join(
                str(residue.get("pdb_residue_number", "")) for residue in hotspot_residues
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
            "pass_backbone_qc": "true" if pass_qc else "false",
            "qc_failure_reasons": ";".join(failure_reasons),
            "qc_notes": "; ".join(notes),
        }
    )
    return base_row


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
    site_label: str,
    site_numbers: set[str],
    hotspot_numbers: set[str],
) -> str:
    pass_rows = [row for row in rows if row.get("pass_backbone_qc") == "true"]
    crop_only = [row for row in rows if row.get("target_site_recovery_status") == "crop_only_contact"]
    top_rows = sorted(
        pass_rows,
        key=lambda row: (
            -_parse_int(row.get("num_hotspot_contacts", 0)),
            -_parse_int(row.get("num_target_site_contacts", 0)),
            float(row.get("peptide_hotspot_min_distance") or 999.0),
            float(row.get("peptide_site_min_distance") or 999.0),
            str(row.get("design_id", "")),
        ),
    )[: args.top_report]
    columns = [
        "design_id",
        "peptide_chain",
        "peptide_length",
        "num_target_site_contacts",
        "num_hotspot_contacts",
        "peptide_site_min_distance",
        "peptide_hotspot_min_distance",
        "macrocycle_geometry_status",
        "pass_backbone_qc",
    ]
    return f"""# FGA RFpeptides Stage 2 Backbone QC

Status: RFpeptides backbone outputs parsed and checked before sequence design.

Important rule: Stage 2 does not treat whole-crop contact as sufficient.
A design must make direct contacts within `contact_cutoff_A` to the intended
`{site_label}` target-site residues and selected hotspots. Near-only states are
reported but fail by default. Designs that contact only a distant part of the
target crop are flagged as `crop_only_contact`.

Output directory:

```text
{output_dir}
```

Parameters:

```text
stage0_root: {args.stage0_root}
stage1_root: {args.stage1_root}
output_subdir: {args.output_subdir}
selected_sites: {args.selected_sites}
contact_cutoff_A: {args.contact_cutoff}
site_near_distance_A: {args.site_near_distance}
hotspot_near_distance_A: {args.hotspot_near_distance}
severe_clash_distance_A: {args.severe_clash_distance}
min_target_contacts: {args.min_target_contacts}
min_site_contacts: {args.min_site_contacts}
min_hotspot_contacts: {args.min_hotspot_contacts}
allow_near_pass: {args.allow_near_pass}
```

Target-site residue numbers:

```text
{",".join(sorted(site_numbers, key=_residue_sort_key))}
```

These are Stage 0 crop sequence positions, not raw residue numbers in an
RFpeptides output PDB. Each design maps them by target-chain residue order;
the mapped output-PDB residue numbers are recorded per row.

Hotspot residue numbers:

```text
{",".join(sorted(hotspot_numbers, key=_residue_sort_key))}
```

## Counts

```text
total_backbones: {len(rows)}
pass_backbone_qc: {len(pass_rows)}
crop_only_contact: {len(crop_only)}
```

## Target-Site Recovery Status

{_status_lines(rows, "target_site_recovery_status")}

## Hotspot Recovery Status

{_status_lines(rows, "hotspot_recovery_status")}

## Macrocycle Geometry Status

{_status_lines(rows, "macrocycle_geometry_status")}

## Top Passing Backbones

{rows_to_markdown(top_rows, columns, "No backbones passed Stage 2 QC.")}
"""


def _write_pymol_review(
    *,
    rows: list[Mapping[str, Any]],
    output_path: Path,
    target_chain: str,
    peptide_chain_field: str,
    site_numbers: set[str],
    hotspot_numbers: set[str],
    top_n: int,
    pymol_path_style: str,
) -> None:
    pass_rows = [row for row in rows if row.get("pass_backbone_qc") == "true"]
    selected = sorted(
        pass_rows,
        key=lambda row: (
            -_parse_int(row.get("num_hotspot_contacts", 0)),
            -_parse_int(row.get("num_target_site_contacts", 0)),
            float(row.get("peptide_hotspot_min_distance") or 999.0),
            str(row.get("design_id", "")),
        ),
    )[:top_n]
    lines = [
        "reinitialize",
        "set retain_order, 1",
        "hide everything, all",
    ]
    for idx, row in enumerate(selected, start=1):
        obj = _safe_token(str(row.get("design_id", f"design_{idx}")))
        pdb_path = _format_pymol_path(str(row.get("rf_pdb", "")), pymol_path_style)
        peptide_chain = str(row.get(peptide_chain_field, ""))
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
                f"color cyan, {obj} and chain {peptide_chain}",
                f"set sphere_scale, 0.45, {obj}_hotspots",
                f"label {obj}_hotspots and name CA, \"%s%s\" % (chain, resi)",
            ]
        )
    if selected:
        lines.append(f"zoom {_safe_token(str(selected[0].get('design_id', 'design_1')))}_site, 14")
    write_markdown(output_path, "\n".join(lines))


def main() -> int:
    parser = argparse.ArgumentParser(description="Collect and QC RFpeptides Stage 1 backbone outputs before sequence design.")
    parser.add_argument("--stage0-root", default="results/rfpeptides_article_route_clean_20260615_fpocket")
    parser.add_argument("--stage1-root", default="results/rfpeptides_article_route_clean_20260615_fpocket_stage1_N1000_no_traj")
    parser.add_argument("--output-root", default="", help="Defaults to --stage1-root.")
    parser.add_argument("--output-subdir", default="03_backbone_qc")
    parser.add_argument("--selected-sites", default="RFpep_Site_2")
    parser.add_argument("--stage0-summary-csv", default="")
    parser.add_argument("--stage1-jobs-csv", default="")
    parser.add_argument("--contact-cutoff", type=float, default=5.0)
    parser.add_argument("--site-near-distance", type=float, default=6.0)
    parser.add_argument("--hotspot-near-distance", type=float, default=8.0)
    parser.add_argument("--severe-clash-distance", type=float, default=1.2)
    parser.add_argument("--min-target-contacts", type=int, default=3)
    parser.add_argument("--min-site-contacts", type=int, default=1)
    parser.add_argument("--min-hotspot-contacts", type=int, default=1)
    parser.add_argument(
        "--allow-near-pass",
        action="store_true",
        help="Compatibility option only. By default, near-only site/hotspot proximity fails Stage 2.",
    )
    parser.add_argument("--top-report", type=int, default=50)
    parser.add_argument("--top-pymol", type=int, default=20)
    parser.add_argument(
        "--pymol-path-style",
        choices=["windows", "wsl", "native"],
        default="windows",
        help="Path style for generated PyMOL review scripts. Use windows for Windows PyMOL.",
    )
    args = parser.parse_args()

    logger = setup_logger("21_collect_rfpeptides_backbones")
    append_run_header(logger, "21_collect_rfpeptides_backbones.py")

    if args.contact_cutoff <= 0:
        raise RuntimeError("--contact-cutoff must be > 0")
    if args.site_near_distance < args.contact_cutoff:
        raise RuntimeError("--site-near-distance must be >= --contact-cutoff")
    if args.hotspot_near_distance < args.contact_cutoff:
        raise RuntimeError("--hotspot-near-distance must be >= --contact-cutoff")
    if args.severe_clash_distance <= 0:
        raise RuntimeError("--severe-clash-distance must be > 0")

    stage0_root = _resolve_mixed_path(args.stage0_root)
    stage1_root = _resolve_mixed_path(args.stage1_root)
    output_root = _resolve_mixed_path(args.output_root) if args.output_root else stage1_root
    output_dir = output_root / args.output_subdir
    stage0_summary_csv = (
        _resolve_mixed_path(args.stage0_summary_csv)
        if args.stage0_summary_csv
        else stage0_root / "00_target_inputs" / "FGA_rfpeptides_stage0_target_inputs_summary.csv"
    )
    stage1_jobs_csv = (
        _resolve_mixed_path(args.stage1_jobs_csv)
        if args.stage1_jobs_csv
        else stage1_root / "01_rfpeptides_jobs" / "FGA_rfpeptides_stage1_jobs.csv"
    )

    stage0 = _stage0_lookup(_read_required_csv(stage0_summary_csv))
    stage1_jobs = _stage1_job_lookup(_read_required_csv(stage1_jobs_csv))
    selected_sites = _split_csv(args.selected_sites)
    if selected_sites != ["RFpep_Site_2"]:
        raise RuntimeError("Current Stage 2 script is scoped to RFpep_Site_2 only.")

    all_rows: list[dict[str, Any]] = []
    last_site_numbers: set[str] = set()
    last_hotspot_numbers: set[str] = set()
    last_target_chain = ""
    for selected_id in selected_sites:
        site_row = stage0.get(selected_id)
        job_row = stage1_jobs.get(selected_id)
        if site_row is None:
            raise RuntimeError(f"Selected site not found in Stage 0 summary: {selected_id}")
        if job_row is None:
            raise RuntimeError(f"Selected site not found in Stage 1 jobs table: {selected_id}")

        site_label = str(site_row.get("site_label") or selected_id)
        site_id = str(site_row.get("site_id", ""))
        target_chain, _, _ = _parse_rf_range(str(site_row.get("rfpeptides_residue_range", "")))
        last_target_chain = target_chain
        mapping_csv = _resolve_mixed_path(str(site_row.get("crop_renumbering_mapping_csv", "")))
        site_numbers, hotspot_numbers = _load_site_mapping(mapping_csv)
        last_site_numbers = site_numbers
        last_hotspot_numbers = hotspot_numbers

        output_prefix = _resolve_mixed_path(str(job_row.get("output_prefix", "")))
        num_designs = _parse_int(job_row.get("num_designs", "0"))
        length_min = _parse_int(job_row.get("length_min", "0"))
        length_max = _parse_int(job_row.get("length_max", "0"))
        if num_designs <= 0:
            raise RuntimeError(f"Invalid num_designs in Stage 1 jobs table for {site_label}: {num_designs}")
        if length_min <= 0 or length_max < length_min:
            raise RuntimeError(f"Invalid length range in Stage 1 jobs table for {site_label}: {length_min}-{length_max}")

        for design_index, pdb_path, trb_path in _find_design_files(output_prefix, num_designs):
            all_rows.append(
                _qc_row_for_design(
                    design_index=design_index,
                    pdb_path=pdb_path,
                    trb_path=trb_path,
                    site_label=site_label,
                    site_id=site_id,
                    target_chain=target_chain,
                    length_min=length_min,
                    length_max=length_max,
                    site_residue_numbers=site_numbers,
                    hotspot_residue_numbers=hotspot_numbers,
                    contact_cutoff=args.contact_cutoff,
                    site_near_distance=args.site_near_distance,
                    hotspot_near_distance=args.hotspot_near_distance,
                    severe_clash_distance=args.severe_clash_distance,
                    min_target_contacts=args.min_target_contacts,
                    min_site_contacts=args.min_site_contacts,
                    min_hotspot_contacts=args.min_hotspot_contacts,
                    allow_near_pass=args.allow_near_pass,
                )
            )

    pass_rows = [row for row in all_rows if row.get("pass_backbone_qc") == "true"]
    write_csv(output_dir / "FGA_rfpeptides_backbones_qc.csv", all_rows, QC_FIELDS)
    write_csv(output_dir / "FGA_rfpeptides_backbones_qc_pass.csv", pass_rows, QC_FIELDS)
    write_markdown(
        output_dir / "FGA_rfpeptides_backbones_qc.md",
        _summary_markdown(
            rows=all_rows,
            args=args,
            output_dir=output_dir,
            site_label="RFpep_Site_2",
            site_numbers=last_site_numbers,
            hotspot_numbers=last_hotspot_numbers,
        ),
    )
    _write_pymol_review(
        rows=all_rows,
        output_path=output_dir / "RFpep_Site_2_stage2_top_pass_review.pml",
        target_chain=last_target_chain,
        peptide_chain_field="peptide_chain",
        site_numbers=last_site_numbers,
        hotspot_numbers=last_hotspot_numbers,
        top_n=args.top_pymol,
        pymol_path_style=args.pymol_path_style,
    )

    logger.info("Parsed Stage 1 backbones: %s", len(all_rows))
    logger.info("Passed Stage 2 backbone QC: %s", len(pass_rows))
    logger.info("Output directory: %s", output_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
