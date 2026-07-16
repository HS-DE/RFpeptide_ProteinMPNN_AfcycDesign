from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path
import re
import shlex
import shutil
import subprocess
from typing import Any, Iterable, Mapping

from common import append_run_header, load_config, read_csv, resolve_path, rows_to_markdown, setup_logger, write_csv, write_markdown
from pdb_utils import ca_coord, centroid, distance, parse_residues


HYDROPHOBIC_RESNAMES = {"ALA", "VAL", "ILE", "LEU", "MET", "PRO", "PHE", "TRP", "TYR"}
AROMATIC_RESNAMES = {"PHE", "TYR", "TRP", "HIS"}
CHARGED_RESNAMES = {"ASP", "GLU", "LYS", "ARG", "HIS"}
STAGE0_ALLOWED_TIERS = {"high", "medium"}
MAX_ACCESSIBLE_SURFACE_AREA = {
    "ALA": 129.0,
    "ARG": 274.0,
    "ASN": 195.0,
    "ASP": 193.0,
    "CYS": 167.0,
    "GLN": 225.0,
    "GLU": 223.0,
    "GLY": 104.0,
    "HIS": 224.0,
    "ILE": 197.0,
    "LEU": 201.0,
    "LYS": 236.0,
    "MET": 224.0,
    "PHE": 240.0,
    "PRO": 159.0,
    "SER": 155.0,
    "THR": 172.0,
    "TRP": 285.0,
    "TYR": 263.0,
    "VAL": 174.0,
}
FPOCKET_FIELDS = [
    "fpocket_status",
    "fpocket_error_summary",
    "nearest_fpocket_pocket_id",
    "nearest_fpocket_distance_to_hotspot",
    "nearest_fpocket_distance_to_site",
    "nearest_fpocket_pocket_score",
    "nearest_fpocket_druggability_score",
    "nearest_fpocket_volume",
    "fpocket_support_status",
    "fpocket_support_reason",
]


def _float_arg(value: str, name: str, low: float | None = None, high: float | None = None) -> float:
    out = float(value)
    if low is not None and out < low:
        raise argparse.ArgumentTypeError(f"{name} must be >= {low}")
    if high is not None and out > high:
        raise argparse.ArgumentTypeError(f"{name} must be <= {high}")
    return out


def _fraction_arg(value: str, name: str) -> float:
    out = float(value)
    if out > 1.0:
        out /= 100.0
    if out < 0.0 or out > 1.0:
        raise argparse.ArgumentTypeError(f"{name} must be a fraction from 0 to 1, or a percent from 0 to 100")
    return out


def _percentile(values: list[float], fraction: float) -> float:
    if not values:
        return 0.0
    ordered = sorted(values)
    idx = min(len(ordered) - 1, max(0, int(round((len(ordered) - 1) * fraction))))
    return ordered[idx]


def _parse_int(value: str, default: int = 0) -> int:
    try:
        return int(str(value).strip())
    except ValueError:
        return default


def _split_csv(value: str) -> list[str]:
    return [item.strip() for item in str(value or "").split(",") if item.strip()]


def _unique_preserve(values: Iterable[str]) -> list[str]:
    seen: set[str] = set()
    out: list[str] = []
    for value in values:
        if value not in seen:
            seen.add(value)
            out.append(value)
    return out


def _residue_label(row: Mapping[str, Any]) -> str:
    return f"{row['chain_id']}{row['pdb_residue_number']}"


def _safe_mean(values: Iterable[float]) -> float:
    vals = list(values)
    return sum(vals) / len(vals) if vals else 0.0


def _row_coord(row: Mapping[str, Any]) -> tuple[float, float, float]:
    return (float(row["x"]), float(row["y"]), float(row["z"]))


def _resname(row: Mapping[str, Any]) -> str:
    return str(row.get("pdb_residue_name", "")).strip().upper()


def _calculate_rsa(resname: str, sasa_abs: float | None) -> float | None:
    if sasa_abs is None:
        return None
    max_asa = MAX_ACCESSIBLE_SURFACE_AREA.get(resname)
    if not max_asa:
        return None
    return sasa_abs / max_asa


def _calculate_sasa_by_residue(pdb_path: str | Path) -> tuple[str, dict[tuple[str, str], float]]:
    try:
        import freesasa  # type: ignore
    except ModuleNotFoundError:
        return "proxy_only", {}

    try:
        structure = freesasa.Structure(str(resolve_path(pdb_path)))
        result = freesasa.calc(structure)
        residue_areas = result.residueAreas()
    except Exception as exc:  # pragma: no cover - depends on optional native package
        return f"proxy_only:freesasa_error:{exc.__class__.__name__}", {}

    lookup: dict[tuple[str, str], float] = {}
    for chain_id, residues in residue_areas.items():
        chain_key = str(chain_id).strip() or "_"
        for residue_number, area in residues.items():
            total = getattr(area, "total", None)
            if total is None:
                continue
            lookup[(chain_key, str(residue_number).strip())] = float(total)
    if not lookup:
        return "proxy_only:freesasa_empty_result", {}
    return "freesasa", lookup


def _lookup_sasa(row: Mapping[str, Any], sasa_lookup: Mapping[tuple[str, str], float]) -> float | None:
    chain_id = str(row.get("chain_id", "")).strip() or "_"
    residue_number = str(row.get("pdb_residue_number", "")).strip()
    if (chain_id, residue_number) in sasa_lookup:
        return sasa_lookup[(chain_id, residue_number)]
    residue_number_no_icode = "".join(ch for ch in residue_number if ch.isdigit() or ch == "-")
    if residue_number_no_icode and (chain_id, residue_number_no_icode) in sasa_lookup:
        return sasa_lookup[(chain_id, residue_number_no_icode)]
    return None


def _is_surface_for_quality(row: Mapping[str, Any], sasa_status: str, rsa_threshold: float) -> bool:
    if sasa_status == "freesasa" and str(row.get("rsa", "")) not in {"", "None"}:
        return float(row["rsa"]) >= rsa_threshold
    value = row.get("is_surface_candidate", False)
    if isinstance(value, bool):
        return value
    return str(value).lower() == "true"


def _quality_status(values: Iterable[str]) -> str:
    unique = _unique_preserve(str(value) for value in values if str(value))
    if not unique:
        return "unknown"
    if len(unique) == 1:
        return unique[0]
    return "mixed:" + ",".join(unique)


def _numeric_values(rows: Iterable[Mapping[str, Any]], field: str) -> list[float]:
    values: list[float] = []
    for row in rows:
        value = row.get(field, "")
        if str(value) in {"", "None"}:
            continue
        values.append(float(value))
    return values


def _short_error_summary(text: str, limit: int = 500) -> str:
    summary = " ".join(str(text or "").split())
    if len(summary) <= limit:
        return summary
    return summary[: limit - 3] + "..."


def _fpocket_command_prefix(fpocket_bin: str) -> list[str]:
    return shlex.split(str(fpocket_bin))


def _resolve_fpocket_command(fpocket_bin: str) -> list[str] | None:
    prefix = _fpocket_command_prefix(fpocket_bin)
    if not prefix:
        return None
    executable = prefix[0]
    exe_path = Path(executable)
    if exe_path.exists():
        return [str(exe_path)] + prefix[1:]
    resolved = shutil.which(executable)
    if resolved is None:
        return None
    return [resolved] + prefix[1:]


def _parse_atom_coords(path: Path) -> list[tuple[float, float, float]]:
    coords: list[tuple[float, float, float]] = []
    if not path.exists():
        return coords
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            try:
                coords.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
            except ValueError:
                parts = line.split()
                if len(parts) >= 9:
                    try:
                        coords.append((float(parts[6]), float(parts[7]), float(parts[8])))
                    except ValueError:
                        continue
    return coords


def _parse_pocket_id(value: str) -> str | None:
    match = re.search(r"pocket[\s_]*(\d+)", value, flags=re.IGNORECASE)
    if not match:
        return None
    return f"pocket{int(match.group(1))}"


def _first_float(value: str) -> float | None:
    match = re.search(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?", value)
    if not match:
        return None
    return float(match.group(0))


def _parse_fpocket_info(info_path: Path) -> dict[str, dict[str, float]]:
    pocket_info: dict[str, dict[str, float]] = defaultdict(dict)
    current_id: str | None = None
    if not info_path.exists():
        return pocket_info
    with info_path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            pocket_id = _parse_pocket_id(line)
            if pocket_id is not None and line.lower().startswith("pocket"):
                current_id = pocket_id
                pocket_info.setdefault(current_id, {})
                continue
            if current_id is None or ":" not in line:
                continue
            key, raw_value = line.split(":", 1)
            key_norm = key.strip().lower()
            number = _first_float(raw_value)
            if number is None:
                continue
            if key_norm in {"score", "pocket score"}:
                pocket_info[current_id]["pocket_score"] = number
            elif "druggability" in key_norm and "score" in key_norm:
                pocket_info[current_id]["druggability_score"] = number
            elif "volume" in key_norm and "volume" not in pocket_info[current_id]:
                pocket_info[current_id]["volume"] = number
    return pocket_info


def _parse_fpocket_pockets(fpocket_dir: Path) -> list[dict[str, Any]]:
    info: dict[str, dict[str, float]] = {}
    for info_path in sorted(fpocket_dir.glob("*_out/*_info.txt")):
        info.update(_parse_fpocket_info(info_path))

    coords_by_pocket: dict[str, list[tuple[float, float, float]]] = defaultdict(list)
    for coord_path in sorted(fpocket_dir.glob("*_out/pockets/pocket*_atm.pdb")):
        pocket_id = _parse_pocket_id(coord_path.name)
        if pocket_id is not None:
            coords_by_pocket[pocket_id].extend(_parse_atom_coords(coord_path))
    for coord_path in sorted(fpocket_dir.glob("*_out/pockets/pocket*_vert.pqr")):
        pocket_id = _parse_pocket_id(coord_path.name)
        if pocket_id is not None and not coords_by_pocket.get(pocket_id):
            coords_by_pocket[pocket_id].extend(_parse_atom_coords(coord_path))

    pocket_ids = sorted(set(info) | set(coords_by_pocket), key=lambda value: _parse_int(value.replace("pocket", ""), 0))
    pockets: list[dict[str, Any]] = []
    for pocket_id in pocket_ids:
        pocket_info = info.get(pocket_id, {})
        pockets.append(
            {
                "pocket_id": pocket_id,
                "coords": coords_by_pocket.get(pocket_id, []),
                "pocket_score": pocket_info.get("pocket_score", ""),
                "druggability_score": pocket_info.get("druggability_score", ""),
                "volume": pocket_info.get("volume", ""),
            }
        )
    return pockets


def _write_fpocket_status(fpocket_dir: Path, status: str, message: str) -> None:
    write_markdown(fpocket_dir / "fpocket_status.txt", f"status: {status}\nmessage: {message}\n")


def _run_fpocket_cross_check(
    *,
    enabled: bool,
    cleaned_pdb: str | Path,
    fpocket_dir: Path,
    fpocket_bin: str,
) -> dict[str, Any]:
    if not enabled:
        return {"status": "not_run", "pockets": [], "error_summary": ""}

    fpocket_dir.mkdir(parents=True, exist_ok=True)
    command_prefix = _resolve_fpocket_command(fpocket_bin)
    if command_prefix is None:
        message = f"fpocket executable not found: {fpocket_bin}"
        _write_fpocket_status(fpocket_dir, "not_available", message)
        return {"status": "not_available", "pockets": [], "error_summary": message}

    source_pdb = resolve_path(cleaned_pdb)
    native_input = fpocket_dir / source_pdb.name
    try:
        shutil.copy2(source_pdb, native_input)
    except OSError as exc:
        message = f"could not copy native cleaned PDB for fpocket: {exc.__class__.__name__}: {exc}"
        _write_fpocket_status(fpocket_dir, "failed", message)
        return {"status": "failed", "pockets": [], "error_summary": _short_error_summary(message)}

    stale_output = fpocket_dir / f"{native_input.stem}_out"
    try:
        if stale_output.exists() and fpocket_dir.resolve() in stale_output.resolve().parents:
            shutil.rmtree(stale_output)
    except OSError:
        pass

    command = command_prefix + ["-f", str(native_input)]
    write_markdown(fpocket_dir / "fpocket_command.txt", " ".join(command) + "\n")
    try:
        completed = subprocess.run(command, cwd=fpocket_dir, capture_output=True, text=True, check=False)
    except OSError as exc:
        message = f"fpocket launch failed: {exc.__class__.__name__}: {exc}"
        _write_fpocket_status(fpocket_dir, "failed", message)
        return {"status": "failed", "pockets": [], "error_summary": _short_error_summary(message)}

    write_markdown(fpocket_dir / "fpocket_stdout.txt", completed.stdout or "")
    write_markdown(fpocket_dir / "fpocket_stderr.txt", completed.stderr or "")
    if completed.returncode != 0:
        message = _short_error_summary(completed.stderr or completed.stdout or f"return code {completed.returncode}")
        _write_fpocket_status(fpocket_dir, "failed", message)
        return {"status": "failed", "pockets": [], "error_summary": message}

    pockets = _parse_fpocket_pockets(fpocket_dir)
    _write_fpocket_status(fpocket_dir, "completed", f"parsed_pockets: {len(pockets)}")
    return {"status": "completed", "pockets": pockets, "error_summary": ""}


def _min_coord_distance(
    left: Iterable[tuple[float, float, float]],
    right: Iterable[tuple[float, float, float]],
) -> float | None:
    left_coords = list(left)
    right_coords = list(right)
    if not left_coords or not right_coords:
        return None
    return min(distance(a, b) for a in left_coords for b in right_coords)


def _blank_fpocket_fields(status: str, error_summary: str = "", support_reason: str = "") -> dict[str, Any]:
    return {
        "fpocket_status": status,
        "fpocket_error_summary": error_summary,
        "nearest_fpocket_pocket_id": "",
        "nearest_fpocket_distance_to_hotspot": "",
        "nearest_fpocket_distance_to_site": "",
        "nearest_fpocket_pocket_score": "",
        "nearest_fpocket_druggability_score": "",
        "nearest_fpocket_volume": "",
        "fpocket_support_status": "not_evaluated" if status != "completed" else "no_nearby_pocket",
        "fpocket_support_reason": support_reason,
    }


def _fpocket_support_reason(
    pocket_id: str,
    hotspot_distance: float | None,
    site_distance: float | None,
    support_status: str,
) -> str:
    hotspot_text = "NA" if hotspot_distance is None else f"{hotspot_distance:.2f}A"
    site_text = "NA" if site_distance is None else f"{site_distance:.2f}A"
    return (
        f"{support_status}: supplementary fpocket pocket/groove evidence from {pocket_id}; "
        f"nearest_hotspot_distance={hotspot_text}; nearest_site_distance={site_text}. "
        "Not used as a final pass/fail gate."
    )


def _annotate_candidates_with_fpocket(
    candidates: list[dict[str, Any]],
    fpocket_result: Mapping[str, Any],
    strong_distance: float,
    moderate_distance: float,
    weak_distance: float,
) -> None:
    status = str(fpocket_result.get("status", "not_run"))
    error_summary = str(fpocket_result.get("error_summary", ""))
    if status != "completed":
        reason = {
            "not_run": "fpocket cross-check was not requested.",
            "not_available": "fpocket executable was not available; supplementary pocket evidence was not evaluated.",
            "failed": f"fpocket failed; supplementary pocket evidence was not evaluated. {error_summary}",
        }.get(status, "fpocket supplementary pocket evidence was not evaluated.")
        for candidate in candidates:
            candidate.update(_blank_fpocket_fields(status, error_summary, reason))
        return

    pockets = [pocket for pocket in fpocket_result.get("pockets", []) if pocket.get("coords")]
    if not pockets:
        for candidate in candidates:
            candidate.update(
                _blank_fpocket_fields(
                    "completed",
                    "",
                    "fpocket completed, but no pocket coordinates were parsed; no supplementary pocket/groove support assigned.",
                )
            )
        return

    for candidate in candidates:
        site_rows = list(candidate.get("site_rows", []))
        site_coords = [_row_coord(row) for row in site_rows]
        hotspot_labels = set(_split_csv(str(candidate.get("initial_hotspot_candidates", ""))))
        hotspot_coords = [_row_coord(row) for row in site_rows if _residue_label(row) in hotspot_labels]
        pocket_metrics: list[dict[str, Any]] = []
        for pocket in pockets:
            pocket_coords = pocket.get("coords", [])
            site_distance = _min_coord_distance(site_coords, pocket_coords)
            hotspot_distance = _min_coord_distance(hotspot_coords, pocket_coords)
            comparable = [value for value in [hotspot_distance, site_distance] if value is not None]
            if not comparable:
                continue
            metric = min(comparable)
            pocket_metrics.append(
                {
                    "pocket": pocket,
                    "site_distance": site_distance,
                    "hotspot_distance": hotspot_distance,
                    "metric": metric,
                }
            )

        if not pocket_metrics:
            candidate.update(
                _blank_fpocket_fields(
                    "completed",
                    "",
                    "fpocket completed, but no pocket/site distance could be calculated; no supplementary pocket/groove support assigned.",
                )
            )
            continue

        nearest = min(pocket_metrics, key=lambda item: float(item["metric"]))
        support_pool = []
        for item in pocket_metrics:
            hotspot_distance = item["hotspot_distance"]
            site_distance = item["site_distance"]
            if hotspot_distance is not None and hotspot_distance <= strong_distance:
                support_pool.append((0, hotspot_distance, item, "strong_support"))
            elif (
                hotspot_distance is not None
                and hotspot_distance <= moderate_distance
                or site_distance is not None
                and site_distance <= strong_distance
            ):
                support_pool.append((1, min(value for value in [hotspot_distance, site_distance] if value is not None), item, "moderate_support"))
            elif site_distance is not None and site_distance <= weak_distance:
                support_pool.append((2, site_distance, item, "weak_support"))

        if support_pool:
            _, _, support_item, support_status = min(support_pool, key=lambda item: (item[0], item[1]))
            support_pocket = support_item["pocket"]
            support_reason = _fpocket_support_reason(
                str(support_pocket["pocket_id"]),
                support_item["hotspot_distance"],
                support_item["site_distance"],
                support_status,
            )
        else:
            support_status = "no_nearby_pocket"
            support_reason = _fpocket_support_reason(
                str(nearest["pocket"]["pocket_id"]),
                nearest["hotspot_distance"],
                nearest["site_distance"],
                support_status,
            )

        nearest_pocket = nearest["pocket"]
        candidate.update(
            {
                "fpocket_status": "completed",
                "fpocket_error_summary": "",
                "nearest_fpocket_pocket_id": nearest_pocket.get("pocket_id", ""),
                "nearest_fpocket_distance_to_hotspot": round(nearest["hotspot_distance"], 3)
                if nearest["hotspot_distance"] is not None
                else "",
                "nearest_fpocket_distance_to_site": round(nearest["site_distance"], 3) if nearest["site_distance"] is not None else "",
                "nearest_fpocket_pocket_score": nearest_pocket.get("pocket_score", ""),
                "nearest_fpocket_druggability_score": nearest_pocket.get("druggability_score", ""),
                "nearest_fpocket_volume": nearest_pocket.get("volume", ""),
                "fpocket_support_status": support_status,
                "fpocket_support_reason": support_reason,
            }
        )


def _pairwise_distance_stats(rows: list[Mapping[str, Any]]) -> tuple[float, float]:
    if len(rows) < 2:
        return 0.0, 0.0
    distances = []
    for idx, left in enumerate(rows):
        for right in rows[idx + 1 :]:
            distances.append(distance(_row_coord(left), _row_coord(right)))
    return _safe_mean(distances), max(distances) if distances else 0.0


def _other_chain_metrics(
    site_rows: list[Mapping[str, Any]],
    all_context_residues: list[Mapping[str, Any]],
) -> tuple[float, int, int]:
    site_chain_ids = {str(row["chain_id"]) for row in site_rows}
    site_coords = [_row_coord(row) for row in site_rows]
    nearest = float("inf")
    within_6: set[str] = set()
    within_10: set[str] = set()
    for residue in all_context_residues:
        if str(residue["chain_id"]) in site_chain_ids:
            continue
        other_coord = residue["coord"]
        min_dist = min((distance(site_coord, other_coord) for site_coord in site_coords), default=float("inf"))
        if min_dist < nearest:
            nearest = min_dist
        key = f"{residue['chain_id']}:{residue['pdb_residue_number']}"
        if min_dist <= 6.0:
            within_6.add(key)
        if min_dist <= 10.0:
            within_10.add(key)
    if nearest == float("inf"):
        nearest = 999.0
    return nearest, len(within_6), len(within_10)


def _macrocycle_span_status(hotspot_mean: float, hotspot_max: float, hotspot_count: int) -> str:
    if hotspot_count < 2:
        return "insufficient_hotspots"
    if hotspot_max < 5.0:
        return "too_compact"
    if hotspot_max <= 18.0 and hotspot_mean <= 14.0:
        return "good_span"
    if hotspot_max <= 24.0:
        return "broad_but_possible"
    return "too_broad"


def _chemical_anchor_status(
    hydrophobic_count: int,
    aromatic_count: int,
    charged_count: int,
    hotspot_resnames: list[str],
) -> str:
    hotspot_aromatic = sum(1 for resname in hotspot_resnames if resname in AROMATIC_RESNAMES)
    hotspot_charged = sum(1 for resname in hotspot_resnames if resname in CHARGED_RESNAMES)
    hotspot_hydrophobic = sum(1 for resname in hotspot_resnames if resname in HYDROPHOBIC_RESNAMES)
    if (hotspot_aromatic >= 1 and hotspot_charged >= 1) or (hotspot_hydrophobic >= 2 and hotspot_charged >= 1):
        return "strong_anchor_mix"
    if aromatic_count >= 1 or charged_count >= 2 or hydrophobic_count >= 4:
        return "moderate_anchor_features"
    return "weak_anchor_features"


def _occlusion_status(nearest_other_chain_distance: float, within_6: int, within_10: int) -> str:
    if nearest_other_chain_distance < 4.0 or within_6 >= 8:
        return "occluded"
    if nearest_other_chain_distance < 6.0 or within_6 >= 4 or within_10 >= 20:
        return "partly_occluded"
    return "open"


def _site_quality_tier_and_reason(
    surface_fraction: float,
    hotspot_surface_fraction: float,
    macrocycle_span_status: str,
    chemical_anchor_status: str,
    occlusion_status: str,
    nearest_other_chain_distance: float,
    other_chain_residues_within_6a: int,
) -> tuple[str, str]:
    hard_reject_reasons = []
    if surface_fraction < 0.35:
        hard_reject_reasons.append(f"low site surface fraction ({surface_fraction:.2f})")
    if hotspot_surface_fraction < 0.50:
        hard_reject_reasons.append(f"low hotspot surface fraction ({hotspot_surface_fraction:.2f})")
    if macrocycle_span_status in {"insufficient_hotspots", "too_broad"}:
        hard_reject_reasons.append(f"macrocycle span status is {macrocycle_span_status}")
    if occlusion_status == "occluded" and nearest_other_chain_distance < 3.0:
        hard_reject_reasons.append("site is strongly occluded by another chain")
    if hard_reject_reasons:
        return "reject", "Rejected because " + "; ".join(hard_reject_reasons) + "."

    if (
        surface_fraction >= 0.70
        and hotspot_surface_fraction >= 0.75
        and macrocycle_span_status == "good_span"
        and chemical_anchor_status != "weak_anchor_features"
        and occlusion_status in {"open", "partly_occluded"}
    ):
        return (
            "high",
            "High quality: exposed site, surface-accessible hotspots, compatible hotspot span, usable chemical anchors, and no severe occlusion.",
        )

    if (
        surface_fraction >= 0.50
        and hotspot_surface_fraction >= 0.50
        and macrocycle_span_status in {"good_span", "broad_but_possible"}
        and occlusion_status != "occluded"
    ):
        reason = (
            "Medium quality: geometry and exposure are usable for Stage 0, but at least one feature is not ideal "
            f"(chemical_anchor_status={chemical_anchor_status}, occlusion_status={occlusion_status}, "
            f"other_chain_residues_within_6A={other_chain_residues_within_6a})."
        )
        return "medium", reason

    return (
        "low",
        "Low quality: the site is not an automatic reject, but exposure, span, anchor chemistry, or occlusion is weaker than preferred.",
    )


def _legacy_patch_lookup(path: Path) -> dict[str, dict[str, set[str]]]:
    rows = read_csv(path)
    lookup: dict[str, dict[str, set[str]]] = {}
    for row in rows:
        patch_id = row.get("patch_id", "")
        if not patch_id:
            continue
        pdb_labels = set(_split_csv(row.get("hotspot_pdb_residues", "")))
        pdb_labels.update(_split_csv(row.get("representative_hotspot_residues", "")))
        uniprot_residues = set(_split_csv(row.get("uniprot_residue_numbers", "")))
        uniprot_residues.update(_split_csv(row.get("representative_hotspot_uniprot_residues", "")))
        lookup[patch_id] = {"pdb_labels": pdb_labels, "uniprot_residues": uniprot_residues}
    return lookup


def _legacy_overlap(site_rows: list[Mapping[str, Any]], legacy: dict[str, dict[str, set[str]]]) -> str:
    if not legacy:
        return "not_evaluated"
    site_labels = {_residue_label(row) for row in site_rows}
    site_uniprots = {str(row["uniprot_residue_number"]) for row in site_rows}
    notes = []
    for patch_id, patch_data in sorted(legacy.items()):
        pdb_overlap = len(site_labels & patch_data["pdb_labels"])
        uni_overlap = len(site_uniprots & patch_data["uniprot_residues"])
        if pdb_overlap or uni_overlap:
            notes.append(f"{patch_id}:pdb={pdb_overlap},uniprot={uni_overlap}")
    return ";".join(notes) if notes else "none"


def _nearest_context_count(
    center: tuple[float, float, float],
    coords: list[tuple[float, float, float]],
    radius: float,
) -> int:
    return sum(1 for coord in coords if distance(center, coord) <= radius)


def _candidate_from_seed(
    seed: Mapping[str, Any],
    rows_by_chain: dict[str, list[dict[str, Any]]],
    all_context_coords: list[tuple[float, float, float]],
    all_context_residues: list[Mapping[str, Any]],
    site_radius: float,
    context_radius: float,
    hotspots_per_site: int,
    legacy: dict[str, dict[str, set[str]]],
    sasa_status: str,
    rsa_surface_threshold: float,
) -> dict[str, Any] | None:
    chain_id = str(seed["chain_id"])
    seed_coord = _row_coord(seed)
    site_rows = [row for row in rows_by_chain[chain_id] if distance(seed_coord, _row_coord(row)) <= site_radius]
    if not site_rows:
        return None

    coords = [_row_coord(row) for row in site_rows]
    center = centroid(coords)
    avg_radius = _safe_mean(distance(center, coord) for coord in coords)
    max_radius = max((distance(center, coord) for coord in coords), default=0.0)
    mean_exposure = _safe_mean(float(row["exposure_proxy"]) for row in site_rows)
    mean_neighbor_count = _safe_mean(float(row["neighbor_count"]) for row in site_rows)
    context_count = _nearest_context_count(center, all_context_coords, context_radius)

    site_size = len(site_rows)
    size_target = max(10, hotspots_per_site * 3)
    size_score = max(0.0, 1.0 - abs(site_size - size_target) / max(size_target, 1))
    compactness_score = max(0.0, 1.0 - (avg_radius / max(site_radius, 1.0)))
    exposure_score = mean_exposure / 100.0
    hotspot_capacity_score = min(1.0, site_size / max(hotspots_per_site, 1))
    score = (0.45 * exposure_score) + (0.25 * size_score) + (0.20 * compactness_score) + (0.10 * hotspot_capacity_score)

    hotspot_rows = sorted(
        site_rows,
        key=lambda row: (
            float(row["neighbor_count"]),
            distance(center, _row_coord(row)),
            _parse_int(row.get("uniprot_residue_number", "0")),
        ),
    )[:hotspots_per_site]
    hotspot_mean_distance, hotspot_max_distance = _pairwise_distance_stats(hotspot_rows)
    n_surface_residues = sum(1 for row in site_rows if _is_surface_for_quality(row, sasa_status, rsa_surface_threshold))
    hotspot_surface_count = sum(1 for row in hotspot_rows if _is_surface_for_quality(row, sasa_status, rsa_surface_threshold))
    surface_fraction = n_surface_residues / max(1, site_size)
    hotspot_surface_fraction = hotspot_surface_count / max(1, len(hotspot_rows))
    site_rsa_values = _numeric_values(site_rows, "rsa")
    hotspot_rsa_values = _numeric_values(hotspot_rows, "rsa")

    resnames = [_resname(row) for row in site_rows]
    hotspot_resnames = [_resname(row) for row in hotspot_rows]
    hydrophobic_count = sum(1 for resname in resnames if resname in HYDROPHOBIC_RESNAMES)
    aromatic_count = sum(1 for resname in resnames if resname in AROMATIC_RESNAMES)
    charged_count = sum(1 for resname in resnames if resname in CHARGED_RESNAMES)
    nearest_other_chain_distance, other_chain_within_6, other_chain_within_10 = _other_chain_metrics(site_rows, all_context_residues)
    span_status = _macrocycle_span_status(hotspot_mean_distance, hotspot_max_distance, len(hotspot_rows))
    anchor_status = _chemical_anchor_status(hydrophobic_count, aromatic_count, charged_count, hotspot_resnames)
    site_occlusion_status = _occlusion_status(nearest_other_chain_distance, other_chain_within_6, other_chain_within_10)
    quality_tier, quality_reason = _site_quality_tier_and_reason(
        surface_fraction,
        hotspot_surface_fraction,
        span_status,
        anchor_status,
        site_occlusion_status,
        nearest_other_chain_distance,
        other_chain_within_6,
    )

    site_rows_sorted = sorted(site_rows, key=lambda row: _parse_int(row.get("uniprot_residue_number", "0")))
    uniprots = [str(row["uniprot_residue_number"]) for row in site_rows_sorted]
    uniprot_nums = [_parse_int(value) for value in uniprots]
    pdb_residues = [str(row["pdb_residue_number"]) for row in site_rows_sorted]
    hotspot_labels = [_residue_label(row) for row in hotspot_rows]

    return {
        "seed_chain_id": chain_id,
        "seed_pdb_residue_number": seed["pdb_residue_number"],
        "seed_uniprot_residue_number": seed["uniprot_residue_number"],
        "site_rows": site_rows_sorted,
        "site_residue_keys": {f"{row['chain_id']}:{row['pdb_residue_number']}" for row in site_rows},
        "original_chain_ids": chain_id,
        "original_pdb_residue_numbers": ",".join(pdb_residues),
        "uniprot_residue_numbers": ",".join(uniprots),
        "uniprot_min": min(uniprot_nums) if uniprot_nums else "",
        "uniprot_max": max(uniprot_nums) if uniprot_nums else "",
        "uniprot_center": round(_safe_mean(float(value) for value in uniprot_nums), 3) if uniprot_nums else "",
        "center_x": round(center[0], 3),
        "center_y": round(center[1], 3),
        "center_z": round(center[2], 3),
        "n_site_residues": site_size,
        "n_surface_residues_in_site": n_surface_residues,
        "surface_residue_fraction": round(surface_fraction, 3),
        "mean_exposure_proxy": round(mean_exposure, 3),
        "mean_neighbor_count": round(mean_neighbor_count, 3),
        "avg_radius": round(avg_radius, 3),
        "max_radius": round(max_radius, 3),
        "native_context_residues_within_radius": context_count,
        "initial_hotspot_candidates": ",".join(hotspot_labels),
        "hotspot_pairwise_distance_mean": round(hotspot_mean_distance, 3),
        "hotspot_pairwise_distance_max": round(hotspot_max_distance, 3),
        "hotspot_surface_fraction": round(hotspot_surface_fraction, 3),
        "mean_site_rsa": round(_safe_mean(site_rsa_values), 4) if site_rsa_values else "",
        "mean_hotspot_rsa": round(_safe_mean(hotspot_rsa_values), 4) if hotspot_rsa_values else "",
        "hydrophobic_residue_count": hydrophobic_count,
        "aromatic_residue_count": aromatic_count,
        "charged_residue_count": charged_count,
        "nearest_other_chain_distance": round(nearest_other_chain_distance, 3),
        "other_chain_residues_within_6A": other_chain_within_6,
        "other_chain_residues_within_10A": other_chain_within_10,
        "macrocycle_span_status": span_status,
        "chemical_anchor_status": anchor_status,
        "occlusion_status": site_occlusion_status,
        "sasa_status": sasa_status,
        "rsa_status": _quality_status(row.get("rsa_status", "") for row in site_rows),
        "site_quality_tier": quality_tier,
        "site_quality_reason": quality_reason,
        "stage0_crop_allowed": "true" if quality_tier in STAGE0_ALLOWED_TIERS else "false",
        "legacy_patch_overlap": _legacy_overlap(site_rows, legacy),
        "candidate_score": round(score, 4),
        "exposure_summary": f"mean_exposure_proxy={mean_exposure:.2f}; mean_neighbor_count={mean_neighbor_count:.2f}",
        "structural_context_summary": (
            f"chain={chain_id}; residues={site_size}; avg_radius={avg_radius:.2f}; "
            f"context_residues_within_{context_radius:.1f}A={context_count}"
        ),
        "macrocycle_accessibility_rationale": (
            "Compact low-neighbor-count FGA surface cluster suitable for RFpeptides hotspot guidance; "
            "requires manual structure review before Stage 0."
        ),
    }


def _deduplicate_candidates(
    candidates: list[dict[str, Any]],
    min_center_distance: float,
    max_overlap_fraction: float,
) -> list[dict[str, Any]]:
    kept: list[dict[str, Any]] = []
    for candidate in sorted(candidates, key=lambda row: float(row["candidate_score"]), reverse=True):
        center = (float(candidate["center_x"]), float(candidate["center_y"]), float(candidate["center_z"]))
        keys = set(candidate["site_residue_keys"])
        duplicate = False
        for old in kept:
            old_center = (float(old["center_x"]), float(old["center_y"]), float(old["center_z"]))
            old_keys = set(old["site_residue_keys"])
            overlap = len(keys & old_keys) / max(1, min(len(keys), len(old_keys)))
            if str(candidate["original_chain_ids"]) == str(old["original_chain_ids"]) and (
                overlap >= max_overlap_fraction or distance(center, old_center) < min_center_distance
            ):
                duplicate = True
                break
        if not duplicate:
            kept.append(candidate)
    return kept


def _center_distance(a: Mapping[str, Any], b: Mapping[str, Any]) -> float:
    return distance(
        (float(a["center_x"]), float(a["center_y"]), float(a["center_z"])),
        (float(b["center_x"]), float(b["center_y"]), float(b["center_z"])),
    )


def _uniprot_overlap_fraction(a: Mapping[str, Any], b: Mapping[str, Any]) -> float:
    a_set = set(_split_csv(str(a.get("uniprot_residue_numbers", ""))))
    b_set = set(_split_csv(str(b.get("uniprot_residue_numbers", ""))))
    return len(a_set & b_set) / max(1, min(len(a_set), len(b_set)))


def _can_propose_candidate(
    candidate: Mapping[str, Any],
    proposed: list[Mapping[str, Any]],
    max_uniprot_overlap_fraction: float,
    min_uniprot_center_distance: float,
    min_center_distance: float,
) -> bool:
    if str(candidate.get("site_quality_tier", "")).lower() not in STAGE0_ALLOWED_TIERS:
        return False
    for old in proposed:
        if _uniprot_overlap_fraction(candidate, old) > max_uniprot_overlap_fraction:
            return False
        if abs(float(candidate["uniprot_center"]) - float(old["uniprot_center"])) < min_uniprot_center_distance:
            return False
        if str(candidate["original_chain_ids"]) == str(old["original_chain_ids"]) and _center_distance(candidate, old) < min_center_distance:
            return False
    return True


def _format_candidate_rows(
    candidates: list[dict[str, Any]],
    propose_sites: int,
    max_uniprot_overlap_fraction: float,
    min_uniprot_center_distance: float,
    min_center_distance: float,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], int]:
    proposed_indices: set[int] = set()
    proposed: list[dict[str, Any]] = []
    for idx, candidate in enumerate(candidates):
        if len(proposed) >= propose_sites:
            break
        if _can_propose_candidate(
            candidate,
            proposed,
            max_uniprot_overlap_fraction,
            min_uniprot_center_distance,
            min_center_distance,
        ):
            proposed_indices.add(idx)
            proposed.append(candidate)

    formatted: list[dict[str, Any]] = []
    deferred: list[dict[str, Any]] = []
    for idx, candidate in enumerate(candidates, start=1):
        row = {key: value for key, value in candidate.items() if key not in {"site_rows", "site_residue_keys"}}
        row["site_id"] = f"RFpep_Candidate_{idx:03d}"
        if (idx - 1) in proposed_indices:
            site_number = sum(1 for proposed_idx in sorted(proposed_indices) if proposed_idx <= (idx - 1))
            row["site_label"] = f"RFpep_Site_{site_number}"
            row["selection_status"] = "proposed_pending_manual_review"
            row["selection_notes"] = (
                "Diversity-filtered computational clean-restart candidate. Review in molecular viewer before Stage 0 crop/hotspot finalization."
            )
            formatted.append(row)
        else:
            row["site_label"] = ""
            row["selection_status"] = "deferred_candidate"
            row["selection_notes"] = "Not selected for first RFpeptides pilot proposal."
            formatted.append(row)
            deferred.append(row)
    return formatted, deferred, len(proposed_indices)


def _selection_markdown(rows: list[dict[str, Any]], args: argparse.Namespace, output_root: Path) -> str:
    proposed = [row for row in rows if row["selection_status"] == "proposed_pending_manual_review"]
    columns = [
        "site_label",
        "site_id",
        "original_chain_ids",
        "initial_hotspot_candidates",
        "n_site_residues",
        "mean_exposure_proxy",
        "avg_radius",
        "candidate_score",
        "legacy_patch_overlap",
    ]
    return f"""# FGA RFpeptides Target-Site Rediscovery

Status: computational first pass; manual review required before Stage 0.

This file was generated for the RFpeptides article-style clean-restart route.
It does not reuse previous Patch_A / Patch_B / Patch_C definitions as design
inputs. Any `legacy_patch_overlap` value is historical context only.

## Parameters

```text
output_root: {output_root}
min_uniprot: {args.min_uniprot}
max_uniprot: {args.max_uniprot}
surface_quantile: {args.surface_quantile}
exposure_radius: {args.exposure_radius}
site_radius: {args.site_radius}
context_radius: {args.context_radius}
min_site_residues: {args.min_site_residues}
max_candidates: {args.max_candidates}
propose_sites: {args.propose_sites}
hotspots_per_site: {args.hotspots_per_site}
rsa_surface_threshold: {args.rsa_surface_threshold}
enable_fpocket: {args.enable_fpocket}
fpocket_bin: {args.fpocket_bin}
fpocket_distance_strong: {args.fpocket_distance_strong}
fpocket_distance_moderate: {args.fpocket_distance_moderate}
fpocket_distance_weak: {args.fpocket_distance_weak}
proposal_max_uniprot_overlap_fraction: {args.proposal_max_uniprot_overlap_fraction}
proposal_min_uniprot_center_distance: {args.proposal_min_uniprot_center_distance}
proposal_min_center_distance: {args.proposal_min_center_distance}
```

## Proposed Target Sites

{rows_to_markdown(proposed, columns, "No target sites were proposed.")}

## Required Manual Review Before Stage 0

- Confirm each proposed target site is accessible in the native fibrinogen
  context.
- Confirm the hotspot residues are compact, exposed, and biologically/design
  defensible.
- Confirm that old Patch_A / Patch_B overlap, if present, is treated only as
  historical convergence and not as the selection reason.
- Confirm whether the proposed `RFpep_Site_1` / `RFpep_Site_2` labels should be
  accepted, renamed, or rejected.
- Do not generate RFpeptides backbones until this review is complete.
"""


def _quality_review_markdown(rows: list[dict[str, Any]], args: argparse.Namespace, output_root: Path) -> str:
    tier_counts = defaultdict(int)
    for row in rows:
        tier_counts[row.get("site_quality_tier", "missing")] += 1

    proposed = [row for row in rows if row["selection_status"] == "proposed_pending_manual_review"]
    allowed = [row for row in rows if row.get("stage0_crop_allowed") == "true"]
    quality_columns = [
        "site_label",
        "site_id",
        "site_quality_tier",
        "stage0_crop_allowed",
        "sasa_status",
        "rsa_status",
        "fpocket_status",
        "fpocket_support_status",
        "nearest_fpocket_distance_to_hotspot",
        "nearest_fpocket_distance_to_site",
        "nearest_fpocket_pocket_score",
        "nearest_fpocket_druggability_score",
        "surface_residue_fraction",
        "hotspot_surface_fraction",
        "mean_site_rsa",
        "mean_hotspot_rsa",
        "hotspot_pairwise_distance_max",
        "macrocycle_span_status",
        "chemical_anchor_status",
        "occlusion_status",
        "nearest_other_chain_distance",
        "other_chain_residues_within_6A",
        "candidate_score",
    ]
    reason_columns = ["site_label", "site_id", "site_quality_tier", "site_quality_reason", "legacy_patch_overlap"]

    tier_lines = "\n".join(f"- {tier}: {tier_counts[tier]}" for tier in ["high", "medium", "low", "reject", "missing"])
    return f"""# FGA RFpeptides Target-Site Quality Review

Status: semi-automatic Stage -1 quality review. This review converts the manual
pre-Stage-0 checklist into measurable site-quality fields.

Important rule: `legacy_patch_overlap` is historical context only and is not
used in site quality scoring.

Important rule: fpocket is supplementary pocket/groove evidence only. It is
not used as a final pass/fail gate and does not replace RSA, geometry,
occlusion, chemical-anchor, or manual structural review.

## Quality Tier Counts

{tier_lines}

## Stage 0 Gate

Only sites with `site_quality_tier` equal to `high` or `medium` are allowed to
enter Stage 0 target crop preparation. Sites marked `low` or `reject` stay in
the discovery table for review but should not be used for RFpeptides target
inputs.

## Quality Parameters

```text
output_root: {output_root}
rsa_surface_threshold: {args.rsa_surface_threshold}
enable_fpocket: {args.enable_fpocket}
fpocket_bin: {args.fpocket_bin}
fpocket_native_context_dir: {output_root / "00_site_discovery" / "fpocket_native_context"}
fpocket_distance_strong: {args.fpocket_distance_strong}
fpocket_distance_moderate: {args.fpocket_distance_moderate}
fpocket_distance_weak: {args.fpocket_distance_weak}
proposal_max_uniprot_overlap_fraction: {args.proposal_max_uniprot_overlap_fraction}
proposal_min_uniprot_center_distance: {args.proposal_min_uniprot_center_distance}
proposal_min_center_distance: {args.proposal_min_center_distance}
```

## Proposed Sites

{rows_to_markdown(proposed, quality_columns, "No high/medium proposed target sites were selected.")}

## All Stage-0-Allowed Sites

{rows_to_markdown(allowed, quality_columns, "No high/medium target sites were found.")}

## Quality Reasons For Proposed Sites

{rows_to_markdown(proposed, reason_columns, "No proposed target sites.")}
"""


def _pymol_residue_selection(chain_id: str, residue_csv: str) -> str:
    residues = _split_csv(residue_csv)
    if not residues:
        return "none"
    residue_terms = "+".join(residues)
    return f"(chain {chain_id} and resi {residue_terms})"


def _as_posix_path(path: str | Path) -> str:
    return str(path).replace("\\", "/")


def _write_pymol_review_scripts(rows: list[dict[str, Any]], cleaned_pdb: str | Path, discovery_dir: Path) -> None:
    pdb_path = _as_posix_path(resolve_path(cleaned_pdb))
    proposed = [row for row in rows if row["selection_status"] == "proposed_pending_manual_review"]
    for row in proposed:
        label = row.get("site_label") or row["site_id"]
        chain_id = str(row["original_chain_ids"]).split(",")[0]
        site_sel = _pymol_residue_selection(chain_id, str(row["original_pdb_residue_numbers"]))
        hotspot_residues = ",".join(label_value[1:] for label_value in _split_csv(str(row["initial_hotspot_candidates"])) if label_value.startswith(chain_id))
        hotspot_sel = _pymol_residue_selection(chain_id, hotspot_residues)
        pml = f"""reinitialize
load {pdb_path}, target
hide everything, all
show cartoon, target
color gray80, target

select {label}_site, target and {site_sel}
select {label}_hotspots, target and {hotspot_sel}
select {label}_other_chain_8A, byres ((target and not chain {chain_id}) within 8 of {label}_site)

show sticks, {label}_site
show spheres, {label}_hotspots
show sticks, {label}_other_chain_8A

color orange, {label}_site
color red, {label}_hotspots
color cyan, {label}_other_chain_8A

set sphere_scale, 0.45, {label}_hotspots
set stick_radius, 0.18, {label}_site
set stick_radius, 0.12, {label}_other_chain_8A

label {label}_hotspots and name CA, "%s%s" % (chain, resi)
zoom {label}_site, 14
orient {label}_site
"""
        write_markdown(discovery_dir / f"{label}_review.pml", pml)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Clean-restart RFpeptides target-site rediscovery for FGA. Does not run RFpeptides.",
    )
    parser.add_argument("--config", default="config/project.yaml")
    parser.add_argument("--output-root", default="results/rfpeptides_article_route_clean_20260609")
    parser.add_argument("--mapping-csv", default="data/annotations/FGA_structure_mapping.csv")
    parser.add_argument("--legacy-patches-csv", default="data/annotations/FGA_epitope_candidates.csv")
    parser.add_argument("--no-legacy-overlap", action="store_true")
    parser.add_argument("--target-chain-ids", default="", help="Optional comma-separated PDB chain IDs to include.")
    parser.add_argument("--min-uniprot", type=int, default=None)
    parser.add_argument("--max-uniprot", type=int, default=None)
    parser.add_argument("--surface-quantile", type=lambda value: _float_arg(value, "surface-quantile", 0.01, 0.95), default=0.35)
    parser.add_argument("--exposure-radius", type=lambda value: _float_arg(value, "exposure-radius", 1.0), default=10.0)
    parser.add_argument("--rsa-surface-threshold", type=lambda value: _fraction_arg(value, "rsa-surface-threshold"), default=0.20)
    parser.add_argument("--site-radius", type=lambda value: _float_arg(value, "site-radius", 1.0), default=12.0)
    parser.add_argument("--context-radius", type=lambda value: _float_arg(value, "context-radius", 1.0), default=14.0)
    parser.add_argument("--min-center-distance", type=lambda value: _float_arg(value, "min-center-distance", 0.0), default=8.0)
    parser.add_argument("--max-overlap-fraction", type=lambda value: _float_arg(value, "max-overlap-fraction", 0.0, 1.0), default=0.55)
    parser.add_argument(
        "--proposal-max-uniprot-overlap-fraction",
        type=lambda value: _float_arg(value, "proposal-max-uniprot-overlap-fraction", 0.0, 1.0),
        default=0.20,
    )
    parser.add_argument(
        "--proposal-min-uniprot-center-distance",
        type=lambda value: _float_arg(value, "proposal-min-uniprot-center-distance", 0.0),
        default=25.0,
    )
    parser.add_argument(
        "--proposal-min-center-distance",
        type=lambda value: _float_arg(value, "proposal-min-center-distance", 0.0),
        default=18.0,
    )
    parser.add_argument("--min-site-residues", type=int, default=6)
    parser.add_argument("--max-candidates", type=int, default=30)
    parser.add_argument("--propose-sites", type=int, default=2)
    parser.add_argument("--hotspots-per-site", type=int, default=4)
    parser.add_argument("--enable-fpocket", action="store_true")
    parser.add_argument("--fpocket-bin", default="fpocket")
    parser.add_argument("--fpocket-distance-strong", type=lambda value: _float_arg(value, "fpocket-distance-strong", 0.0), default=6.0)
    parser.add_argument("--fpocket-distance-moderate", type=lambda value: _float_arg(value, "fpocket-distance-moderate", 0.0), default=10.0)
    parser.add_argument("--fpocket-distance-weak", type=lambda value: _float_arg(value, "fpocket-distance-weak", 0.0), default=12.0)
    args = parser.parse_args()

    logger = setup_logger("18_discover_rfpeptides_target_sites")
    append_run_header(logger, "18_discover_rfpeptides_target_sites.py")
    config = load_config(args.config)

    if args.min_uniprot is None:
        args.min_uniprot = int(config.get("target_regions", {}).get("main_chain", {}).get("start", 36))
    if args.max_uniprot is None:
        args.max_uniprot = int(config.get("target_regions", {}).get("main_chain", {}).get("end", 866))
    if args.fpocket_distance_moderate < args.fpocket_distance_strong:
        raise RuntimeError("--fpocket-distance-moderate must be >= --fpocket-distance-strong")
    if args.fpocket_distance_weak < args.fpocket_distance_moderate:
        raise RuntimeError("--fpocket-distance-weak must be >= --fpocket-distance-moderate")

    output_root = resolve_path(args.output_root)
    discovery_dir = output_root / "00_site_discovery"
    mapping_rows = read_csv(args.mapping_csv)
    if not mapping_rows:
        raise RuntimeError(f"Missing or empty mapping CSV: {args.mapping_csv}")

    cleaned_pdb = config["structures"]["cleaned_pdb_file"]
    chains = parse_residues(cleaned_pdb)
    residue_lookup = {}
    all_context_coords: list[tuple[float, float, float]] = []
    all_context_residues: list[dict[str, Any]] = []
    for chain_id, residues in chains.items():
        for residue in residues:
            coord = ca_coord(residue)
            if coord is None:
                continue
            residue_lookup[(chain_id, residue["pdb_residue_number"])] = residue
            all_context_coords.append(coord)
            all_context_residues.append(
                {
                    "chain_id": chain_id,
                    "pdb_residue_number": residue["pdb_residue_number"],
                    "pdb_residue_name": residue["pdb_residue_name"],
                    "coord": coord,
                }
            )

    sasa_status, sasa_lookup = _calculate_sasa_by_residue(cleaned_pdb)

    target_chain_ids = set(_split_csv(args.target_chain_ids))
    mapped: list[dict[str, Any]] = []
    for row in mapping_rows:
        chain_id = row.get("chain_id", "")
        if target_chain_ids and chain_id not in target_chain_ids:
            continue
        uni = _parse_int(row.get("uniprot_residue_number", "0"))
        if uni < args.min_uniprot or uni > args.max_uniprot:
            continue
        residue = residue_lookup.get((chain_id, row.get("pdb_residue_number", "")))
        coord = ca_coord(residue) if residue else None
        if coord is None:
            continue
        item = dict(row)
        item["x"], item["y"], item["z"] = coord
        sasa = _lookup_sasa(item, sasa_lookup)
        rsa = _calculate_rsa(_resname(item), sasa)
        item["sasa_abs"] = round(sasa, 3) if sasa is not None else ""
        item["sasa_status"] = sasa_status
        item["rsa"] = round(rsa, 4) if rsa is not None else ""
        if rsa is not None:
            item["rsa_status"] = "freesasa_rsa"
        elif sasa_status == "freesasa":
            item["rsa_status"] = "rsa_unavailable"
        else:
            item["rsa_status"] = "proxy_only"
        mapped.append(item)

    if not mapped:
        raise RuntimeError("No mapped FGA residues with coordinates passed the target filters.")

    for idx, row in enumerate(mapped):
        coord = _row_coord(row)
        row["neighbor_count"] = sum(1 for other in all_context_coords if distance(coord, other) <= args.exposure_radius) - 1
    threshold = _percentile([float(row["neighbor_count"]) for row in mapped], args.surface_quantile)
    cap = max(threshold * 2.0, 1.0)
    for row in mapped:
        row["is_surface_candidate"] = float(row["neighbor_count"]) <= threshold
        row["exposure_proxy"] = round(max(0.0, 100.0 * (1.0 - min(float(row["neighbor_count"]), cap) / cap)), 3)
        row["is_surface_for_quality"] = _is_surface_for_quality(row, sasa_status, args.rsa_surface_threshold)

    surface_rows = [row for row in mapped if row["is_surface_candidate"]]
    if len(surface_rows) < args.min_site_residues:
        raise RuntimeError("Too few surface candidate residues; adjust --surface-quantile or input filters.")

    rows_by_chain: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in mapped:
        rows_by_chain[str(row["chain_id"])].append(row)

    legacy = {} if args.no_legacy_overlap else _legacy_patch_lookup(resolve_path(args.legacy_patches_csv))
    raw_candidates: list[dict[str, Any]] = []
    for seed in sorted(surface_rows, key=lambda row: (float(row["neighbor_count"]), _parse_int(row["uniprot_residue_number"]))):
        candidate = _candidate_from_seed(
            seed,
            rows_by_chain,
            all_context_coords,
            all_context_residues,
            args.site_radius,
            args.context_radius,
            args.hotspots_per_site,
            legacy,
            sasa_status,
            args.rsa_surface_threshold,
        )
        if candidate is None or int(candidate["n_site_residues"]) < args.min_site_residues:
            continue
        raw_candidates.append(candidate)

    candidates = _deduplicate_candidates(raw_candidates, args.min_center_distance, args.max_overlap_fraction)
    candidates = candidates[: args.max_candidates]
    fpocket_result = _run_fpocket_cross_check(
        enabled=args.enable_fpocket,
        cleaned_pdb=cleaned_pdb,
        fpocket_dir=discovery_dir / "fpocket_native_context",
        fpocket_bin=args.fpocket_bin,
    )
    _annotate_candidates_with_fpocket(
        candidates,
        fpocket_result,
        args.fpocket_distance_strong,
        args.fpocket_distance_moderate,
        args.fpocket_distance_weak,
    )
    candidate_rows, deferred_rows, proposed_count = _format_candidate_rows(
        candidates,
        args.propose_sites,
        args.proposal_max_uniprot_overlap_fraction,
        args.proposal_min_uniprot_center_distance,
        args.proposal_min_center_distance,
    )

    candidate_fields = [
        "site_id",
        "site_label",
        "selection_status",
        "selection_notes",
        "original_chain_ids",
        "original_pdb_residue_numbers",
        "uniprot_residue_numbers",
        "uniprot_min",
        "uniprot_max",
        "uniprot_center",
        "center_x",
        "center_y",
        "center_z",
        "n_site_residues",
        "n_surface_residues_in_site",
        "surface_residue_fraction",
        "mean_exposure_proxy",
        "mean_neighbor_count",
        "avg_radius",
        "max_radius",
        "native_context_residues_within_radius",
        "exposure_summary",
        "structural_context_summary",
        "macrocycle_accessibility_rationale",
        "initial_hotspot_candidates",
        "hotspot_pairwise_distance_mean",
        "hotspot_pairwise_distance_max",
        "hotspot_surface_fraction",
        "mean_site_rsa",
        "mean_hotspot_rsa",
        "hydrophobic_residue_count",
        "aromatic_residue_count",
        "charged_residue_count",
        "nearest_other_chain_distance",
        "other_chain_residues_within_6A",
        "other_chain_residues_within_10A",
        "macrocycle_span_status",
        "chemical_anchor_status",
        "occlusion_status",
        "sasa_status",
        "rsa_status",
        *FPOCKET_FIELDS,
        "site_quality_tier",
        "site_quality_reason",
        "stage0_crop_allowed",
        "candidate_score",
        "seed_chain_id",
        "seed_pdb_residue_number",
        "seed_uniprot_residue_number",
        "legacy_patch_overlap",
    ]
    write_csv(discovery_dir / "FGA_rfpeptides_target_site_candidates.csv", candidate_rows, candidate_fields)
    write_csv(discovery_dir / "FGA_rfpeptides_rejected_or_deferred_sites.csv", deferred_rows, candidate_fields)
    write_markdown(
        discovery_dir / "FGA_rfpeptides_target_site_selection.md",
        _selection_markdown(candidate_rows, args, output_root),
    )
    write_markdown(
        discovery_dir / "FGA_rfpeptides_target_site_quality_review.md",
        _quality_review_markdown(candidate_rows, args, output_root),
    )
    _write_pymol_review_scripts(candidate_rows, cleaned_pdb, discovery_dir)

    inventory_fields = [
        "pdb_id",
        "chain_id",
        "pdb_residue_number",
        "pdb_residue_name",
        "uniprot_id",
        "uniprot_residue_number",
        "uniprot_residue_name",
        "neighbor_count",
        "exposure_proxy",
        "sasa_abs",
        "sasa_status",
        "rsa",
        "rsa_status",
        "is_surface_candidate",
        "is_surface_for_quality",
        "x",
        "y",
        "z",
    ]
    write_csv(discovery_dir / "FGA_rfpeptides_visible_residue_inventory.csv", mapped, inventory_fields)

    logger.info("Wrote target-site candidates: %s", len(candidate_rows))
    logger.info("Proposed pending manual review: %s", proposed_count)
    logger.info("SASA status: %s", sasa_status)
    logger.info("RSA surface threshold: %.3f", args.rsa_surface_threshold)
    logger.info("fpocket status: %s", fpocket_result.get("status", "not_run"))
    if proposed_count < args.propose_sites:
        logger.warning("Only %s diverse proposed site(s) found for requested propose-sites=%s.", proposed_count, args.propose_sites)
    logger.info("Output root: %s", output_root)
    logger.info("No RFpeptides generation was run.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
