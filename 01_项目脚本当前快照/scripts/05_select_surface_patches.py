from __future__ import annotations

import argparse
from collections import defaultdict
from typing import Callable

from common import append_run_header, import_available, load_config, read_csv, setup_logger, write_csv
from pdb_utils import ca_coord, centroid, distance, parse_residues, write_patch_pdb


def _percentile(values: list[float], fraction: float) -> float:
    if not values:
        return 0.0
    values = sorted(values)
    idx = min(len(values) - 1, max(0, int(round((len(values) - 1) * fraction))))
    return values[idx]


def _choose_patch(rows: list[dict], predicate: Callable[[dict], bool], limit: int) -> tuple[list[dict], bool]:
    eligible = [row for row in rows if predicate(row)]
    used_fallback = False
    if not eligible:
        eligible = rows
        used_fallback = True
    eligible = sorted(eligible, key=lambda row: (float(row["neighbor_count"]), int(row["uniprot_residue_number"])))
    return eligible[:limit], used_fallback


def _unique_preserve_order(values: list[str]) -> list[str]:
    seen = set()
    out = []
    for value in values:
        if value not in seen:
            seen.add(value)
            out.append(value)
    return out


def _unique_ints(values: list[str]) -> list[int]:
    seen = set()
    out: list[int] = []
    for value in values:
        item = int(value)
        if item not in seen:
            seen.add(item)
            out.append(item)
    return out


def _projection_lookup(mapping_rows: list[dict]) -> dict[str, dict[int, str]]:
    lookup: dict[str, dict[int, str]] = defaultdict(dict)
    for row in mapping_rows:
        lookup[row["chain_id"]][int(row["uniprot_residue_number"])] = row["pdb_residue_number"]
    return lookup


def _project_hotspots(chain_id: str, uniprot_residues: list[int], lookup: dict[str, dict[int, str]]) -> str:
    labels = []
    for uni_residue in sorted(set(uniprot_residues)):
        pdb_residue = lookup.get(chain_id, {}).get(uni_residue)
        if pdb_residue:
            labels.append(f"{chain_id}{pdb_residue}")
    return ",".join(_unique_preserve_order(labels))


def _representative_chain(uniprot_residues: list[int], lookup: dict[str, dict[int, str]]) -> str:
    """Choose one PDB chain by projected UniProt coverage, then prefer canonical chain order."""

    preferred_order = {chain_id: idx for idx, chain_id in enumerate(["A", "D", "G", "J"])}
    best_chain = ""
    best_key = (-1, -999)
    for chain_id in sorted(lookup):
        coverage = sum(1 for residue in set(uniprot_residues) if residue in lookup[chain_id])
        key = (coverage, -preferred_order.get(chain_id, 999))
        if key > best_key:
            best_key = key
            best_chain = chain_id
    return best_chain


def _chain_qualified_hotspots(chosen: list[dict]) -> str:
    labels = [
        f"{row['chain_id']}{row['pdb_residue_number']}"
        for row in sorted(chosen, key=lambda row: (row["chain_id"], int(row["pdb_residue_number"])))
    ]
    return ",".join(_unique_preserve_order(labels))


def main() -> int:
    parser = argparse.ArgumentParser(description="Select exposed FGA surface patches from native fibrinogen.")
    parser.add_argument("--config", default="config/project.yaml")
    args = parser.parse_args()

    logger = setup_logger("05_select_surface_patches")
    append_run_header(logger, "05_select_surface_patches.py")
    config = load_config(args.config)

    mapping_rows = read_csv("data/annotations/FGA_structure_mapping.csv")
    if not mapping_rows:
        raise RuntimeError("Missing FGA_structure_mapping.csv; cannot select patches.")
    projection = _projection_lookup(mapping_rows)
    chains = parse_residues(config["structures"]["cleaned_pdb_file"])
    residue_lookup = {}
    for chain_id, residues in chains.items():
        for residue in residues:
            residue_lookup[(chain_id, residue["pdb_residue_number"])] = residue

    if import_available("freesasa"):
        logger.info("freesasa detected; still writing neighbor-count proxy for this pipeline.")
    else:
        logger.warning("freesasa is unavailable; using residue neighbor count as exposure proxy.")

    mapped = []
    for row in mapping_rows:
        residue = residue_lookup.get((row["chain_id"], row["pdb_residue_number"]))
        coord = ca_coord(residue) if residue else None
        if coord is None:
            continue
        item = dict(row)
        item["x"], item["y"], item["z"] = coord
        mapped.append(item)
    if not mapped:
        raise RuntimeError("No mapped residues have usable CA coordinates; cannot select patches.")

    coords = [(float(row["x"]), float(row["y"]), float(row["z"])) for row in mapped]
    for i, row in enumerate(mapped):
        neighbor_count = sum(1 for j, coord in enumerate(coords) if i != j and distance(coords[i], coord) <= 10.0)
        row["neighbor_count"] = neighbor_count
        row["sasa_proxy"] = round(max(0.0, 100.0 * (1.0 - min(neighbor_count, 40) / 40.0)), 3)
    threshold = _percentile([float(row["neighbor_count"]) for row in mapped], 0.35)
    for row in mapped:
        row["is_surface_proxy"] = float(row["neighbor_count"]) <= threshold

    surface_rows = [row for row in mapped if row["is_surface_proxy"]]
    if len(surface_rows) < 10:
        surface_rows = mapped

    patch_specs = [
        (
            "Patch_A",
            "native_visible_surface",
            lambda row: True,
            "high",
            "medium",
            "Visible exposed FGA surface in 3GHG; requires real design and prediction confirmation.",
        ),
        (
            "Patch_B",
            "N_terminal_main_chain_surface",
            lambda row: 36 <= int(row["uniprot_residue_number"]) <= 200,
            "medium",
            "medium",
            "Visible exposed FGA 36-200 region, avoiding signal peptide and fibrinopeptide A.",
        ),
        (
            "Patch_C",
            "C_terminal_alphaC_exploratory",
            lambda row: int(row["uniprot_residue_number"]) >= 600,
            "exploratory",
            "high",
            "FGA C-terminal / alphaC exploratory region; high risk when not visible in 3GHG.",
        ),
    ]

    epitope_rows = []
    patch_residues = defaultdict(set)
    for patch_id, patch_type, predicate, priority, risk, note in patch_specs:
        chosen, used_fallback = _choose_patch(surface_rows, predicate, 18)
        selected_uniprots = _unique_ints([row["uniprot_residue_number"] for row in chosen])
        rep_chain = _representative_chain(selected_uniprots, projection)
        rep_hotspots = _project_hotspots(rep_chain, selected_uniprots, projection) if rep_chain else ""
        chosen_coords = [(float(row["x"]), float(row["y"]), float(row["z"])) for row in chosen]
        center = centroid(chosen_coords)

        for row in chosen:
            patch_residues[row["chain_id"]].add(row["pdb_residue_number"])

        if used_fallback:
            logger.warning("%s did not have region-specific visible surface residues; reused global visible surface rows.", patch_id)

        epitope_rows.append(
            {
                "patch_id": patch_id,
                "patch_type": patch_type,
                "chain_id": ",".join(sorted({row["chain_id"] for row in chosen})),
                "hotspot_pdb_residues": _chain_qualified_hotspots(chosen),
                "representative_chain": rep_chain,
                "representative_hotspot_residues": rep_hotspots,
                "representative_hotspot_uniprot_residues": ",".join(str(residue) for residue in sorted(set(selected_uniprots))),
                "uniprot_residue_numbers": ",".join(str(row["uniprot_residue_number"]) for row in chosen),
                "pdb_residue_numbers": ",".join(str(row["pdb_residue_number"]) for row in chosen),
                "center_x": round(center[0], 3),
                "center_y": round(center[1], 3),
                "center_z": round(center[2], 3),
                "n_surface_residues": len(chosen),
                "mean_sasa": round(sum(float(row["sasa_proxy"]) for row in chosen) / len(chosen), 3) if chosen else 0.0,
                "priority": priority,
                "risk_level": risk,
                "selection_status": "fallback_no_region_specific_surface_residues" if used_fallback else "region_specific_surface_residues",
                "note": note,
            }
        )
        logger.info("%s selected residues: %s", patch_id, len(chosen))

    surface_fields = [
        "pdb_id",
        "chain_id",
        "pdb_residue_number",
        "pdb_residue_name",
        "uniprot_id",
        "uniprot_residue_number",
        "uniprot_residue_name",
        "neighbor_count",
        "sasa_proxy",
        "is_surface_proxy",
        "x",
        "y",
        "z",
    ]
    write_csv("data/annotations/FGA_surface_residues.csv", mapped, surface_fields)

    epitope_fields = [
        "patch_id",
        "patch_type",
        "chain_id",
        "hotspot_pdb_residues",
        "representative_chain",
        "representative_hotspot_residues",
        "representative_hotspot_uniprot_residues",
        "uniprot_residue_numbers",
        "pdb_residue_numbers",
        "center_x",
        "center_y",
        "center_z",
        "n_surface_residues",
        "mean_sasa",
        "priority",
        "risk_level",
        "selection_status",
        "note",
    ]
    write_csv("data/annotations/FGA_epitope_candidates.csv", epitope_rows, epitope_fields)
    kept_atoms = write_patch_pdb(config["structures"]["cleaned_pdb_file"], "data/structures/prepared/FGA_target_patches.pdb", patch_residues)
    logger.info("Wrote FGA_surface_residues.csv / FGA_epitope_candidates.csv / FGA_target_patches.pdb")
    logger.info("Patch PDB atom count: %s", kept_atoms)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
