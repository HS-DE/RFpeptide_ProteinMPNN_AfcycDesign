from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any, Iterable, Mapping

from common import append_run_header, load_config, read_csv, resolve_path, rows_to_markdown, setup_logger, write_csv, write_markdown
from pdb_utils import ca_coord, distance, parse_residues


STAGE0_ALLOWED_TIERS = {"high", "medium"}

MAPPING_FIELDS = [
    "site_label",
    "site_id",
    "rfpeptides_chain",
    "rfpeptides_residue_number",
    "rfpeptides_residue_label",
    "rfpeptides_residue_name",
    "original_pdb_id",
    "original_chain_id",
    "original_pdb_residue_number",
    "original_insertion_code",
    "original_pdb_residue_name",
    "uniprot_accession",
    "uniprot_residue_number",
    "uniprot_residue_name",
    "sasa_abs",
    "rsa",
    "rsa_status",
    "is_target_site_residue",
    "is_selected_hotspot",
    "selection_or_exclusion_note",
]

SUMMARY_FIELDS = [
    "site_label",
    "site_id",
    "site_quality_tier",
    "target_pdb",
    "crop_renumbering_mapping_csv",
    "site_residues_csv",
    "hotspots_txt",
    "hotspot_selection_rationale_md",
    "rfpeptides_chain",
    "rfpeptides_residue_range",
    "rfpeptides_hotspots",
    "original_chain_id",
    "original_crop_residue_range",
    "original_site_residue_count",
    "crop_residue_count",
    "crop_atom_count",
    "crop_context_radius",
    "candidate_score",
    "legacy_patch_overlap",
    "status",
    "notes",
]


def _split_csv(value: str) -> list[str]:
    return [item.strip() for item in str(value or "").split(",") if item.strip()]


def _truthy(value: Any) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes", "pass"}


def _safe_token(value: str) -> str:
    keep = []
    for ch in str(value):
        if ch.isalnum() or ch in {"_", "-", "."}:
            keep.append(ch)
        else:
            keep.append("_")
    return "".join(keep).strip("_") or "site"


def _residue_sort_key(value: str) -> tuple[int, str]:
    digits = "".join(ch for ch in value if ch.isdigit() or ch == "-")
    suffix = "".join(ch for ch in value if not (ch.isdigit() or ch == "-"))
    try:
        number = int(digits)
    except ValueError:
        number = 0
    return number, suffix


def _read_required_csv(path: Path) -> list[dict[str, str]]:
    rows = read_csv(path)
    if not rows:
        raise RuntimeError(f"Missing or empty CSV: {path}")
    return rows


def _candidate_lookup(rows: Iterable[Mapping[str, str]]) -> dict[str, dict[str, str]]:
    lookup: dict[str, dict[str, str]] = {}
    for row in rows:
        item = dict(row)
        for key in [item.get("site_label", ""), item.get("site_id", "")]:
            if key:
                lookup[key] = item
    return lookup


def _row_by_residue(rows: Iterable[Mapping[str, str]]) -> dict[tuple[str, str], dict[str, str]]:
    lookup: dict[tuple[str, str], dict[str, str]] = {}
    for row in rows:
        chain_id = str(row.get("chain_id", "")).strip() or "_"
        residue_number = str(row.get("pdb_residue_number", "")).strip()
        if chain_id and residue_number:
            lookup[(chain_id, residue_number)] = dict(row)
    return lookup


def _parse_hotspot_residues(labels_csv: str, chain_id: str, max_count: int) -> list[str]:
    residues: list[str] = []
    for label in _split_csv(labels_csv):
        if not label.startswith(chain_id):
            continue
        residue_number = label[len(chain_id) :].strip()
        if residue_number:
            residues.append(residue_number)
    return residues[:max_count]


def _validate_selected_candidate(row: Mapping[str, str], selected_id: str) -> None:
    tier = str(row.get("site_quality_tier", "")).strip().lower()
    if tier not in STAGE0_ALLOWED_TIERS:
        raise RuntimeError(f"{selected_id} has site_quality_tier={tier}; only high/medium sites may enter Stage 0.")
    if not _truthy(row.get("stage0_crop_allowed", "")):
        raise RuntimeError(f"{selected_id} has stage0_crop_allowed={row.get('stage0_crop_allowed', '')}; refusing Stage 0.")


def _site_chain_and_residues(row: Mapping[str, str]) -> tuple[str, list[str]]:
    chain_ids = _split_csv(str(row.get("original_chain_ids", "")))
    if len(chain_ids) != 1:
        raise RuntimeError(f"{row.get('site_label') or row.get('site_id')} has non-single-chain site: {chain_ids}")
    residues = _split_csv(str(row.get("original_pdb_residue_numbers", "")))
    if not residues:
        raise RuntimeError(f"{row.get('site_label') or row.get('site_id')} has no original_pdb_residue_numbers")
    return chain_ids[0], residues


def _crop_chain_residues(
    chain_residues: list[Mapping[str, Any]],
    site_residue_numbers: list[str],
    crop_context_radius: float,
    max_crop_residues: int,
) -> list[Mapping[str, Any]]:
    index_by_residue = {str(residue["pdb_residue_number"]): idx for idx, residue in enumerate(chain_residues)}
    missing = [residue for residue in site_residue_numbers if residue not in index_by_residue]
    if missing:
        raise RuntimeError(f"Site residue(s) missing from source PDB chain: {','.join(missing)}")

    site_indices = [index_by_residue[residue] for residue in site_residue_numbers]
    site_coords = [ca_coord(chain_residues[idx]) for idx in site_indices]
    site_coords = [coord for coord in site_coords if coord is not None]
    if not site_coords:
        raise RuntimeError("Selected site residues have no usable coordinates.")

    include_indices = set(site_indices)
    for idx, residue in enumerate(chain_residues):
        coord = ca_coord(residue)
        if coord is None:
            continue
        if any(distance(coord, site_coord) <= crop_context_radius for site_coord in site_coords):
            include_indices.add(idx)

    start = min(include_indices)
    end = max(include_indices)
    crop = chain_residues[start : end + 1]
    if len(crop) > max_crop_residues:
        raise RuntimeError(
            f"Crop would contain {len(crop)} residues, above --max-crop-residues={max_crop_residues}. "
            "Reduce --crop-context-radius or raise the limit after review."
        )
    return crop


def _rewrite_atom_line(line: str, atom_serial: int, target_chain_id: str, new_residue_number: int) -> str:
    raw = line.rstrip("\n")
    tail = raw[27:] if len(raw) > 27 else ""
    return f"{raw[:6]}{atom_serial:5d}{raw[11:21]}{target_chain_id}{new_residue_number:4d} {tail}".rstrip() + "\n"


def _write_cropped_pdb(
    source_pdb: Path,
    output_pdb: Path,
    source_chain_id: str,
    old_to_new_residue: Mapping[str, int],
    target_chain_id: str,
) -> int:
    output_pdb.parent.mkdir(parents=True, exist_ok=True)
    atom_serial = 1
    kept_atoms = 0
    with source_pdb.open("r", encoding="utf-8", errors="replace") as inp, output_pdb.open(
        "w", encoding="utf-8", newline="\n"
    ) as out:
        for line in inp:
            if not line.startswith("ATOM  "):
                continue
            chain_id = line[21].strip() or "_"
            if chain_id != source_chain_id:
                continue
            residue_number = line[22:26].strip() + line[26].strip()
            new_residue_number = old_to_new_residue.get(residue_number)
            if new_residue_number is None:
                continue
            out.write(_rewrite_atom_line(line, atom_serial, target_chain_id, new_residue_number))
            atom_serial += 1
            kept_atoms += 1
        out.write(f"TER   {atom_serial:5d}      {target_chain_id}{max(old_to_new_residue.values()):4d}\n")
        out.write("END\n")
    if kept_atoms == 0:
        raise RuntimeError(f"No atoms were written to {output_pdb}")
    return kept_atoms


def _mapping_row(
    *,
    site_label: str,
    site_id: str,
    rf_chain: str,
    rf_residue_number: int,
    residue: Mapping[str, Any],
    mapping_lookup: Mapping[tuple[str, str], Mapping[str, str]],
    inventory_lookup: Mapping[tuple[str, str], Mapping[str, str]],
    site_residues: set[str],
    hotspot_residues: set[str],
) -> dict[str, Any]:
    chain_id = str(residue["chain_id"])
    residue_number = str(residue["pdb_residue_number"])
    annotation = mapping_lookup.get((chain_id, residue_number), {})
    inventory = inventory_lookup.get((chain_id, residue_number), {})
    is_site = residue_number in site_residues
    is_hotspot = residue_number in hotspot_residues
    if is_hotspot:
        note = "selected_stage0_hotspot"
    elif is_site:
        note = "target_site_residue_not_hotspot_keep_for_contact_recovery"
    else:
        note = "crop_context_residue"
    return {
        "site_label": site_label,
        "site_id": site_id,
        "rfpeptides_chain": rf_chain,
        "rfpeptides_residue_number": rf_residue_number,
        "rfpeptides_residue_label": f"{rf_chain}{rf_residue_number}",
        "rfpeptides_residue_name": residue.get("pdb_residue_name", ""),
        "original_pdb_id": annotation.get("pdb_id", "3GHG"),
        "original_chain_id": chain_id,
        "original_pdb_residue_number": residue_number,
        "original_insertion_code": residue.get("pdb_icode", ""),
        "original_pdb_residue_name": residue.get("pdb_residue_name", ""),
        "uniprot_accession": annotation.get("uniprot_id", ""),
        "uniprot_residue_number": annotation.get("uniprot_residue_number", ""),
        "uniprot_residue_name": annotation.get("uniprot_residue_name", ""),
        "sasa_abs": inventory.get("sasa_abs", ""),
        "rsa": inventory.get("rsa", ""),
        "rsa_status": inventory.get("rsa_status", ""),
        "is_target_site_residue": "true" if is_site else "false",
        "is_selected_hotspot": "true" if is_hotspot else "false",
        "selection_or_exclusion_note": note,
    }


def _hotspot_text(site_label: str, hotspot_rows: list[Mapping[str, Any]], original_hotspots: list[str], source_chain_id: str) -> str:
    rf_hotspots = [str(row["rfpeptides_residue_label"]) for row in hotspot_rows]
    rf_python_list = "[" + ",".join(f"'{hotspot}'" for hotspot in rf_hotspots) + "]"
    original_labels = [f"{source_chain_id}{residue}" for residue in original_hotspots]
    return f"""# {site_label} RFpeptides hotspot residues

RFpeptides hotspot CSV:
{",".join(rf_hotspots)}

RFpeptides ppi.hotspot_res Python list:
{rf_python_list}

Original 3GHG hotspot residues:
{",".join(original_labels)}
"""


def _rationale_markdown(
    *,
    site_row: Mapping[str, str],
    target_pdb: Path,
    mapping_csv: Path,
    site_residue_csv: Path,
    hotspot_rows: list[Mapping[str, Any]],
    crop_rows: list[Mapping[str, Any]],
    source_chain_id: str,
    crop_context_radius: float,
    target_chain_id: str,
) -> str:
    site_label = str(site_row.get("site_label") or site_row.get("site_id"))
    rf_hotspots = [str(row["rfpeptides_residue_label"]) for row in hotspot_rows]
    crop_range = f"{target_chain_id}1-{target_chain_id}{len(crop_rows)}"
    original_range = (
        f"{source_chain_id}{crop_rows[0]['original_pdb_residue_number']}-"
        f"{source_chain_id}{crop_rows[-1]['original_pdb_residue_number']}"
    )
    hotspot_columns = [
        "rfpeptides_residue_label",
        "original_chain_id",
        "original_pdb_residue_number",
        "original_pdb_residue_name",
        "uniprot_residue_number",
        "sasa_abs",
        "rsa",
        "rsa_status",
    ]
    site_columns = [
        "site_label",
        "site_id",
        "site_quality_tier",
        "surface_residue_fraction",
        "hotspot_surface_fraction",
        "mean_site_rsa",
        "mean_hotspot_rsa",
        "macrocycle_span_status",
        "chemical_anchor_status",
        "occlusion_status",
        "candidate_score",
    ]
    return f"""# {site_label} Stage 0 Target Input

Status: target crop and hotspot mapping prepared for manual review. No
RFpeptides backbone generation was run.

## Source Site

{rows_to_markdown([site_row], site_columns, "No source site metadata.")}

Quality reason:

```text
{site_row.get("site_quality_reason", "")}
```

Legacy patch overlap is historical context only:

```text
{site_row.get("legacy_patch_overlap", "")}
```

## Crop

```text
target_pdb: {target_pdb}
mapping_csv: {mapping_csv}
site_residue_csv: {site_residue_csv}
original_chain: {source_chain_id}
original_crop_range: {original_range}
rfpeptides_chain: {target_chain_id}
rfpeptides_crop_range: {crop_range}
crop_context_radius_A: {crop_context_radius:g}
crop_residue_count: {len(crop_rows)}
```

The crop first includes same-chain residues within the requested CA distance of
the selected target site, then fills the interval between the earliest and
latest included residue to keep one continuous RFpeptides target chain.

## Hotspots

Use these RFpeptides-renumbered hotspot residues after manual review:

```text
{",".join(rf_hotspots)}
```

Command-list form:

```text
{ "[" + ",".join(f"'{hotspot}'" for hotspot in rf_hotspots) + "]" }
```

{rows_to_markdown(hotspot_rows, hotspot_columns, "No hotspots mapped.")}

Rationale:

- Hotspots are a deliberately small subset of the broader target site, not the
  whole site definition.
- Each hotspot is inside the selected Stage -1 target site and has RSA/SASA
  evidence exported above when FreeSASA was available.
- Non-hotspot site residues remain in `{site_residue_csv.name}` for downstream
  contact-recovery checks.
- The mapping table records how RFpeptides target-chain numbering maps back to
  original 3GHG chain/residue IDs and UniProt positions.
"""


def _summary_markdown(rows: list[Mapping[str, Any]], output_dir: Path) -> str:
    columns = [
        "site_label",
        "site_quality_tier",
        "target_pdb",
        "rfpeptides_residue_range",
        "rfpeptides_hotspots",
        "crop_residue_count",
        "crop_atom_count",
        "status",
    ]
    return f"""# FGA RFpeptides Stage 0 Target Inputs

Status: target PDBs, crop-renumbering maps, and hotspot rationales prepared.
No RFpeptides backbone generation was run.

Output directory:

```text
{output_dir}
```

{rows_to_markdown(rows, columns, "No Stage 0 target inputs were prepared.")}
"""


def main() -> int:
    parser = argparse.ArgumentParser(description="Prepare RFpeptides Stage 0 target PDBs and hotspot maps for selected FGA sites.")
    parser.add_argument("--config", default="config/project.yaml")
    parser.add_argument("--input-root", default="results/rfpeptides_article_route_clean_20260612")
    parser.add_argument("--output-root", default="results/rfpeptides_article_route_clean_20260612")
    parser.add_argument("--selected-sites", required=True, help="Comma-separated site labels or candidate IDs from Stage -1.")
    parser.add_argument("--crop-context-radius", type=float, default=10.0)
    parser.add_argument("--hotspots-per-site", type=int, default=4)
    parser.add_argument("--target-chain-id", default="A")
    parser.add_argument("--max-crop-residues", type=int, default=250)
    parser.add_argument("--candidate-csv", default="")
    parser.add_argument("--mapping-csv", default="data/annotations/FGA_structure_mapping.csv")
    parser.add_argument("--inventory-csv", default="")
    args = parser.parse_args()

    logger = setup_logger("19_prepare_rfpeptides_article_inputs")
    append_run_header(logger, "19_prepare_rfpeptides_article_inputs.py")
    config = load_config(args.config)

    selected_sites = _split_csv(args.selected_sites)
    if not selected_sites:
        raise RuntimeError("--selected-sites must contain at least one site label or candidate ID.")
    if args.crop_context_radius <= 0:
        raise RuntimeError("--crop-context-radius must be > 0")
    if args.hotspots_per_site <= 0:
        raise RuntimeError("--hotspots-per-site must be > 0")

    input_root = resolve_path(args.input_root)
    output_root = resolve_path(args.output_root)
    discovery_dir = input_root / "00_site_discovery"
    output_dir = output_root / "00_target_inputs"
    candidate_csv = resolve_path(args.candidate_csv) if args.candidate_csv else discovery_dir / "FGA_rfpeptides_target_site_candidates.csv"
    inventory_csv = resolve_path(args.inventory_csv) if args.inventory_csv else discovery_dir / "FGA_rfpeptides_visible_residue_inventory.csv"

    candidate_rows = _read_required_csv(candidate_csv)
    candidates = _candidate_lookup(candidate_rows)
    mapping_lookup = _row_by_residue(_read_required_csv(resolve_path(args.mapping_csv)))
    inventory_lookup = _row_by_residue(read_csv(inventory_csv))

    cleaned_pdb = resolve_path(config["structures"]["cleaned_pdb_file"])
    chains = parse_residues(cleaned_pdb)
    if not chains:
        raise RuntimeError(f"No residues parsed from cleaned PDB: {cleaned_pdb}")

    summary_rows: list[dict[str, Any]] = []
    for selected_id in selected_sites:
        site_row = candidates.get(selected_id)
        if site_row is None:
            raise RuntimeError(f"Selected site not found in {candidate_csv}: {selected_id}")
        _validate_selected_candidate(site_row, selected_id)

        site_label = site_row.get("site_label") or site_row.get("site_id") or selected_id
        site_id = site_row.get("site_id", "")
        safe_site = _safe_token(site_label)
        source_chain_id, site_residue_numbers = _site_chain_and_residues(site_row)
        source_chain_residues = chains.get(source_chain_id)
        if not source_chain_residues:
            raise RuntimeError(f"Chain {source_chain_id} not found in {cleaned_pdb}")

        hotspot_residue_numbers = _parse_hotspot_residues(
            str(site_row.get("initial_hotspot_candidates", "")),
            source_chain_id,
            args.hotspots_per_site,
        )
        if len(hotspot_residue_numbers) < args.hotspots_per_site:
            raise RuntimeError(
                f"{site_label} has only {len(hotspot_residue_numbers)} hotspot(s), "
                f"fewer than --hotspots-per-site={args.hotspots_per_site}."
            )

        crop_residues = _crop_chain_residues(
            source_chain_residues,
            site_residue_numbers,
            args.crop_context_radius,
            args.max_crop_residues,
        )
        old_to_new = {str(residue["pdb_residue_number"]): idx for idx, residue in enumerate(crop_residues, start=1)}
        missing_hotspots = [residue for residue in hotspot_residue_numbers if residue not in old_to_new]
        if missing_hotspots:
            raise RuntimeError(f"{site_label} hotspot(s) missing from crop: {','.join(missing_hotspots)}")

        target_pdb = output_dir / f"{safe_site}_target.pdb"
        mapping_csv = output_dir / f"{safe_site}_crop_renumbering_mapping.csv"
        site_residue_csv = output_dir / f"{safe_site}_site_residues.csv"
        hotspots_txt = output_dir / f"{safe_site}_hotspots.txt"
        rationale_md = output_dir / f"{safe_site}_hotspot_selection_rationale.md"

        crop_atom_count = _write_cropped_pdb(cleaned_pdb, target_pdb, source_chain_id, old_to_new, args.target_chain_id)

        site_residue_set = set(site_residue_numbers)
        hotspot_residue_set = set(hotspot_residue_numbers)
        mapping_rows = [
            _mapping_row(
                site_label=site_label,
                site_id=site_id,
                rf_chain=args.target_chain_id,
                rf_residue_number=old_to_new[str(residue["pdb_residue_number"])],
                residue=residue,
                mapping_lookup=mapping_lookup,
                inventory_lookup=inventory_lookup,
                site_residues=site_residue_set,
                hotspot_residues=hotspot_residue_set,
            )
            for residue in crop_residues
        ]
        mapping_by_original_residue = {str(row["original_pdb_residue_number"]): row for row in mapping_rows}
        hotspot_rows = [mapping_by_original_residue[residue] for residue in hotspot_residue_numbers]
        site_rows = [row for row in mapping_rows if row["is_target_site_residue"] == "true"]

        write_csv(mapping_csv, mapping_rows, MAPPING_FIELDS)
        write_csv(site_residue_csv, site_rows, MAPPING_FIELDS)
        write_markdown(hotspots_txt, _hotspot_text(site_label, hotspot_rows, hotspot_residue_numbers, source_chain_id))
        write_markdown(
            rationale_md,
            _rationale_markdown(
                site_row=site_row,
                target_pdb=target_pdb,
                mapping_csv=mapping_csv,
                site_residue_csv=site_residue_csv,
                hotspot_rows=hotspot_rows,
                crop_rows=mapping_rows,
                source_chain_id=source_chain_id,
                crop_context_radius=args.crop_context_radius,
                target_chain_id=args.target_chain_id,
            ),
        )

        rf_range = f"{args.target_chain_id}1-{args.target_chain_id}{len(crop_residues)}"
        original_range = (
            f"{source_chain_id}{crop_residues[0]['pdb_residue_number']}-"
            f"{source_chain_id}{crop_residues[-1]['pdb_residue_number']}"
        )
        summary_rows.append(
            {
                "site_label": site_label,
                "site_id": site_id,
                "site_quality_tier": site_row.get("site_quality_tier", ""),
                "target_pdb": target_pdb,
                "crop_renumbering_mapping_csv": mapping_csv,
                "site_residues_csv": site_residue_csv,
                "hotspots_txt": hotspots_txt,
                "hotspot_selection_rationale_md": rationale_md,
                "rfpeptides_chain": args.target_chain_id,
                "rfpeptides_residue_range": rf_range,
                "rfpeptides_hotspots": ",".join(str(row["rfpeptides_residue_label"]) for row in hotspot_rows),
                "original_chain_id": source_chain_id,
                "original_crop_residue_range": original_range,
                "original_site_residue_count": len(site_residue_numbers),
                "crop_residue_count": len(crop_residues),
                "crop_atom_count": crop_atom_count,
                "crop_context_radius": args.crop_context_radius,
                "candidate_score": site_row.get("candidate_score", ""),
                "legacy_patch_overlap": site_row.get("legacy_patch_overlap", ""),
                "status": "prepared_pending_manual_review",
                "notes": "Stage 0 target input only; RFpeptides generation not run.",
            }
        )
        logger.info(
            "Prepared %s: crop residues=%s, atoms=%s, hotspots=%s",
            site_label,
            len(crop_residues),
            crop_atom_count,
            ",".join(str(row["rfpeptides_residue_label"]) for row in hotspot_rows),
        )

    write_csv(output_dir / "FGA_rfpeptides_stage0_target_inputs_summary.csv", summary_rows, SUMMARY_FIELDS)
    write_markdown(output_dir / "FGA_rfpeptides_stage0_target_inputs_summary.md", _summary_markdown(summary_rows, output_dir))

    logger.info("Prepared Stage 0 target input site(s): %s", len(summary_rows))
    logger.info("Output directory: %s", output_dir)
    logger.info("No RFpeptides generation was run.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
