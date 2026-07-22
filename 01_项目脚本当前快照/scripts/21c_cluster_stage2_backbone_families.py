from __future__ import annotations

import argparse
import math
import os
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np

from common import append_run_header, read_csv, resolve_path, rows_to_markdown, setup_logger, write_csv, write_markdown
from pdb_utils import ca_coord, parse_residues


FAMILY_FIELDS = [
    "backbone_family_id",
    "backbone_family_method",
    "backbone_family_rmsd_threshold_A",
    "backbone_family_size",
    "backbone_family_batches",
    "family_representative_global_backbone_id",
    "family_member_quality_rank",
    "target_alignment_ca_rmsd_A",
    "peptide_ca_rmsd_to_family_representative_A",
    "cyclic_shift_to_family_representative",
    "stage2_5_selected",
    "stage2_5_selection_rank",
    "stage2_5_selection_reason",
]

FAMILY_SUMMARY_FIELDS = [
    "backbone_family_id",
    "peptide_length",
    "family_size",
    "family_batches",
    "representative_global_backbone_id",
    "representative_global_backbone_label",
    "representative_batch",
    "representative_pdb",
    "num_target_contacts",
    "num_target_site_contacts",
    "num_hotspot_contacts",
    "peptide_site_min_distance",
    "peptide_hotspot_min_distance",
    "macrocycle_terminal_cn_distance",
    "selected_for_stage3_screen",
    "stage2_5_selection_rank",
]


def _resolve_mixed_path(value: str | Path) -> Path:
    text = str(value).strip().replace("\\", "/")
    if os.name == "nt" and text.startswith("/mnt/") and len(text) > 7 and text[6] == "/":
        return Path(f"{text[5].upper()}:/{text[7:]}")
    if os.name != "nt" and len(text) >= 3 and text[1] == ":" and text[2] == "/":
        return Path(f"/mnt/{text[0].lower()}{text[2:]}")
    path = Path(text)
    return path if path.is_absolute() else resolve_path(path)


def _as_float(value: Any, default: float) -> float:
    try:
        parsed = float(str(value).strip())
    except ValueError:
        return default
    return parsed if math.isfinite(parsed) else default


def _as_int(value: Any, default: int = 0) -> int:
    try:
        return int(float(str(value).strip()))
    except ValueError:
        return default


def _quality_key(row: Mapping[str, Any]) -> tuple[Any, ...]:
    return (
        -_as_int(row.get("num_hotspot_contacts", 0)),
        -_as_int(row.get("num_target_site_contacts", 0)),
        _as_float(row.get("peptide_hotspot_min_distance", ""), 999.0),
        _as_float(row.get("peptide_site_min_distance", ""), 999.0),
        abs(_as_float(row.get("macrocycle_terminal_cn_distance", ""), 99.0) - 1.33),
        -_as_int(row.get("num_target_contacts", 0)),
        str(row.get("global_backbone_id", "")),
    )


def _ca_coords(residues: Sequence[Mapping[str, Any]], label: str) -> np.ndarray:
    points = []
    for residue in residues:
        coord = ca_coord(residue)
        if coord is None:
            raise RuntimeError(f"Missing CA coordinate in {label}")
        points.append(coord)
    if len(points) < 3:
        raise RuntimeError(f"Too few CA coordinates in {label}: {len(points)}")
    return np.asarray(points, dtype=float)


def _kabsch(mobile: np.ndarray, reference: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    if mobile.shape != reference.shape or mobile.shape[0] < 3:
        raise RuntimeError(f"Cannot align CA arrays with shapes {mobile.shape} and {reference.shape}")
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


def _cyclic_shift_min_rmsd(mobile: np.ndarray, reference: np.ndarray) -> tuple[float, int]:
    if mobile.shape != reference.shape:
        raise RuntimeError(f"Peptide CA shape mismatch: {mobile.shape} != {reference.shape}")
    shifted = np.stack([np.roll(reference, shift, axis=0) for shift in range(reference.shape[0])])
    squared = np.sum((shifted - mobile[None, :, :]) ** 2, axis=2)
    values = np.sqrt(np.mean(squared, axis=1))
    best_shift = int(np.argmin(values))
    return float(values[best_shift]), best_shift


def _load_aligned_peptide(
    row: Mapping[str, str],
    reference_target_ca: np.ndarray,
) -> tuple[np.ndarray, float]:
    pdb_path = _resolve_mixed_path(str(row.get("rf_pdb", "")))
    if not pdb_path.is_file():
        raise RuntimeError(f"Missing backbone PDB: {pdb_path}")
    chains = parse_residues(pdb_path)
    target_chain = str(row.get("target_chain", "")).strip()
    peptide_chain = str(row.get("peptide_chain", "")).strip()
    if target_chain not in chains or peptide_chain not in chains:
        raise RuntimeError(
            f"Expected target/peptide chains {target_chain}/{peptide_chain} were not found in {pdb_path}"
        )
    target_ca = _ca_coords(chains[target_chain], f"target chain {target_chain} in {pdb_path}")
    peptide_ca = _ca_coords(chains[peptide_chain], f"peptide chain {peptide_chain} in {pdb_path}")
    expected_length = _as_int(row.get("peptide_length", 0))
    if peptide_ca.shape[0] != expected_length:
        raise RuntimeError(
            f"Peptide length mismatch for {row.get('global_backbone_id')}: "
            f"CSV={expected_length}, PDB={peptide_ca.shape[0]}"
        )
    rotation, translation, target_rmsd = _kabsch(target_ca, reference_target_ca)
    return peptide_ca @ rotation.T + translation, target_rmsd


def _cluster_rows(
    rows: list[dict[str, Any]],
    threshold: float,
) -> list[dict[str, Any]]:
    families: list[dict[str, Any]] = []
    by_length: dict[int, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        by_length[_as_int(row.get("peptide_length", 0))].append(row)

    for length in sorted(by_length):
        length_rows = sorted(by_length[length], key=_quality_key)
        length_families: list[dict[str, Any]] = []
        for row in length_rows:
            coords = row.pop("_aligned_peptide_ca")
            centroid = coords.mean(axis=0)
            closest_family: dict[str, Any] | None = None
            closest_rmsd = float("inf")
            closest_shift = 0
            for family in length_families:
                if float(np.linalg.norm(centroid - family["centroid"])) > threshold:
                    continue
                rmsd, shift = _cyclic_shift_min_rmsd(coords, family["representative_coords"])
                if rmsd < closest_rmsd:
                    closest_family = family
                    closest_rmsd = rmsd
                    closest_shift = shift
            if closest_family is None or closest_rmsd > threshold:
                family = {
                    "length": length,
                    "representative": row,
                    "representative_coords": coords,
                    "centroid": centroid,
                    "members": [(row, 0.0, 0)],
                }
                length_families.append(family)
            else:
                closest_family["members"].append((row, closest_rmsd, closest_shift))

        for family_index, family in enumerate(length_families, start=1):
            family["family_id"] = f"S25_L{length:02d}_F{family_index:04d}"
            families.append(family)
    return families


def _select_family_representatives(
    families: Sequence[dict[str, Any]],
    max_selected: int,
    max_per_batch_length: int,
) -> list[dict[str, Any]]:
    queues: dict[tuple[int, str], list[dict[str, Any]]] = defaultdict(list)
    for family in families:
        representative = family["representative"]
        key = (int(family["length"]), str(representative["source_batch_label"]))
        queues[key].append(family)
    for queue in queues.values():
        queue.sort(key=lambda family: _quality_key(family["representative"]))

    lengths = sorted({key[0] for key in queues})
    batches = sorted({key[1] for key in queues})
    selected: list[dict[str, Any]] = []
    for slot in range(max_per_batch_length):
        for length in lengths:
            for batch in batches:
                queue = queues.get((length, batch), [])
                if len(queue) <= slot:
                    continue
                selected.append(queue[slot])
                if len(selected) >= max_selected:
                    return selected
    return selected


def _median(values: Iterable[float]) -> float | str:
    finite = [value for value in values if math.isfinite(value)]
    return round(float(np.median(finite)), 3) if finite else ""


def _summary_markdown(
    *,
    rows: Sequence[Mapping[str, Any]],
    families: Sequence[Mapping[str, Any]],
    selected: Sequence[Mapping[str, Any]],
    args: argparse.Namespace,
    output_dir: Path,
) -> str:
    length_rows = []
    for length in sorted({_as_int(row.get("peptide_length", 0)) for row in rows}):
        members = [row for row in rows if _as_int(row.get("peptide_length", 0)) == length]
        length_families = [family for family in families if int(family["length"]) == length]
        length_selected = [
            family for family in selected if int(family["length"]) == length
        ]
        length_rows.append(
            {
                "length": length,
                "stage2_pass": len(members),
                "families": len(length_families),
                "selected": len(length_selected),
            }
        )
    batch_rows = []
    for batch in sorted({str(row.get("source_batch_label", "")) for row in rows}):
        batch_rows.append(
            {
                "batch": batch,
                "stage2_pass": sum(1 for row in rows if row.get("source_batch_label") == batch),
                "selected": sum(
                    1 for family in selected if family["representative"].get("source_batch_label") == batch
                ),
            }
        )
    target_rmsds = [_as_float(row.get("target_alignment_ca_rmsd_A", ""), float("nan")) for row in rows]
    return f"""# Stage 2.5 Backbone Family Clustering And Diversity Selection

Input consists only of Stage 2 QC-pass backbones. No Stage 2 result is changed.

## Method

1. Align each output target chain to the Stage 0 RFpep_Site_2 target using all
   target CA atoms.
2. Compare peptide CA coordinates in that common target frame.
3. Minimize RMSD over cyclic residue-index shifts, without independently
   superposing the peptide.
4. Cluster only peptides of the same length with a deterministic greedy-leader
   algorithm.
5. Select at most one representative from each family, balanced by peptide
   length and source batch.

The {args.family_rmsd_threshold:g} A family threshold is an operational
diversity threshold, not a biological pass/fail standard.

```text
stage2_pass_backbones: {len(rows)}
backbone_families: {len(families)}
selected_for_stage3_screen: {len(selected)}
family_rmsd_threshold_A: {args.family_rmsd_threshold:g}
max_selected: {args.max_selected}
max_per_batch_length: {args.max_per_batch_length}
median_target_alignment_CA_RMSD_A: {_median(target_rmsds)}
```

## Length Coverage

{rows_to_markdown(length_rows, ['length', 'stage2_pass', 'families', 'selected'], 'No length groups.')}

## Batch Coverage

{rows_to_markdown(batch_rows, ['batch', 'stage2_pass', 'selected'], 'No batch groups.')}

## Downstream Rule

Downstream preparation must identify structures by `global_backbone_id`, not
the batch-local `design_id`. These are screening backbones, not final peptide
candidates.

Output directory:

```text
{output_dir}
```
"""


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Cluster Stage 2 QC-pass cyclic backbones and select a batch/length-balanced diversity panel."
    )
    parser.add_argument(
        "--manifest-pass-csv",
        required=True,
    )
    parser.add_argument(
        "--stage0-target-pdb",
        default=(
            "results/rfpeptides_article_route_clean_20260615_fpocket/"
            "00_target_inputs/RFpep_Site_2_target.pdb"
        ),
    )
    parser.add_argument(
        "--output-root",
        required=True,
    )
    parser.add_argument("--family-rmsd-threshold", type=float, default=4.0)
    parser.add_argument("--max-selected", type=int, default=260)
    parser.add_argument("--max-per-batch-length", type=int, default=4)
    args = parser.parse_args()

    logger = setup_logger("21c_cluster_stage2_backbone_families")
    append_run_header(logger, "21c_cluster_stage2_backbone_families.py")

    if args.family_rmsd_threshold <= 0:
        raise RuntimeError("--family-rmsd-threshold must be > 0")
    if args.max_selected <= 0 or args.max_per_batch_length <= 0:
        raise RuntimeError("--max-selected and --max-per-batch-length must be > 0")

    manifest_path = _resolve_mixed_path(args.manifest_pass_csv)
    rows = [dict(row) for row in read_csv(manifest_path)]
    if not rows:
        raise RuntimeError(f"Missing or empty Stage 2.5 pass manifest: {manifest_path}")
    if any(str(row.get("pass_backbone_qc", "")).lower() != "true" for row in rows):
        raise RuntimeError("The Stage 2.5 pass manifest contains a row that did not pass Stage 2 QC")
    global_ids = [str(row.get("global_backbone_id", "")).strip() for row in rows]
    if not all(global_ids) or len(set(global_ids)) != len(global_ids):
        raise RuntimeError("Stage 2.5 pass manifest has blank or duplicate global_backbone_id values")

    stage0_target_pdb = _resolve_mixed_path(args.stage0_target_pdb)
    stage0_chains = parse_residues(stage0_target_pdb)
    if "A" not in stage0_chains:
        raise RuntimeError(f"Stage 0 target PDB lacks chain A: {stage0_target_pdb}")
    reference_target_ca = _ca_coords(stage0_chains["A"], f"Stage 0 target {stage0_target_pdb}")

    for index, row in enumerate(rows, start=1):
        peptide_ca, target_rmsd = _load_aligned_peptide(row, reference_target_ca)
        row["_aligned_peptide_ca"] = peptide_ca
        row["target_alignment_ca_rmsd_A"] = round(target_rmsd, 3)
        if index % 250 == 0 or index == len(rows):
            logger.info("Loaded and target-aligned backbones: %s/%s", index, len(rows))

    families = _cluster_rows(rows, args.family_rmsd_threshold)
    selected_families = _select_family_representatives(
        families,
        args.max_selected,
        args.max_per_batch_length,
    )
    selected_family_ids = {str(family["family_id"]) for family in selected_families}
    selected_ranks = {
        str(family["family_id"]): rank for rank, family in enumerate(selected_families, start=1)
    }

    family_summary: list[dict[str, Any]] = []
    for family in families:
        family_id = str(family["family_id"])
        representative = family["representative"]
        members = list(family["members"])
        batches = sorted({str(member[0]["source_batch_label"]) for member in members})
        for member_rank, (row, rmsd, shift) in enumerate(members, start=1):
            is_representative = row["global_backbone_id"] == representative["global_backbone_id"]
            is_selected = is_representative and family_id in selected_family_ids
            if is_selected:
                selection_reason = (
                    "selected_family_representative;batch_length_balanced;"
                    "stage2_site_hotspot_macrocycle_clash_qc_pass"
                )
            elif is_representative:
                selection_reason = "not_selected_batch_length_quota"
            else:
                selection_reason = "not_selected_nonrepresentative_same_backbone_family"
            row.update(
                {
                    "backbone_family_id": family_id,
                    "backbone_family_method": "target_aligned_cyclic_shift_minimized_peptide_CA_greedy_leader",
                    "backbone_family_rmsd_threshold_A": args.family_rmsd_threshold,
                    "backbone_family_size": len(members),
                    "backbone_family_batches": ",".join(batches),
                    "family_representative_global_backbone_id": representative["global_backbone_id"],
                    "family_member_quality_rank": member_rank,
                    "peptide_ca_rmsd_to_family_representative_A": round(float(rmsd), 3),
                    "cyclic_shift_to_family_representative": shift,
                    "stage2_5_selected": "true" if is_selected else "false",
                    "stage2_5_selection_rank": selected_ranks.get(family_id, "") if is_selected else "",
                    "stage2_5_selection_reason": selection_reason,
                }
            )
        family_summary.append(
            {
                "backbone_family_id": family_id,
                "peptide_length": family["length"],
                "family_size": len(members),
                "family_batches": ",".join(batches),
                "representative_global_backbone_id": representative["global_backbone_id"],
                "representative_global_backbone_label": representative["global_backbone_label"],
                "representative_batch": representative["source_batch_label"],
                "representative_pdb": representative["rf_pdb"],
                "num_target_contacts": representative.get("num_target_contacts", ""),
                "num_target_site_contacts": representative.get("num_target_site_contacts", ""),
                "num_hotspot_contacts": representative.get("num_hotspot_contacts", ""),
                "peptide_site_min_distance": representative.get("peptide_site_min_distance", ""),
                "peptide_hotspot_min_distance": representative.get("peptide_hotspot_min_distance", ""),
                "macrocycle_terminal_cn_distance": representative.get("macrocycle_terminal_cn_distance", ""),
                "selected_for_stage3_screen": "true" if family_id in selected_family_ids else "false",
                "stage2_5_selection_rank": selected_ranks.get(family_id, ""),
            }
        )

    selected_rows = sorted(
        [row for row in rows if row.get("stage2_5_selected") == "true"],
        key=lambda row: _as_int(row.get("stage2_5_selection_rank", 0)),
    )
    output_dir = _resolve_mixed_path(args.output_root) / "04_backbone_diversity"
    original_fields = [field for field in rows[0].keys() if not field.startswith("_")]
    output_fields = FAMILY_FIELDS + [field for field in original_fields if field not in FAMILY_FIELDS]
    write_csv(output_dir / "FGA_rfpeptides_stage2_5_backbone_family_members.csv", rows, output_fields)
    write_csv(output_dir / "FGA_rfpeptides_stage2_5_backbone_family_summary.csv", family_summary, FAMILY_SUMMARY_FIELDS)
    write_csv(output_dir / "FGA_rfpeptides_stage2_5_selected_backbones.csv", selected_rows, output_fields)
    write_markdown(
        output_dir / "FGA_rfpeptides_stage2_5_backbone_diversity.md",
        _summary_markdown(
            rows=rows,
            families=families,
            selected=selected_families,
            args=args,
            output_dir=output_dir,
        ),
    )

    logger.info("Stage 2 pass backbones clustered: %s", len(rows))
    logger.info("Backbone families: %s", len(families))
    logger.info("Selected family representatives: %s", len(selected_rows))
    logger.info("Output directory: %s", output_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
