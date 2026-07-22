from __future__ import annotations

import argparse
import hashlib
import os
import re
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable, Mapping

from common import append_run_header, read_csv, resolve_path, rows_to_markdown, setup_logger, write_csv, write_markdown


IDENTITY_FIELDS = [
    "global_backbone_id",
    "global_backbone_label",
    "source_batch_label",
    "source_stage2_root",
    "source_local_design_id",
    "source_route_run_id",
    "source_route_batch_id",
    "pdb_sha256",
    "atom_coordinate_sha256",
    "local_design_id_occurrences",
    "local_design_id_batches",
    "exact_pdb_duplicate_count",
    "exact_pdb_duplicate_global_ids",
    "coordinate_duplicate_count",
    "coordinate_duplicate_global_ids",
]

DUPLICATE_AUDIT_FIELDS = [
    "duplicate_type",
    "duplicate_key",
    "occurrence_count",
    "batch_count",
    "batches",
    "global_backbone_ids",
    "all_pass_stage2_qc",
]


def _resolve_mixed_path(value: str | Path) -> Path:
    text = str(value).strip().replace("\\", "/")
    if os.name == "nt" and text.startswith("/mnt/") and len(text) > 7 and text[6] == "/":
        return Path(f"{text[5].upper()}:/{text[7:]}")
    if os.name != "nt" and len(text) >= 3 and text[1] == ":" and text[2] == "/":
        return Path(f"/mnt/{text[0].lower()}{text[2:]}")
    path = Path(text)
    return path if path.is_absolute() else resolve_path(path)


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _atom_coordinate_sha256(path: Path) -> str:
    """Hash atom identity and coordinates while ignoring PDB headers and B factors."""

    digest = hashlib.sha256()
    atom_count = 0
    with path.open("r", encoding="ascii", errors="strict") as handle:
        for line in handle:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            if len(line) < 54:
                raise RuntimeError(f"Short atom record in {path}: {line.rstrip()}")
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError as exc:
                raise RuntimeError(f"Invalid atom coordinates in {path}: {line.rstrip()}") from exc
            canonical = "|".join(
                [
                    line[0:6].strip(),
                    line[12:16].strip(),
                    line[16:17].strip(),
                    line[17:20].strip(),
                    line[21:22].strip() or "_",
                    line[22:26].strip(),
                    line[26:27].strip(),
                    f"{x:.3f}",
                    f"{y:.3f}",
                    f"{z:.3f}",
                    line[76:78].strip() if len(line) >= 78 else "",
                ]
            )
            digest.update(canonical.encode("ascii"))
            digest.update(b"\n")
            atom_count += 1
    if atom_count == 0:
        raise RuntimeError(f"No ATOM/HETATM records found in backbone PDB: {path}")
    return digest.hexdigest()


def _single_value(rows: Iterable[Mapping[str, str]], field: str, label: str) -> str:
    values = {str(row.get(field, "")).strip() for row in rows if str(row.get(field, "")).strip()}
    if len(values) != 1:
        raise RuntimeError(f"Expected exactly one non-empty {field} in {label}; found {sorted(values)}")
    return next(iter(values))


def _batch_label(route_batch_id: str, index: int) -> str:
    match = re.search(r"batch[_-]?(\d+)", route_batch_id, flags=re.IGNORECASE)
    return f"batch{int(match.group(1)):02d}" if match else f"batch{index:02d}"


def _joined_ids(rows: Iterable[Mapping[str, Any]]) -> str:
    return ",".join(sorted(str(row["global_backbone_id"]) for row in rows))


def _duplicate_rows(
    groups: Mapping[str, list[dict[str, Any]]],
    duplicate_type: str,
) -> list[dict[str, Any]]:
    output: list[dict[str, Any]] = []
    for key, rows in sorted(groups.items()):
        if len(rows) < 2:
            continue
        batches = sorted({str(row["source_batch_label"]) for row in rows})
        output.append(
            {
                "duplicate_type": duplicate_type,
                "duplicate_key": key,
                "occurrence_count": len(rows),
                "batch_count": len(batches),
                "batches": ",".join(batches),
                "global_backbone_ids": _joined_ids(rows),
                "all_pass_stage2_qc": "true"
                if all(str(row.get("pass_backbone_qc", "")).lower() == "true" for row in rows)
                else "false",
            }
        )
    return output


def _summary_markdown(
    *,
    rows: list[dict[str, Any]],
    pass_rows: list[dict[str, Any]],
    source_counts: Mapping[str, int],
    local_groups: Mapping[str, list[dict[str, Any]]],
    pdb_groups: Mapping[str, list[dict[str, Any]]],
    coordinate_groups: Mapping[str, list[dict[str, Any]]],
    output_dir: Path,
) -> str:
    repeated_local = {key: value for key, value in local_groups.items() if len(value) > 1}
    repeated_pdb = {key: value for key, value in pdb_groups.items() if len(value) > 1}
    repeated_coordinates = {key: value for key, value in coordinate_groups.items() if len(value) > 1}
    source_rows = [
        {
            "batch": batch,
            "all_backbones": count,
            "stage2_pass": sum(
                1
                for row in rows
                if row["source_batch_label"] == batch and row.get("pass_backbone_qc") == "true"
            ),
        }
        for batch, count in sorted(source_counts.items())
    ]
    return f"""# Stage 2.5 Global Backbone Identity Audit

This manifest merges Stage 2 outputs without changing any Stage 2 QC result.

Identity rules:

```text
global_backbone_id = route_run_id + "__" + source_local_design_id
global_backbone_label = batch label + "__" + source_local_design_id
```

`global_backbone_id` is the machine key for downstream stages. The local
`design_id` is retained only as source provenance and must not be used as a
cross-batch primary key.

Two content hashes are recorded:

- `pdb_sha256`: exact source-file bytes.
- `atom_coordinate_sha256`: ordered atom identity and XYZ coordinates, ignoring
  headers, occupancy, and B factors.

## Counts

```text
merged_backbones: {len(rows)}
stage2_qc_pass: {len(pass_rows)}
unique_global_backbone_ids: {len({row['global_backbone_id'] for row in rows})}
unique_local_design_ids: {len(local_groups)}
local_design_ids_repeated_across_batches: {len(repeated_local)}
unique_exact_pdb_hashes: {len(pdb_groups)}
exact_pdb_duplicate_groups: {len(repeated_pdb)}
unique_atom_coordinate_hashes: {len(coordinate_groups)}
coordinate_duplicate_groups: {len(repeated_coordinates)}
```

## Batch Counts

{rows_to_markdown(source_rows, ['batch', 'all_backbones', 'stage2_pass'], 'No batches were read.')}

## Interpretation

Repeated local IDs are expected when each batch starts numbering at zero. They
are identity collisions, not evidence that the structures are identical.
Exact-file and coordinate duplicate groups are reported separately in the
duplicate audit CSV.

Output directory:

```text
{output_dir}
```
"""


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Merge Stage 2 backbone QC tables and assign globally unique backbone identities."
    )
    parser.add_argument("--stage2-roots", nargs="+", required=True)
    parser.add_argument("--output-root", required=True)
    args = parser.parse_args()

    logger = setup_logger("21b_build_stage2_global_backbone_manifest")
    append_run_header(logger, "21b_build_stage2_global_backbone_manifest.py")

    output_root = _resolve_mixed_path(args.output_root)
    output_dir = output_root / "03_backbone_qc_merged"
    merged_rows: list[dict[str, Any]] = []
    original_fields: list[str] = []
    seen_roots: set[Path] = set()
    seen_route_run_ids: set[str] = set()
    seen_route_batch_ids: set[str] = set()
    source_counts: dict[str, int] = {}

    for root_index, root_value in enumerate(args.stage2_roots, start=1):
        stage2_root = _resolve_mixed_path(root_value).resolve()
        if stage2_root in seen_roots:
            raise RuntimeError(f"Stage 2 root was supplied more than once: {stage2_root}")
        seen_roots.add(stage2_root)
        qc_csv = stage2_root / "03_backbone_qc" / "FGA_rfpeptides_backbones_qc.csv"
        rows = read_csv(qc_csv)
        if not rows:
            raise RuntimeError(f"Missing or empty Stage 2 QC table: {qc_csv}")
        if not original_fields:
            original_fields = list(rows[0].keys())
        elif list(rows[0].keys()) != original_fields:
            raise RuntimeError(f"Stage 2 QC columns differ from the first input: {qc_csv}")

        route_run_id = _single_value(rows, "route_run_id", str(qc_csv))
        route_batch_id = _single_value(rows, "route_batch_id", str(qc_csv))
        if route_run_id in seen_route_run_ids:
            raise RuntimeError(f"Duplicate route_run_id across Stage 2 roots: {route_run_id}")
        if route_batch_id in seen_route_batch_ids:
            raise RuntimeError(f"Duplicate route_batch_id across Stage 2 roots: {route_batch_id}")
        seen_route_run_ids.add(route_run_id)
        seen_route_batch_ids.add(route_batch_id)
        batch_label = _batch_label(route_batch_id, root_index)

        local_ids = [str(row.get("design_id", "")).strip() for row in rows]
        if not all(local_ids):
            raise RuntimeError(f"Blank design_id found in {qc_csv}")
        if len(set(local_ids)) != len(local_ids):
            raise RuntimeError(f"Duplicate design_id values within one Stage 2 root: {qc_csv}")

        source_counts[batch_label] = len(rows)
        for row in rows:
            local_id = str(row["design_id"]).strip()
            pdb_path = _resolve_mixed_path(str(row.get("rf_pdb", "")))
            if not pdb_path.is_file():
                raise RuntimeError(f"Missing backbone PDB for {batch_label}/{local_id}: {pdb_path}")
            merged_rows.append(
                {
                    "global_backbone_id": f"{route_run_id}__{local_id}",
                    "global_backbone_label": f"{batch_label}__{local_id}",
                    "source_batch_label": batch_label,
                    "source_stage2_root": stage2_root,
                    "source_local_design_id": local_id,
                    "source_route_run_id": route_run_id,
                    "source_route_batch_id": route_batch_id,
                    "pdb_sha256": _sha256_file(pdb_path),
                    "atom_coordinate_sha256": _atom_coordinate_sha256(pdb_path),
                    **dict(row),
                }
            )

    global_ids = [str(row["global_backbone_id"]) for row in merged_rows]
    if len(set(global_ids)) != len(global_ids):
        collisions = sorted(key for key in set(global_ids) if global_ids.count(key) > 1)
        raise RuntimeError(f"global_backbone_id collision(s): {collisions[:10]}")

    local_groups: dict[str, list[dict[str, Any]]] = defaultdict(list)
    pdb_groups: dict[str, list[dict[str, Any]]] = defaultdict(list)
    coordinate_groups: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in merged_rows:
        local_groups[str(row["source_local_design_id"])].append(row)
        pdb_groups[str(row["pdb_sha256"])].append(row)
        coordinate_groups[str(row["atom_coordinate_sha256"])].append(row)

    for row in merged_rows:
        local_peers = local_groups[str(row["source_local_design_id"])]
        pdb_peers = pdb_groups[str(row["pdb_sha256"])]
        coordinate_peers = coordinate_groups[str(row["atom_coordinate_sha256"])]
        row.update(
            {
                "local_design_id_occurrences": len(local_peers),
                "local_design_id_batches": ",".join(
                    sorted({str(peer["source_batch_label"]) for peer in local_peers})
                ),
                "exact_pdb_duplicate_count": len(pdb_peers),
                "exact_pdb_duplicate_global_ids": _joined_ids(pdb_peers) if len(pdb_peers) > 1 else "",
                "coordinate_duplicate_count": len(coordinate_peers),
                "coordinate_duplicate_global_ids": _joined_ids(coordinate_peers)
                if len(coordinate_peers) > 1
                else "",
            }
        )

    pass_rows = [row for row in merged_rows if str(row.get("pass_backbone_qc", "")).lower() == "true"]
    duplicate_audit = (
        _duplicate_rows(local_groups, "local_design_id_reused_across_batches")
        + _duplicate_rows(pdb_groups, "exact_pdb_sha256_duplicate")
        + _duplicate_rows(coordinate_groups, "atom_coordinate_sha256_duplicate")
    )
    output_fields = IDENTITY_FIELDS + [field for field in original_fields if field not in IDENTITY_FIELDS]
    write_csv(output_dir / "FGA_rfpeptides_stage2_global_backbone_manifest.csv", merged_rows, output_fields)
    write_csv(output_dir / "FGA_rfpeptides_stage2_global_backbone_manifest_pass.csv", pass_rows, output_fields)
    write_csv(output_dir / "FGA_rfpeptides_stage2_duplicate_audit.csv", duplicate_audit, DUPLICATE_AUDIT_FIELDS)
    write_markdown(
        output_dir / "FGA_rfpeptides_stage2_global_backbone_manifest.md",
        _summary_markdown(
            rows=merged_rows,
            pass_rows=pass_rows,
            source_counts=source_counts,
            local_groups=local_groups,
            pdb_groups=pdb_groups,
            coordinate_groups=coordinate_groups,
            output_dir=output_dir,
        ),
    )

    logger.info("Merged Stage 2 backbones: %s", len(merged_rows))
    logger.info("Stage 2 QC pass rows: %s", len(pass_rows))
    logger.info("Unique global backbone IDs: %s", len(set(global_ids)))
    logger.info("Unique exact PDB hashes: %s", len(pdb_groups))
    logger.info("Unique atom-coordinate hashes: %s", len(coordinate_groups))
    logger.info("Output directory: %s", output_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
