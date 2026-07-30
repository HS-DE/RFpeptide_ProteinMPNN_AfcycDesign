from __future__ import annotations

import argparse
import copy
import hashlib
import os
import re
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable, Mapping

from common import (
    ROUTE_PROVENANCE_FIELDS,
    SOURCE_ROUTE_PROVENANCE_FIELDS,
    assert_active_route_path,
    append_run_header,
    canonical_json_sha256,
    load_route_manifest,
    read_csv,
    resolve_path,
    route_provenance_fields,
    rows_to_markdown,
    setup_logger,
    validate_route_project_config,
    validate_row_route_provenance,
    write_csv,
    write_markdown,
    write_route_manifest,
)


IDENTITY_FIELDS = [
    "global_backbone_id",
    "global_backbone_label",
    "source_batch_label",
    "source_stage2_root",
    "source_local_design_id",
    "source_route_run_id",
    "source_route_batch_id",
    *SOURCE_ROUTE_PROVENANCE_FIELDS,
    "pdb_sha256",
    "atom_coordinate_sha256",
    "local_design_id_occurrences",
    "local_design_id_batches",
    "exact_pdb_duplicate_count",
    "exact_pdb_duplicate_global_ids",
    "coordinate_duplicate_count",
    "coordinate_duplicate_global_ids",
]

MANIFEST_COMPATIBILITY_FIELDS = [
    "route_name",
    "route_protocol_version",
    "hotspot_mapping_version",
    "production_route",
    "cyclization",
    "site_labels",
    "protocol_peptide_length_min",
    "protocol_peptide_length_max",
    "run_peptide_length_min",
    "run_peptide_length_max",
    "project_config_sha256",
    "effective_project_config_sha256",
    "stage0_sites",
    "stage1_protocol",
    "rfpeptides_runtime",
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


def _manifest_compatibility_identity(manifest: Mapping[str, Any]) -> dict[str, Any]:
    return {field: copy.deepcopy(manifest.get(field)) for field in MANIFEST_COMPATIBILITY_FIELDS}


def _validate_compatible_manifest(
    manifest: Mapping[str, Any],
    expected_identity: Mapping[str, Any],
    label: str,
) -> None:
    observed = _manifest_compatibility_identity(manifest)
    for field in MANIFEST_COMPATIBILITY_FIELDS:
        if observed[field] != expected_identity[field]:
            raise RuntimeError(
                f"{label} is incompatible with the first Stage 2 route manifest for {field}: "
                f"observed={observed[field]!r}, expected={expected_identity[field]!r}"
            )


def _validate_all_row_provenance(
    rows: list[dict[str, str]],
    expected: Mapping[str, str],
    label: str,
) -> None:
    # The common validator confirms the referenced manifest bytes; the loop then
    # checks every row without repeatedly hashing the same small JSON file.
    validate_row_route_provenance(rows[0], expected, f"{label} row 1")
    expected_manifest = assert_active_route_path(expected["route_manifest_path"], f"{label} route manifest")
    for index, row in enumerate(rows, start=1):
        row_label = f"{label} row {index}"
        for field in ROUTE_PROVENANCE_FIELDS:
            observed = str(row.get(field, "")).strip()
            if field == "route_manifest_path":
                observed_path = assert_active_route_path(observed, f"{row_label} route manifest")
                if observed_path.resolve() != expected_manifest.resolve():
                    raise RuntimeError(
                        f"{row_label} route manifest path mismatch: "
                        f"observed={observed_path}, expected={expected_manifest}"
                    )
            elif observed != str(expected[field]):
                raise RuntimeError(
                    f"{row_label} route provenance mismatch for {field}: "
                    f"observed={observed!r}, expected={expected[field]!r}"
                )


def _aggregate_manifest_payload(
    source_manifests: list[dict[str, Any]],
    merged_count: int,
) -> dict[str, Any]:
    first = source_manifests[0]["manifest"]
    payload = copy.deepcopy(
        {
            key: value
            for key, value in first.items()
            if key not in {"batch_id", "created_at", "num_designs_requested", "run_id"}
        }
    )
    source_records = sorted(
        [
            {
                "run_id": str(item["manifest"]["run_id"]),
                "batch_id": str(item["manifest"]["batch_id"]),
                "manifest_path": str(item["path"]),
                "manifest_sha256": str(item["sha256"]),
            }
            for item in source_manifests
        ],
        key=lambda record: (record["batch_id"], record["run_id"]),
    )
    source_identity_sha256 = canonical_json_sha256(
        [
            {
                "run_id": record["run_id"],
                "batch_id": record["batch_id"],
                "manifest_sha256": record["manifest_sha256"],
            }
            for record in source_records
        ]
    )
    payload.update(
        {
            "batch_id": f"stage2_5_merge_{len(source_records)}b_{source_identity_sha256[:12]}",
            "num_designs_requested": merged_count,
            "source_route_manifests": source_records,
            "stage2_5_operation": "global_identity_merge_and_backbone_family_selection",
            "stage2_5_source_count": len(source_records),
            "stage2_5_source_identity_sha256": source_identity_sha256,
        }
    )
    return payload


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
    parser.add_argument("--project-config", required=True)
    parser.add_argument("--output-root", required=True)
    args = parser.parse_args()

    logger = setup_logger("21b_build_stage2_global_backbone_manifest")
    append_run_header(logger, "21b_build_stage2_global_backbone_manifest.py")

    output_root = assert_active_route_path(
        _resolve_mixed_path(args.output_root),
        "Stage 21b aggregate output root",
        must_exist=False,
    )
    output_dir = output_root / "03_backbone_qc_merged"
    merged_rows: list[dict[str, Any]] = []
    original_fields: list[str] = []
    seen_roots: set[Path] = set()
    seen_route_run_ids: set[str] = set()
    seen_route_batch_ids: set[str] = set()
    source_counts: dict[str, int] = {}
    source_manifests: list[dict[str, Any]] = []
    expected_manifest_identity: dict[str, Any] | None = None

    for root_index, root_value in enumerate(args.stage2_roots, start=1):
        stage2_root = assert_active_route_path(
            _resolve_mixed_path(root_value),
            f"Stage 21b source Stage 2 root {root_index}",
        )
        if stage2_root in seen_roots:
            raise RuntimeError(f"Stage 2 root was supplied more than once: {stage2_root}")
        seen_roots.add(stage2_root)
        manifest_path, route_manifest, route_manifest_sha256 = load_route_manifest(stage2_root)
        validate_route_project_config(args.project_config, route_manifest)
        source_route_provenance = route_provenance_fields(
            manifest_path,
            route_manifest,
            route_manifest_sha256,
        )
        if expected_manifest_identity is None:
            expected_manifest_identity = _manifest_compatibility_identity(route_manifest)
        else:
            _validate_compatible_manifest(
                route_manifest,
                expected_manifest_identity,
                f"Stage 21b source manifest {manifest_path}",
            )

        qc_csv = assert_active_route_path(
            stage2_root / "03_backbone_qc" / "FGA_rfpeptides_backbones_qc.csv",
            f"Stage 21b source Stage 2 QC table {root_index}",
        )
        rows = read_csv(qc_csv)
        if not rows:
            raise RuntimeError(f"Missing or empty Stage 2 QC table: {qc_csv}")
        _validate_all_row_provenance(
            rows,
            source_route_provenance,
            f"Stage 21b {qc_csv}",
        )
        if len(rows) != int(route_manifest["num_designs_requested"]):
            raise RuntimeError(
                f"Stage 2 row count does not match route manifest for {qc_csv}: "
                f"rows={len(rows)}, requested={route_manifest['num_designs_requested']}"
            )
        if not original_fields:
            original_fields = list(rows[0].keys())
        elif list(rows[0].keys()) != original_fields:
            raise RuntimeError(f"Stage 2 QC columns differ from the first input: {qc_csv}")

        route_run_id = _single_value(rows, "route_run_id", str(qc_csv))
        route_batch_id = _single_value(rows, "route_batch_id", str(qc_csv))
        if route_run_id != str(route_manifest["run_id"]) or route_batch_id != str(route_manifest["batch_id"]):
            raise RuntimeError(f"Stage 2 row route identity does not match manifest: {qc_csv}")
        if route_run_id in seen_route_run_ids:
            raise RuntimeError(f"Duplicate route_run_id across Stage 2 roots: {route_run_id}")
        if route_batch_id in seen_route_batch_ids:
            raise RuntimeError(f"Duplicate route_batch_id across Stage 2 roots: {route_batch_id}")
        seen_route_run_ids.add(route_run_id)
        seen_route_batch_ids.add(route_batch_id)
        batch_label = _batch_label(route_batch_id, root_index)
        source_manifests.append(
            {
                "path": manifest_path,
                "manifest": route_manifest,
                "sha256": route_manifest_sha256,
            }
        )

        local_ids = [str(row.get("design_id", "")).strip() for row in rows]
        if not all(local_ids):
            raise RuntimeError(f"Blank design_id found in {qc_csv}")
        if len(set(local_ids)) != len(local_ids):
            raise RuntimeError(f"Duplicate design_id values within one Stage 2 root: {qc_csv}")

        source_counts[batch_label] = len(rows)
        for row in rows:
            local_id = str(row["design_id"]).strip()
            pdb_path = assert_active_route_path(
                _resolve_mixed_path(str(row.get("rf_pdb", ""))),
                f"Stage 21b source backbone PDB for {batch_label}/{local_id}",
            )
            if not pdb_path.is_file():
                raise RuntimeError(f"Missing backbone PDB for {batch_label}/{local_id}: {pdb_path}")
            merged_rows.append(
                {
                    **dict(row),
                    "global_backbone_id": f"{route_run_id}__{local_id}",
                    "global_backbone_label": f"{batch_label}__{local_id}",
                    "source_batch_label": batch_label,
                    "source_stage2_root": stage2_root,
                    "source_local_design_id": local_id,
                    "source_route_run_id": route_run_id,
                    "source_route_batch_id": route_batch_id,
                    "source_run_id": route_run_id,
                    "source_batch_id": route_batch_id,
                    "source_route_manifest": manifest_path,
                    "source_route_manifest_sha256": route_manifest_sha256,
                    "pdb_sha256": _sha256_file(pdb_path),
                    "atom_coordinate_sha256": _atom_coordinate_sha256(pdb_path),
                }
            )

    aggregate_payload = _aggregate_manifest_payload(source_manifests, len(merged_rows))
    existing_manifest_path = output_root / "route_manifest.json"
    if existing_manifest_path.exists():
        _, existing_manifest, _ = load_route_manifest(output_root)
        stable_identity_fields = [
            "batch_id",
            "num_designs_requested",
            "stage2_5_operation",
            "stage2_5_source_count",
            "stage2_5_source_identity_sha256",
        ]
        for field in stable_identity_fields:
            if existing_manifest.get(field) != aggregate_payload.get(field):
                raise RuntimeError(
                    f"Stage 21b output root already contains a different aggregate route manifest for {field}"
                )
        aggregate_payload["source_route_manifests"] = existing_manifest["source_route_manifests"]
        aggregate_payload["created_at"] = existing_manifest["created_at"]
    aggregate_manifest_path, aggregate_manifest, aggregate_manifest_sha256 = write_route_manifest(
        output_root,
        aggregate_payload,
    )
    aggregate_route_provenance = route_provenance_fields(
        aggregate_manifest_path,
        aggregate_manifest,
        aggregate_manifest_sha256,
    )
    for row in merged_rows:
        row.update(aggregate_route_provenance)

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
