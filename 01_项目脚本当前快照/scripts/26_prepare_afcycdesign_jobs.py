from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import shutil
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable, Mapping

from common import (
    ROUTE_PROVENANCE_FIELDS,
    SOURCE_ROUTE_PROVENANCE_FIELDS,
    assert_active_route_path,
    canonical_json_sha256,
    append_run_header,
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


COLABDESIGN_GAMMA_COMMIT = "5ab4efaba2321a6c3c314b82d2fff8e0241f5c2d"
PROTOCOL_VERSION = "stage5A_v2_single_sequence_mlm015"

AA3_TO_AA1 = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LE": "L",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "MSE": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
}
LEGAL_AA = set("ACDEFGHIKLMNPQRSTVWY")

MANIFEST_FIELDS = [
    "stage5_candidate_id",
    "peptide_sequence_hash",
    "target_sequence_hash",
    "protocol_hash",
    "protocol_version",
    "selection_order",
    "batch",
    "backbone_id",
    "source_stage4_design_id",
    "peptide_sequence",
    "peptide_length",
    "site_label",
    "site_id",
    "design_pdb",
    "staged_design_pdb",
    "target_pdb",
    "target_chain",
    "peptide_chain",
    "complex_rosetta_total_score",
    "separated_rosetta_total_score",
    "ddg_proxy_no_repack",
    "ddg_proxy_per_peptide_residue",
    "ddg_proxy_per_target_contact",
    "ddg_proxy_per_site_contact",
    "num_target_contacts",
    "num_target_site_contacts",
    "num_hotspot_contacts",
    "peptide_hotspot_min_distance",
    "macrocycle_terminal_cn_distance",
    "clash_status",
    "detached_or_collapsed_flag",
    "sequence_liability_notes",
    "stage4_priority_rank",
    "stage4_priority_class",
    "selection_role",
    "selection_reason",
    "prediction_protocol",
    "validation_test_type",
    "uses_design_pose_as_model_input",
    "uses_template",
    "uses_initial_guess",
    "cyclic_topology_encoding",
    "seeds_per_candidate",
    "models_per_seed",
    "recycles",
    "requested_recycles",
    "forward_passes",
    "target_msa_mode",
    "peptide_msa_mode",
    "use_mlm",
    "mlm_replace_fraction",
    "prediction_job_count",
    "status",
] + ROUTE_PROVENANCE_FIELDS + SOURCE_ROUTE_PROVENANCE_FIELDS

JOB_FIELDS = [
    "stage5_job_id",
    "stage5_candidate_id",
    "peptide_sequence_hash",
    "protocol_hash",
    "batch",
    "backbone_id",
    "seed",
    "target_sequence_length",
    "target_sequence_hash",
    "peptide_sequence",
    "peptide_length",
    "cyclic_sequence_notation",
    "input_fasta",
    "job_spec_json",
    "design_pose_for_posthoc_comparison",
    "prediction_output_dir",
    "run_script",
    "prediction_protocol",
    "validation_test_type",
    "uses_template",
    "uses_initial_guess",
    "target_msa_mode",
    "peptide_msa_mode",
    "use_mlm",
    "mlm_replace_fraction",
    "requested_recycles",
    "forward_passes",
    "status",
] + ROUTE_PROVENANCE_FIELDS + SOURCE_ROUTE_PROVENANCE_FIELDS


def _safe_token(value: str) -> str:
    keep = []
    for ch in str(value):
        keep.append(ch if ch.isalnum() or ch in {"_", "-", "."} else "_")
    return "".join(keep).strip("_") or "item"


def _sha1_short(value: str, length: int = 8) -> str:
    return hashlib.sha1(value.encode("utf-8")).hexdigest()[:length]


def _sequence_hash(sequence: str) -> str:
    return _sha1_short(sequence.strip().upper())


def _protocol_hash(payload: Mapping[str, Any]) -> str:
    canonical = json.dumps(dict(payload), sort_keys=True, separators=(",", ":"))
    return _sha1_short(canonical)


def _validate_sequence(sequence: str, label: str, expected_length: int | None = None) -> str:
    normalized = str(sequence).strip().upper()
    if not normalized:
        raise RuntimeError(f"{label} is empty")
    invalid = sorted(set(normalized) - LEGAL_AA)
    if invalid:
        raise RuntimeError(f"{label} contains unsupported amino acids: {','.join(invalid)}")
    if expected_length is not None and len(normalized) != expected_length:
        raise RuntimeError(f"{label} length is {len(normalized)}, expected {expected_length}")
    return normalized


def _load_site_indices(mapping_csv: Path, expected_target_length: int) -> tuple[list[int], list[int], str]:
    rows = _read_required_csv(mapping_csv)
    ordered = sorted(rows, key=lambda row: _parse_int(row.get("rfpeptides_residue_number", "")))
    if len(ordered) != expected_target_length:
        raise RuntimeError(
            f"Stage 0 mapping has {len(ordered)} target rows, expected {expected_target_length}: {mapping_csv}"
        )
    observed_numbers = [_parse_int(row.get("rfpeptides_residue_number", "")) for row in ordered]
    if observed_numbers != list(range(1, expected_target_length + 1)):
        raise RuntimeError("Stage 0 mapping target numbering is not contiguous from 1 to the expected target length")
    mapping_sequence = "".join(
        AA3_TO_AA1.get(str(row.get("rfpeptides_residue_name", "")).strip().upper(), "X") for row in ordered
    )
    mapping_sequence = _validate_sequence(mapping_sequence, "Stage 0 mapping target sequence", expected_target_length)
    site_indices = [
        _parse_int(row["rfpeptides_residue_number"])
        for row in ordered
        if _truthy(row.get("is_target_site_residue", ""))
    ]
    hotspot_indices = [
        _parse_int(row["rfpeptides_residue_number"])
        for row in ordered
        if _truthy(row.get("is_selected_hotspot", ""))
    ]
    if not site_indices or not hotspot_indices:
        raise RuntimeError("Stage 0 mapping does not define both Site_2 and hotspot residues")
    return site_indices, hotspot_indices, mapping_sequence


def _resolve_mixed_path(value: str | Path) -> Path:
    text = str(value).strip().replace("\\", "/")
    if os.name == "nt" and text.startswith("/mnt/") and len(text) > 7 and text[6] == "/":
        return Path(f"{text[5].upper()}:/{text[7:]}")
    if os.name != "nt" and len(text) >= 3 and text[1] == ":" and text[2] == "/":
        return Path(f"/mnt/{text[0].lower()}{text[2:]}")
    path = Path(text)
    return path if path.is_absolute() else resolve_path(path)


def _to_wsl_path(path: str | Path) -> str:
    text = str(path).replace("\\", "/")
    if len(text) >= 3 and text[1] == ":" and text[2] == "/":
        return f"/mnt/{text[0].lower()}{text[2:]}"
    return text


def _quote_bash(value: str | Path) -> str:
    return "'" + str(value).replace("'", "'\"'\"'") + "'"


def _write_text_lf(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(text.rstrip() + "\n")


def _read_required_csv(path: Path) -> list[dict[str, str]]:
    rows = read_csv(path)
    if not rows:
        raise RuntimeError(f"Missing or empty CSV: {path}")
    return rows


def _parse_float(value: Any) -> float | None:
    try:
        out = float(str(value).strip())
    except ValueError:
        return None
    return None if math.isnan(out) or math.isinf(out) else out


def _parse_int(value: Any, default: int = 0) -> int:
    try:
        return int(float(str(value).strip()))
    except ValueError:
        return default


def _truthy(value: Any) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes", "pass"}


def _normalized(value: float | None, denominator: int) -> float | str:
    if value is None or denominator <= 0:
        return ""
    return round(value / denominator, 4)


def _warning_count(value: Any) -> int:
    text = str(value or "").strip()
    if not text or text.lower() == "none":
        return 0
    return len([item for item in text.split(";") if item.strip()])


def _candidate_sort_key(row: Mapping[str, Any]) -> tuple[Any, ...]:
    rank = _parse_int(row.get("stage4_priority_rank", ""), 999999)
    warning_count = _warning_count(row.get("sequence_liability_notes", ""))
    ddg = _parse_float(row.get("ddg_proxy_no_repack", ""))
    ddg_sort = ddg if ddg is not None else 999999.0
    return (
        warning_count,
        rank,
        ddg_sort,
        -_parse_int(row.get("num_hotspot_contacts", "")),
        -_parse_int(row.get("num_target_site_contacts", "")),
        str(row.get("stage4_design_id", "")),
    )


def _eligible_stage4_rows(rows: Iterable[Mapping[str, str]], batch: str) -> list[dict[str, Any]]:
    eligible: list[dict[str, Any]] = []
    for source in rows:
        row: dict[str, Any] = dict(source)
        row["batch"] = batch
        if not _truthy(row.get("pass_stage4_qc", "")):
            continue
        if row.get("macrocycle_geometry_status") != "pass_head_to_tail_macrocycle":
            continue
        if row.get("target_site_recovery_status") != "site_contact_pass":
            continue
        if row.get("hotspot_recovery_status") != "hotspot_contact_pass":
            continue
        if row.get("clash_status") != "pass_no_severe_clash":
            continue
        if row.get("detached_or_collapsed_flag") != "pass_basic_pose_geometry":
            continue
        sequence = str(row.get("peptide_sequence", "")).strip()
        row["peptide_length"] = _parse_int(row.get("peptide_length", ""), len(sequence)) or len(sequence)
        ddg = _parse_float(row.get("ddg_proxy_no_repack", ""))
        row["ddg_proxy_per_peptide_residue"] = _normalized(ddg, int(row["peptide_length"]))
        row["ddg_proxy_per_target_contact"] = _normalized(ddg, _parse_int(row.get("num_target_contacts", "")))
        row["ddg_proxy_per_site_contact"] = _normalized(ddg, _parse_int(row.get("num_target_site_contacts", "")))
        eligible.append(row)
    return eligible


def _source_route_compatibility_sha256(manifest: Mapping[str, Any]) -> str:
    stage0_sites = []
    for source_site in manifest["stage0_sites"]:
        site = dict(source_site)
        site.pop("target_pdb", None)
        site.pop("mapping_csv", None)
        stage0_sites.append(site)
    return canonical_json_sha256(
        {
            "route_name": manifest["route_name"],
            "route_protocol_version": manifest["route_protocol_version"],
            "hotspot_mapping_version": manifest["hotspot_mapping_version"],
            "cyclization": manifest["cyclization"],
            "project_config_sha256": manifest["project_config_sha256"],
            "effective_project_config_sha256": manifest["effective_project_config_sha256"],
            "site_labels": manifest["site_labels"],
            "stage0_sites": stage0_sites,
        }
    )


def _validate_source_route_set(records: Iterable[Mapping[str, Any]]) -> str:
    rows = list(records)
    run_ids = [str(row["manifest"]["run_id"]) for row in rows]
    manifest_hashes = [str(row["manifest_sha256"]) for row in rows]
    if len(run_ids) != len(set(run_ids)):
        raise RuntimeError("Stage 26 source routes contain a duplicate run_id")
    if len(manifest_hashes) != len(set(manifest_hashes)):
        raise RuntimeError("Stage 26 source routes contain a duplicate manifest SHA-256")
    compatibility = {_source_route_compatibility_sha256(row["manifest"]) for row in rows}
    if len(compatibility) != 1:
        raise RuntimeError("Stage 26 source routes are not provenance-compatible")
    return next(iter(compatibility))


def _select_candidates(rows: list[dict[str, Any]], candidate_count: int) -> list[dict[str, Any]]:
    by_batch_backbone: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        by_batch_backbone[(str(row["batch"]), str(row["backbone_id"]))].append(row)

    representatives = [sorted(group, key=_candidate_sort_key)[0] for group in by_batch_backbone.values()]
    reps_by_batch: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in representatives:
        reps_by_batch[str(row["batch"])].append(row)
    for batch in reps_by_batch:
        reps_by_batch[batch].sort(key=_candidate_sort_key)

    batches = sorted(reps_by_batch)
    selected: list[dict[str, Any]] = []
    minimum_per_batch = 2 if candidate_count >= 4 and len(batches) >= 2 else 1
    for batch in batches:
        selected.extend(reps_by_batch[batch][:minimum_per_batch])

    selected_ids = {str(row["stage4_design_id"]) for row in selected}
    remaining_reps = sorted(
        [row for row in representatives if str(row["stage4_design_id"]) not in selected_ids],
        key=_candidate_sort_key,
    )
    for row in remaining_reps:
        if len(selected) >= candidate_count:
            break
        selected.append(row)
        selected_ids.add(str(row["stage4_design_id"]))

    if len(selected) < candidate_count:
        selected_backbone_counts: dict[tuple[str, str], int] = defaultdict(int)
        for row in selected:
            selected_backbone_counts[(str(row["batch"]), str(row["backbone_id"]))] += 1
        remaining_rows = sorted(
            [row for row in rows if str(row["stage4_design_id"]) not in selected_ids],
            key=_candidate_sort_key,
        )
        for row in remaining_rows:
            key = (str(row["batch"]), str(row["backbone_id"]))
            if selected_backbone_counts[key] >= 2:
                continue
            selected.append(row)
            selected_ids.add(str(row["stage4_design_id"]))
            selected_backbone_counts[key] += 1
            if len(selected) >= candidate_count:
                break

    if len(selected) < candidate_count:
        raise RuntimeError(f"Only {len(selected)} eligible Stage 4 candidates were available for requested count={candidate_count}.")

    selected = sorted(selected, key=lambda row: (str(row["batch"]), _candidate_sort_key(row)))
    backbone_counts = defaultdict(int)
    for row in selected:
        backbone_counts[(str(row["batch"]), str(row["backbone_id"]))] += 1
    for row in selected:
        key = (str(row["batch"]), str(row["backbone_id"]))
        row["selection_role"] = "backbone_family_representative" if backbone_counts[key] == 1 else "sequence_sensitivity_test"
    return selected


def _pdb_chain_sequence(pdb_path: Path, chain_id: str) -> str:
    seen: set[tuple[str, str]] = set()
    sequence: list[str] = []
    with pdb_path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith("ATOM  "):
                continue
            chain = line[21].strip() or "_"
            if chain != chain_id:
                continue
            residue_key = (line[22:26].strip(), line[26].strip())
            if residue_key in seen:
                continue
            seen.add(residue_key)
            resname = line[17:20].strip().upper()
            aa = AA3_TO_AA1.get(resname)
            if aa is None:
                raise RuntimeError(f"Unsupported residue {resname} in {pdb_path}, chain {chain_id}, residue {residue_key}.")
            sequence.append(aa)
    if not sequence:
        raise RuntimeError(f"No residues found for chain {chain_id} in {pdb_path}")
    return "".join(sequence)


def _selection_reason(row: Mapping[str, Any]) -> str:
    batch = str(row["batch"])
    backbone = str(row["backbone_id"])
    role = str(row["selection_role"])
    ddg = row.get("ddg_proxy_no_repack", "")
    site_contacts = row.get("num_target_site_contacts", "")
    hotspot_contacts = row.get("num_hotspot_contacts", "")
    if role == "sequence_sensitivity_test":
        return f"Second sequence for {backbone}; included only to test sequence sensitivity on the same backbone."
    return (
        f"{batch} representative for distinct backbone {backbone}; Stage 4 hard QC passed; "
        f"ddg_proxy_no_repack={ddg}, site_contacts={site_contacts}, hotspot_contacts={hotspot_contacts}."
    )


def _protocol_audit_markdown(
    *,
    output_dir: Path,
    target_pdb: Path,
    target_sequence: str,
    candidate_count: int,
    seeds_per_candidate: int,
    models_per_seed: int,
    recycles: int,
    protocol_hash: str,
    use_mlm: bool,
    mlm_replace_fraction: float,
    afcycdesign_python: str,
    colabdesign_source: str,
    python_overlay: str,
    af_params: str,
) -> str:
    return f"""# Stage 5 AfCycDesign Protocol Audit

Status: method-level support confirmed; local prediction environment requires
the generated preflight check before any prediction is allowed.

## Topology Support

AfCycDesign supports head-to-tail cyclic peptides by replacing the linear
relative-position offsets within the peptide chain with cyclic offsets. For a
target-peptide complex, the cyclic offset is applied only to the peptide chain;
target and inter-chain offsets remain standard.

The planned sequence notation is:

```text
TARGET_SEQUENCE:(CYCLIC_PEPTIDE_SEQUENCE)
```

The runner additionally calls `add_cyclic_offset(model, [1])`, where chain
index 1 is the peptide.

Important limitation: this is a model-input positional encoding, not an
explicit chemical bond record or a Rosetta-style covalent topology declaration.
The predicted PDB may not contain a `CONECT` record for the terminal peptide
bond. Therefore every output must still pass terminal C-N distance and cyclic
topology QC.

## Validation Mode

Primary prepared mode:

```text
validation_test_type: independent_recovery
template_mode: none
use_initial_guess: false
target_msa_mode: single_sequence
peptide_msa_mode: single_sequence
msa_rows_input: 1
use_mlm: {str(use_mlm).lower()}
mlm_replace_fraction: {mlm_replace_fraction:g}
target_structure_used_as_model_input: false
design_peptide_pose_used_as_model_input: false
```

The Stage 4 design PDB is staged only for post-prediction target alignment,
peptide RMSD, and contact-recovery analysis. The prediction runner does not
load it.

If a future protocol uses the design PDB as a template or initial guess, that
run must be labeled `initial_guess_pose_stability_test`, not independent
recovery.

Ordinary linear AF2/RF2 without a verified cyclic-offset implementation is not
an accepted substitute for this Stage 5 protocol.

## Planned Compute

```text
candidate_count: {candidate_count}
seeds_per_candidate: {seeds_per_candidate}
independent_seed_jobs: {candidate_count * seeds_per_candidate}
models_per_seed: {models_per_seed}
requested_recycles: {recycles}
forward_passes: {recycles + 1}
protocol_version: {PROTOCOL_VERSION}
protocol_hash: {protocol_hash}
template_mode: none
use_initial_guess: false
target_pdb_for_sequence_and_posthoc_reference: {target_pdb}
target_sequence_length: {len(target_sequence)}
```

Each seed job evaluates all {models_per_seed} AlphaFold model parameter sets.
Thus the first round contains {candidate_count * seeds_per_candidate} seed-level
predictions and {candidate_count * seeds_per_candidate * models_per_seed}
model evaluations.

`num_msa=512` is model capacity only. This protocol supplies one A3M row and
must not be described as using 512 homolog sequences. The fixed gamma source
stores `model.aux["plddt"]` as a 0-1 fraction and multiplies it by 100 only
when writing PDB B-factors; Stage 5 metrics report both scales explicitly.

## Local Preflight

Expected Python:

```text
{afcycdesign_python}
```

Pinned gamma source checkout:

```text
{colabdesign_source}
```

Stage 5-only Python dependency overlay:

```text
{python_overlay}
```

Expected AlphaFold parameters:

```text
{af_params}
```

The currently deployed `ColabDesign-cyclic-binder` environment contains a
cyclic binder design implementation but was found to lack the gamma prediction
modules `colabdesign.af.contrib.predict` and
`colabdesign.af.contrib.cyclic`. Run `jobs/check_stage5_afcycdesign_protocol.sh`
after deploying the pinned gamma protocol. Predictions remain blocked until
that preflight passes.

Pinned ColabDesign gamma commit:

```text
{COLABDESIGN_GAMMA_COMMIT}
```

## Sources

- [AfCycDesign paper](https://www.nature.com/articles/s41467-025-59940-7)
- [Official ColabDesign gamma prediction notebook](https://github.com/sokrypton/ColabDesign/blob/gamma/af/examples/predict.ipynb)
- [RFpeptides paper](https://www.nature.com/articles/s41589-025-01929-w)
- [Pinned ColabDesign commit](https://github.com/sokrypton/ColabDesign/commit/{COLABDESIGN_GAMMA_COMMIT})
"""


def _manifest_markdown(
    *,
    rows: list[Mapping[str, Any]],
    eligible_counts: Mapping[tuple[str, str], int],
    output_dir: Path,
    protocol_audit: Path,
    jobs_csv: Path,
) -> str:
    columns = [
        "selection_order",
        "batch",
        "backbone_id",
        "peptide_sequence",
        "peptide_length",
        "ddg_proxy_no_repack",
        "ddg_proxy_per_peptide_residue",
        "ddg_proxy_per_target_contact",
        "ddg_proxy_per_site_contact",
        "num_target_site_contacts",
        "num_hotspot_contacts",
        "selection_role",
    ]
    count_lines = "\n".join(
        f"- {batch} / {backbone}: {eligible_counts[(batch, backbone)]} Stage 4 hard-QC candidate(s)"
        for batch, backbone in sorted(eligible_counts)
    )
    return f"""# FGA RFpeptides Stage 5 Candidate Manifest

Status: candidate merge, protocol audit, inputs, and gated jobs prepared. No
structure prediction was run.

Output directory:

```text
{output_dir}
```

## Selection Rules

- Every supplied source batch is represented when enough eligible candidates
  are available.
- Distinct backbone families are preferred over near-duplicate sequences.
- A selected backbone has one representative by default and at most two
  sequences. A second sequence is allowed only as a sequence-sensitivity test.
- Raw `ddg_proxy_no_repack` is not treated as a length-independent binding
  energy. Per-residue and per-contact normalizations are reported alongside it.
- All selected rows passed Stage 4 contact, macrocycle, severe-clash, and
  detached/collapsed gates.

## Eligible Backbone Families

{count_lines}

Backbone coverage and candidate counts are derived from the supplied source
route manifests and the current Stage 4 hard-QC tables; no batch or backbone is
hard-coded into this preparation step.

## Selected Candidates

{rows_to_markdown(rows, columns, "No candidates selected.")}

## Prediction Classification

The prepared primary route is sequence/target-based AfCycDesign independent
recovery. It uses no template and no initial guess. Design PDBs are retained
only for post-hoc target alignment and recovery metrics.

See the full protocol audit:

```text
{protocol_audit}
```

Seed-level job table:

```text
{jobs_csv}
```

## Required Stage 5 Result Metrics

- Target-aligned peptide backbone RMSD.
- RFpep_Site_2 contact recovery.
- Hotspot contact count and hotspot minimum distance.
- Off-site contacts and `same_target_site_flag`.
- Macrocycle terminal C-N distance and cyclic topology status.
- Severe clash status.
- Peptide pLDDT.
- Interface PAE or model-provided interface confidence.
- Recovery success count across five seeds per candidate.

These candidates are validation inputs, not final peptide candidates.
"""


def _fetch_protocol_script(colabdesign_source: str) -> str:
    return f"""#!/bin/bash
set -euo pipefail

SOURCE_DIR="${{COLABDESIGN_GAMMA_SOURCE:-{colabdesign_source}}}"
PINNED_COMMIT="{COLABDESIGN_GAMMA_COMMIT}"
COMMIT_MARKER="$SOURCE_DIR/.stage5_colabdesign_commit"

if [[ -f "$COMMIT_MARKER" ]]; then
  current_commit="$(tr -d '[:space:]' < "$COMMIT_MARKER")"
  if [[ "$current_commit" != "$PINNED_COMMIT" ]]; then
    echo "FAIL: existing gamma source is at $current_commit, expected $PINNED_COMMIT" >&2
    exit 2
  fi
  if [[ -f "$SOURCE_DIR/colabdesign/af/contrib/predict.py" ]] && [[ -f "$SOURCE_DIR/colabdesign/af/contrib/cyclic.py" ]]; then
    echo "PASS: pinned ColabDesign gamma source already exists."
    exit 0
  fi
  echo "FAIL: commit marker exists but required gamma prediction files are missing." >&2
  exit 2
fi

if [[ -e "$SOURCE_DIR" ]]; then
  echo "FAIL: source path exists without a valid Stage 5 commit marker: $SOURCE_DIR" >&2
  exit 2
fi

mkdir -p "$(dirname "$SOURCE_DIR")"
TMP_DIR="$(mktemp -d)"
trap 'rm -rf "$TMP_DIR"' EXIT
curl -fL --retry 5 --retry-delay 2 \\
  "https://codeload.github.com/sokrypton/ColabDesign/tar.gz/$PINNED_COMMIT" \\
  -o "$TMP_DIR/colabdesign.tar.gz"
tar -xzf "$TMP_DIR/colabdesign.tar.gz" -C "$TMP_DIR"
EXTRACTED_DIR="$(find "$TMP_DIR" -mindepth 1 -maxdepth 1 -type d -name 'ColabDesign-*' | head -n 1)"
if [[ -z "$EXTRACTED_DIR" ]]; then
  echo "FAIL: downloaded archive did not contain the expected ColabDesign source directory." >&2
  exit 2
fi
mv "$EXTRACTED_DIR" "$SOURCE_DIR"
printf '%s\n' "$PINNED_COMMIT" > "$COMMIT_MARKER"
echo "Prepared pinned ColabDesign gamma source at $SOURCE_DIR"
"""


def _python_overlay_script(afcycdesign_python: str, python_overlay: str) -> str:
    return f"""#!/bin/bash
set -euo pipefail

PYTHON_BIN="${{AFCYCDESIGN_PYTHON:-{afcycdesign_python}}}"
OVERLAY_DIR="${{AFCYCDESIGN_PYTHON_OVERLAY:-{python_overlay}}}"

if [[ ! -x "$PYTHON_BIN" ]]; then
  echo "FAIL: AfCycDesign Python is missing or not executable." >&2
  exit 2
fi

mkdir -p "$OVERLAY_DIR"
if PYTHONPATH="$OVERLAY_DIR" "$PYTHON_BIN" -c 'import IPython' >/dev/null 2>&1; then
  echo "PASS: Stage 5 Python overlay already provides IPython."
  exit 0
fi

"$PYTHON_BIN" -m pip install --upgrade --target "$OVERLAY_DIR" "IPython==8.37.0"
PYTHONPATH="$OVERLAY_DIR" "$PYTHON_BIN" -c 'import IPython; print("Prepared Stage 5 IPython overlay:", IPython.__version__)'
"""


def _preflight_script(
    afcycdesign_python: str,
    colabdesign_source: str,
    python_overlay: str,
    af_params: str,
) -> str:
    return f"""#!/bin/bash
set -euo pipefail

PYTHON_BIN="${{AFCYCDESIGN_PYTHON:-{afcycdesign_python}}}"
SOURCE_DIR="${{COLABDESIGN_GAMMA_SOURCE:-{colabdesign_source}}}"
OVERLAY_DIR="${{AFCYCDESIGN_PYTHON_OVERLAY:-{python_overlay}}}"
AF_PARAMS_DIR="${{AF_PARAMS:-{af_params}}}"

echo "[Stage5 preflight] python: $PYTHON_BIN"
echo "[Stage5 preflight] source: $SOURCE_DIR"
echo "[Stage5 preflight] overlay: $OVERLAY_DIR"
echo "[Stage5 preflight] params: $AF_PARAMS_DIR"

if [[ ! -x "$PYTHON_BIN" ]]; then
  echo "FAIL: AfCycDesign Python is missing or not executable." >&2
  exit 2
fi
if [[ ! -d "$AF_PARAMS_DIR" ]]; then
  echo "FAIL: AlphaFold parameter directory is missing." >&2
  exit 2
fi
if [[ ! -d "$OVERLAY_DIR" ]]; then
  echo "FAIL: Stage 5 Python overlay is missing. Run prepare_stage5_afcycdesign_python_overlay.sh first." >&2
  exit 2
fi
if [[ ! -f "$SOURCE_DIR/colabdesign/af/contrib/predict.py" ]] || [[ ! -f "$SOURCE_DIR/colabdesign/af/contrib/cyclic.py" ]]; then
  echo "FAIL: pinned ColabDesign gamma prediction source is missing. Run prepare_pinned_colabdesign_gamma_source.sh first." >&2
  exit 2
fi

COMMIT_MARKER="$SOURCE_DIR/.stage5_colabdesign_commit"
if [[ ! -f "$COMMIT_MARKER" ]]; then
  echo "FAIL: Stage 5 ColabDesign commit marker is missing." >&2
  exit 2
fi
current_commit="$(tr -d '[:space:]' < "$COMMIT_MARKER")"
if [[ "$current_commit" != "{COLABDESIGN_GAMMA_COMMIT}" ]]; then
  echo "FAIL: ColabDesign source commit is $current_commit, expected {COLABDESIGN_GAMMA_COMMIT}." >&2
  exit 2
fi

export PYTHONPATH="$SOURCE_DIR:$OVERLAY_DIR${{PYTHONPATH:+:$PYTHONPATH}}"

"$PYTHON_BIN" - <<'PY'
import importlib.util
import sys

required = [
    "colabdesign",
    "colabdesign.af.contrib.predict",
    "colabdesign.af.contrib.cyclic",
]
missing = [name for name in required if importlib.util.find_spec(name) is None]
if missing:
    print("FAIL: missing modules: " + ", ".join(missing), file=sys.stderr)
    raise SystemExit(3)

from colabdesign.af.contrib.cyclic import add_cyclic_offset
print("PASS: gamma prediction modules and add_cyclic_offset are importable.")
print("Protocol requirement: template_mode=none; use_initial_guess=false; cyclic chain index=1.")
PY

echo "PASS: Stage 5 AfCycDesign protocol preflight completed."
"""


def _job_script(
    *,
    project_root: Path,
    runner: Path,
    job_spec: Path,
    afcycdesign_python: str,
    colabdesign_source: str,
    python_overlay: str,
    af_params: str,
) -> str:
    return f"""#!/bin/bash
set -euo pipefail

if [[ "${{RUN_STAGE5_PREDICTIONS:-NO}}" != "YES" ]]; then
  echo "Stage 5 prediction is review-gated. Set RUN_STAGE5_PREDICTIONS=YES only after inspecting the manifest and protocol audit." >&2
  exit 3
fi

PYTHON_BIN="${{AFCYCDESIGN_PYTHON:-{afcycdesign_python}}}"
SOURCE_DIR="${{COLABDESIGN_GAMMA_SOURCE:-{colabdesign_source}}}"
OVERLAY_DIR="${{AFCYCDESIGN_PYTHON_OVERLAY:-{python_overlay}}}"
AF_PARAMS_DIR="${{AF_PARAMS:-{af_params}}}"
export COLABDESIGN_GAMMA_SOURCE="$SOURCE_DIR"
export PYTHONPATH="$SOURCE_DIR:$OVERLAY_DIR${{PYTHONPATH:+:$PYTHONPATH}}"

cd {_quote_bash(_to_wsl_path(project_root))}
"$PYTHON_BIN" {_quote_bash(_to_wsl_path(runner))} \\
  --job-spec {_quote_bash(_to_wsl_path(job_spec))} \\
  --af-params "$AF_PARAMS_DIR"
"""


def _master_script(job_scripts: list[Path], preflight: Path) -> str:
    lines = [
        "#!/bin/bash",
        "set -euo pipefail",
        "",
        "if [[ \"${RUN_STAGE5_PREDICTIONS:-NO}\" != \"YES\" ]]; then",
        '  echo "Stage 5 prediction is review-gated. No prediction has been started." >&2',
        '  echo "After review: export RUN_STAGE5_PREDICTIONS=YES" >&2',
        "  exit 3",
        "fi",
        "",
        f"bash {_quote_bash(_to_wsl_path(preflight))}",
        "",
    ]
    for script in job_scripts:
        lines.append(f"bash {_quote_bash(_to_wsl_path(script))}")
    return "\n".join(lines) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description="Prepare merged Stage 5 AfCycDesign independent-recovery jobs.")
    parser.add_argument(
        "--source-run-root",
        action="append",
        required=True,
        help="Upstream run root containing 06_rosetta_scoring. Repeat for each source run.",
    )
    parser.add_argument("--stage0-root", required=True)
    parser.add_argument("--output-root", required=True)
    parser.add_argument("--project-config", required=True)
    parser.add_argument("--candidate-count", type=int, default=5)
    parser.add_argument("--seeds-per-candidate", type=int, default=5)
    parser.add_argument("--models-per-seed", type=int, default=5)
    parser.add_argument("--recycles", type=int, default=6)
    parser.add_argument("--expected-target-length", type=int, default=86)
    parser.add_argument(
        "--afcycdesign-python",
        default="$HOME/fga_model_envs/colabdesign-py310/.pixi/envs/default/bin/python",
    )
    parser.add_argument(
        "--colabdesign-source",
        default="$HOME/fga_model_envs/sources/ColabDesign-gamma-stage5",
    )
    parser.add_argument(
        "--python-overlay",
        default="$HOME/fga_model_envs/stage5_afcycdesign_python_overlay",
    )
    parser.add_argument("--af-params", default="$HOME/fga_model_envs/af_params")
    args = parser.parse_args()

    logger = setup_logger("26_prepare_afcycdesign_jobs")
    append_run_header(logger, "26_prepare_afcycdesign_jobs.py")

    if args.candidate_count < 2:
        raise RuntimeError("--candidate-count must be >= 2")
    if args.seeds_per_candidate < 1:
        raise RuntimeError("--seeds-per-candidate must be >= 1")
    if args.models_per_seed != 5:
        raise RuntimeError("Stage 5 independent recovery requires all 5 model parameter sets; use --models-per-seed=5.")
    if args.recycles < 1:
        raise RuntimeError("--recycles must be >= 1")
    if args.expected_target_length < 1:
        raise RuntimeError("--expected-target-length must be >= 1")

    project_root = resolve_path(".")
    source_roots = [_resolve_mixed_path(value) for value in args.source_run_root]
    if len(source_roots) != len({str(path.resolve()) for path in source_roots}):
        raise RuntimeError("--source-run-root contains duplicate source roots")
    all_rows: list[dict[str, Any]] = []
    source_manifests: list[dict[str, Any]] = []
    source_route_records: list[dict[str, Any]] = []
    loaded_source_manifests: list[dict[str, Any]] = []
    for root in source_roots:
        assert_active_route_path(root, "Stage 26 source run root")
        manifest_path, source_manifest, source_manifest_sha256 = load_route_manifest(root)
        validate_route_project_config(args.project_config, source_manifest)
        source_run_id = str(source_manifest["run_id"])
        batch = str(source_manifest["batch_id"])
        source_route_records.append(
            {"manifest": source_manifest, "manifest_sha256": source_manifest_sha256}
        )
        source_provenance = route_provenance_fields(manifest_path, source_manifest, source_manifest_sha256)
        stage4_csv = root / "06_rosetta_scoring" / "FGA_rfpeptides_stage4_rosetta_interface_scores_pass.csv"
        assert_active_route_path(stage4_csv, f"Stage 26 Stage 4 pass CSV for {batch}")
        source_rows = _read_required_csv(stage4_csv)
        for row in source_rows:
            validate_row_route_provenance(row, source_provenance, f"Stage 26 Stage 4 row in {batch}")
        eligible_rows = _eligible_stage4_rows(source_rows, batch)
        for row in eligible_rows:
            row.update(
                {
                    "source_run_id": source_run_id,
                    "source_batch_id": batch,
                    "source_route_manifest": str(manifest_path),
                    "source_route_manifest_sha256": source_manifest_sha256,
                }
            )
        all_rows.extend(eligible_rows)
        source_manifests.append(
            {
                "run_id": source_run_id,
                "batch_id": batch,
                "manifest_path": str(manifest_path),
                "manifest_sha256": source_manifest_sha256,
            }
        )
        loaded_source_manifests.append(source_manifest)
    _validate_source_route_set(source_route_records)
    if not all_rows:
        raise RuntimeError("No eligible Stage 4 hard-QC rows were found.")

    selected = _select_candidates(all_rows, args.candidate_count)
    output_root = _resolve_mixed_path(args.output_root)
    assert_active_route_path(output_root, "Stage 26 output root", must_exist=False)
    output_dir = output_root / "07_structure_validation"
    inputs_dir = output_dir / "inputs"
    design_dir = inputs_dir / "design_poses"
    fasta_dir = inputs_dir / "fastas"
    specs_dir = inputs_dir / "job_specs"
    jobs_dir = output_dir / "jobs"
    target_dir = inputs_dir / "target"

    stage0_root = _resolve_mixed_path(args.stage0_root)
    assert_active_route_path(stage0_root, "Stage 26 Stage 0 root")
    source_target_pdb = stage0_root / "00_target_inputs" / "RFpep_Site_2_target.pdb"
    assert_active_route_path(source_target_pdb, "Stage 26 Stage 0 target PDB")
    if not source_target_pdb.exists():
        raise RuntimeError(f"Missing Stage 0 target PDB: {source_target_pdb}")
    target_pdb = target_dir / source_target_pdb.name
    target_pdb.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source_target_pdb, target_pdb)
    target_sequence = _validate_sequence(
        _pdb_chain_sequence(target_pdb, "A"),
        "Stage 0 target PDB chain A sequence",
        args.expected_target_length,
    )
    mapping_csv = stage0_root / "00_target_inputs" / "RFpep_Site_2_crop_renumbering_mapping.csv"
    assert_active_route_path(mapping_csv, "Stage 26 Stage 0 mapping CSV")
    site2_indices, hotspot_indices, mapping_target_sequence = _load_site_indices(
        mapping_csv,
        args.expected_target_length,
    )
    if mapping_target_sequence != target_sequence:
        raise RuntimeError("Stage 0 mapping sequence does not match Stage 0 target PDB chain A sequence")
    target_sequence_hash = _sequence_hash(target_sequence)
    use_mlm = True
    mlm_replace_fraction = 0.15
    protocol_payload = {
        "protocol_version": PROTOCOL_VERSION,
        "colabdesign_commit": COLABDESIGN_GAMMA_COMMIT,
        "model_type": "alphafold2_multimer_v3",
        "target_msa_mode": "single_sequence",
        "peptide_msa_mode": "single_sequence",
        "use_mlm": use_mlm,
        "mlm_replace_fraction": mlm_replace_fraction,
        "cyclic_chain_index": 1,
        "template_mode": "none",
        "use_initial_guess": False,
        "models_per_seed": args.models_per_seed,
        "requested_recycles": args.recycles,
    }
    protocol_hash = _protocol_hash(protocol_payload)
    base_manifest = loaded_source_manifests[0]
    merged_source_identity = canonical_json_sha256(sorted(source_manifests, key=lambda item: item["run_id"]))
    merged_manifest_path, merged_manifest, merged_manifest_sha256 = write_route_manifest(
        output_root,
        {
            "batch_id": f"merged_{merged_source_identity[:12]}",
            "site_labels": list(base_manifest["site_labels"]),
            "protocol_peptide_length_min": int(base_manifest["protocol_peptide_length_min"]),
            "protocol_peptide_length_max": int(base_manifest["protocol_peptide_length_max"]),
            "run_peptide_length_min": min(int(item["run_peptide_length_min"]) for item in loaded_source_manifests),
            "run_peptide_length_max": max(int(item["run_peptide_length_max"]) for item in loaded_source_manifests),
            "num_designs_requested": sum(int(item["num_designs_requested"]) for item in loaded_source_manifests),
            "project_config": base_manifest["project_config"],
            "project_config_sha256": base_manifest["project_config_sha256"],
            "effective_project_config_sha256": base_manifest["effective_project_config_sha256"],
            "stage0_sites": list(base_manifest["stage0_sites"]),
            "source_route_manifests": sorted(source_manifests, key=lambda item: item["run_id"]),
            "merged_source_identity_sha256": merged_source_identity,
            "stage5_protocol": protocol_payload,
        },
    )
    route_provenance = route_provenance_fields(
        merged_manifest_path,
        merged_manifest,
        merged_manifest_sha256,
    )
    _write_text_lf(target_dir / "RFpep_Site_2_target.fasta", f">RFpep_Site_2_target\n{target_sequence}")

    runner = project_root / "scripts" / "external" / "run_afcycdesign_independent_recovery.py"
    if not runner.exists():
        raise RuntimeError(f"Missing Stage 5 prediction runner: {runner}")

    manifest_rows: list[dict[str, Any]] = []
    job_rows: list[dict[str, Any]] = []
    job_scripts: list[Path] = []
    selected.sort(key=lambda row: (str(row["batch"]), _candidate_sort_key(row)))
    for selection_order, row in enumerate(selected, start=1):
        source_pdb = _resolve_mixed_path(str(row.get("scored_pdb", "")))
        assert_active_route_path(source_pdb, "Stage 26 Stage 4 scored PDB")
        if not source_pdb.exists():
            raise RuntimeError(f"Missing Stage 4 scored PDB: {source_pdb}")
        sequence = _validate_sequence(str(row.get("peptide_sequence", "")), "Stage 4 peptide_sequence")
        recorded_length = _parse_int(row.get("peptide_length", ""), len(sequence))
        if recorded_length != len(sequence):
            raise RuntimeError(
                f"Stage 4 peptide_length={recorded_length} does not match peptide_sequence length={len(sequence)} "
                f"for {row.get('stage4_design_id', '')}"
            )
        design_target_sequence = _validate_sequence(
            _pdb_chain_sequence(source_pdb, "A"),
            f"Stage 4 design target chain A sequence for {source_pdb.name}",
            args.expected_target_length,
        )
        design_peptide_sequence = _validate_sequence(
            _pdb_chain_sequence(source_pdb, "B"),
            f"Stage 4 design peptide chain B sequence for {source_pdb.name}",
            len(sequence),
        )
        if design_target_sequence != target_sequence:
            raise RuntimeError(f"Stage 4 design target sequence does not match Stage 0 target sequence: {source_pdb}")
        if design_peptide_sequence != sequence:
            raise RuntimeError(
                f"Stage 4 CSV peptide_sequence does not match design PDB chain B sequence: {source_pdb}"
            )
        peptide_sequence_hash = _sequence_hash(sequence)
        candidate_id = (
            f"S5V2_{selection_order:02d}_{row['batch']}_{_safe_token(str(row['backbone_id']))}_"
            f"seq{peptide_sequence_hash}"
        )
        staged_pdb = design_dir / f"{candidate_id}_design_pose.pdb"
        staged_pdb.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source_pdb, staged_pdb)

        cyclic_notation = f"{target_sequence}:({sequence})"
        fasta = fasta_dir / f"{candidate_id}.fasta"
        _write_text_lf(
            fasta,
            f">{candidate_id}|mode=independent_recovery|cyclic_chain=2|template=none|initial_guess=false\n{cyclic_notation}",
        )

        selection_reason = _selection_reason(row)
        manifest_rows.append(
            {
                "stage5_candidate_id": candidate_id,
                "peptide_sequence_hash": peptide_sequence_hash,
                "target_sequence_hash": target_sequence_hash,
                "protocol_hash": protocol_hash,
                "protocol_version": PROTOCOL_VERSION,
                "selection_order": selection_order,
                "batch": row["batch"],
                "backbone_id": row.get("backbone_id", ""),
                "source_stage4_design_id": row.get("stage4_design_id", ""),
                "peptide_sequence": sequence,
                "peptide_length": len(sequence),
                "site_label": row.get("site_label", ""),
                "site_id": row.get("site_id", ""),
                "design_pdb": source_pdb,
                "staged_design_pdb": staged_pdb,
                "target_pdb": target_pdb,
                "target_chain": "A",
                "peptide_chain": "B",
                "complex_rosetta_total_score": row.get("complex_rosetta_total_score", ""),
                "separated_rosetta_total_score": row.get("separated_rosetta_total_score", ""),
                "ddg_proxy_no_repack": row.get("ddg_proxy_no_repack", ""),
                "ddg_proxy_per_peptide_residue": row.get("ddg_proxy_per_peptide_residue", ""),
                "ddg_proxy_per_target_contact": row.get("ddg_proxy_per_target_contact", ""),
                "ddg_proxy_per_site_contact": row.get("ddg_proxy_per_site_contact", ""),
                "num_target_contacts": row.get("num_target_contacts", ""),
                "num_target_site_contacts": row.get("num_target_site_contacts", ""),
                "num_hotspot_contacts": row.get("num_hotspot_contacts", ""),
                "peptide_hotspot_min_distance": row.get("peptide_hotspot_min_distance", ""),
                "macrocycle_terminal_cn_distance": row.get("macrocycle_terminal_cn_distance", ""),
                "clash_status": row.get("clash_status", ""),
                "detached_or_collapsed_flag": row.get("detached_or_collapsed_flag", ""),
                "sequence_liability_notes": row.get("sequence_liability_notes", ""),
                "stage4_priority_rank": row.get("stage4_priority_rank", ""),
                "stage4_priority_class": row.get("stage4_priority_class", ""),
                "selection_role": row.get("selection_role", ""),
                "selection_reason": selection_reason,
                "prediction_protocol": "afcycdesign_gamma_cyclic_offset_no_template",
                "validation_test_type": "independent_recovery",
                "uses_design_pose_as_model_input": "false",
                "uses_template": "false",
                "uses_initial_guess": "false",
                "cyclic_topology_encoding": "peptide_chain_relative_position_cyclic_offset",
                "seeds_per_candidate": args.seeds_per_candidate,
                "models_per_seed": args.models_per_seed,
                "recycles": args.recycles,
                "requested_recycles": args.recycles,
                "forward_passes": args.recycles + 1,
                "target_msa_mode": "single_sequence",
                "peptide_msa_mode": "single_sequence",
                "use_mlm": str(use_mlm).lower(),
                "mlm_replace_fraction": mlm_replace_fraction,
                "prediction_job_count": args.seeds_per_candidate,
                "status": "prepared_pending_protocol_preflight_and_manual_review",
                **{key: row[key] for key in SOURCE_ROUTE_PROVENANCE_FIELDS},
                **route_provenance,
            }
        )

        for seed in range(args.seeds_per_candidate):
            job_id = f"{candidate_id}_prot{protocol_hash}_seed{seed:02d}"
            output_path = output_dir / "predictions" / candidate_id / f"protocol_{protocol_hash}" / f"seed_{seed:02d}"
            job_spec = specs_dir / f"{job_id}.json"
            spec = {
                "stage5_job_id": job_id,
                "stage5_candidate_id": candidate_id,
                "peptide_sequence_hash": peptide_sequence_hash,
                "target_sequence_hash": target_sequence_hash,
                "protocol_hash": protocol_hash,
                "protocol_version": PROTOCOL_VERSION,
                "batch": row["batch"],
                "backbone_id": row.get("backbone_id", ""),
                "target_sequence": target_sequence,
                "target_sequence_length": len(target_sequence),
                "peptide_sequence": sequence,
                "peptide_sequence_length": len(sequence),
                "site2_target_indices_1based": site2_indices,
                "hotspot_target_indices_1based": hotspot_indices,
                "cyclic_chain_index": 1,
                "seed": seed,
                "models": "all",
                "models_per_seed": args.models_per_seed,
                "recycles": args.recycles,
                "requested_recycles": args.recycles,
                "forward_passes": args.recycles + 1,
                "target_msa_mode": "single_sequence",
                "peptide_msa_mode": "single_sequence",
                "msa_rows_input": 1,
                "use_mlm": use_mlm,
                "mlm_replace_fraction": mlm_replace_fraction,
                "template_mode": "none",
                "use_initial_guess": False,
                "validation_test_type": "independent_recovery",
                "design_pose_for_posthoc_comparison": _to_wsl_path(staged_pdb),
                "prediction_output_dir": _to_wsl_path(output_path),
                "colabdesign_commit": COLABDESIGN_GAMMA_COMMIT,
                **{key: row[key] for key in SOURCE_ROUTE_PROVENANCE_FIELDS},
                **route_provenance,
            }
            _write_text_lf(job_spec, json.dumps(spec, indent=2, sort_keys=True))

            run_script = jobs_dir / f"run_{job_id}.sh"
            _write_text_lf(
                run_script,
                _job_script(
                    project_root=project_root,
                    runner=runner,
                    job_spec=job_spec,
                    afcycdesign_python=args.afcycdesign_python,
                    colabdesign_source=args.colabdesign_source,
                    python_overlay=args.python_overlay,
                    af_params=args.af_params,
                ),
            )
            job_scripts.append(run_script)
            job_rows.append(
                {
                    "stage5_job_id": job_id,
                    "stage5_candidate_id": candidate_id,
                    "peptide_sequence_hash": peptide_sequence_hash,
                    "protocol_hash": protocol_hash,
                    "batch": row["batch"],
                    "backbone_id": row.get("backbone_id", ""),
                    "seed": seed,
                    "target_sequence_length": len(target_sequence),
                    "target_sequence_hash": target_sequence_hash,
                    "peptide_sequence": sequence,
                    "peptide_length": len(sequence),
                    "cyclic_sequence_notation": cyclic_notation,
                    "input_fasta": fasta,
                    "job_spec_json": job_spec,
                    "design_pose_for_posthoc_comparison": staged_pdb,
                    "prediction_output_dir": output_path,
                    "run_script": run_script,
                    "prediction_protocol": "afcycdesign_gamma_cyclic_offset_no_template",
                    "validation_test_type": "independent_recovery",
                    "uses_template": "false",
                    "uses_initial_guess": "false",
                    "target_msa_mode": "single_sequence",
                    "peptide_msa_mode": "single_sequence",
                    "use_mlm": str(use_mlm).lower(),
                    "mlm_replace_fraction": mlm_replace_fraction,
                    "requested_recycles": args.recycles,
                    "forward_passes": args.recycles + 1,
                    "status": "prepared_not_run",
                    **{key: row[key] for key in SOURCE_ROUTE_PROVENANCE_FIELDS},
                    **route_provenance,
                }
            )

    manifest_csv = output_dir / "FGA_rfpeptides_stage5_candidate_manifest.csv"
    jobs_csv = inputs_dir / "FGA_rfpeptides_stage5_prediction_jobs.csv"
    protocol_audit = output_dir / "FGA_rfpeptides_stage5_protocol_audit.md"
    fetch_protocol = jobs_dir / "prepare_pinned_colabdesign_gamma_source.sh"
    prepare_overlay = jobs_dir / "prepare_stage5_afcycdesign_python_overlay.sh"
    preflight = jobs_dir / "check_stage5_afcycdesign_protocol.sh"
    master_script = jobs_dir / "run_stage5_afcycdesign_independent_recovery_all.sh"
    write_csv(manifest_csv, manifest_rows, MANIFEST_FIELDS)
    write_csv(jobs_csv, job_rows, JOB_FIELDS)
    write_markdown(
        protocol_audit,
        _protocol_audit_markdown(
            output_dir=output_dir,
            target_pdb=target_pdb,
            target_sequence=target_sequence,
            candidate_count=len(manifest_rows),
            seeds_per_candidate=args.seeds_per_candidate,
            models_per_seed=args.models_per_seed,
            recycles=args.recycles,
            protocol_hash=protocol_hash,
            use_mlm=use_mlm,
            mlm_replace_fraction=mlm_replace_fraction,
            afcycdesign_python=args.afcycdesign_python,
            colabdesign_source=args.colabdesign_source,
            python_overlay=args.python_overlay,
            af_params=args.af_params,
        ),
    )
    _write_text_lf(fetch_protocol, _fetch_protocol_script(args.colabdesign_source))
    _write_text_lf(
        prepare_overlay,
        _python_overlay_script(args.afcycdesign_python, args.python_overlay),
    )
    _write_text_lf(
        preflight,
        _preflight_script(
            args.afcycdesign_python,
            args.colabdesign_source,
            args.python_overlay,
            args.af_params,
        ),
    )
    _write_text_lf(master_script, _master_script(job_scripts, preflight))

    eligible_counts: dict[tuple[str, str], int] = defaultdict(int)
    for row in all_rows:
        eligible_counts[(str(row["batch"]), str(row["backbone_id"]))] += 1
    write_markdown(
        output_dir / "FGA_rfpeptides_stage5_candidate_manifest.md",
        _manifest_markdown(
            rows=manifest_rows,
            eligible_counts=eligible_counts,
            output_dir=output_dir,
            protocol_audit=protocol_audit,
            jobs_csv=jobs_csv,
        ),
    )

    logger.info("Eligible Stage 4 hard-QC rows: %s", len(all_rows))
    logger.info("Selected Stage 5 candidates: %s", len(manifest_rows))
    logger.info("Prepared seed-level jobs: %s", len(job_rows))
    logger.info("Output directory: %s", output_dir)
    logger.info("No Stage 5 structure prediction was run.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
