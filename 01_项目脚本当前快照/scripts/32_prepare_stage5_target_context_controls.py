from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shlex
from pathlib import Path
from typing import Any, Iterable, Mapping

from common import (
    ROUTE_PROVENANCE_FIELDS,
    assert_active_route_path,
    append_run_header,
    canonical_json_sha256,
    load_route_manifest,
    read_csv,
    resolve_path,
    route_provenance_fields,
    setup_logger,
    sha256_file,
    validate_route_project_config,
    write_csv,
    write_markdown,
    write_route_manifest,
)


COLABDESIGN_GAMMA_COMMIT = "5ab4efaba2321a6c3c314b82d2fff8e0241f5c2d"
PROTOCOL_VERSION = "stage5_target_context_controls_v1"
MODEL_NAMES = [f"model_{index}_multimer_v3" for index in range(1, 6)]
LEGAL_AA = set("ACDEFGHIKLMNPQRSTVWY")
AA3_TO_AA1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "MSE": "M", "PHE": "F",
    "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y",
    "VAL": "V",
}

CONTROL_FIELDS = [
    "context_id",
    "context_label",
    "context_description",
    "context_pdb",
    "context_pdb_sha256",
    "context_mapping_csv",
    "context_mapping_csv_sha256",
    "template_mode",
    "template_sequence_masked",
    "template_sidechains_masked",
    "template_interchain_features_masked",
    "chain_ids",
    "chain_lengths",
    "total_residues",
    "fga_residue_range",
    "partner_segments",
    "seeds_planned",
    "models_per_seed",
    "model_predictions_planned",
    "protocol_hash",
    "status",
    "notes",
] + ROUTE_PROVENANCE_FIELDS

JOB_FIELDS = [
    "stage5_context_job_id",
    "context_id",
    "protocol_hash",
    "seed",
    "requested_recycles",
    "forward_passes",
    "models_per_seed",
    "context_pdb",
    "context_pdb_sha256",
    "job_spec_json",
    "prediction_output_dir",
    "run_script",
    "template_mode",
    "template_sequence_masked",
    "template_sidechains_masked",
    "template_interchain_features_masked",
    "target_msa_mode",
    "peptide_included",
    "use_initial_guess",
    "use_dropout",
    "status",
] + ROUTE_PROVENANCE_FIELDS


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


def _write_text(path: Path, value: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(value.rstrip() + "\n", encoding="utf-8", newline="\n")


def _sha1_text(value: str) -> str:
    return hashlib.sha1(value.encode("utf-8")).hexdigest()


def _truthy(value: Any) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes"}


def _parse_source_pdb(path: Path) -> dict[str, dict[int, dict[str, Any]]]:
    chains: dict[str, dict[int, dict[str, Any]]] = {}
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            if not raw_line.startswith("ATOM  "):
                continue
            line = raw_line.rstrip("\r\n")
            if len(line) < 54:
                continue
            altloc = line[16]
            if altloc not in {" ", "A"}:
                continue
            chain = line[21].strip() or "_"
            insertion = line[26].strip()
            if insertion:
                raise RuntimeError(f"Insertion code is not supported in context extraction: {chain}{line[22:27].strip()}")
            residue_number = int(line[22:26])
            resname = line[17:20].strip().upper()
            if resname not in AA3_TO_AA1:
                raise RuntimeError(f"Unsupported residue {resname} in {path}")
            residue = chains.setdefault(chain, {}).setdefault(
                residue_number,
                {"resname": resname, "lines": []},
            )
            if residue["resname"] != resname:
                raise RuntimeError(f"Conflicting residue names for {chain}{residue_number}")
            if altloc == "A":
                chars = list(line)
                chars[16] = " "
                line = "".join(chars)
            residue["lines"].append(line)
    return chains


def _extract_context(
    *,
    source: Mapping[str, Mapping[int, Mapping[str, Any]]],
    segments: Iterable[tuple[str, int, int, str, str]],
    output_pdb: Path,
    mapping_csv: Path,
    site_original_numbers: set[int],
    hotspot_original_numbers: set[int],
) -> tuple[list[str], list[str], list[int]]:
    output_lines: list[str] = []
    mapping_rows: list[dict[str, Any]] = []
    chain_ids: list[str] = []
    sequences: list[str] = []
    lengths: list[int] = []
    serial = 1
    for original_chain, start, end, output_chain, role in segments:
        if output_chain in chain_ids:
            raise RuntimeError(f"Duplicate output chain in context: {output_chain}")
        chain_ids.append(output_chain)
        residues = source.get(original_chain, {})
        missing = [number for number in range(start, end + 1) if number not in residues]
        if missing:
            raise RuntimeError(
                f"Context segment {original_chain}{start}-{end} is incomplete; missing {missing[:10]}"
            )
        sequence: list[str] = []
        for output_number, original_number in enumerate(range(start, end + 1), start=1):
            residue = residues[original_number]
            resname = str(residue["resname"])
            aa = AA3_TO_AA1[resname]
            sequence.append(aa)
            has_ca = False
            for original_line in residue["lines"]:
                atom_name = original_line[12:16].strip()
                has_ca = has_ca or atom_name == "CA"
                chars = list(original_line.ljust(80))
                chars[6:11] = list(f"{serial:5d}")
                chars[21] = output_chain
                chars[22:26] = list(f"{output_number:4d}")
                chars[26] = " "
                output_lines.append("".join(chars).rstrip())
                serial += 1
            if not has_ca:
                raise RuntimeError(f"Context residue {original_chain}{original_number} lacks CA")
            is_fga = role == "fga"
            mapping_rows.append(
                {
                    "context_chain": output_chain,
                    "context_residue_number": output_number,
                    "context_residue_name": resname,
                    "context_residue_aa": aa,
                    "context_role": role,
                    "original_pdb_id": "3GHG",
                    "original_chain_id": original_chain,
                    "original_pdb_residue_number": original_number,
                    "is_fga_chain": str(is_fga).lower(),
                    "is_target_site_residue": str(is_fga and original_number in site_original_numbers).lower(),
                    "is_selected_hotspot": str(is_fga and original_number in hotspot_original_numbers).lower(),
                }
            )
        sequences.append("".join(sequence))
        lengths.append(len(sequence))
        output_lines.append(
            f"TER   {serial:5d}      {str(residues[end]['resname']):>3s} {output_chain}{len(sequence):4d}"
        )
        serial += 1
    output_lines.append("END")
    _write_text(output_pdb, "\n".join(output_lines))
    fields = [
        "context_chain", "context_residue_number", "context_residue_name", "context_residue_aa",
        "context_role", "original_pdb_id", "original_chain_id", "original_pdb_residue_number",
        "is_fga_chain", "is_target_site_residue", "is_selected_hotspot",
    ]
    mapping_csv.parent.mkdir(parents=True, exist_ok=True)
    with mapping_csv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(mapping_rows)
    return chain_ids, sequences, lengths


def _stage0_original_positions(stage0_root: Path) -> tuple[Path, Path, set[int], set[int]]:
    target_pdb = assert_active_route_path(
        stage0_root / "00_target_inputs" / "RFpep_Site_2_target.pdb",
        "Stage 32 Stage 0 target PDB",
    )
    mapping_csv = assert_active_route_path(
        stage0_root / "00_target_inputs" / "RFpep_Site_2_crop_renumbering_mapping.csv",
        "Stage 32 Stage 0 mapping CSV",
    )
    rows = read_csv(mapping_csv)
    if len(rows) != 86:
        raise RuntimeError(f"Stage 0 Site_2 mapping must contain 86 rows, observed {len(rows)}")
    rows.sort(key=lambda row: int(row["rfpeptides_residue_number"]))
    if {str(row["original_chain_id"]) for row in rows} != {"G"}:
        raise RuntimeError("Stage 0 Site_2 crop must map entirely to original 3GHG chain G")
    original_numbers = [int(row["original_pdb_residue_number"]) for row in rows]
    if original_numbers != list(range(115, 201)):
        raise RuntimeError("Stage 0 Site_2 crop must map exactly to original 3GHG G115-G200")
    site = {int(row["original_pdb_residue_number"]) for row in rows if _truthy(row.get("is_target_site_residue"))}
    hotspots = {int(row["original_pdb_residue_number"]) for row in rows if _truthy(row.get("is_selected_hotspot"))}
    if not site or hotspots != {196, 198, 199, 200}:
        raise RuntimeError(f"Unexpected Stage 0 Site_2/hotspot mapping: site={sorted(site)}, hotspots={sorted(hotspots)}")
    return target_pdb, mapping_csv, site, hotspots


def _job_script(
    *,
    project_root: Path,
    runner: Path,
    job_spec: Path,
    python_bin: str,
    source_dir: str,
    overlay_dir: str,
    af_params: str,
) -> str:
    return f'''#!/bin/bash
set -euo pipefail

if [[ "${{RUN_STAGE5_CONTEXT_CONTROLS:-NO}}" != "YES" ]]; then
  echo "Target-context controls are review-gated. Set RUN_STAGE5_CONTEXT_CONTROLS=YES after inspection." >&2
  exit 3
fi

PYTHON_BIN="${{AFCYCDESIGN_PYTHON:-{python_bin}}}"
SOURCE_DIR="${{COLABDESIGN_GAMMA_SOURCE:-{source_dir}}}"
OVERLAY_DIR="${{AFCYCDESIGN_PYTHON_OVERLAY:-{overlay_dir}}}"
AF_PARAMS_DIR="${{AF_PARAMS:-{af_params}}}"
export COLABDESIGN_GAMMA_SOURCE="$SOURCE_DIR"
export PYTHONPATH="$SOURCE_DIR:$OVERLAY_DIR${{PYTHONPATH:+:$PYTHONPATH}}"

cd {shlex.quote(_to_wsl_path(project_root))}
"$PYTHON_BIN" {shlex.quote(_to_wsl_path(runner))} \
  --job-spec {shlex.quote(_to_wsl_path(job_spec))} \
  --af-params "$AF_PARAMS_DIR"
'''


def _preflight_script(
    *,
    project_root: Path,
    runner: Path,
    specs: Iterable[Path],
    python_bin: str,
    source_dir: str,
    overlay_dir: str,
    af_params: str,
) -> str:
    commands = "\n".join(
        f'"$PYTHON_BIN" {shlex.quote(_to_wsl_path(runner))} --job-spec '
        f'{shlex.quote(_to_wsl_path(spec))} --af-params "$AF_PARAMS_DIR" --preflight-only'
        for spec in specs
    )
    return f'''#!/bin/bash
set -euo pipefail
PYTHON_BIN="${{AFCYCDESIGN_PYTHON:-{python_bin}}}"
SOURCE_DIR="${{COLABDESIGN_GAMMA_SOURCE:-{source_dir}}}"
OVERLAY_DIR="${{AFCYCDESIGN_PYTHON_OVERLAY:-{overlay_dir}}}"
AF_PARAMS_DIR="${{AF_PARAMS:-{af_params}}}"
export COLABDESIGN_GAMMA_SOURCE="$SOURCE_DIR"
export PYTHONPATH="$SOURCE_DIR:$OVERLAY_DIR${{PYTHONPATH:+:$PYTHONPATH}}"
cd {shlex.quote(_to_wsl_path(project_root))}
{commands}
echo "PASS: all Stage 5 target-context control preflights completed."
'''


def _master_script(preflight: Path, scripts: Iterable[Path]) -> str:
    commands = "\n".join(f"bash {shlex.quote(_to_wsl_path(path))}" for path in scripts)
    return f'''#!/bin/bash
set -euo pipefail
bash {shlex.quote(_to_wsl_path(preflight))}
{commands}
'''


def _plan_markdown(controls: list[Mapping[str, Any]], output_dir: Path) -> str:
    rows = "\n".join(
        f"- `{row['context_id']}`: {row['context_description']} "
        f"({row['total_residues']} aa; template mode `{row['template_mode']}`)."
        for row in controls
    )
    return f"""# Stage 5 Expanded Target Context Controls

Purpose: determine whether target recovery improves when the 86-aa Site_2 crop
is replaced by a larger experimentally observed 3GHG context and when target
template information is retained.

These are target-only diagnostics. No peptide sequence, peptide template,
cyclic offset, or initial guess is used. They do not validate or reject any
peptide candidate.

## Controls

{rows}

Each control uses one seed and all five AlphaFold multimer-v3 parameter sets by
default. Predictions are deterministic (`dropout=false`) so this first screen
compares context definitions, not stochastic seed diversity.

The largest context is not full-length FGA. It is the complete experimentally
visible FGA chain G segment G27-G200 plus local native H118-H184 and I68-I127
segments from 3GHG.

Output root:

```text
{output_dir}
```

Prediction scripts remain gated until `RUN_STAGE5_CONTEXT_CONTROLS=YES` is set.
"""


def main() -> int:
    parser = argparse.ArgumentParser(description="Prepare small expanded-target AfCycDesign recovery controls.")
    parser.add_argument("--source-run-root", required=True)
    parser.add_argument("--stage0-root", required=True)
    parser.add_argument("--output-root", required=True)
    parser.add_argument("--project-config", required=True)
    parser.add_argument("--seeds-per-context", type=int, default=1)
    parser.add_argument("--models-per-seed", type=int, default=5)
    parser.add_argument("--recycles", type=int, default=6)
    parser.add_argument(
        "--afcycdesign-python",
        default="/home/luomi/fga_model_envs/colabdesign-py310/.pixi/envs/default/bin/python",
    )
    parser.add_argument(
        "--colabdesign-source",
        default="/home/luomi/fga_model_envs/sources/ColabDesign-gamma-stage5",
    )
    parser.add_argument(
        "--python-overlay",
        default="/home/luomi/fga_model_envs/stage5_afcycdesign_python_overlay",
    )
    parser.add_argument("--af-params", default="/home/luomi/fga_model_envs/af_params")
    args = parser.parse_args()

    logger = setup_logger("32_prepare_stage5_target_context_controls")
    append_run_header(logger, "32_prepare_stage5_target_context_controls.py")
    if args.seeds_per_context < 1:
        raise RuntimeError("--seeds-per-context must be positive")
    if args.models_per_seed != 5:
        raise RuntimeError("This diagnostic requires all five AlphaFold multimer-v3 model sets")
    if args.recycles < 1:
        raise RuntimeError("--recycles must be positive")

    source_root = assert_active_route_path(_resolve_mixed_path(args.source_run_root), "Stage 32 source run root")
    stage0_root = assert_active_route_path(_resolve_mixed_path(args.stage0_root), "Stage 32 Stage 0 root")
    output_root = assert_active_route_path(
        _resolve_mixed_path(args.output_root), "Stage 32 output root", must_exist=False
    )
    source_manifest_path, source_manifest, source_manifest_sha256 = load_route_manifest(source_root)
    config = validate_route_project_config(args.project_config, source_manifest)
    stage0_target, stage0_mapping, site_original, hotspot_original = _stage0_original_positions(stage0_root)
    source_site = source_manifest["stage0_sites"][0]
    if sha256_file(stage0_target) != str(source_site["target_pdb_sha256"]):
        raise RuntimeError("Supplied Stage 0 target PDB differs from the source route manifest")
    if sha256_file(stage0_mapping) != str(source_site["mapping_csv_sha256"]):
        raise RuntimeError("Supplied Stage 0 mapping differs from the source route manifest")

    cleaned_pdb = assert_active_route_path(
        _resolve_mixed_path(config["structures"]["cleaned_pdb_file"]),
        "Stage 32 cleaned 3GHG PDB",
    )
    source_pdb = _parse_source_pdb(cleaned_pdb)
    project_root = Path(__file__).resolve().parents[1]
    output_dir = output_root / "07_structure_validation_target_context_controls"
    contexts_dir = output_dir / "inputs" / "contexts"
    mappings_dir = output_dir / "inputs" / "mappings"
    specs_dir = output_dir / "inputs" / "job_specs"
    jobs_dir = output_dir / "jobs"
    predictions_dir = output_dir / "predictions"
    runner = project_root / "scripts" / "external" / "run_afcycdesign_target_context_control.py"
    if not runner.is_file():
        raise RuntimeError(f"Missing Stage 5 target-context runner: {runner}")

    definitions = [
        {
            "context_id": "C0_crop86_masked_template",
            "description": "Original G115-G200 86-aa crop with Stage 5B-style template masking",
            "segments": [("G", 115, 200, "A", "fga")],
            "mask": (True, True, True),
        },
        {
            "context_id": "C1_crop86_full_template",
            "description": "Original G115-G200 86-aa crop with target template sequence and sidechains retained",
            "segments": [("G", 115, 200, "A", "fga")],
            "mask": (False, False, False),
        },
        {
            "context_id": "C2_fga_visible174_full_template",
            "description": "Complete experimentally visible FGA chain G27-G200 with full template information",
            "segments": [("G", 27, 200, "A", "fga")],
            "mask": (False, False, False),
        },
        {
            "context_id": "C3_native_GHI301_full_template",
            "description": "Visible FGA G27-G200 plus local native FGB H118-H184 and FGG I68-I127 context",
            "segments": [
                ("G", 27, 200, "A", "fga"),
                ("H", 118, 184, "B", "native_partner_fgb"),
                ("I", 68, 127, "C", "native_partner_fgg"),
            ],
            "mask": (False, False, False),
        },
    ]

    context_records: list[dict[str, Any]] = []
    prepared: list[dict[str, Any]] = []
    for definition in definitions:
        context_id = str(definition["context_id"])
        context_pdb = contexts_dir / f"{context_id}.pdb"
        context_mapping = mappings_dir / f"{context_id}_mapping.csv"
        chain_ids, chain_sequences, chain_lengths = _extract_context(
            source=source_pdb,
            segments=definition["segments"],
            output_pdb=context_pdb,
            mapping_csv=context_mapping,
            site_original_numbers=site_original,
            hotspot_original_numbers=hotspot_original,
        )
        if any(set(sequence) - LEGAL_AA for sequence in chain_sequences):
            raise RuntimeError(f"Illegal amino acid in extracted context {context_id}")
        mapping_rows = read_csv(context_mapping)
        site_indices = [
            int(row["context_residue_number"])
            for row in mapping_rows
            if row["context_chain"] == "A" and _truthy(row["is_target_site_residue"])
        ]
        hotspot_indices = [
            int(row["context_residue_number"])
            for row in mapping_rows
            if row["context_chain"] == "A" and _truthy(row["is_selected_hotspot"])
        ]
        if not site_indices or len(hotspot_indices) != 4:
            raise RuntimeError(f"Extracted context {context_id} lost Site_2/hotspot mapping")
        rm_seq, rm_sc, rm_ic = definition["mask"]
        template_mode = "target_context_masked" if any(definition["mask"]) else "target_context_full"
        payload = {
            "protocol_version": PROTOCOL_VERSION,
            "context_id": context_id,
            "context_pdb_sha256": sha256_file(context_pdb),
            "context_mapping_csv_sha256": sha256_file(context_mapping),
            "chain_ids": chain_ids,
            "chain_lengths": chain_lengths,
            "chain_sequences": chain_sequences,
            "site2_fga_indices_1based": site_indices,
            "hotspot_fga_indices_1based": hotspot_indices,
            "template_mode": template_mode,
            "template_sequence_masked": rm_seq,
            "template_sidechains_masked": rm_sc,
            "template_interchain_features_masked": rm_ic,
            "target_msa_mode": "single_sequence",
            "peptide_included": False,
            "use_initial_guess": False,
            "use_dropout": False,
            "requested_recycles": args.recycles,
            "models_per_seed": args.models_per_seed,
            "colabdesign_commit": COLABDESIGN_GAMMA_COMMIT,
        }
        protocol_hash = canonical_json_sha256(payload)[:12]
        prepared.append(
            {
                "definition": definition,
                "context_pdb": context_pdb,
                "context_mapping": context_mapping,
                "chain_ids": chain_ids,
                "chain_sequences": chain_sequences,
                "chain_lengths": chain_lengths,
                "site_indices": site_indices,
                "hotspot_indices": hotspot_indices,
                "template_mode": template_mode,
                "rm_seq": rm_seq,
                "rm_sc": rm_sc,
                "rm_ic": rm_ic,
                "protocol_hash": protocol_hash,
            }
        )

    route_manifest_path, route_manifest, route_manifest_sha256 = write_route_manifest(
        output_root,
        {
            "batch_id": f"stage5_target_context_controls_{canonical_json_sha256([item['protocol_hash'] for item in prepared])[:12]}",
            "site_labels": list(source_manifest["site_labels"]),
            "protocol_peptide_length_min": source_manifest["protocol_peptide_length_min"],
            "protocol_peptide_length_max": source_manifest["protocol_peptide_length_max"],
            "run_peptide_length_min": source_manifest["run_peptide_length_min"],
            "run_peptide_length_max": source_manifest["run_peptide_length_max"],
            "num_designs_requested": source_manifest["num_designs_requested"],
            "project_config": source_manifest["project_config"],
            "project_config_sha256": source_manifest["project_config_sha256"],
            "effective_project_config_sha256": source_manifest["effective_project_config_sha256"],
            "stage0_sites": source_manifest["stage0_sites"],
            "source_route_manifests": [
                {
                    "run_id": source_manifest["run_id"],
                    "batch_id": source_manifest["batch_id"],
                    "manifest_path": str(source_manifest_path),
                    "manifest_sha256": source_manifest_sha256,
                }
            ],
            "stage5_target_context_control_protocol": {
                "protocol_version": PROTOCOL_VERSION,
                "contexts": [item["protocol_hash"] for item in prepared],
                "seeds_per_context": args.seeds_per_context,
                "models_per_seed": args.models_per_seed,
                "requested_recycles": args.recycles,
                "peptide_included": False,
            },
        },
    )
    provenance = route_provenance_fields(route_manifest_path, route_manifest, route_manifest_sha256)
    jobs: list[dict[str, Any]] = []
    specs: list[Path] = []
    run_scripts: list[Path] = []
    for item in prepared:
        definition = item["definition"]
        context_id = str(definition["context_id"])
        segments = definition["segments"]
        partner_segments = ",".join(
            f"{chain}{start}-{end}" for chain, start, end, _, role in segments if role != "fga"
        ) or "none"
        context_records.append(
            {
                "context_id": context_id,
                "context_label": context_id,
                "context_description": definition["description"],
                "context_pdb": item["context_pdb"],
                "context_pdb_sha256": sha256_file(item["context_pdb"]),
                "context_mapping_csv": item["context_mapping"],
                "context_mapping_csv_sha256": sha256_file(item["context_mapping"]),
                "template_mode": item["template_mode"],
                "template_sequence_masked": str(item["rm_seq"]).lower(),
                "template_sidechains_masked": str(item["rm_sc"]).lower(),
                "template_interchain_features_masked": str(item["rm_ic"]).lower(),
                "chain_ids": ",".join(item["chain_ids"]),
                "chain_lengths": ",".join(str(value) for value in item["chain_lengths"]),
                "total_residues": sum(item["chain_lengths"]),
                "fga_residue_range": f"G{segments[0][1]}-G{segments[0][2]}",
                "partner_segments": partner_segments,
                "seeds_planned": args.seeds_per_context,
                "models_per_seed": args.models_per_seed,
                "model_predictions_planned": args.seeds_per_context * args.models_per_seed,
                "protocol_hash": item["protocol_hash"],
                "status": "prepared_not_run",
                "notes": "target-only diagnostic; no peptide and no initial guess",
                **provenance,
            }
        )
        for seed in range(args.seeds_per_context):
            job_id = f"S5CTX_{context_id}_prot{item['protocol_hash']}_seed{seed:02d}"
            spec_path = specs_dir / f"{job_id}.json"
            prediction_dir = predictions_dir / context_id / f"protocol_{item['protocol_hash']}" / f"seed_{seed:02d}"
            spec = {
                "stage5_context_job_id": job_id,
                "context_id": context_id,
                "context_label": context_id,
                "context_description": definition["description"],
                "protocol_hash": item["protocol_hash"],
                "protocol_version": PROTOCOL_VERSION,
                "context_pdb": _to_wsl_path(item["context_pdb"]),
                "context_pdb_sha256": sha256_file(item["context_pdb"]),
                "context_mapping_csv": _to_wsl_path(item["context_mapping"]),
                "context_mapping_csv_sha256": sha256_file(item["context_mapping"]),
                "source_cleaned_pdb": _to_wsl_path(cleaned_pdb),
                "source_cleaned_pdb_sha256": sha256_file(cleaned_pdb),
                "chain_ids": item["chain_ids"],
                "chain_lengths": item["chain_lengths"],
                "chain_sequences": item["chain_sequences"],
                "fga_chain_index": 0,
                "site2_fga_indices_1based": item["site_indices"],
                "hotspot_fga_indices_1based": item["hotspot_indices"],
                "template_mode": item["template_mode"],
                "use_templates": True,
                "use_batch_as_template": False,
                "template_sequence_masked": item["rm_seq"],
                "template_sidechains_masked": item["rm_sc"],
                "template_interchain_features_masked": item["rm_ic"],
                "target_msa_mode": "single_sequence",
                "msa_rows_expected": 1,
                "peptide_included": False,
                "use_initial_guess": False,
                "use_dropout": False,
                "cyclic_topology_encoding": "not_applicable_no_peptide",
                "validation_test_type": "target_context_template_recovery_control",
                "seed": seed,
                "requested_recycles": args.recycles,
                "forward_passes": args.recycles + 1,
                "models_per_seed": args.models_per_seed,
                "model_names_expected": MODEL_NAMES,
                "colabdesign_commit": COLABDESIGN_GAMMA_COMMIT,
                "prediction_output_dir": _to_wsl_path(prediction_dir),
                **provenance,
            }
            _write_text(spec_path, json.dumps(spec, indent=2, sort_keys=True))
            run_script = jobs_dir / f"run_{job_id}.sh"
            _write_text(
                run_script,
                _job_script(
                    project_root=project_root,
                    runner=runner,
                    job_spec=spec_path,
                    python_bin=args.afcycdesign_python,
                    source_dir=args.colabdesign_source,
                    overlay_dir=args.python_overlay,
                    af_params=args.af_params,
                ),
            )
            specs.append(spec_path)
            run_scripts.append(run_script)
            jobs.append(
                {
                    "stage5_context_job_id": job_id,
                    "context_id": context_id,
                    "protocol_hash": item["protocol_hash"],
                    "seed": seed,
                    "requested_recycles": args.recycles,
                    "forward_passes": args.recycles + 1,
                    "models_per_seed": args.models_per_seed,
                    "context_pdb": item["context_pdb"],
                    "context_pdb_sha256": sha256_file(item["context_pdb"]),
                    "job_spec_json": spec_path,
                    "prediction_output_dir": prediction_dir,
                    "run_script": run_script,
                    "template_mode": item["template_mode"],
                    "template_sequence_masked": str(item["rm_seq"]).lower(),
                    "template_sidechains_masked": str(item["rm_sc"]).lower(),
                    "template_interchain_features_masked": str(item["rm_ic"]).lower(),
                    "target_msa_mode": "single_sequence",
                    "peptide_included": "false",
                    "use_initial_guess": "false",
                    "use_dropout": "false",
                    "status": "prepared_not_run",
                    **provenance,
                }
            )

    preflight = jobs_dir / "check_stage5_target_context_controls.sh"
    master = jobs_dir / "run_stage5_target_context_controls_all.sh"
    _write_text(
        preflight,
        _preflight_script(
            project_root=project_root,
            runner=runner,
            specs=specs,
            python_bin=args.afcycdesign_python,
            source_dir=args.colabdesign_source,
            overlay_dir=args.python_overlay,
            af_params=args.af_params,
        ),
    )
    _write_text(master, _master_script(preflight, run_scripts))
    write_csv(output_dir / "FGA_rfpeptides_stage5_target_context_control_manifest.csv", context_records, CONTROL_FIELDS)
    write_csv(output_dir / "inputs" / "FGA_rfpeptides_stage5_target_context_control_jobs.csv", jobs, JOB_FIELDS)
    write_markdown(
        output_dir / "FGA_rfpeptides_stage5_target_context_control_plan.md",
        _plan_markdown(context_records, output_dir),
    )
    logger.info("Expanded target contexts prepared: %s", len(context_records))
    logger.info("Seed jobs prepared: %s", len(jobs))
    logger.info("Planned model predictions: %s", len(jobs) * args.models_per_seed)
    logger.info("Output directory: %s", output_dir)
    logger.info("No target-context prediction was run.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
