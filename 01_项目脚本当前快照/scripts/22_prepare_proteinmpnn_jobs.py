from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
from typing import Any, Iterable, Mapping

from common import append_run_header, read_csv, resolve_path, rows_to_markdown, setup_logger, write_csv, write_markdown


JOB_FIELDS = [
    "stage3_job_id",
    "run_group_id",
    "protocol_identity_sha256",
    "stage3_mode",
    "preferred_route",
    "design_id",
    "backbone_id",
    "site_label",
    "site_id",
    "source_backbone_pdb",
    "source_backbone_pdb_sha256",
    "stage2_pass_csv",
    "stage2_pass_csv_sha256",
    "input_tag",
    "input_tags",
    "input_pdb",
    "input_pdb_sha256",
    "input_pdbs",
    "input_pdb_sha256s",
    "input_pdb_dir",
    "runlist",
    "output_pdb_dir",
    "dl_binder_design_root",
    "dl_interface_design_script",
    "conda_env",
    "seqs_per_backbone",
    "relax_cycles",
    "proteinmpnn_temperature",
    "omit_aas",
    "command",
    "run_script",
    "log_file",
    "status",
    "notes",
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


def _to_wsl_path(path: str | Path) -> str:
    text = str(path).replace("\\", "/")
    if len(text) >= 3 and text[1] == ":" and text[2] == "/":
        return f"/mnt/{text[0].lower()}{text[2:]}"
    return text


def _resolve_mixed_path(value: str | Path) -> Path:
    text = str(value).strip().replace("\\", "/")
    if not text:
        return resolve_path(text)
    if os.name == "nt" and text.startswith("/mnt/") and len(text) > 7 and text[6] == "/":
        return Path(f"{text[5].upper()}:/{text[7:]}")
    if os.name != "nt" and len(text) >= 3 and text[1] == ":" and text[2] == "/":
        return Path(f"/mnt/{text[0].lower()}{text[2:]}")
    path = Path(text)
    if path.is_absolute():
        return path
    return resolve_path(path)


def _quote_bash(value: str) -> str:
    return "'" + str(value).replace("'", "'\"'\"'") + "'"


def _source_path_arg(value: str) -> str:
    text = str(value).strip()
    if text.startswith("~/"):
        return text
    return _quote_bash(_to_wsl_path(text))


def _read_required_csv(path: Path) -> list[dict[str, str]]:
    rows = read_csv(path)
    if not rows:
        raise RuntimeError(f"Missing or empty CSV: {path}")
    return rows


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _canonical_sha256(value: Mapping[str, Any]) -> str:
    payload = json.dumps(value, ensure_ascii=True, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _strict_design_lookup(rows: Iterable[Mapping[str, str]]) -> dict[str, dict[str, str]]:
    lookup: dict[str, dict[str, str]] = {}
    for row in rows:
        item = dict(row)
        design_id = str(item.get("design_id", "")).strip()
        if not design_id:
            raise RuntimeError("Stage 2 pass row is missing design_id")
        if "," in design_id:
            raise RuntimeError(f"Stage 2 pass row has aggregate design_id: {design_id}")
        if design_id in lookup:
            raise RuntimeError(f"Duplicate Stage 2 design_id: {design_id}")
        lookup[design_id] = item
    return lookup


def _chain_id_from_pdb_line(line: str, source_pdb: Path) -> str:
    if len(line) <= 21:
        raise RuntimeError(f"Coordinate record is too short to contain a chain ID in {source_pdb}: {line!r}")
    return line[21].strip() or "_"


def _pdb_record_name(line: str) -> str:
    return line[:6].strip().upper()


def _validate_normalized_pdb(
    output_pdb: Path,
    expected_by_chain: Mapping[str, list[str]],
    ordered_chains: list[str],
) -> None:
    lines = output_pdb.read_text(encoding="utf-8").splitlines()
    if not lines or lines[-1] != "END" or sum(line == "END" for line in lines) != 1:
        raise RuntimeError(f"Normalized PDB must end with exactly one END: {output_pdb}")
    first_coordinate = next(
        (idx for idx, line in enumerate(lines) if _pdb_record_name(line) in {"ATOM", "HETATM"}),
        None,
    )
    if first_coordinate is None:
        raise RuntimeError(f"Normalized PDB contains no ATOM/HETATM records: {output_pdb}")
    if any(_pdb_record_name(line) == "TER" for line in lines[:first_coordinate]):
        raise RuntimeError(f"Normalized PDB contains TER before the first coordinate: {output_pdb}")

    actual_by_chain: dict[str, list[str]] = {}
    ter_count = 0
    for line in lines:
        record = _pdb_record_name(line)
        if record in {"ATOM", "HETATM", "ANISOU", "SIGATM", "SIGUIJ"}:
            chain_id = _chain_id_from_pdb_line(line, output_pdb)
            actual_by_chain.setdefault(chain_id, []).append(line)
        elif record == "TER":
            ter_count += 1
    if actual_by_chain != dict(expected_by_chain):
        raise RuntimeError(f"Coordinate records changed during PDB normalization: {output_pdb}")
    if ter_count != len(ordered_chains):
        raise RuntimeError(
            f"Normalized PDB TER count {ter_count} does not match coordinate chain count {len(ordered_chains)}: {output_pdb}"
        )
    actual_order: list[str] = []
    for line in lines:
        if _pdb_record_name(line) not in {"ATOM", "HETATM"}:
            continue
        chain_id = _chain_id_from_pdb_line(line, output_pdb)
        if not actual_order or actual_order[-1] != chain_id:
            actual_order.append(chain_id)
    if actual_order != ordered_chains:
        raise RuntimeError(f"Normalized PDB chain order mismatch: expected={ordered_chains}, actual={actual_order}")


def _copy_or_reorder_pdb(source_pdb: Path, output_pdb: Path, peptide_chain: str, target_chain: str) -> str:
    if not peptide_chain or not target_chain:
        raise RuntimeError("Peptide and target chain IDs must both be provided")
    if peptide_chain == target_chain:
        raise RuntimeError(f"Peptide and target chain IDs must differ: {peptide_chain}")
    output_pdb.parent.mkdir(parents=True, exist_ok=True)
    lines_by_chain: dict[str, list[str]] = {}
    chain_order: list[str] = []
    header_lines: list[str] = []
    footer_lines: list[str] = []
    seen_coordinate = False
    coordinate_records = {"ATOM", "HETATM"}
    associated_records = {"ANISOU", "SIGATM", "SIGUIJ"}
    footer_records = {"CONECT", "MASTER"}
    with source_pdb.open("r", encoding="utf-8", errors="strict") as handle:
        for line in handle:
            text = line.rstrip("\r\n")
            record = _pdb_record_name(text)
            if record in {"MODEL", "ENDMDL"}:
                raise RuntimeError(f"Multi-model PDB input is not supported for Stage 3 preparation: {source_pdb}")
            if record in {"TER", "END"}:
                continue
            if record in coordinate_records | associated_records:
                chain_id = _chain_id_from_pdb_line(text, source_pdb)
                if record in coordinate_records and chain_id not in chain_order:
                    chain_order.append(chain_id)
                if record in associated_records and chain_id not in lines_by_chain:
                    raise RuntimeError(
                        f"Coordinate-associated record appears before its chain coordinates in {source_pdb}: {text!r}"
                    )
                lines_by_chain.setdefault(chain_id, []).append(text)
                seen_coordinate = True
            elif record in footer_records or seen_coordinate:
                footer_lines.append(text)
            else:
                header_lines.append(text)

    missing = [chain for chain in [peptide_chain, target_chain] if chain not in lines_by_chain]
    if missing:
        raise RuntimeError(f"Missing required PDB coordinate chain(s) {missing} in {source_pdb}; found {chain_order}")

    ordered_chains = [peptide_chain, target_chain]
    ordered_chains.extend(chain for chain in chain_order if chain not in set(ordered_chains))
    with output_pdb.open("w", encoding="utf-8", newline="\n") as out:
        for line in header_lines:
            out.write(line + "\n")
        for chain in ordered_chains:
            for line in lines_by_chain.get(chain, []):
                out.write(line + "\n")
            if lines_by_chain.get(chain):
                out.write("TER\n")
        for line in footer_lines:
            out.write(line + "\n")
        out.write("END\n")
    expected_by_chain = {chain: lines_by_chain[chain] for chain in ordered_chains}
    _validate_normalized_pdb(output_pdb, expected_by_chain, ordered_chains)
    return "normalized_peptide_chain_first"


def _write_text_lf(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(text.rstrip() + "\n")


def _command_lines(
    *,
    dl_interface_design_script: Path,
    input_pdb_dir: Path,
    output_pdb_dir: Path,
    runlist: Path,
    checkpoint: Path,
    relax_cycles: int,
    seqs_per_struct: int,
    temperature: float,
    omit_aas: str,
    extra_args: list[str],
) -> list[str]:
    lines = [
        f"python {_quote_bash(_to_wsl_path(dl_interface_design_script))} \\",
        f"-pdbdir {_quote_bash(_to_wsl_path(input_pdb_dir))} \\",
        f"-outpdbdir {_quote_bash(_to_wsl_path(output_pdb_dir))} \\",
        f"-runlist {_quote_bash(_to_wsl_path(runlist))} \\",
        f"-checkpoint_name {_quote_bash(_to_wsl_path(checkpoint))} \\",
        f"-relax_cycles {relax_cycles} \\",
        f"-seqs_per_struct {seqs_per_struct} \\",
        f"-temperature {temperature:g} \\",
        f"-omit_AAs {_quote_bash(omit_aas)}" + (" \\" if extra_args else ""),
    ]
    for idx, extra_arg in enumerate(extra_args):
        suffix = " \\" if idx < len(extra_args) - 1 else ""
        lines.append(f"{extra_arg}{suffix}")
    return lines


def _command_one_line(lines: list[str]) -> str:
    parts = []
    for line in lines:
        stripped = line.strip()
        if stripped.endswith("\\"):
            stripped = stripped[:-1].strip()
        parts.append(stripped)
    return " ".join(parts)


def _run_script_text(
    *,
    stage3_mode: str,
    conda_setup: str,
    conda_env: str,
    work_dir: Path,
    checkpoint: Path,
    output_pdb_dir: Path,
    log_file: Path,
    command_multiline: str,
) -> str:
    return "\n".join(
        [
            "#!/bin/bash",
            "set -eo pipefail",
            "",
            f"mkdir -p {_quote_bash(_to_wsl_path(work_dir))}",
            f"mkdir -p {_quote_bash(_to_wsl_path(output_pdb_dir))}",
            f"mkdir -p {_quote_bash(_to_wsl_path(log_file.parent))}",
            f"mkdir -p {_quote_bash(_to_wsl_path(checkpoint.parent))}",
            f"cd {_quote_bash(_to_wsl_path(work_dir))}",
            f"rm -f {_quote_bash(_to_wsl_path(checkpoint))}",
            f"source {_source_path_arg(conda_setup)}",
            f"conda activate {_quote_bash(conda_env)}",
            "set -u",
            "",
            f"echo '[{stage3_mode}] starting Stage 3'",
            "{",
            command_multiline,
            f"}} 2>&1 | tee {_quote_bash(_to_wsl_path(log_file))}",
            f"echo '[{stage3_mode}] finished Stage 3'",
            "",
        ]
    )


def _summary_markdown(job_rows: list[Mapping[str, Any]], run_script: Path) -> str:
    columns = [
        "stage3_job_id",
        "stage3_mode",
        "design_id",
        "seqs_per_backbone",
        "relax_cycles",
        "proteinmpnn_temperature",
        "omit_aas",
        "output_pdb_dir",
        "status",
    ]
    return f"""# FGA RFpeptides Stage 3 ProteinMPNN Jobs

Status: ProteinMPNN inputs and run script prepared. No sequence design was run
by this preparation script.

Run script:

```text
{run_script}
```

Important:

- The input complex keeps the RFpeptides peptide chain first and the FGA target
  crop second, matching `dl_interface_design.py` expectations.
- FGA target chain remains visible/fixed; ProteinMPNN designs the peptide chain.
- If `relax_cycles > 0`, `dl_interface_design.py` disallows
  `seqs_per_struct > 1`, so this preparation script duplicates the selected
  backbone into independent input tags and runs one FastRelax sequence per tag.
- If `relax_cycles == 0`, this preparation script keeps one input copy per
  backbone and asks ProteinMPNN to generate multiple sequences from that fixed
  RFpeptides backbone.
- Stage 3 outputs are sequence-designed intermediate structures, not final
  peptide candidates.

{rows_to_markdown(job_rows, columns, "No Stage 3 jobs were prepared.")}
"""


def main() -> int:
    parser = argparse.ArgumentParser(description="Prepare Stage 3 ProteinMPNN jobs for RFpeptides backbones.")
    parser.add_argument("--stage2-root", required=True)
    parser.add_argument("--output-root", default="", help="Defaults to --stage2-root.")
    parser.add_argument("--stage2-pass-csv", default="")
    parser.add_argument("--selected-backbones", required=True)
    parser.add_argument("--dl-binder-design-root", default="C:/SH/peptide_str/dl_binder_design")
    parser.add_argument("--seqs-per-backbone", type=int, default=8)
    parser.add_argument("--relax-cycles", type=int, default=0)
    parser.add_argument("--temperature", type=float, default=0.10)
    parser.add_argument("--omit-aas", default="CX")
    parser.add_argument("--conda-setup", default="~/fga_model_envs/miniforge3/etc/profile.d/conda.sh")
    parser.add_argument("--conda-env", default="proteinmpnn_binder_design")
    parser.add_argument("--run-script-name", default="", help="Defaults to a mode-specific Stage 3B script name.")
    parser.add_argument("--extra-arg", action="append", default=[], help="Additional dl_interface_design.py argument. Repeatable.")
    args = parser.parse_args()

    logger = setup_logger("22_prepare_proteinmpnn_jobs")
    append_run_header(logger, "22_prepare_proteinmpnn_jobs.py")

    if args.seqs_per_backbone <= 0:
        raise RuntimeError("--seqs-per-backbone must be > 0")
    if args.relax_cycles < 0:
        raise RuntimeError("--relax-cycles must be >= 0")
    if args.temperature < 0:
        raise RuntimeError("--temperature must be >= 0")
    if not args.omit_aas:
        raise RuntimeError("--omit-aas must not be empty")

    stage2_root = _resolve_mixed_path(args.stage2_root)
    output_root = _resolve_mixed_path(args.output_root) if args.output_root else stage2_root
    stage2_pass_csv = (
        _resolve_mixed_path(args.stage2_pass_csv)
        if args.stage2_pass_csv
        else stage2_root / "03_backbone_qc" / "FGA_rfpeptides_backbones_qc_pass.csv"
    )
    dl_binder_design_root = _resolve_mixed_path(args.dl_binder_design_root)
    dl_interface_design_script = dl_binder_design_root / "mpnn_fr" / "dl_interface_design.py"
    if not dl_interface_design_script.exists():
        raise RuntimeError(f"Missing dl_interface_design.py: {dl_interface_design_script}")
    default_checkpoint = dl_binder_design_root / "mpnn_fr" / "ProteinMPNN" / "vanilla_model_weights" / "v_48_020.pt"
    if not default_checkpoint.exists():
        raise RuntimeError(f"Missing ProteinMPNN checkpoint: {default_checkpoint}")

    pass_rows = _read_required_csv(stage2_pass_csv)
    pass_lookup = _strict_design_lookup(pass_rows)
    selected = _split_csv(args.selected_backbones)
    if not selected:
        raise RuntimeError("--selected-backbones must not be empty")
    if len(selected) != len(set(selected)):
        raise RuntimeError("--selected-backbones contains duplicate backbone IDs")

    selected_rows: list[dict[str, Any]] = []
    source_pdb_hashes: dict[str, str] = {}
    for backbone_id in selected:
        row = pass_lookup.get(backbone_id)
        if row is None:
            raise RuntimeError(f"Selected backbone not found in Stage 2 pass CSV: {backbone_id}")
        if str(row.get("pass_backbone_qc", "")).strip().lower() != "true":
            raise RuntimeError(f"Selected backbone is not marked pass_backbone_qc=true: {backbone_id}")
        source_pdb = _resolve_mixed_path(str(row.get("rf_pdb", "")))
        if not source_pdb.exists():
            raise RuntimeError(f"Missing source backbone PDB for {backbone_id}: {source_pdb}")
        peptide_chain = str(row.get("peptide_chain", "")).strip()
        target_chain = str(row.get("target_chain", "")).strip()
        if not peptide_chain or not target_chain:
            raise RuntimeError(f"Missing peptide/target chain metadata for {backbone_id}")
        if peptide_chain == target_chain:
            raise RuntimeError(f"Peptide and target chains are identical for {backbone_id}: {peptide_chain}")
        item = dict(row)
        item["_resolved_source_pdb"] = source_pdb
        item["_peptide_chain"] = peptide_chain
        item["_target_chain"] = target_chain
        selected_rows.append(item)
        source_pdb_hashes[backbone_id] = _sha256_file(source_pdb)

    stage3_mode = "proteinmpnn_only" if args.relax_cycles == 0 else "proteinmpnn_fastrelax"
    stage2_pass_csv_sha256 = _sha256_file(stage2_pass_csv)
    protocol_identity = {
        "selected_backbones": sorted(selected),
        "source_backbone_pdb_sha256": {key: source_pdb_hashes[key] for key in sorted(source_pdb_hashes)},
        "stage2_pass_csv_sha256": stage2_pass_csv_sha256,
        "stage3_mode": stage3_mode,
        "seqs_per_backbone": args.seqs_per_backbone,
        "relax_cycles": args.relax_cycles,
        "temperature": args.temperature,
        "omit_aas": args.omit_aas,
        "extra_args": list(args.extra_arg),
        "dl_interface_design_script_sha256": _sha256_file(dl_interface_design_script),
        "proteinmpnn_checkpoint_sha256": _sha256_file(default_checkpoint),
    }
    protocol_identity_sha256 = _canonical_sha256(protocol_identity)
    run_group_id = f"stage3_{len(selected)}bp_{protocol_identity_sha256[:12]}"
    tag_mode = "mpnnonly" if stage3_mode == "proteinmpnn_only" else "mpnnfr"
    input_subdir = "pdbs_proteinmpnn_only" if stage3_mode == "proteinmpnn_only" else "pdbs"
    output_subdir = "proteinmpnn_only_pdbs" if stage3_mode == "proteinmpnn_only" else "fastrelax_pdbs"
    input_dir = output_root / "04_proteinmpnn_inputs" / input_subdir / run_group_id
    jobs_dir = output_root / "04_proteinmpnn_inputs"
    work_dir = jobs_dir / "work" / run_group_id
    output_pdb_dir = output_root / "05_proteinmpnn_sequences" / output_subdir / run_group_id
    logs_dir = output_root / "logs"
    checkpoint = jobs_dir / "checkpoints" / f"{run_group_id}.checkpoint"
    runlist = jobs_dir / f"FGA_rfpeptides_{run_group_id}_runlist.txt"
    if args.run_script_name:
        requested_script = Path(args.run_script_name)
        suffix = requested_script.suffix or ".sh"
        run_script_name = f"{_safe_token(requested_script.stem)}_{run_group_id}{suffix}"
    else:
        run_script_name = f"run_{run_group_id}.sh"
    run_script = jobs_dir / run_script_name
    log_file = logs_dir / f"{run_group_id}.log"
    seqs_per_struct = 1 if args.relax_cycles > 0 else args.seqs_per_backbone
    copies_per_backbone = args.seqs_per_backbone if args.relax_cycles > 0 else 1

    runlist_tags: list[str] = []
    job_rows: list[dict[str, Any]] = []
    prepared_inputs: dict[str, list[tuple[str, Path, str]]] = {}
    for row in selected_rows:
        backbone_id = str(row["design_id"])
        source_pdb = Path(row["_resolved_source_pdb"])
        peptide_chain = str(row["_peptide_chain"])
        target_chain = str(row["_target_chain"])
        prepared_inputs[backbone_id] = []
        for copy_idx in range(1, copies_per_backbone + 1):
            tag = _safe_token(f"{backbone_id}_{tag_mode}_{copy_idx:03d}")
            input_pdb = input_dir / f"{tag}.pdb"
            copy_note = _copy_or_reorder_pdb(source_pdb, input_pdb, peptide_chain, target_chain)
            prepared_inputs[backbone_id].append((tag, input_pdb, copy_note))
            runlist_tags.append(tag)

    _write_text_lf(runlist, "\n".join(runlist_tags))
    command_lines = _command_lines(
        dl_interface_design_script=dl_interface_design_script,
        input_pdb_dir=input_dir,
        output_pdb_dir=output_pdb_dir,
        runlist=runlist,
        checkpoint=checkpoint,
        relax_cycles=args.relax_cycles,
        seqs_per_struct=seqs_per_struct,
        temperature=args.temperature,
        omit_aas=args.omit_aas,
        extra_args=args.extra_arg,
    )
    command_multiline = "\n".join(command_lines)

    for row in selected_rows:
        backbone_id = str(row["design_id"])
        input_items = prepared_inputs[backbone_id]
        tags = [item[0] for item in input_items]
        input_pdbs = [item[1] for item in input_items]
        copy_notes = [f"{item[0]}:{item[2]}" for item in input_items]
        job_id = _safe_token(f"{backbone_id}_{stage3_mode}_{protocol_identity_sha256[:12]}")
        job_rows.append(
            {
            "stage3_job_id": job_id,
            "run_group_id": run_group_id,
            "protocol_identity_sha256": protocol_identity_sha256,
            "stage3_mode": stage3_mode,
            "preferred_route": "true" if stage3_mode == "proteinmpnn_only" else "false",
            "design_id": backbone_id,
            "backbone_id": backbone_id,
            "site_label": row.get("site_label", ""),
            "site_id": row.get("site_id", ""),
            "source_backbone_pdb": row["_resolved_source_pdb"],
            "source_backbone_pdb_sha256": source_pdb_hashes[backbone_id],
            "stage2_pass_csv": stage2_pass_csv,
            "stage2_pass_csv_sha256": stage2_pass_csv_sha256,
            "input_tag": tags[0],
            "input_tags": ",".join(tags),
            "input_pdb": input_pdbs[0],
            "input_pdb_sha256": _sha256_file(input_pdbs[0]),
            "input_pdbs": ",".join(str(path) for path in input_pdbs),
            "input_pdb_sha256s": ",".join(_sha256_file(path) for path in input_pdbs),
            "input_pdb_dir": input_dir,
            "runlist": runlist,
            "output_pdb_dir": output_pdb_dir,
            "dl_binder_design_root": dl_binder_design_root,
            "dl_interface_design_script": dl_interface_design_script,
            "conda_env": args.conda_env,
            "seqs_per_backbone": args.seqs_per_backbone,
            "relax_cycles": args.relax_cycles,
            "proteinmpnn_temperature": args.temperature,
            "omit_aas": args.omit_aas,
            "command": _command_one_line(command_lines),
            "run_script": run_script,
            "log_file": log_file,
            "status": "pending_manual_execution",
            "notes": f"Stage 3 {stage3_mode} command only; sequence design not run by this script. "
            + "Input copy notes: "
            + ";".join(copy_notes),
            }
        )

    _write_text_lf(
        run_script,
        _run_script_text(
            stage3_mode=stage3_mode,
            conda_setup=args.conda_setup,
            conda_env=args.conda_env,
            work_dir=work_dir,
            checkpoint=checkpoint,
            output_pdb_dir=output_pdb_dir,
            log_file=log_file,
            command_multiline=command_multiline,
        ),
    )
    jobs_csv = jobs_dir / f"FGA_rfpeptides_{run_group_id}_jobs.csv"
    jobs_md = jobs_dir / f"FGA_rfpeptides_{run_group_id}_jobs.md"
    write_csv(jobs_csv, job_rows, JOB_FIELDS)
    write_markdown(jobs_md, _summary_markdown(job_rows, run_script))

    logger.info("Prepared Stage 3 %s jobs: %s", stage3_mode, len(job_rows))
    logger.info("Input PDB copies: %s", len(runlist_tags))
    logger.info("Run script: %s", run_script)
    logger.info("Jobs table: %s", jobs_csv)
    logger.info("No ProteinMPNN generation was run.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
