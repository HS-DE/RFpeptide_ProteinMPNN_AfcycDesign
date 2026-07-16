from __future__ import annotations

import argparse
import os
import shutil
from pathlib import Path
from typing import Any, Iterable, Mapping

from common import append_run_header, read_csv, resolve_path, rows_to_markdown, setup_logger, write_csv, write_markdown


JOB_FIELDS = [
    "stage3_job_id",
    "stage3_mode",
    "design_id",
    "site_label",
    "site_id",
    "source_backbone_pdb",
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


def _lookup_rows(rows: Iterable[Mapping[str, str]]) -> dict[str, dict[str, str]]:
    lookup: dict[str, dict[str, str]] = {}
    for row in rows:
        item = dict(row)
        for key in [item.get("design_id", ""), item.get("site_label", ""), item.get("site_id", "")]:
            if key:
                lookup[key] = item
    return lookup


def _first_atom_chain_order(pdb_path: Path) -> list[str]:
    seen: set[str] = set()
    order: list[str] = []
    with pdb_path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith("ATOM  "):
                continue
            chain_id = line[21].strip() or "_"
            if chain_id not in seen:
                seen.add(chain_id)
                order.append(chain_id)
    return order


def _copy_or_reorder_pdb(source_pdb: Path, output_pdb: Path, peptide_chain: str, target_chain: str) -> str:
    output_pdb.parent.mkdir(parents=True, exist_ok=True)
    chain_order = _first_atom_chain_order(source_pdb)
    if chain_order and chain_order[0] == peptide_chain:
        shutil.copy2(source_pdb, output_pdb)
        return "copied_original_chain_order"

    lines_by_chain: dict[str, list[str]] = {}
    other_lines: list[str] = []
    with source_pdb.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("ATOM  "):
                chain_id = line[21].strip() or "_"
                lines_by_chain.setdefault(chain_id, []).append(line.rstrip("\n"))
            elif line.startswith("END"):
                continue
            else:
                other_lines.append(line.rstrip("\n"))

    ordered_chains = [peptide_chain, target_chain]
    ordered_chains.extend(chain for chain in chain_order if chain not in set(ordered_chains))
    with output_pdb.open("w", encoding="utf-8", newline="\n") as out:
        for line in other_lines:
            out.write(line + "\n")
        for chain in ordered_chains:
            for line in lines_by_chain.get(chain, []):
                out.write(line + "\n")
            if lines_by_chain.get(chain):
                out.write("TER\n")
        out.write("END\n")
    return "reordered_peptide_chain_first"


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
    parser.add_argument("--stage2-root", default="results/rfpeptides_article_route_clean_20260615_fpocket_stage1_N1000_no_traj")
    parser.add_argument("--output-root", default="", help="Defaults to --stage2-root.")
    parser.add_argument("--stage2-pass-csv", default="")
    parser.add_argument("--selected-backbones", default="RFpep_Site_2_0007")
    parser.add_argument("--dl-binder-design-root", default="C:/SH/peptide_str/dl_binder_design")
    parser.add_argument("--seqs-per-backbone", type=int, default=8)
    parser.add_argument("--relax-cycles", type=int, default=1)
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
    pass_lookup = _lookup_rows(pass_rows)
    selected = _split_csv(args.selected_backbones)
    if not selected:
        raise RuntimeError("--selected-backbones must not be empty")

    stage3_mode = "proteinmpnn_only" if args.relax_cycles == 0 else "proteinmpnn_fastrelax"
    tag_mode = "mpnnonly" if stage3_mode == "proteinmpnn_only" else "mpnnfr"
    input_subdir = "pdbs_proteinmpnn_only" if stage3_mode == "proteinmpnn_only" else "pdbs"
    output_subdir = "proteinmpnn_only_pdbs" if stage3_mode == "proteinmpnn_only" else "fastrelax_pdbs"
    input_dir = output_root / "04_proteinmpnn_inputs" / input_subdir
    jobs_dir = output_root / "04_proteinmpnn_inputs"
    work_dir = jobs_dir / "work"
    output_pdb_dir = output_root / "05_proteinmpnn_sequences" / output_subdir
    logs_dir = output_root / "logs"
    checkpoint = jobs_dir / "checkpoints" / f"{stage3_mode}_stage3B.checkpoint"
    runlist = jobs_dir / f"FGA_rfpeptides_stage3B_{stage3_mode}_runlist.txt"
    run_script_name = args.run_script_name or f"run_stage3B_{stage3_mode}_site2_0007.sh"
    run_script = jobs_dir / run_script_name
    log_file = logs_dir / f"{stage3_mode}_stage3B_RFpep_Site_2_0007.log"
    seqs_per_struct = 1 if args.relax_cycles > 0 else args.seqs_per_backbone
    copies_per_backbone = args.seqs_per_backbone if args.relax_cycles > 0 else 1

    runlist_tags: list[str] = []
    job_rows: list[dict[str, Any]] = []
    copy_notes: list[str] = []
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

        for copy_idx in range(1, copies_per_backbone + 1):
            tag = _safe_token(f"{backbone_id}_{tag_mode}_{copy_idx:03d}")
            input_pdb = input_dir / f"{tag}.pdb"
            copy_note = _copy_or_reorder_pdb(source_pdb, input_pdb, peptide_chain, target_chain)
            copy_notes.append(f"{tag}:{copy_note}")
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

    job_id = _safe_token(f"RFpep_Site_2_0007_N{args.seqs_per_backbone}_{stage3_mode}_FR{args.relax_cycles}")
    job_rows.append(
        {
            "stage3_job_id": job_id,
            "stage3_mode": stage3_mode,
            "design_id": ",".join(selected),
            "site_label": pass_lookup[selected[0]].get("site_label", ""),
            "site_id": pass_lookup[selected[0]].get("site_id", ""),
            "source_backbone_pdb": pass_lookup[selected[0]].get("rf_pdb", ""),
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
    jobs_csv = jobs_dir / f"FGA_rfpeptides_stage3B_{stage3_mode}_jobs.csv"
    jobs_md = jobs_dir / f"FGA_rfpeptides_stage3B_{stage3_mode}_jobs.md"
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
