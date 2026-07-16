from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any, Iterable, Mapping

from common import append_run_header, read_csv, resolve_path, rows_to_markdown, setup_logger, write_csv, write_markdown


JOB_FIELDS = [
    "rfpeptides_job_id",
    "site_label",
    "site_id",
    "site_quality_tier",
    "stage0_target_pdb",
    "stage0_mapping_csv",
    "stage0_hotspots_txt",
    "rfpeptides_root",
    "working_directory",
    "output_prefix",
    "output_dir",
    "log_file",
    "config_name",
    "num_designs",
    "length_min",
    "length_max",
    "contigmap_contigs",
    "hotspot_res",
    "cyclic",
    "cyc_chains",
    "diffuser_T",
    "use_potentials",
    "command",
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
    return "".join(keep).strip("_") or "job"


def _to_wsl_path(path: str | Path) -> str:
    text = str(path).replace("\\", "/")
    if len(text) >= 3 and text[1] == ":" and text[2] == "/":
        drive = text[0].lower()
        return f"/mnt/{drive}{text[2:]}"
    return text


def _quote_bash(value: str) -> str:
    return "'" + value.replace("'", "'\"'\"'") + "'"


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


def _stage0_lookup(rows: Iterable[Mapping[str, str]]) -> dict[str, dict[str, str]]:
    lookup: dict[str, dict[str, str]] = {}
    for row in rows:
        item = dict(row)
        for key in [item.get("site_label", ""), item.get("site_id", "")]:
            if key:
                lookup[key] = item
    return lookup


def _parse_rf_range(value: str) -> tuple[str, int, int]:
    text = str(value).strip()
    if "-" not in text:
        raise RuntimeError(f"Invalid rfpeptides_residue_range: {value}")
    left, right = text.split("-", 1)
    left_chain = left[0]
    right_chain = right[0]
    if left_chain != right_chain:
        raise RuntimeError(f"RFpeptides residue range crosses chains: {value}")
    try:
        start = int(left[1:])
        end = int(right[1:])
    except ValueError as exc:
        raise RuntimeError(f"Invalid rfpeptides_residue_range: {value}") from exc
    if start < 1 or end < start:
        raise RuntimeError(f"Invalid rfpeptides_residue_range: {value}")
    return left_chain, start, end


def _hotspot_list_arg(hotspots: list[str]) -> str:
    return "[" + ",".join(f"'{hotspot}'" for hotspot in hotspots) + "]"


def _hotspot_shell_arg(hotspots: list[str]) -> str:
    # README binder guidance accepts the whole Hydra override as one quoted
    # argument, e.g. 'ppi.hotspot_res=[A30,A33,A34]'.
    return "'ppi.hotspot_res=[" + ",".join(hotspots) + "]'"


def _write_text_lf(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(text.rstrip() + "\n")


def _build_command_lines(
    *,
    rfpeptides_root: Path,
    config_name: str,
    output_prefix: Path,
    num_designs: int,
    contig: str,
    input_pdb: Path,
    diffuser_t: int,
    cyc_chains: str,
    hotspots: list[str],
    extra_overrides: list[str],
) -> list[str]:
    lines = [
        f"./scripts/run_inference.py --config-name {config_name} \\",
        f"inference.output_prefix={_to_wsl_path(output_prefix)} \\",
        f"inference.num_designs={num_designs} \\",
        f"'contigmap.contigs=[{contig}]' \\",
        f"inference.input_pdb={_to_wsl_path(input_pdb)} \\",
        "inference.cyclic=True \\",
        f"diffuser.T={diffuser_t} \\",
        f"inference.cyc_chains='{cyc_chains}' \\",
        _hotspot_shell_arg(hotspots) + (" \\" if extra_overrides else ""),
    ]
    for idx, override in enumerate(extra_overrides):
        suffix = " \\" if idx < len(extra_overrides) - 1 else ""
        lines.append(f"{override}{suffix}")
    return lines


def _command_one_line(command_lines: list[str]) -> str:
    parts: list[str] = []
    for line in command_lines:
        stripped = line.strip()
        if stripped.endswith("\\"):
            stripped = stripped[:-1].strip()
        parts.append(stripped)
    return " ".join(parts)


def _run_script_text(
    *,
    job_rows: list[Mapping[str, Any]],
    conda_setup: str,
    conda_env: str,
    rfpeptides_root: Path,
) -> str:
    lines = [
        "#!/bin/bash",
        "set -eo pipefail",
        "",
        f"cd {_quote_bash(_to_wsl_path(rfpeptides_root))}",
        f"source {_source_path_arg(conda_setup)}",
        f"conda activate {_quote_bash(conda_env)}",
        "set -u",
        "",
        "python - <<'PY'",
        "from types import SimpleNamespace",
        "from rfdiffusion.inference.utils import get_idx0_hotspots",
        "mappings = {",
        "    'receptor_con_ref_pdb_idx': [('A', 82), ('A', 84)],",
        "    'receptor_con_hal_idx0': [81, 83],",
        "}",
        "observed = get_idx0_hotspots(mappings, SimpleNamespace(hotspot_res=['A82', 'A84']), binderlen=17)",
        "expected = [98, 100]",
        "if list(observed or []) != expected:",
        "    raise RuntimeError(f'RFpeptides hotspot indexing preflight failed: observed={observed}, expected={expected}')",
        "print('[RFpeptides preflight] hotspot indices include the binder-length offset: PASS')",
        "PY",
        "",
        "mkdir -p " + " ".join(_quote_bash(_to_wsl_path(row["output_dir"])) for row in job_rows),
        "mkdir -p " + " ".join(_quote_bash(_to_wsl_path(Path(str(row["log_file"])).parent)) for row in job_rows),
        "",
    ]
    for row in job_rows:
        lines.extend(
            [
                f"echo '[RFpeptides] starting {row['rfpeptides_job_id']}'",
                "{",
                str(row["command_multiline"]),
                f"}} 2>&1 | tee {_quote_bash(_to_wsl_path(row['log_file']))}",
                f"echo '[RFpeptides] finished {row['rfpeptides_job_id']}'",
                "",
            ]
        )
    return "\n".join(lines).rstrip() + "\n"


def _summary_markdown(job_rows: list[Mapping[str, Any]], run_script: Path) -> str:
    columns = [
        "rfpeptides_job_id",
        "site_label",
        "num_designs",
        "contigmap_contigs",
        "hotspot_res",
        "output_prefix",
        "status",
    ]
    return f"""# FGA RFpeptides Stage 1 Jobs

Status: RFpeptides command table and run script prepared. No backbone
generation was run by this preparation script.

Run script:

```text
{run_script}
```

{rows_to_markdown(job_rows, columns, "No RFpeptides jobs were prepared.")}

Important:

- This is a small pilot in design count, not a reduced-constraint run.
- The command keeps the RFpeptides binder requirements: target PDB, target
  contig, cyclic generation, cyclic chain, diffuser timesteps, and hotspot
  residues.
- RFpep_Site_3 and RFpep_Site_4 remain deferred and are not included in this
  Stage 1 pilot job table.
"""


def main() -> int:
    parser = argparse.ArgumentParser(description="Prepare RFpeptides Stage 1 job commands for FGA article-route pilot.")
    parser.add_argument("--input-root", default="results/rfpeptides_article_route_clean_20260612")
    parser.add_argument("--output-root", default="results/rfpeptides_article_route_clean_20260612")
    parser.add_argument("--selected-sites", default="RFpep_Site_2")
    parser.add_argument("--stage0-summary-csv", default="")
    parser.add_argument("--rfpeptides-root", default="C:/SH/peptide_str/rfd_macro")
    parser.add_argument("--num-designs", type=int, default=10)
    parser.add_argument("--length-min", type=int, default=12)
    parser.add_argument("--length-max", type=int, default=18)
    parser.add_argument("--diffuser-t", type=int, default=50)
    parser.add_argument("--cyc-chains", default="a")
    parser.add_argument("--config-name", default="base")
    parser.add_argument("--conda-setup", default="~/fga_model_envs/miniforge3/etc/profile.d/conda.sh")
    parser.add_argument("--conda-env", default="SE3nv")
    parser.add_argument("--run-script-name", default="run_rfpeptides_stage1_site2.sh")
    parser.add_argument("--extra-override", action="append", default=[], help="Additional Hydra override. Repeatable.")
    parser.add_argument("--allow-potentials", action="store_true", help="Record that potentials are intentionally enabled via extra overrides.")
    args = parser.parse_args()

    logger = setup_logger("20_make_rfpeptides_article_jobs")
    append_run_header(logger, "20_make_rfpeptides_article_jobs.py")

    if args.num_designs <= 0:
        raise RuntimeError("--num-designs must be > 0")
    if args.length_min <= 0 or args.length_max < args.length_min:
        raise RuntimeError("--length-min/--length-max define an invalid length range")
    if args.diffuser_t <= 0:
        raise RuntimeError("--diffuser-t must be > 0")
    if not args.cyc_chains:
        raise RuntimeError("--cyc-chains must not be empty")

    input_root = resolve_path(args.input_root)
    output_root = resolve_path(args.output_root)
    jobs_dir = output_root / "01_rfpeptides_jobs"
    backbones_root = output_root / "02_rfpeptides_backbones"
    logs_dir = output_root / "logs"
    rfpeptides_root = resolve_path(args.rfpeptides_root)
    if not (rfpeptides_root / "scripts" / "run_inference.py").exists():
        raise RuntimeError(f"Missing RFpeptides run_inference.py under --rfpeptides-root: {rfpeptides_root}")

    stage0_summary_csv = (
        resolve_path(args.stage0_summary_csv)
        if args.stage0_summary_csv
        else input_root / "00_target_inputs" / "FGA_rfpeptides_stage0_target_inputs_summary.csv"
    )
    stage0_rows = _read_required_csv(stage0_summary_csv)
    stage0 = _stage0_lookup(stage0_rows)

    selected_sites = _split_csv(args.selected_sites)
    if not selected_sites:
        raise RuntimeError("--selected-sites must contain at least one Stage 0 site label or candidate ID")

    if not args.allow_potentials:
        potential_overrides = [override for override in args.extra_override if override.startswith("potentials.")]
        if potential_overrides:
            raise RuntimeError(
                "Potential overrides were provided but --allow-potentials was not set. "
                "The first FGA RFpeptides pilot should keep potentials off unless explicitly documented."
            )

    job_rows: list[dict[str, Any]] = []
    for selected_id in selected_sites:
        row = stage0.get(selected_id)
        if row is None:
            raise RuntimeError(f"Selected Stage 0 site not found in {stage0_summary_csv}: {selected_id}")
        site_label = row.get("site_label") or selected_id
        if site_label != "RFpep_Site_2":
            raise RuntimeError(
                f"{site_label} is not in the current Stage 1 pilot scope. "
                "Current decision is to generate RFpep_Site_2 only."
            )
        target_pdb = Path(str(row.get("target_pdb", "")))
        if not target_pdb.exists():
            raise RuntimeError(f"Missing Stage 0 target PDB for {site_label}: {target_pdb}")
        hotspots = _split_csv(str(row.get("rfpeptides_hotspots", "")))
        if not 3 <= len(hotspots) <= 6:
            raise RuntimeError(f"{site_label} has {len(hotspots)} hotspot(s); expected 3-6 for binder guidance.")
        chain_id, start, end = _parse_rf_range(str(row.get("rfpeptides_residue_range", "")))
        if any(not hotspot.startswith(chain_id) for hotspot in hotspots):
            raise RuntimeError(f"{site_label} hotspot chain IDs do not match target range chain {chain_id}: {hotspots}")

        length_tag = f"L{args.length_min}_{args.length_max}"
        job_id = _safe_token(f"{site_label}_{length_tag}_N{args.num_designs}")
        output_prefix = backbones_root / site_label / f"{site_label}_{length_tag}"
        output_dir = output_prefix.parent
        log_file = logs_dir / f"rfpeptides_stage1_{job_id}.log"
        contig = f"{args.length_min}-{args.length_max} {chain_id}{start}-{end}/0"
        command_lines = _build_command_lines(
            rfpeptides_root=rfpeptides_root,
            config_name=args.config_name,
            output_prefix=output_prefix,
            num_designs=args.num_designs,
            contig=contig,
            input_pdb=target_pdb,
            diffuser_t=args.diffuser_t,
            cyc_chains=args.cyc_chains,
            hotspots=hotspots,
            extra_overrides=args.extra_override,
        )
        command_multiline = "\n".join(command_lines)
        job_rows.append(
            {
                "rfpeptides_job_id": job_id,
                "site_label": site_label,
                "site_id": row.get("site_id", ""),
                "site_quality_tier": row.get("site_quality_tier", ""),
                "stage0_target_pdb": target_pdb,
                "stage0_mapping_csv": row.get("crop_renumbering_mapping_csv", ""),
                "stage0_hotspots_txt": row.get("hotspots_txt", ""),
                "rfpeptides_root": rfpeptides_root,
                "working_directory": rfpeptides_root,
                "output_prefix": output_prefix,
                "output_dir": output_dir,
                "log_file": log_file,
                "config_name": args.config_name,
                "num_designs": args.num_designs,
                "length_min": args.length_min,
                "length_max": args.length_max,
                "contigmap_contigs": f"[{contig}]",
                "hotspot_res": _hotspot_list_arg(hotspots),
                "cyclic": "true",
                "cyc_chains": args.cyc_chains,
                "diffuser_T": args.diffuser_t,
                "use_potentials": "true" if args.allow_potentials else "false",
                "command": _command_one_line(command_lines),
                "command_multiline": command_multiline,
                "status": "pending_manual_execution",
                "notes": (
                    "Stage 1 command only; generation is not run by this preparation script. "
                    "The generated runner aborts unless hotspot indices include the binder-length offset."
                ),
            }
        )

    run_script = jobs_dir / args.run_script_name
    run_script_text = _run_script_text(
        job_rows=job_rows,
        conda_setup=args.conda_setup,
        conda_env=args.conda_env,
        rfpeptides_root=rfpeptides_root,
    )
    _write_text_lf(run_script, run_script_text)
    write_csv(jobs_dir / "FGA_rfpeptides_stage1_jobs.csv", job_rows, JOB_FIELDS)
    write_markdown(jobs_dir / "FGA_rfpeptides_stage1_jobs.md", _summary_markdown(job_rows, run_script))

    logger.info("Prepared RFpeptides Stage 1 jobs: %s", len(job_rows))
    logger.info("Run script: %s", run_script)
    logger.info("No RFpeptides generation was run.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
