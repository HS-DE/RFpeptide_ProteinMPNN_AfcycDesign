from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess
from typing import Any, Iterable, Mapping

from common import append_run_header, read_csv, resolve_path, rows_to_markdown, setup_logger, write_csv, write_markdown


JOB_FIELDS = [
    "rfpeptides_job_id",
    "site_label",
    "site_id",
    "site_quality_tier",
    "stage0_target_pdb",
    "stage0_target_pdb_sha256",
    "stage0_mapping_csv",
    "stage0_mapping_csv_sha256",
    "stage0_hotspots_txt",
    "stage0_hotspots_normalized",
    "stage0_hotspots_sha256",
    "rfpeptides_root",
    "runtime_git_commit",
    "runtime_inference_utils_sha256",
    "runtime_model_runners_sha256",
    "runtime_run_inference_sha256",
    "runtime_util_sha256",
    "working_directory",
    "output_prefix",
    "output_dir",
    "log_file",
    "complete_log_file",
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


def _normalize_hotspots(values: Iterable[str]) -> list[str]:
    parsed: list[tuple[str, int]] = []
    for value in values:
        text = str(value).strip()
        match = re.fullmatch(r"([A-Za-z])(\-?\d+)", text)
        if match is None:
            raise RuntimeError(f"Invalid hotspot residue label: {value}")
        parsed.append((match.group(1), int(match.group(2))))
    if len(set(parsed)) != len(parsed):
        raise RuntimeError(f"Duplicate hotspot residue labels are not allowed: {values}")
    return [f"{chain}{number}" for chain, number in sorted(parsed)]


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


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _sha256_text(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def _git_head(path: Path) -> str:
    completed = subprocess.run(
        ["git", "-C", str(path), "rev-parse", "HEAD"],
        capture_output=True,
        text=True,
        check=False,
    )
    if completed.returncode != 0:
        raise RuntimeError(f"Could not resolve RFpeptides runtime Git commit under {path}: {completed.stderr.strip()}")
    return completed.stdout.strip()


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
        f"python ./scripts/run_inference.py --config-name {config_name} \\",
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
    runtime_git_commit: str,
    runtime_hashes: Mapping[str, str],
    complete_log_file: Path,
) -> str:
    runtime_root = _to_wsl_path(rfpeptides_root)
    preflight_jobs = [
        {
            "job_id": str(row["rfpeptides_job_id"]),
            "target_pdb": _to_wsl_path(row["stage0_target_pdb"]),
            "mapping_csv": _to_wsl_path(row["stage0_mapping_csv"]),
            "target_pdb_sha256": str(row["stage0_target_pdb_sha256"]),
            "mapping_csv_sha256": str(row["stage0_mapping_csv_sha256"]),
            "hotspots": list(row["hotspots_list"]),
            "hotspots_normalized": str(row["stage0_hotspots_normalized"]),
            "hotspots_sha256": str(row["stage0_hotspots_sha256"]),
            "contig": str(row["contig_raw"]),
        }
        for row in job_rows
    ]
    preflight_jobs_json = json.dumps(preflight_jobs, sort_keys=True)
    lines = [
        "#!/bin/bash",
        "set -eo pipefail",
        "",
        f"cd {_quote_bash(_to_wsl_path(rfpeptides_root))}",
        f"source {_source_path_arg(conda_setup)}",
        f"conda activate {_quote_bash(conda_env)}",
        "set -u",
        "mkdir -p " + " ".join(_quote_bash(_to_wsl_path(row["output_dir"])) for row in job_rows),
        "mkdir -p " + " ".join(_quote_bash(_to_wsl_path(Path(str(row["log_file"])).parent)) for row in job_rows),
        f"exec > >(tee {_quote_bash(_to_wsl_path(complete_log_file))}) 2>&1",
        f"export RFPEPTIDES_RUNTIME_ROOT={_quote_bash(runtime_root)}",
        "export RFPEPTIDES_REQUIRE_PROVENANCE_CLOSURE='1'",
        "export PYTHONPATH=\"${RFPEPTIDES_RUNTIME_ROOT}${PYTHONPATH:+:${PYTHONPATH}}\"",
        f"export EXPECTED_RFPEPTIDES_COMMIT={_quote_bash(runtime_git_commit)}",
        f"export EXPECTED_INFERENCE_UTILS_SHA256={_quote_bash(runtime_hashes['inference_utils'])}",
        f"export EXPECTED_MODEL_RUNNERS_SHA256={_quote_bash(runtime_hashes['model_runners'])}",
        f"export EXPECTED_RUN_INFERENCE_SHA256={_quote_bash(runtime_hashes['run_inference'])}",
        f"export EXPECTED_UTIL_SHA256={_quote_bash(runtime_hashes['util'])}",
        "",
        "python - <<'PY'",
        "import csv",
        "import hashlib",
        "import json",
        "import os",
        "from pathlib import Path",
        "import subprocess",
        "from types import SimpleNamespace",
        "import rfdiffusion.util as rfd_util",
        "from rfdiffusion.contigs import ContigMap",
        "import rfdiffusion.inference.model_runners as model_runners",
        "import rfdiffusion.inference.utils as inference_utils",
        "from rfdiffusion.inference.utils import get_idx0_hotspots, process_target",
        f"job_specs = json.loads({preflight_jobs_json!r})",
        "runtime_root = Path(os.environ['RFPEPTIDES_RUNTIME_ROOT']).resolve()",
        "expected_paths = {",
        "    'inference_utils': runtime_root / 'rfdiffusion/inference/utils.py',",
        "    'model_runners': runtime_root / 'rfdiffusion/inference/model_runners.py',",
        "    'run_inference': runtime_root / 'scripts/run_inference.py',",
        "    'util': runtime_root / 'rfdiffusion/util.py',",
        "}",
        "loaded_paths = {",
        "    'inference_utils': Path(inference_utils.__file__).resolve(),",
        "    'model_runners': Path(model_runners.__file__).resolve(),",
        "    'util': Path(rfd_util.__file__).resolve(),",
        "}",
        "for key, loaded_path in loaded_paths.items():",
        "    if loaded_path != expected_paths[key]:",
        "        raise RuntimeError(f'RFpeptides runtime split-brain: {key} loaded from {loaded_path}, expected {expected_paths[key]}')",
        "expected_hashes = {",
        "    'inference_utils': os.environ['EXPECTED_INFERENCE_UTILS_SHA256'],",
        "    'model_runners': os.environ['EXPECTED_MODEL_RUNNERS_SHA256'],",
        "    'run_inference': os.environ['EXPECTED_RUN_INFERENCE_SHA256'],",
        "    'util': os.environ['EXPECTED_UTIL_SHA256'],",
        "}",
        "for key, path in expected_paths.items():",
        "    observed_hash = hashlib.sha256(path.read_bytes()).hexdigest()",
        "    if observed_hash != expected_hashes[key]:",
        "        raise RuntimeError(f'RFpeptides runtime file changed: {path}; observed={observed_hash}, expected={expected_hashes[key]}')",
        "observed_commit = subprocess.check_output(['git', '-C', str(runtime_root), 'rev-parse', 'HEAD'], text=True).strip()",
        "if observed_commit != os.environ['EXPECTED_RFPEPTIDES_COMMIT']:",
        "    raise RuntimeError(f'RFpeptides runtime commit changed: observed={observed_commit}, expected={os.environ[\"EXPECTED_RFPEPTIDES_COMMIT\"]}')",
        "for job in job_specs:",
        "    target_pdb = Path(job['target_pdb']).resolve()",
        "    mapping_csv = Path(job['mapping_csv']).resolve()",
        "    if hashlib.sha256(target_pdb.read_bytes()).hexdigest() != job['target_pdb_sha256']:",
        "        raise RuntimeError(f\"{job['job_id']}: Stage 0 target PDB hash mismatch\")",
        "    if hashlib.sha256(mapping_csv.read_bytes()).hexdigest() != job['mapping_csv_sha256']:",
        "        raise RuntimeError(f\"{job['job_id']}: Stage 0 mapping CSV hash mismatch\")",
        "    if hashlib.sha256(job['hotspots_normalized'].encode('utf-8')).hexdigest() != job['hotspots_sha256']:",
        "        raise RuntimeError(f\"{job['job_id']}: normalized hotspot hash mismatch\")",
        "    with mapping_csv.open('r', encoding='utf-8-sig', newline='') as handle:",
        "        mapping_rows = list(csv.DictReader(handle))",
        "    mapping_rows.sort(key=lambda row: int(row['rfpeptides_residue_number']))",
        "    if not mapping_rows:",
        "        raise RuntimeError(f\"{job['job_id']}: Stage 0 mapping CSV is empty\")",
        "    expected_refs = [",
        "        (str(row['rfpeptides_chain']), int(row['rfpeptides_residue_number']))",
        "        for row in mapping_rows",
        "    ]",
        "    if len(set(expected_refs)) != len(expected_refs):",
        "        raise RuntimeError(f\"{job['job_id']}: duplicate Stage 0 mapping residues\")",
        "    mapping_hotspots = sorted(",
        "        f\"{row['rfpeptides_chain']}{int(row['rfpeptides_residue_number'])}\"",
        "        for row in mapping_rows",
        "        if str(row.get('is_selected_hotspot', '')).strip().lower() == 'true'",
        "    )",
        "    if mapping_hotspots != sorted(job['hotspots']):",
        "        raise RuntimeError(",
        "            f\"{job['job_id']}: Stage 0 mapping hotspots differ from the job: \"",
        "            f\"mapping={mapping_hotspots} job={sorted(job['hotspots'])}\"",
        "        )",
        "    target_feats = process_target(str(target_pdb), parse_hetatom=True, center=False)",
        "    target_pdb_idx = [(str(chain), int(number)) for chain, number in target_feats['pdb_idx']]",
        "    if target_pdb_idx != expected_refs:",
        "        raise RuntimeError(",
        "            f\"{job['job_id']}: real Stage 0 target PDB residue order differs from mapping CSV\"",
        "        )",
        "    contig_map = ContigMap(target_feats, contigs=[job['contig']])",
        "    mappings = contig_map.get_mappings()",
        "    binderlen = len(contig_map.inpaint)",
        "    receptor_refs = [(str(chain), int(number)) for chain, number in mappings['receptor_con_ref_pdb_idx']]",
        "    receptor_hal = [int(value) for value in mappings['receptor_con_hal_idx0']]",
        "    complex_refs = [(str(chain), int(number)) for chain, number in mappings['complex_con_ref_pdb_idx']]",
        "    complex_hal = [int(value) for value in mappings['complex_con_hal_idx0']]",
        "    expected_receptor_hal = list(range(len(expected_refs)))",
        "    expected_complex_hal = list(range(binderlen, binderlen + len(expected_refs)))",
        "    if receptor_refs != expected_refs or complex_refs != expected_refs:",
        "        raise RuntimeError(f\"{job['job_id']}: real ContigMap target references differ from Stage 0 mapping\")",
        "    if receptor_hal != expected_receptor_hal:",
        "        raise RuntimeError(f\"{job['job_id']}: receptor-local ContigMap indices are not contiguous from zero\")",
        "    if complex_hal != expected_complex_hal:",
        "        raise RuntimeError(f\"{job['job_id']}: complex-global ContigMap indices do not include binderlen\")",
        "    helper_indices = sorted(",
        "        int(value)",
        "        for value in (get_idx0_hotspots(mappings, SimpleNamespace(hotspot_res=job['hotspots']), binderlen) or [])",
        "    )",
        "    hotspot_set = set(job['hotspots'])",
        "    expected_indices = sorted(",
        "        binderlen + index",
        "        for index, (chain, number) in enumerate(expected_refs)",
        "        if f'{chain}{number}' in hotspot_set",
        "    )",
        "    contig_indices = sorted(",
        "        complex_hal[index]",
        "        for index, (chain, number) in enumerate(complex_refs)",
        "        if f'{chain}{number}' in hotspot_set",
        "    )",
        "    if len(expected_indices) != len(hotspot_set):",
        "        raise RuntimeError(f\"{job['job_id']}: not all hotspots occur in the real Stage 0 mapping\")",
        "    if helper_indices != expected_indices or contig_indices != expected_indices:",
        "        raise RuntimeError(",
        "            f\"{job['job_id']}: real ContigMap hotspot provenance mismatch: \"",
        "            f\"helper={helper_indices} contig={contig_indices} expected={expected_indices}\"",
        "        )",
        "    print('[RFpeptides preflight] real_contigmap=' + json.dumps({",
        "        'job_id': job['job_id'],",
        "        'binderlen': binderlen,",
        "        'receptor_con_ref_pdb_idx': receptor_refs,",
        "        'receptor_con_hal_idx0': receptor_hal,",
        "        'complex_con_ref_pdb_idx': complex_refs,",
        "        'complex_con_hal_idx0': complex_hal,",
        "        'expected_hotspot_complex_global_indices': expected_indices,",
        "    }, separators=(',', ':')))",
        "print(f'[RFpeptides preflight] runtime={runtime_root}')",
        "print(f'[RFpeptides preflight] commit={observed_commit}')",
        "print('[RFpeptides preflight] module identity, hashes, Stage 0 inputs, and real ContigMap: PASS')",
        "PY",
        "",
    ]
    for row in job_rows:
        lines.extend(
            [
                f"echo '[RFpeptides] starting {row['rfpeptides_job_id']}'",
                f"export RFPEPTIDES_JOB_ID={_quote_bash(str(row['rfpeptides_job_id']))}",
                f"export STAGE0_TARGET_PDB_SHA256={_quote_bash(str(row['stage0_target_pdb_sha256']))}",
                f"export STAGE0_MAPPING_CSV={_quote_bash(_to_wsl_path(row['stage0_mapping_csv']))}",
                f"export STAGE0_MAPPING_CSV_SHA256={_quote_bash(str(row['stage0_mapping_csv_sha256']))}",
                f"export STAGE0_HOTSPOTS_NORMALIZED={_quote_bash(str(row['stage0_hotspots_normalized']))}",
                f"export STAGE0_HOTSPOTS_SHA256={_quote_bash(str(row['stage0_hotspots_sha256']))}",
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

- The runner pins one RFpeptides runtime root, Git commit, and four source-file
  hashes. It aborts if preflight and inference would import different code.
- The Stage 0 target PDB, residue-mapping CSV, and normalized hotspot list are
  SHA256-locked in the job table and rechecked immediately before inference.
- Preflight constructs a real ContigMap from the locked Stage 0 PDB, requested
  contig, and hotspots. Runtime then hard-checks helper, ContigMap-derived, and
  actual model hotspot-tensor indices before any output can be accepted.
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
    parser.add_argument("--rfpeptides-root", default="/home/luomi/fga_model_envs/rfpeptides/RFdiffusion")
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
    runtime_files = {
        "inference_utils": rfpeptides_root / "rfdiffusion" / "inference" / "utils.py",
        "model_runners": rfpeptides_root / "rfdiffusion" / "inference" / "model_runners.py",
        "run_inference": rfpeptides_root / "scripts" / "run_inference.py",
        "util": rfpeptides_root / "rfdiffusion" / "util.py",
    }
    missing_runtime_files = [str(path) for path in runtime_files.values() if not path.exists()]
    if missing_runtime_files:
        raise RuntimeError(f"Missing RFpeptides runtime files: {missing_runtime_files}")
    runtime_git_commit = _git_head(rfpeptides_root)
    runtime_hashes = {key: _sha256(path) for key, path in runtime_files.items()}

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

    complete_log_file = logs_dir / f"{Path(args.run_script_name).stem}.complete.log"
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
        mapping_csv = Path(str(row.get("crop_renumbering_mapping_csv", "")))
        if not mapping_csv.exists():
            raise RuntimeError(f"Missing Stage 0 mapping CSV for {site_label}: {mapping_csv}")
        hotspots = _normalize_hotspots(_split_csv(str(row.get("rfpeptides_hotspots", ""))))
        if not 3 <= len(hotspots) <= 6:
            raise RuntimeError(f"{site_label} has {len(hotspots)} hotspot(s); expected 3-6 for binder guidance.")
        chain_id, start, end = _parse_rf_range(str(row.get("rfpeptides_residue_range", "")))
        if any(not hotspot.startswith(chain_id) for hotspot in hotspots):
            raise RuntimeError(f"{site_label} hotspot chain IDs do not match target range chain {chain_id}: {hotspots}")
        mapping_hotspots = _normalize_hotspots(
            f"{mapping_row.get('rfpeptides_chain', '')}{mapping_row.get('rfpeptides_residue_number', '')}"
            for mapping_row in _read_required_csv(mapping_csv)
            if str(mapping_row.get("is_selected_hotspot", "")).strip().lower() == "true"
        )
        if mapping_hotspots != hotspots:
            raise RuntimeError(
                f"{site_label} hotspot list differs between Stage 0 summary and mapping CSV: "
                f"summary={hotspots}, mapping={mapping_hotspots}"
            )
        hotspots_normalized = ",".join(hotspots)
        target_pdb_sha256 = _sha256(target_pdb)
        mapping_csv_sha256 = _sha256(mapping_csv)
        hotspots_sha256 = _sha256_text(hotspots_normalized)

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
                "stage0_target_pdb_sha256": target_pdb_sha256,
                "stage0_mapping_csv": mapping_csv,
                "stage0_mapping_csv_sha256": mapping_csv_sha256,
                "stage0_hotspots_txt": row.get("hotspots_txt", ""),
                "stage0_hotspots_normalized": hotspots_normalized,
                "stage0_hotspots_sha256": hotspots_sha256,
                "rfpeptides_root": rfpeptides_root,
                "runtime_git_commit": runtime_git_commit,
                "runtime_inference_utils_sha256": runtime_hashes["inference_utils"],
                "runtime_model_runners_sha256": runtime_hashes["model_runners"],
                "runtime_run_inference_sha256": runtime_hashes["run_inference"],
                "runtime_util_sha256": runtime_hashes["util"],
                "working_directory": rfpeptides_root,
                "output_prefix": output_prefix,
                "output_dir": output_dir,
                "log_file": log_file,
                "complete_log_file": complete_log_file,
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
                "contig_raw": contig,
                "hotspots_list": hotspots,
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
        runtime_git_commit=runtime_git_commit,
        runtime_hashes=runtime_hashes,
        complete_log_file=complete_log_file,
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
