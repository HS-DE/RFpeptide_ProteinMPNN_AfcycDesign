from __future__ import annotations

import argparse
import hashlib
import json
import os
import shlex
from pathlib import Path
from typing import Any, Iterable, Mapping

from common import assert_active_route_path, append_run_header, read_csv, resolve_path, rows_to_markdown, setup_logger, write_csv, write_markdown
from pdb_utils import parse_residues, residue_sequence


COLABDESIGN_GAMMA_COMMIT = "5ab4efaba2321a6c3c314b82d2fff8e0241f5c2d"
PROTOCOL_VERSION = "stage5_target_only_controls_v1"
LEGAL_AA = set("ACDEFGHIKLMNPQRSTVWY")
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

CONTROL_FIELDS = [
    "control_group",
    "control_label",
    "target_msa_mode",
    "target_a3m",
    "target_a3m_sha1",
    "msa_rows",
    "use_mlm",
    "mlm_replace_fraction",
    "seeds_planned",
    "models_per_seed",
    "jobs_planned",
    "jobs_runnable",
    "model_evaluations_planned",
    "template_mode",
    "use_initial_guess",
    "peptide_included",
    "validation_test_type",
    "protocol_hash",
    "status",
    "notes",
]

JOB_FIELDS = [
    "stage5_control_job_id",
    "control_group",
    "protocol_hash",
    "target_msa_mode",
    "target_a3m",
    "target_a3m_sha1",
    "msa_rows_expected",
    "use_mlm",
    "mlm_replace_fraction",
    "seed",
    "requested_recycles",
    "forward_passes",
    "models_per_seed",
    "job_spec_json",
    "prediction_output_dir",
    "run_script",
    "template_mode",
    "use_initial_guess",
    "peptide_included",
    "status",
]


def _sha1_text(value: str) -> str:
    return hashlib.sha1(value.encode("utf-8")).hexdigest()


def _sha1_file(path: Path) -> str:
    digest = hashlib.sha1()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _protocol_hash(payload: Mapping[str, Any]) -> str:
    canonical = json.dumps(dict(payload), sort_keys=True, separators=(",", ":"))
    return _sha1_text(canonical)[:8]


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


def _write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text.rstrip() + "\n", encoding="utf-8", newline="\n")


def _validate_sequence(value: str, label: str, expected_length: int | None = None) -> str:
    sequence = str(value).strip().upper()
    if not sequence:
        raise RuntimeError(f"{label} is empty")
    invalid = sorted(set(sequence) - LEGAL_AA)
    if invalid:
        raise RuntimeError(f"{label} contains unsupported amino acids: {','.join(invalid)}")
    if expected_length is not None and len(sequence) != expected_length:
        raise RuntimeError(f"{label} length is {len(sequence)}, expected {expected_length}")
    return sequence


def _truthy(value: Any) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes"}


def _read_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header = ""
    chunks: list[str] = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header:
                    records.append((header, "".join(chunks)))
                header = line[1:].strip() or f"row_{len(records) + 1}"
                chunks = []
            elif not line.startswith("#"):
                if not header:
                    raise RuntimeError(f"Sequence encountered before FASTA/A3M header: {path}")
                chunks.append(line)
    if header:
        records.append((header, "".join(chunks)))
    if not records:
        raise RuntimeError(f"Missing FASTA/A3M records: {path}")
    return records


def _write_a3m(path: Path, records: Iterable[tuple[str, str]]) -> None:
    lines: list[str] = []
    for header, sequence in records:
        lines.extend([f">{header}", sequence])
    _write_text(path, "\n".join(lines))


def _a3m_match_tokens(sequence: str, expected_match_states: int) -> list[str]:
    tokens: list[str] = []
    for char in sequence.replace(".", "-"):
        if char.islower():
            if not tokens:
                raise RuntimeError("A3M row has an insertion before its first match-state column")
            tokens[-1] += char
        else:
            tokens.append(char)
    if len(tokens) != expected_match_states:
        raise RuntimeError(
            f"A3M row has {len(tokens)} match-state columns, expected {expected_match_states}"
        )
    return tokens


def _crop_full_target_a3m(
    *,
    source_a3m: Path,
    full_sequence: str,
    crop_start: int,
    crop_end: int,
    expected_crop_sequence: str,
    max_rows: int,
) -> list[tuple[str, str]]:
    records = _read_fasta(source_a3m)
    query_tokens = _a3m_match_tokens(records[0][1], len(full_sequence))
    query_match = "".join(token[0] for token in query_tokens).replace("-", "").upper()
    if query_match != full_sequence or any(token[0] == "-" for token in query_tokens):
        raise RuntimeError(
            "Real homolog A3M query must match data/input/FGA_full_length_1_866.fasta exactly; "
            "domain-only or differently numbered A3M inputs require an explicit mapping workflow"
        )
    cropped: list[tuple[str, str]] = []
    seen_sequences: set[str] = set()
    for header, sequence in records:
        tokens = _a3m_match_tokens(sequence, len(full_sequence))
        selected: list[str] = []
        for position in range(crop_start - 1, crop_end):
            token = tokens[position]
            selected.append(token if position < crop_end - 1 else token[0])
        cropped_sequence = "".join(selected)
        match_chars = "".join(ch for ch in cropped_sequence if not ch.islower())
        if set(match_chars) == {"-"}:
            continue
        if cropped_sequence in seen_sequences:
            continue
        seen_sequences.add(cropped_sequence)
        cropped.append((header, cropped_sequence))
        if len(cropped) >= max_rows:
            break
    if not cropped:
        raise RuntimeError("No target MSA rows remained after crop")
    query_crop = "".join(ch for ch in cropped[0][1] if not ch.islower() and ch != "-").upper()
    if query_crop != expected_crop_sequence:
        raise RuntimeError("Cropped homolog A3M query does not match the Stage 0 target sequence")
    if len(cropped) < 2 or not any(
        "".join(ch for ch in sequence if not ch.islower()).upper() != expected_crop_sequence
        for _, sequence in cropped[1:]
    ):
        raise RuntimeError("Homolog MSA control requires at least one real non-query homolog row")
    return cropped


def _stage0_context(
    stage0_root: Path,
    expected_target_length: int,
) -> tuple[Path, str, list[int], list[int], int, int]:
    target_pdb = stage0_root / "00_target_inputs" / "RFpep_Site_2_target.pdb"
    mapping_csv = stage0_root / "00_target_inputs" / "RFpep_Site_2_crop_renumbering_mapping.csv"
    assert_active_route_path(target_pdb, "Stage 28 Stage 0 target PDB")
    assert_active_route_path(mapping_csv, "Stage 28 Stage 0 mapping CSV")
    if not target_pdb.is_file():
        raise RuntimeError(f"Missing Stage 0 target PDB: {target_pdb}")
    chains = parse_residues(target_pdb)
    if "A" not in chains or set(chains) != {"A"}:
        raise RuntimeError("Stage 0 target PDB must contain only target chain A")
    target_sequence = _validate_sequence(
        residue_sequence(chains["A"]),
        "Stage 0 target chain A",
        expected_target_length,
    )
    rows = read_csv(mapping_csv)
    if len(rows) != expected_target_length:
        raise RuntimeError(f"Stage 0 mapping has {len(rows)} rows, expected {expected_target_length}")
    rows.sort(key=lambda row: int(row["rfpeptides_residue_number"]))
    expected_numbers = list(range(1, expected_target_length + 1))
    observed_numbers = [int(row["rfpeptides_residue_number"]) for row in rows]
    if observed_numbers != expected_numbers:
        raise RuntimeError("Stage 0 crop numbering must be contiguous from 1")
    mapping_sequence = "".join(
        AA3_TO_AA1.get(str(row.get("rfpeptides_residue_name", "")).upper(), "X") for row in rows
    )
    mapping_sequence = _validate_sequence(mapping_sequence, "Stage 0 mapping sequence", expected_target_length)
    if mapping_sequence != target_sequence:
        raise RuntimeError("Stage 0 mapping sequence does not match target PDB chain A")
    uniprot_numbers = [int(row["uniprot_residue_number"]) for row in rows]
    if uniprot_numbers != list(range(min(uniprot_numbers), max(uniprot_numbers) + 1)):
        raise RuntimeError("Stage 0 crop does not map to a contiguous UniProt interval")
    site_indices = [int(row["rfpeptides_residue_number"]) for row in rows if _truthy(row.get("is_target_site_residue"))]
    hotspot_indices = [int(row["rfpeptides_residue_number"]) for row in rows if _truthy(row.get("is_selected_hotspot"))]
    if not site_indices or not hotspot_indices:
        raise RuntimeError("Stage 0 mapping lacks Site_2 or hotspot indices")
    return target_pdb, target_sequence, site_indices, hotspot_indices, min(uniprot_numbers), max(uniprot_numbers)


def _job_script(
    *,
    project_root: Path,
    runner: Path,
    job_spec: Path,
    python_bin: str,
    source_dir: str,
    overlay_dir: str,
    af_params: str,
    runnable: bool,
) -> str:
    if not runnable:
        return """#!/bin/bash
set -euo pipefail
echo "BLOCKED: a verified real full-target FGA homolog A3M is required. Regenerate jobs with --full-target-a3m." >&2
exit 4
"""
    return f"""#!/bin/bash
set -euo pipefail

if [[ "${{RUN_STAGE5_TARGET_CONTROLS:-NO}}" != "YES" ]]; then
  echo "Target-only controls are review-gated. Set RUN_STAGE5_TARGET_CONTROLS=YES after inspection." >&2
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
"""


def _preflight_script(python_bin: str, source_dir: str, overlay_dir: str, af_params: str) -> str:
    return f"""#!/bin/bash
set -euo pipefail
PYTHON_BIN="${{AFCYCDESIGN_PYTHON:-{python_bin}}}"
SOURCE_DIR="${{COLABDESIGN_GAMMA_SOURCE:-{source_dir}}}"
OVERLAY_DIR="${{AFCYCDESIGN_PYTHON_OVERLAY:-{overlay_dir}}}"
AF_PARAMS_DIR="${{AF_PARAMS:-{af_params}}}"
COMMIT_MARKER="$SOURCE_DIR/.stage5_colabdesign_commit"

[[ -x "$PYTHON_BIN" ]] || {{ echo "FAIL: AfCycDesign Python missing" >&2; exit 2; }}
[[ -d "$AF_PARAMS_DIR" ]] || {{ echo "FAIL: AlphaFold parameters missing" >&2; exit 2; }}
[[ -f "$COMMIT_MARKER" ]] || {{ echo "FAIL: commit marker missing" >&2; exit 2; }}
[[ "$(tr -d '[:space:]' < "$COMMIT_MARKER")" == "{COLABDESIGN_GAMMA_COMMIT}" ]] || {{ echo "FAIL: commit mismatch" >&2; exit 2; }}
export PYTHONPATH="$SOURCE_DIR:$OVERLAY_DIR${{PYTHONPATH:+:$PYTHONPATH}}"
"$PYTHON_BIN" - <<'PY'
import pathlib
import colabdesign
from colabdesign import mk_af_model
from colabdesign.af.contrib import predict
print("PASS: target-only control modules importable from", pathlib.Path(colabdesign.__file__).resolve())
print("Protocol: target only; template=none; initial_guess=false; no peptide; no cyclic offset.")
PY
"""


def _master_script(preflight: Path, runnable_scripts: list[Path]) -> str:
    lines = [
        "#!/bin/bash",
        "set -euo pipefail",
        "",
        'if [[ "${RUN_STAGE5_TARGET_CONTROLS:-NO}" != "YES" ]]; then',
        '  echo "No target-control prediction started. Set RUN_STAGE5_TARGET_CONTROLS=YES after review." >&2',
        "  exit 3",
        "fi",
        "",
        f"bash {shlex.quote(_to_wsl_path(preflight))}",
        "",
    ]
    lines.extend(f"bash {shlex.quote(_to_wsl_path(script))}" for script in runnable_scripts)
    return "\n".join(lines) + "\n"


def _report(
    *,
    controls: list[Mapping[str, Any]],
    jobs: list[Mapping[str, Any]],
    output_dir: Path,
    homolog_source: str,
    crop_start: int,
    crop_end: int,
) -> str:
    runnable = [row for row in jobs if row["status"] == "prepared_not_run"]
    blocked = [row for row in jobs if str(row["status"]).startswith("blocked")]
    columns = [
        "control_group",
        "target_msa_mode",
        "msa_rows",
        "use_mlm",
        "jobs_planned",
        "jobs_runnable",
        "model_evaluations_planned",
        "status",
        "notes",
    ]
    return f"""# Stage 5 Target-Only Recovery Controls

These controls isolate target recovery before any Stage 5B decision. They are
sequence-based target-only predictions: no peptide, no target template, and no
initial guess. No cyclic offset is applied because no peptide is present.

## Four Protocol Groups

{rows_to_markdown(controls, columns, "No controls were prepared.")}

## Counts

```text
planned_groups: 4
planned_seed_jobs: {len(jobs)}
runnable_seed_jobs: {len(runnable)}
blocked_seed_jobs: {len(blocked)}
planned_model_evaluations: {sum(int(row['models_per_seed']) for row in jobs)}
output_directory: {output_dir}
```

Each seed evaluates all five AlphaFold parameter sets. With the default three
seeds, the complete four-group comparison is 12 seed jobs and 60 model
evaluations. `requested_recycles=6` means one initial forward plus six recycle
for seven forward passes.

## Real Homolog MSA Rule

```text
full_target_a3m_source: {homolog_source or 'not_provided'}
crop_uniprot_columns: {crop_start}-{crop_end}
```

Groups C and D become runnable only after a real, unpaired homolog A3M whose
query exactly matches `data/input/FGA_full_length_1_866.fasta` is provided.
The preparation script crops its match-state columns to the 86-residue Stage 0
target interval. It never constructs a fake paired target-peptide MSA.

No prediction was started by this preparation script.
"""


def main() -> int:
    parser = argparse.ArgumentParser(description="Prepare Stage 5 target-only MLM/MSA recovery controls.")
    parser.add_argument("--stage0-root", required=True)
    parser.add_argument("--output-root", required=True)
    parser.add_argument("--full-target-fasta", default="data/input/FGA_full_length_1_866.fasta")
    parser.add_argument("--full-target-a3m", default="")
    parser.add_argument("--max-homolog-msa-rows", type=int, default=512)
    parser.add_argument("--expected-target-length", type=int, default=86)
    parser.add_argument("--seeds-per-group", type=int, default=3)
    parser.add_argument("--models-per-seed", type=int, default=5)
    parser.add_argument("--recycles", type=int, default=6)
    parser.add_argument(
        "--afcycdesign-python",
        default="$HOME/fga_model_envs/colabdesign-py310/.pixi/envs/default/bin/python",
    )
    parser.add_argument("--colabdesign-source", default="$HOME/fga_model_envs/sources/ColabDesign-gamma-stage5")
    parser.add_argument("--python-overlay", default="$HOME/fga_model_envs/stage5_afcycdesign_python_overlay")
    parser.add_argument("--af-params", default="$HOME/fga_model_envs/af_params")
    args = parser.parse_args()

    logger = setup_logger("28_prepare_stage5_target_controls")
    append_run_header(logger, "28_prepare_stage5_target_controls.py")
    if args.seeds_per_group < 2:
        raise RuntimeError("Use at least two seeds per target-only control group")
    if args.models_per_seed != 5:
        raise RuntimeError("Target controls require all five model parameter sets")
    if args.recycles < 1 or args.max_homolog_msa_rows < 2:
        raise RuntimeError("Invalid recycle or homolog MSA row setting")

    project_root = resolve_path(".")
    stage0_root = _resolve_mixed_path(args.stage0_root)
    output_root = _resolve_mixed_path(args.output_root)
    assert_active_route_path(stage0_root, "Stage 28 Stage 0 root")
    assert_active_route_path(output_root, "Stage 28 output root", must_exist=False)
    output_dir = output_root / "07_structure_validation_target_controls"
    inputs_dir = output_dir / "inputs"
    msa_dir = inputs_dir / "msas"
    specs_dir = inputs_dir / "job_specs"
    jobs_dir = output_dir / "jobs"
    prediction_root = output_dir / "predictions"
    runner = project_root / "scripts" / "external" / "run_afcycdesign_target_recovery_control.py"
    if not runner.is_file():
        raise RuntimeError(f"Missing target-only runner: {runner}")

    target_pdb, target_sequence, site_indices, hotspot_indices, crop_start, crop_end = _stage0_context(
        stage0_root,
        args.expected_target_length,
    )
    full_fasta = _resolve_mixed_path(args.full_target_fasta)
    assert_active_route_path(full_fasta, "Stage 28 full target FASTA")
    full_records = _read_fasta(full_fasta)
    if len(full_records) != 1:
        raise RuntimeError("Full FGA FASTA must contain exactly one query sequence")
    full_sequence = _validate_sequence(full_records[0][1], "full FGA sequence")
    if full_sequence[crop_start - 1 : crop_end] != target_sequence:
        raise RuntimeError("Full FGA sequence interval does not match the Stage 0 target sequence")

    single_a3m = msa_dir / "FGA_Site2_target_single_sequence.a3m"
    _write_a3m(single_a3m, [("FGA_Site2_target_single_sequence", target_sequence)])
    homolog_a3m = msa_dir / f"FGA_full_homolog_MSA_crop_{crop_start}_{crop_end}.a3m"
    homolog_ready = False
    homolog_rows = 0
    homolog_source = ""
    if args.full_target_a3m:
        source_a3m = _resolve_mixed_path(args.full_target_a3m)
        assert_active_route_path(source_a3m, "Stage 28 homolog MSA A3M")
        if not source_a3m.is_file():
            raise RuntimeError(f"Provided full-target A3M does not exist: {source_a3m}")
        cropped_records = _crop_full_target_a3m(
            source_a3m=source_a3m,
            full_sequence=full_sequence,
            crop_start=crop_start,
            crop_end=crop_end,
            expected_crop_sequence=target_sequence,
            max_rows=args.max_homolog_msa_rows,
        )
        _write_a3m(homolog_a3m, cropped_records)
        homolog_ready = True
        homolog_rows = len(cropped_records)
        homolog_source = str(source_a3m)

    controls_config = [
        ("A", "single_sequence_mlm015", "single_sequence", single_a3m, True, 0.15, True),
        ("B", "single_sequence_no_mlm", "single_sequence", single_a3m, False, 0.0, True),
        ("C", "homolog_msa_mlm015", "unpaired_homolog_msa", homolog_a3m, True, 0.15, homolog_ready),
        ("D", "homolog_msa_no_mlm", "unpaired_homolog_msa", homolog_a3m, False, 0.0, homolog_ready),
    ]
    controls: list[dict[str, Any]] = []
    jobs: list[dict[str, Any]] = []
    runnable_scripts: list[Path] = []
    target_hash = _sha1_text(target_sequence)[:8]
    for group, label, msa_mode, a3m_path, use_mlm, replace_fraction, runnable in controls_config:
        msa_rows = 1 if msa_mode == "single_sequence" else homolog_rows
        a3m_sha1 = _sha1_file(a3m_path) if runnable else ""
        protocol_payload = {
            "protocol_version": PROTOCOL_VERSION,
            "control_group": group,
            "target_sequence_hash": target_hash,
            "target_a3m_sha1": a3m_sha1,
            "target_msa_mode": msa_mode,
            "use_mlm": use_mlm,
            "mlm_replace_fraction": replace_fraction,
            "model_type": "alphafold2_multimer_v3",
            "requested_recycles": args.recycles,
            "models_per_seed": args.models_per_seed,
            "template_mode": "none",
            "use_initial_guess": False,
            "peptide_included": False,
            "colabdesign_commit": COLABDESIGN_GAMMA_COMMIT,
        }
        protocol_hash = _protocol_hash(protocol_payload)
        controls.append(
            {
                "control_group": group,
                "control_label": label,
                "target_msa_mode": msa_mode,
                "target_a3m": a3m_path if runnable else "",
                "target_a3m_sha1": a3m_sha1,
                "msa_rows": msa_rows,
                "use_mlm": str(use_mlm).lower(),
                "mlm_replace_fraction": replace_fraction,
                "seeds_planned": args.seeds_per_group,
                "models_per_seed": args.models_per_seed,
                "jobs_planned": args.seeds_per_group,
                "jobs_runnable": args.seeds_per_group if runnable else 0,
                "model_evaluations_planned": args.seeds_per_group * args.models_per_seed,
                "template_mode": "none",
                "use_initial_guess": "false",
                "peptide_included": "false",
                "validation_test_type": "sequence_based_target_recovery_control",
                "protocol_hash": protocol_hash,
                "status": "prepared_not_run" if runnable else "blocked_missing_real_homolog_msa",
                "notes": (
                    "single target sequence; one A3M row"
                    if msa_mode == "single_sequence"
                    else (
                        f"real unpaired FGA homolog MSA cropped to {crop_start}-{crop_end}"
                        if runnable
                        else "requires verified real full-target FGA homolog A3M; no synthetic MSA generated"
                    )
                ),
            }
        )
        for seed in range(args.seeds_per_group):
            job_id = f"S5TC_{group}_{label}_prot{protocol_hash}_seed{seed:02d}"
            spec_path = specs_dir / f"{job_id}.json"
            prediction_dir = prediction_root / label / f"protocol_{protocol_hash}" / f"seed_{seed:02d}"
            spec = {
                "stage5_control_job_id": job_id,
                "control_group": group,
                "control_label": label,
                "protocol_hash": protocol_hash,
                "protocol_version": PROTOCOL_VERSION,
                "target_sequence": target_sequence,
                "target_sequence_length": len(target_sequence),
                "target_sequence_hash": target_hash,
                "target_reference_pdb": _to_wsl_path(target_pdb),
                "target_a3m": _to_wsl_path(a3m_path),
                "target_a3m_sha1": a3m_sha1,
                "target_msa_mode": msa_mode,
                "msa_rows_expected": msa_rows,
                "site2_target_indices_1based": site_indices,
                "hotspot_target_indices_1based": hotspot_indices,
                "seed": seed,
                "models_per_seed": args.models_per_seed,
                "requested_recycles": args.recycles,
                "forward_passes": args.recycles + 1,
                "use_mlm": use_mlm,
                "mlm_replace_fraction": replace_fraction,
                "template_mode": "none",
                "use_initial_guess": False,
                "peptide_included": False,
                "validation_test_type": "sequence_based_target_recovery_control",
                "colabdesign_commit": COLABDESIGN_GAMMA_COMMIT,
                "prediction_output_dir": _to_wsl_path(prediction_dir),
                "status": "prepared_not_run" if runnable else "blocked_missing_real_homolog_msa",
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
                    runnable=runnable,
                ),
            )
            if runnable:
                runnable_scripts.append(run_script)
            jobs.append(
                {
                    "stage5_control_job_id": job_id,
                    "control_group": group,
                    "protocol_hash": protocol_hash,
                    "target_msa_mode": msa_mode,
                    "target_a3m": a3m_path if runnable else "",
                    "target_a3m_sha1": a3m_sha1,
                    "msa_rows_expected": msa_rows,
                    "use_mlm": str(use_mlm).lower(),
                    "mlm_replace_fraction": replace_fraction,
                    "seed": seed,
                    "requested_recycles": args.recycles,
                    "forward_passes": args.recycles + 1,
                    "models_per_seed": args.models_per_seed,
                    "job_spec_json": spec_path,
                    "prediction_output_dir": prediction_dir,
                    "run_script": run_script,
                    "template_mode": "none",
                    "use_initial_guess": "false",
                    "peptide_included": "false",
                    "status": "prepared_not_run" if runnable else "blocked_missing_real_homolog_msa",
                }
            )

    preflight = jobs_dir / "check_stage5_target_controls_protocol.sh"
    master = jobs_dir / "run_stage5_target_controls_runnable.sh"
    _write_text(
        preflight,
        _preflight_script(args.afcycdesign_python, args.colabdesign_source, args.python_overlay, args.af_params),
    )
    _write_text(master, _master_script(preflight, runnable_scripts))
    write_csv(output_dir / "FGA_rfpeptides_stage5_target_control_manifest.csv", controls, CONTROL_FIELDS)
    write_csv(inputs_dir / "FGA_rfpeptides_stage5_target_control_jobs.csv", jobs, JOB_FIELDS)
    write_markdown(
        output_dir / "FGA_rfpeptides_stage5_target_control_plan.md",
        _report(
            controls=controls,
            jobs=jobs,
            output_dir=output_dir,
            homolog_source=homolog_source,
            crop_start=crop_start,
            crop_end=crop_end,
        ),
    )
    write_markdown(
        inputs_dir / "REAL_HOMOLOG_MSA_REQUIRED.md",
        f"""# Real FGA Homolog MSA Requirement

Groups C and D require a real unpaired homolog A3M for the complete 866-aa
FGA query in `{full_fasta}`. The query row must exactly match that FASTA.

Do not create a paired target-peptide MSA. The peptide is absent from these
target-only controls and has no natural homolog MSA.

After obtaining the A3M manually, regenerate the controls with:

```bash
python scripts/28_prepare_stage5_target_controls.py \\
  --full-target-a3m /path/to/FGA_full_length_homologs.a3m
```

The script verifies the full query, crops UniProt columns {crop_start}-{crop_end},
checks that the crop equals the 86-aa Stage 0 target, and requires at least one
non-query homolog before enabling C/D.
""",
    )

    logger.info("Target-only control groups planned: %s", len(controls))
    logger.info("Seed jobs planned: %s", len(jobs))
    logger.info("Seed jobs runnable now: %s", len(runnable_scripts))
    logger.info("Seed jobs blocked pending real homolog MSA: %s", len(jobs) - len(runnable_scripts))
    logger.info("Output directory: %s", output_dir)
    logger.info("No target-control or complex prediction was run.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
