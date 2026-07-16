from __future__ import annotations

import argparse
import json
from collections import OrderedDict
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping

from common import append_run_header, load_config, project_root, read_csv, read_fasta_sequence, setup_logger, write_csv


JOB_FIELDS = [
    "colabfold_job_id",
    "peptide_id",
    "core_sequence",
    "patch_id",
    "seed",
    "input_fasta",
    "output_dir",
    "msa_mode",
    "model_type",
    "num_models",
    "num_recycle",
    "target_length",
    "peptide_length",
    "hotspot_target_positions",
    "command",
    "status",
    "notes",
]


def _truthy(value: Any) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes", "pass"}


def _safe_job_token(value: str) -> str:
    keep = []
    for ch in str(value):
        if ch.isalnum() or ch in {"_", "-", "."}:
            keep.append(ch)
        else:
            keep.append("_")
    return "".join(keep).strip("_") or "candidate"


def _unique_candidates(rows: Iterable[Mapping[str, str]]) -> List[Dict[str, str]]:
    out: "OrderedDict[str, Dict[str, str]]" = OrderedDict()
    for row in rows:
        peptide_id = row.get("peptide_id") or row.get("raw_id") or row.get("boltz_job_id") or ""
        core = row.get("core_sequence", "")
        if not peptide_id or not core:
            continue
        if row.get("sequence_filter_pass") and not _truthy(row.get("sequence_filter_pass")):
            continue
        if peptide_id not in out:
            out[peptide_id] = {
                "peptide_id": peptide_id,
                "core_sequence": core,
                "patch_id": row.get("patch_id", ""),
                "hotspot_target_positions": row.get("hotspot_target_positions", ""),
                "source_notes": row.get("notes", ""),
            }
    return list(out.values())


def _write_multimer_fasta(path: Path, job_id: str, target_seq: str, peptide_seq: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    # ColabFold complex input uses colon-separated chains in one FASTA record.
    path.write_text(f">{job_id}\n{target_seq}:{peptide_seq}\n", encoding="utf-8", newline="\n")


def main() -> int:
    parser = argparse.ArgumentParser(description="Prepare ColabFold/AlphaFold-Multimer jobs for FGA cyclic peptide candidates.")
    parser.add_argument("--config", default="config/project.yaml")
    parser.add_argument(
        "--candidate-csv",
        default="results/boltz_predictions_raw3ghg_template_msa_top30/boltz_jobs.csv",
        help="Candidate/job CSV. Defaults to the Boltz top30 MSA job table and deduplicates peptide_id.",
    )
    parser.add_argument("--output-root", default="results/colabfold_predictions_top30_seed1")
    parser.add_argument("--target-fasta", default="data/input/FGA_chain_36_866.fasta")
    parser.add_argument("--seeds", type=int, default=1)
    parser.add_argument("--msa-mode", default="single_sequence")
    parser.add_argument("--model-type", default="alphafold2_multimer_v3")
    parser.add_argument("--num-models", type=int, default=1)
    parser.add_argument("--num-recycle", type=int, default=0)
    args = parser.parse_args()

    logger = setup_logger("16_prepare_colabfold_prediction_jobs")
    append_run_header(logger, "16_prepare_colabfold_prediction_jobs.py")
    load_config(args.config)

    target_seq = read_fasta_sequence(args.target_fasta)
    target_len = len(target_seq)
    output_root = project_root() / args.output_root
    input_dir = output_root / "inputs" / "fasta"
    rows = _unique_candidates(read_csv(args.candidate_csv))
    if not rows:
        write_csv(output_root / "colabfold_jobs.csv", [], JOB_FIELDS)
        logger.warning("No candidates found in %s", args.candidate_csv)
        return 0

    job_rows: List[Dict[str, Any]] = []
    for idx, cand in enumerate(rows, start=1):
        peptide_id = cand["peptide_id"]
        peptide_seq = cand["core_sequence"]
        for seed in range(1, args.seeds + 1):
            job_id = f"{cand.get('patch_id', 'Patch')}_afm_{idx:05d}_seed{seed}"
            fasta_path = input_dir / f"{_safe_job_token(job_id)}.fasta"
            out_dir = output_root / "outputs" / _safe_job_token(job_id)
            _write_multimer_fasta(fasta_path, job_id, target_seq, peptide_seq)
            rel_fasta = fasta_path.relative_to(project_root()).as_posix()
            rel_out_dir = out_dir.relative_to(project_root()).as_posix()
            cmd = (
                "colabfold_batch "
                f"--msa-mode {args.msa_mode} "
                f"--model-type {args.model_type} "
                f"--num-models {args.num_models} "
                f"--num-recycle {args.num_recycle} "
                f"--num-seeds 1 --random-seed {seed} "
                f"{rel_fasta} {rel_out_dir}"
            )
            notes = {
                "source_candidate_csv": args.candidate_csv,
                "source_candidate_id": peptide_id,
                "source_notes": cand.get("source_notes", ""),
                "colabfold_disulfide_constraint": "not_supported_by_input; check SG-SG geometry post hoc",
                "colabfold_patch_constraint": "not_supported_by_input; check patch contacts post hoc",
            }
            job_rows.append(
                {
                    "colabfold_job_id": job_id,
                    "peptide_id": peptide_id,
                    "core_sequence": peptide_seq,
                    "patch_id": cand.get("patch_id", ""),
                    "seed": seed,
                    "input_fasta": rel_fasta,
                    "output_dir": rel_out_dir,
                    "msa_mode": args.msa_mode,
                    "model_type": args.model_type,
                    "num_models": args.num_models,
                    "num_recycle": args.num_recycle,
                    "target_length": target_len,
                    "peptide_length": len(peptide_seq),
                    "hotspot_target_positions": cand.get("hotspot_target_positions", ""),
                    "command": cmd,
                    "status": "pending_external_model",
                    "notes": json.dumps(notes, ensure_ascii=False, sort_keys=True),
                }
            )

    write_csv(output_root / "colabfold_jobs.csv", job_rows, JOB_FIELDS)
    logger.info("Wrote ColabFold jobs: %s", len(job_rows))
    logger.info("Unique peptides: %s", len(rows))
    logger.info("Output root: %s", output_root)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
