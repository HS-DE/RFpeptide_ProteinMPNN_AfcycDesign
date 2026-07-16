from __future__ import annotations

import argparse

from common import append_run_header, executable_available, load_config, project_root, read_csv, read_fasta_sequence, resolve_path, setup_logger, write_csv


JOB_FIELDS = ["prediction_job_id", "peptide_id", "core_sequence", "target_pdb", "patch_id", "seed", "engine", "input_fasta", "output_dir", "command", "status"]


def _truthy(value: str) -> bool:
    return str(value).lower() in {"true", "1", "yes", "pass"}


def main() -> int:
    parser = argparse.ArgumentParser(description="为 hard-filtered candidates 生成复合物预测 job。")
    parser.add_argument("--config", default="config/project.yaml")
    args = parser.parse_args()

    logger = setup_logger("09_prepare_complex_prediction_jobs")
    append_run_header(logger, "09_prepare_complex_prediction_jobs.py")
    config = load_config(args.config)

    # =====================
    # Step 1. 读取通过硬过滤的真实候选
    # =====================
    candidates = [r for r in read_csv("results/filtered/FGA_hard_filtered_candidates.csv") if _truthy(r.get("sequence_filter_pass", ""))]
    if not candidates:
        write_csv("results/complex_predictions/complex_prediction_jobs.csv", [], JOB_FIELDS)
        logger.warning("没有通过硬过滤的真实候选；输出空复合物预测 job schema。")
        return 0

    target_seq = read_fasta_sequence("data/input/FGA_chain_36_866.fasta")
    input_dir = project_root() / "results/complex_predictions/inputs"
    input_dir.mkdir(parents=True, exist_ok=True)
    colabfold_available = executable_available("colabfold_batch")
    boltz_available = executable_available("boltz")

    # =====================
    # Step 2. 为每个候选和 seed 生成任务
    # =====================
    rows = []
    engines = config["complex_prediction"]["prediction_engines"]
    seeds = range(1, int(config["complex_prediction"]["seeds_per_candidate"]) + 1)
    target_pdb = resolve_path(config["structures"]["cleaned_pdb_file"])
    for cand in candidates:
        peptide_id = cand.get("raw_id", "")
        for engine in engines:
            normalized_engine = "boltz" if "boltz" in engine else "colabfold"
            available = boltz_available if normalized_engine == "boltz" else colabfold_available
            for seed in seeds:
                job_id = f"{peptide_id}_{normalized_engine}_seed{seed}"
                fasta_path = input_dir / f"{job_id}.fasta"
                fasta_path.write_text(
                    f">FGA_chain_36_866__{peptide_id}\n{target_seq}:{cand['core_sequence']}\n",
                    encoding="utf-8",
                    newline="\n",
                )
                out_dir = f"results/complex_predictions/{job_id}"
                if normalized_engine == "colabfold":
                    command = f"colabfold_batch --num-seeds 1 --random-seed {seed} {fasta_path.as_posix()} {out_dir}"
                else:
                    command = f"boltz predict {fasta_path.as_posix()} --out_dir {out_dir} --seed {seed}"
                rows.append(
                    {
                        "prediction_job_id": job_id,
                        "peptide_id": peptide_id,
                        "core_sequence": cand["core_sequence"],
                        "target_pdb": str(target_pdb),
                        "patch_id": cand.get("patch_id", ""),
                        "seed": seed,
                        "engine": normalized_engine,
                        "input_fasta": str(fasta_path.relative_to(project_root())),
                        "output_dir": out_dir,
                        "command": command,
                        "status": "ready" if available else "pending_external_model",
                    }
                )

    # =====================
    # Step 3. 输出预测 job 表
    # =====================
    write_csv("results/complex_predictions/complex_prediction_jobs.csv", rows, JOB_FIELDS)
    logger.info("输出复合物预测 job 数: %s", len(rows))
    if not colabfold_available:
        logger.warning("未检测到 colabfold_batch；ColabFold job 仅作为待运行任务。")
    if not boltz_available:
        logger.warning("未检测到 boltz；Boltz-2 job 仅作为待运行任务。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
