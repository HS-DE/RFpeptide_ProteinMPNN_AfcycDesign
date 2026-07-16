from __future__ import annotations

import argparse

from common import append_run_header, load_config, read_csv, setup_logger, write_csv
from sequence_filters import filter_candidate_sequence


FILTER_FIELDS = [
    "raw_id",
    "method",
    "job_id",
    "patch_id",
    "target_region",
    "core_sequence",
    "core_length",
    "raw_score",
    "source_file",
    "notes",
    "starts_with_cys",
    "ends_with_cys",
    "internal_cys_count",
    "core_length_pass",
    "net_charge",
    "charge_pass",
    "hydrophobic_run_max",
    "hydrophobicity_pass",
    "w_count",
    "m_count",
    "low_complexity_flag",
    "sequence_filter_pass",
    "filter_notes",
]


def main() -> int:
    parser = argparse.ArgumentParser(description="对真实 raw candidates 做 Cys-Cys 环肽硬过滤。")
    parser.add_argument("--config", default="config/project.yaml")
    args = parser.parse_args()

    logger = setup_logger("08_filter_sequences")
    append_run_header(logger, "08_filter_sequences.py")
    config = load_config(args.config)

    # =====================
    # Step 1. 读取 raw candidates
    # =====================
    raw_rows = read_csv("results/raw_designs/FGA_raw_candidates.csv")
    if not raw_rows:
        write_csv("results/filtered/FGA_hard_filtered_candidates.csv", [], FILTER_FIELDS)
        logger.warning("没有真实 raw candidates；输出空 hard-filter schema。")
        return 0

    # =====================
    # Step 2. 执行硬过滤
    # =====================
    out_rows = []
    pass_count = 0
    for row in raw_rows:
        result = filter_candidate_sequence(row.get("core_sequence", ""), config)
        merged = dict(row)
        merged.update(result)
        merged["core_length"] = len(result["core_sequence"])
        out_rows.append(merged)
        if result["sequence_filter_pass"]:
            pass_count += 1

    # =====================
    # Step 3. 输出过滤结果
    # =====================
    write_csv("results/filtered/FGA_hard_filtered_candidates.csv", out_rows, FILTER_FIELDS)
    logger.info("输入候选数: %s；硬过滤通过数: %s", len(raw_rows), pass_count)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
