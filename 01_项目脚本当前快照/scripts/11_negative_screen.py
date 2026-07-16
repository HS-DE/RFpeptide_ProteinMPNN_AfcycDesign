from __future__ import annotations

import argparse

from common import append_run_header, load_config, read_csv, setup_logger, write_csv


NEGATIVE_FIELDS = [
    "peptide_id",
    "negative_target",
    "negative_score",
    "negative_interface_pae",
    "negative_contacts",
    "non_specific_risk",
    "negative_screen_pass",
    "notes",
]


def main() -> int:
    parser = argparse.ArgumentParser(description="汇总真实负筛选结果。")
    parser.add_argument("--config", default="config/project.yaml")
    args = parser.parse_args()

    logger = setup_logger("11_negative_screen")
    append_run_header(logger, "11_negative_screen.py")
    load_config(args.config)

    # =====================
    # Step 1. 读取外部负筛选结果
    # =====================
    manual_rows = read_csv("results/filtered/manual_negative_screen_summary.csv")
    if manual_rows:
        rows = []
        for row in manual_rows:
            rows.append(
                {
                    "peptide_id": row.get("peptide_id", ""),
                    "negative_target": row.get("negative_target", ""),
                    "negative_score": row.get("negative_score", ""),
                    "negative_interface_pae": row.get("negative_interface_pae", ""),
                    "negative_contacts": row.get("negative_contacts", ""),
                    "non_specific_risk": row.get("non_specific_risk", ""),
                    "negative_screen_pass": row.get("negative_screen_pass", ""),
                    "notes": row.get("notes", "real externally supplied negative screen result"),
                }
            )
        write_csv("results/filtered/FGA_negative_screen_summary.csv", rows, NEGATIVE_FIELDS)
        logger.info("输出真实负筛选结果行数: %s", len(rows))
        return 0

    # =====================
    # Step 2. 未运行负筛选时输出空表
    # =====================
    write_csv("results/filtered/FGA_negative_screen_summary.csv", [], NEGATIVE_FIELDS)
    logger.warning("未找到 manual_negative_screen_summary.csv；负筛选未运行，最终排序不会生成可合成 top10。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
