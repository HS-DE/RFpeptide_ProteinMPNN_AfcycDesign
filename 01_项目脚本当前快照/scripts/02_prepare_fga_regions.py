from __future__ import annotations

import argparse

from common import append_run_header, load_config, read_fasta_sequence, setup_logger, write_csv, write_fasta
from region_utils import build_region_rows, fasta_header


def main() -> int:
    parser = argparse.ArgumentParser(description="生成 FGA full-length/extracellular/main-chain 区域 FASTA。")
    parser.add_argument("--config", default="config/project.yaml")
    args = parser.parse_args()

    logger = setup_logger("02_prepare_fga_regions")
    append_run_header(logger, "02_prepare_fga_regions.py")
    config = load_config(args.config)

    # =====================
    # Step 1. 读取 full-length FASTA
    # =====================
    sequence = read_fasta_sequence("data/input/FGA_full_length_1_866.fasta")
    expected_end = int(config["target_regions"]["full_length"]["end"])
    if len(sequence) < expected_end:
        raise ValueError(f"FGA full-length 序列长度 {len(sequence)} 小于配置终点 {expected_end}")
    logger.info("输入 full-length 长度: %s", len(sequence))

    # =====================
    # Step 2. 生成区域序列
    # =====================
    rows = build_region_rows(sequence, config)
    for row in rows:
        if row["region_name"] == "full_length_1_866":
            path = "data/input/FGA_full_length_1_866.fasta"
        elif row["region_name"] == "extracellular_20_866":
            path = "data/input/FGA_extracellular_20_866.fasta"
        else:
            path = "data/input/FGA_chain_36_866.fasta"
        write_fasta(path, fasta_header(row["region_name"], config), row["sequence"])
        logger.info("输出区域 FASTA: %s length=%s", path, row["length"])

    # =====================
    # Step 3. 输出区域注释表
    # =====================
    fields = ["region_name", "start", "end", "length", "use_for_design", "priority", "note", "sequence"]
    write_csv("data/annotations/FGA_regions.csv", rows, fields)
    logger.info("输出文件: data/annotations/FGA_regions.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
