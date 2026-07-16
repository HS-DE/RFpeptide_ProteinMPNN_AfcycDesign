from __future__ import annotations

import argparse

from common import append_run_header, load_config, read_csv, rows_to_markdown, setup_logger, write_markdown


def main() -> int:
    parser = argparse.ArgumentParser(description="导出中文 FGA 环肽设计报告。")
    parser.add_argument("--config", default="config/project.yaml")
    args = parser.parse_args()

    logger = setup_logger("13_export_final_report")
    append_run_header(logger, "13_export_final_report.py")
    config = load_config(args.config)

    # =====================
    # Step 1. 汇总当前输出状态
    # =====================
    regions = read_csv("data/annotations/FGA_regions.csv")
    patches = read_csv("data/annotations/FGA_epitope_candidates.csv")
    design_jobs = read_csv("results/raw_designs/design_jobs.csv")
    raw = read_csv("results/raw_designs/FGA_raw_candidates.csv")
    scored = read_csv("results/filtered/FGA_scored_candidates.csv")
    negative = read_csv("results/filtered/FGA_negative_screen_summary.csv")
    top10 = read_csv("results/final/FGA_top10_synthesis_priority.csv")
    has_real_final = bool(top10)
    model_status = (
        "当前已经存在完成排序的真实候选。"
        if has_real_final
        else "当前尚未运行真实生成/复合物预测模型，或尚未完成负筛选，因此没有可用于合成的真实候选序列。"
    )

    # =====================
    # Step 2. 生成报告正文
    # =====================
    top10_md = rows_to_markdown(
        top10,
        ["peptide_id", "target_patch", "core_sequence", "final_synthesis_format", "final_score", "priority"],
        "当前没有通过完整模型生成、复合物预测评分和负筛选的 top10 候选。",
    )
    patch_md = rows_to_markdown(
        patches,
        ["patch_id", "patch_type", "chain_id", "n_surface_residues", "priority", "risk_level"],
        "当前尚未生成 patch 候选。",
    )
    region_md = rows_to_markdown(
        regions,
        ["region_name", "start", "end", "length", "use_for_design", "priority"],
        "当前尚未生成 FGA 区域注释。",
    )
    report = f"""
# FGA Cys-Cys 环肽计算设计报告

## 1. 项目目的

本项目针对 human FGA / fibrinogen alpha chain / UniProt P02671，建立 Cys-Cys 二硫键环肽计算设计流程。候选肽统一采用 `Biotin-PEG4-GSG-[core_sequence]-NH2` 格式，其中 core sequence 必须以 Cys 开头和结尾，内部不含额外 Cys。

## 2. 输入数据说明

输入 Excel 为 `data/input/高丰度蛋白信息.xlsx`，按 `Gene == FGA` 或 `UniprotID == P02671` 提取 FGA full-length 序列。当前区域输出如下：

{region_md}

## 3. 为什么不用 full-length precursor 直接设计

FGA full-length precursor 包含 signal peptide 和 fibrinopeptide A。1-19 aa 为 signal peptide，不应作为真实血浆靶点；20-35 aa 为 fibrinopeptide A，存在被切除风险。因此第一版主设计对象优先使用 `FGA_chain_36_866`，并保留 full-length 与 20-866 区域用于记录和对照。

## 4. 为什么优先考虑 native fibrinogen 中暴露的 FGA 区域

FGA 在血浆中主要作为 fibrinogen 复合体的一部分存在，不宜只把孤立 FGA 单链作为靶标。设计流程优先使用 native human fibrinogen experimental structure 中可见并暴露的 FGA 区域。

## 5. 结构来源说明

优先结构为 RCSB PDB `3GHG` human fibrinogen X-ray structure。`03_prepare_structures.py` 会下载或使用本地 `data/structures/raw/3GHG.pdb`，并输出清理后的 `data/structures/prepared/fibrinogen_3GHG_clean.pdb`。AlphaFold P02671 单链结构仅作为补充，不作为唯一设计依据。

## 6. Patch 选择方法

`04_map_fga_structure.py` 通过 PDB 链序列和 FGA UniProt 序列局部比对识别 FGA 链，不固定假设 chain ID。`05_select_surface_patches.py` 优先支持 freesasa；若不可用，则使用 residue neighbor count 作为表面暴露近似指标，并输出 Patch_A/Patch_B/Patch_C：

{patch_md}

## 7. 环肽设计规则

本阶段只做方案 A：Cys-Cys disulfide cyclic peptide。核心肽长度限制为 10-18 aa，主力长度为 12/14/16 aa，最终合成格式自动生成为 `Biotin-PEG4-GSG-{{core_sequence}}-NH2`。不生成 lactam、click、head-to-tail、linear peptide 或 miniprotein 方案。

## 8. 序列过滤规则

硬过滤要求 core sequence 以 C 开头、以 C 结尾、内部不含额外 Cys、长度 10-18 aa、净电荷 -3 到 +3、连续强疏水残基不超过 4、W 和 M 各不超过 1，并排除明显低复杂度或 poly-basic/poly-acidic 序列。

## 9. 复合物预测评分规则

真实候选必须经过复合物预测或 docking 评分后才能进入最终排序。评分字段包括 ipTM/pTM、interface PAE、peptide pLDDT、interface contacts、peptide-target minimum distance、patch consistency、multi-seed pose consistency 和 Cys-Cys geometry。当前真实 raw candidates 数: {len(raw)}；真实 scored candidates 数: {len(scored)}。

## 10. 负筛选规则

负筛选目标包括 ALB、APOA1、TF、A2M、C3、IGG_FC。负筛选只用于降低明显非特异 sticky peptide 风险，不能证明绝对特异。当前负筛选结果行数: {len(negative)}。

## 11. Top10 候选表

{top10_md}

## 12. 风险和限制

1. FGA 在血浆中主要位于 fibrinogen 复合体中，因此单链预测存在偏差。
2. 3GHG 是实验结构，但不一定覆盖 FGA 全部柔性区域。
3. FGA alphaC 等区域可能柔性较高，相关候选应标记为 exploratory/high risk。
4. Cys-Cys 二硫键环肽在计算结构中需要检查首尾 Cys 几何合理性。
5. Biotin-PEG4-GSG linker 未必完整进入结构模型，因此最终磁珠固定后的表现仍需实验验证。
6. 所有候选均为计算候选，不能宣称已验证结合。
7. 负筛选只能降低明显非特异风险，不能证明绝对特异。

## 13. 下一步建议

先在 GPU/Colab 环境运行 `results/raw_designs/design_jobs.csv` 中的 ColabDesign-cyclic-binder 或 RFdiffusion/RFpeptides 任务；把真实模型输出放回指定目录后，再运行 raw candidate 收集、硬过滤、复合物预测 job 生成、ColabFold/Boltz-2 预测、评分、负筛选和最终排序。

## 当前状态声明

{model_status}

这些序列若未来出现于 final 表，也只能解释为计算生成的候选环肽 binder；尚未经过实验验证，不能宣称具有已验证结合能力，后续仍需要合成和实验验证。
"""

    # =====================
    # Step 3. 写出报告
    # =====================
    write_markdown("results/final/FGA_design_report.md", report)
    logger.info("输出报告: results/final/FGA_design_report.md")
    logger.info("design job 数: %s", len(design_jobs))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
