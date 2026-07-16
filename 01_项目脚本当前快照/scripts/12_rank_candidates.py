from __future__ import annotations

import argparse
from collections import defaultdict

from common import append_run_header, load_config, read_csv, setup_logger, write_csv
from ranking import FINAL_COLUMNS, rank_rows


def _truthy(value: str) -> bool:
    return str(value).lower() in {"true", "1", "yes", "pass"}


def main() -> int:
    parser = argparse.ArgumentParser(description="综合评分并导出 top50/top10。")
    parser.add_argument("--config", default="config/project.yaml")
    args = parser.parse_args()

    logger = setup_logger("12_rank_candidates")
    append_run_header(logger, "12_rank_candidates.py")
    config = load_config(args.config)

    # =====================
    # Step 1. 读取评分、过滤和负筛选结果
    # =====================
    scored = [r for r in read_csv("results/filtered/FGA_scored_candidates.csv") if _truthy(r.get("complex_score_pass", ""))]
    filtered = {r.get("raw_id", ""): r for r in read_csv("results/filtered/FGA_hard_filtered_candidates.csv") if _truthy(r.get("sequence_filter_pass", ""))}
    negative_rows = read_csv("results/filtered/FGA_negative_screen_summary.csv")
    design_jobs = read_csv("results/raw_designs/design_jobs.csv")
    patch_hotspots = {}
    for job in design_jobs:
        patch_hotspots.setdefault(job.get("patch_id", ""), job.get("hotspot_residues", ""))

    if not scored:
        write_csv("results/final/FGA_top50_candidates.csv", [], FINAL_COLUMNS)
        write_csv("results/final/FGA_top10_synthesis_priority.csv", [], FINAL_COLUMNS)
        logger.warning("没有通过真实复合物评分的候选；输出空 final schema。")
        return 0
    if not negative_rows and config["negative_screen"].get("enabled", True):
        write_csv("results/final/FGA_top50_candidates.csv", [], FINAL_COLUMNS)
        write_csv("results/final/FGA_top10_synthesis_priority.csv", [], FINAL_COLUMNS)
        logger.warning("负筛选未运行；为避免未经完整评分的最终候选，输出空 final schema。")
        return 0

    negative_by_peptide = defaultdict(list)
    for row in negative_rows:
        negative_by_peptide[row.get("peptide_id", "")].append(row)

    # =====================
    # Step 2. 合并字段并计算综合排序
    # =====================
    rows = []
    prefix = config["peptide_design"]["final_format_prefix"]
    suffix = config["peptide_design"]["final_format_suffix"]
    for row in scored:
        peptide_id = row.get("peptide_id", "")
        filt = filtered.get(peptide_id, {})
        negs = negative_by_peptide.get(peptide_id, [])
        if not negs and config["negative_screen"].get("enabled", True):
            negative_flag = "missing"
        else:
            negative_fail = any(str(n.get("negative_screen_pass", "")).lower() != "pass" for n in negs)
            negative_flag = "fail" if negative_fail else "pass"
        core = row.get("core_sequence") or filt.get("core_sequence", "")
        priority = "exclude" if negative_flag != "pass" else "ranked"
        notes = "real model-scored candidate"
        if negative_flag == "fail":
            notes = "negative screen fail; excluded"
        elif negative_flag == "missing":
            notes = "negative screen missing; excluded"
        final_row = {
            "peptide_id": peptide_id,
            "target": config["project"]["target_gene"],
            "target_uniprot": config["project"]["target_uniprot"],
            "target_patch": row.get("patch_id") or filt.get("patch_id", ""),
            "uniprot_region": filt.get("target_region", "FGA_chain_36_866"),
            "hotspot_residues": patch_hotspots.get(row.get("patch_id") or filt.get("patch_id", ""), ""),
            "core_sequence": core,
            "core_length": len(core),
            "final_synthesis_format": f"{prefix}{core}{suffix}" if core else "",
            "cyclization": config["peptide_design"]["cyclization"],
            "mean_iptm": row.get("mean_iptm", ""),
            "best_iptm": row.get("best_iptm", ""),
            "mean_interface_pae": row.get("mean_interface_pae", ""),
            "best_interface_pae": row.get("best_interface_pae", ""),
            "peptide_plddt": row.get("mean_peptide_plddt", ""),
            "interface_contacts": row.get("interface_contacts", ""),
            "pose_consistency_rmsd": row.get("pose_consistency_rmsd", ""),
            "cys_cys_geometry": row.get("cys_cys_geometry", ""),
            "net_charge": filt.get("net_charge", ""),
            "hydrophobicity_flag": "pass" if _truthy(filt.get("hydrophobicity_pass", "")) else "fail",
            "sequence_filter_pass": filt.get("sequence_filter_pass", ""),
            "negative_screen_flag": negative_flag,
            "patch_consistency_flag": row.get("patch_consistency_flag", ""),
            "mean_peptide_plddt": row.get("mean_peptide_plddt", ""),
            "priority": priority,
            "notes": notes,
        }
        rows.append(final_row)

    ranked = [r for r in rank_rows(rows) if r.get("priority") != "exclude" and r.get("negative_screen_flag") == "pass"]
    top50 = ranked[: int(config["ranking"]["top_n_candidates"])]
    top10 = [r for r in ranked if r.get("cys_cys_geometry") == "pass"][: int(config["ranking"]["top_n_synthesis_priority"])]

    # =====================
    # Step 3. 输出 final 表
    # =====================
    write_csv("results/final/FGA_top50_candidates.csv", top50, FINAL_COLUMNS)
    write_csv("results/final/FGA_top10_synthesis_priority.csv", top10, FINAL_COLUMNS)
    logger.info("输出 top50=%s top10=%s", len(top50), len(top10))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
