from __future__ import annotations

import argparse
from collections import defaultdict

from common import append_run_header, load_config, read_csv, setup_logger, write_csv
from ranking import normalized_complex_score


SUMMARY_FIELDS = [
    "peptide_id",
    "core_sequence",
    "patch_id",
    "n_seeds",
    "best_seed",
    "mean_iptm",
    "best_iptm",
    "mean_interface_pae",
    "best_interface_pae",
    "mean_peptide_plddt",
    "interface_contacts",
    "pose_consistency_rmsd",
    "patch_consistency_flag",
    "cys_cys_geometry",
    "complex_score",
    "complex_score_pass",
    "notes",
]


def _to_float(value, default=0.0) -> float:
    try:
        if value in ("", None):
            return default
        return float(value)
    except (TypeError, ValueError):
        return default


def _truthy(value: str) -> bool:
    return str(value).lower() in {"true", "1", "yes", "pass"}


def main() -> int:
    parser = argparse.ArgumentParser(description="解析真实复合物预测输出并计算评分。")
    parser.add_argument("--config", default="config/project.yaml")
    args = parser.parse_args()

    logger = setup_logger("10_score_complex_predictions")
    append_run_header(logger, "10_score_complex_predictions.py")
    config = load_config(args.config)

    # =====================
    # Step 1. 读取用户放回的真实评分表
    # =====================
    manual_rows = read_csv("results/complex_predictions/manual_complex_prediction_summary.csv")
    if not manual_rows:
        write_csv("results/complex_predictions/FGA_complex_prediction_summary.csv", [], SUMMARY_FIELDS)
        write_csv("results/filtered/FGA_scored_candidates.csv", [], SUMMARY_FIELDS)
        logger.warning("未找到真实复合物预测评分文件 manual_complex_prediction_summary.csv；输出空评分 schema。")
        return 0

    # =====================
    # Step 2. 按阈值计算 pass/fail 和 complex_score
    # =====================
    thresholds = config["scoring_thresholds"]
    filtered_lookup = {r.get("raw_id", ""): r for r in read_csv("results/filtered/FGA_hard_filtered_candidates.csv") if _truthy(r.get("sequence_filter_pass", ""))}
    rows = []
    skipped_unknown = 0
    for row in manual_rows:
        peptide_id = row.get("peptide_id") or row.get("raw_id", "")
        source = filtered_lookup.get(peptide_id)
        if not source:
            skipped_unknown += 1
            logger.warning("Skipping complex score for unknown or unfiltered peptide_id: %s", peptide_id)
            continue
        merged = dict(row)
        merged.setdefault("peptide_id", peptide_id)
        merged.setdefault("core_sequence", source.get("core_sequence", ""))
        merged.setdefault("patch_id", source.get("patch_id", ""))
        merged.setdefault("n_seeds", row.get("n_seeds", ""))
        merged.setdefault("best_seed", row.get("best_seed", ""))
        merged.setdefault("pose_consistency_rmsd", row.get("pose_consistency_rmsd", ""))
        merged.setdefault("patch_consistency_flag", row.get("patch_consistency_flag", ""))
        merged.setdefault("cys_cys_geometry", row.get("cys_cys_geometry", ""))
        score = normalized_complex_score(merged)
        pass_flag = (
            _to_float(merged.get("best_interface_pae"), 999) <= float(thresholds["max_interface_pae"])
            and _to_float(merged.get("mean_peptide_plddt"), 0) >= float(thresholds["min_peptide_plddt"])
            and _to_float(merged.get("interface_contacts"), 0) >= float(thresholds["min_interface_contacts"])
            and str(merged.get("patch_consistency_flag", "")).lower() == "pass"
            and str(merged.get("cys_cys_geometry", "")).lower() == "pass"
        )
        merged["complex_score"] = round(score, 4)
        merged["complex_score_pass"] = pass_flag
        merged["notes"] = merged.get("notes", "real externally supplied complex prediction score")
        rows.append(merged)

    # =====================
    # Step 3. 输出评分结果
    # =====================
    write_csv("results/complex_predictions/FGA_complex_prediction_summary.csv", rows, SUMMARY_FIELDS)
    write_csv("results/filtered/FGA_scored_candidates.csv", rows, SUMMARY_FIELDS)
    if skipped_unknown:
        logger.warning("Skipped complex score rows without matching hard-filtered candidate: %s", skipped_unknown)
    logger.info("解析真实复合物预测评分行数: %s", len(rows))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
