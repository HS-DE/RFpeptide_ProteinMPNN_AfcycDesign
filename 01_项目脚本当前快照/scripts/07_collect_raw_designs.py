from __future__ import annotations

import argparse
import json
from pathlib import Path

from common import append_run_header, clean_sequence, load_config, project_root, read_csv, setup_logger, write_csv


RAW_FIELDS = ["raw_id", "method", "job_id", "patch_id", "target_region", "core_sequence", "core_length", "raw_score", "source_file", "notes"]
SKIPPED_FIELDS = RAW_FIELDS + ["skip_reason"]


def _is_candidate_csv(path: Path) -> bool:
    name = path.name.lower()
    if ".demo" in name:
        return False
    if "rejected" in name:
        return False
    return path.suffix.lower() == ".csv"


def _candidate_sources() -> list[Path]:
    root = project_root()
    paths: list[Path] = []
    for rel in [
        "results/raw_designs/manual_candidates.csv",
        "results/raw_designs/colabdesign_outputs",
        "results/raw_designs/rfdiffusion_outputs",
    ]:
        path = root / rel
        if path.is_file() and path.suffix.lower() == ".csv":
            paths.append(path)
        elif path.is_dir():
            paths.extend(sorted(p for p in path.rglob("*.csv") if _is_candidate_csv(p)))
    return paths


def _notes_dict(value: str) -> dict:
    if not value:
        return {}
    try:
        parsed = json.loads(value)
    except (TypeError, json.JSONDecodeError):
        return {}
    return parsed if isinstance(parsed, dict) else {}


def _truthy(value) -> bool:
    return str(value).lower() in {"true", "1", "yes", "pass"}


def _falsey(value) -> bool:
    return str(value).lower() in {"false", "0", "no", "fail"}


def _unoptimized_colabdesign_reason(record: dict, core: str) -> str:
    method = record.get("method", "")
    if "colabdesign" not in method.lower():
        return ""

    notes = _notes_dict(record.get("notes", ""))
    init_seq = clean_sequence(str(notes.get("init_sequence", "")))
    changed = notes.get("final_sequence_changed", "")
    if _truthy(changed):
        return ""
    if _falsey(changed):
        return "final_sequence_changed=false"
    if init_seq and init_seq == core:
        return "core_sequence_equals_init_sequence_without_change_flag"
    return ""


def main() -> int:
    parser = argparse.ArgumentParser(description="收集真实模型生成的 raw candidate peptide。")
    parser.add_argument("--config", default="config/project.yaml")
    args = parser.parse_args()

    logger = setup_logger("07_collect_raw_designs")
    append_run_header(logger, "07_collect_raw_designs.py")
    load_config(args.config)

    # =====================
    # Step 1. 搜索真实模型输出
    # =====================
    rows = []
    skipped_rows = []
    seen = set()
    duplicate_count = 0
    unoptimized_count = 0
    for src in _candidate_sources():
        rel_src = str(src.relative_to(project_root()))
        records = read_csv(src)
        logger.info("读取候选来源: %s rows=%s", rel_src, len(records))
        for idx, record in enumerate(records, start=1):
            core = clean_sequence(record.get("core_sequence") or record.get("sequence") or record.get("peptide") or "")
            if not core:
                logger.warning("跳过缺少 core_sequence 的记录: %s row=%s", rel_src, idx)
                continue
            unoptimized_reason = _unoptimized_colabdesign_reason(record, core)
            if unoptimized_reason:
                unoptimized_count += 1
                skipped = {
                    "raw_id": record.get("raw_id") or f"skipped_{unoptimized_count:06d}",
                    "method": record.get("method", "external_model_output"),
                    "job_id": record.get("job_id", ""),
                    "patch_id": record.get("patch_id", ""),
                    "target_region": record.get("target_region", "FGA_chain_36_866"),
                    "core_sequence": core,
                    "core_length": len(core),
                    "raw_score": record.get("raw_score", ""),
                    "source_file": rel_src,
                    "notes": record.get("notes", ""),
                    "skip_reason": unoptimized_reason,
                }
                skipped_rows.append(skipped)
                continue
            patch_id = record.get("patch_id", "")
            dedup_key = (patch_id, core)
            if dedup_key in seen:
                duplicate_count += 1
                continue
            seen.add(dedup_key)
            raw_id = record.get("raw_id") or f"raw_{len(rows) + 1:06d}"
            rows.append(
                {
                    "raw_id": raw_id,
                    "method": record.get("method", "external_model_output"),
                    "job_id": record.get("job_id", ""),
                    "patch_id": patch_id,
                    "target_region": record.get("target_region", "FGA_chain_36_866"),
                    "core_sequence": core,
                    "core_length": len(core),
                    "raw_score": record.get("raw_score", ""),
                    "source_file": rel_src,
                    "notes": record.get("notes", "real model/user-provided candidate"),
                }
            )

    # =====================
    # Step 2. 输出 raw candidate 表或空 demo 表
    # =====================
    if rows:
        write_csv("results/raw_designs/FGA_raw_candidates.csv", rows, RAW_FIELDS)
        write_csv("results/raw_designs/FGA_unoptimized_candidates_skipped.csv", skipped_rows, SKIPPED_FIELDS)
        if duplicate_count:
            logger.warning("Skipped duplicate raw candidates: %s", duplicate_count)
        if unoptimized_count:
            logger.warning("Skipped unoptimized ColabDesign candidates: %s", unoptimized_count)
        logger.info("输出真实 raw candidates: %s", len(rows))
    else:
        write_csv("results/raw_designs/FGA_raw_candidates.csv", [], RAW_FIELDS)
        write_csv("results/raw_designs/FGA_unoptimized_candidates_skipped.csv", skipped_rows, SKIPPED_FIELDS)
        if unoptimized_count:
            logger.warning("Skipped unoptimized ColabDesign candidates: %s", unoptimized_count)
        demo_rows = [
            {
                "raw_id": "DEMO_DO_NOT_USE_0001",
                "method": "demo_placeholder_not_model_output",
                "job_id": "demo",
                "patch_id": "Patch_A",
                "target_region": "FGA_chain_36_866",
                "core_sequence": "DEMO_SEQUENCE_NOT_A_PEPTIDE",
                "core_length": 27,
                "raw_score": "",
                "source_file": "generated_demo_schema_only",
                "notes": "demo placeholder; not a peptide candidate; never used by pipeline final ranking",
            }
        ]
        write_csv("results/raw_designs/FGA_raw_candidates.demo.csv", demo_rows, RAW_FIELDS)
        logger.warning("未找到真实模型输出；只生成明确标记的 demo/toy 文件，不进入最终候选。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
