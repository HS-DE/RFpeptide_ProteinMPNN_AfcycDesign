from __future__ import annotations

import argparse
import sys

from common import (
    append_run_header,
    executable_available,
    import_available,
    load_config,
    project_root,
    resolve_path,
    setup_logger,
    write_csv,
)


def main() -> int:
    parser = argparse.ArgumentParser(description="检查 FGA 环肽设计流程环境。")
    parser.add_argument("--config", default="config/project.yaml")
    parser.add_argument("--strict-required-deps", action="store_true", help="缺少必需依赖时返回失败。")
    args = parser.parse_args()

    logger = setup_logger("00_check_environment")
    append_run_header(logger, "00_check_environment.py")

    # =====================
    # Step 1. 读取配置并检查目录
    # =====================
    config = load_config(args.config)
    root = project_root()
    required_dirs = [
        "config",
        "data/input",
        "data/structures/raw",
        "data/structures/prepared",
        "data/annotations",
        "scripts",
        "tests",
        "results/raw_designs",
        "results/complex_predictions",
        "results/filtered",
        "results/final",
        "logs",
        "env",
    ]
    for rel in required_dirs:
        path = root / rel
        status = "ok" if path.exists() else "missing"
        logger.info("目录检查: %s -> %s", rel, status)

    # =====================
    # Step 2. 检查 Python 和依赖
    # =====================
    rows = []
    python_ok = sys.version_info >= (3, 10)
    rows.append({"name": "python>=3.10", "type": "required", "available": python_ok, "note": sys.version.split()[0]})
    required_modules = {
        "pandas": "pandas",
        "numpy": "numpy",
        "biopython": "Bio",
        "pyyaml": "yaml",
        "openpyxl": "openpyxl",
        "scipy": "scipy",
        "scikit-learn": "sklearn",
        "pytest": "pytest",
    }
    for label, module in required_modules.items():
        available = import_available(module)
        rows.append({"name": label, "type": "required", "available": available, "note": module})
        if available:
            logger.info("必需依赖可用: %s", label)
        else:
            logger.warning("缺少必需依赖: %s；建议 conda env create -f env/environment.yml", label)

    optional_checks = {
        "freesasa": import_available("freesasa"),
        "mdtraj": import_available("mdtraj"),
        "pymol": import_available("pymol"),
        "colabfold_batch": executable_available("colabfold_batch"),
        "boltz": executable_available("boltz"),
    }
    for label, available in optional_checks.items():
        rows.append({"name": label, "type": "optional", "available": available, "note": "warning only"})
        if not available:
            logger.warning("可选工具不可用: %s；基础流程不中断", label)

    # =====================
    # Step 3. 检查输入文件
    # =====================
    excel_file = resolve_path(config["input"]["excel_file"])
    pdb_file = resolve_path(config["structures"]["primary_pdb_file"])
    rows.append({"name": "input_excel", "type": "required_file", "available": excel_file.exists(), "note": str(excel_file)})
    rows.append({"name": "3GHG_pdb", "type": "structure_file", "available": pdb_file.exists(), "note": str(pdb_file)})
    if not excel_file.exists():
        logger.error("缺少输入 Excel: %s", excel_file)
    else:
        logger.info("输入 Excel 存在: %s", excel_file)
    if not pdb_file.exists():
        logger.warning("3GHG.pdb 尚不存在；03_prepare_structures.py 会尝试从 RCSB 下载")
    else:
        logger.info("3GHG.pdb 存在: %s", pdb_file)

    # =====================
    # Step 4. 输出环境检查表
    # =====================
    write_csv("logs/environment_check.csv", rows, ["name", "type", "available", "note"])
    missing_required = [r["name"] for r in rows if r["type"] == "required" and not r["available"]]
    if args.strict_required_deps and missing_required:
        logger.error("严格模式下缺少必需依赖: %s", ", ".join(missing_required))
        return 2
    logger.info("环境检查完成；缺少依赖已记录，基础准备流程会尽量继续。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
