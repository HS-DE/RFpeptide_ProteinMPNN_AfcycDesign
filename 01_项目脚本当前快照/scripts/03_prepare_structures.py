from __future__ import annotations

import argparse
import urllib.request

from common import append_run_header, load_config, resolve_path, setup_logger
from pdb_utils import clean_pdb


def main() -> int:
    parser = argparse.ArgumentParser(description="下载/清理 3GHG native fibrinogen PDB。")
    parser.add_argument("--config", default="config/project.yaml")
    args = parser.parse_args()

    logger = setup_logger("03_prepare_structures")
    append_run_header(logger, "03_prepare_structures.py")
    config = load_config(args.config)

    # =====================
    # Step 1. 检查或下载 3GHG
    # =====================
    pdb_id = config["structures"]["primary_pdb"]
    raw_pdb = resolve_path(config["structures"]["primary_pdb_file"])
    raw_pdb.parent.mkdir(parents=True, exist_ok=True)
    if not raw_pdb.exists():
        url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        logger.warning("未找到 %s，尝试从 RCSB 下载: %s", raw_pdb, url)
        try:
            with urllib.request.urlopen(url, timeout=60) as response:
                raw_pdb.write_bytes(response.read())
            logger.info("下载完成: %s", raw_pdb)
        except Exception as exc:
            raise RuntimeError(
                f"无法下载 {pdb_id}.pdb。请手动下载并放入 {raw_pdb} 后重试。原始错误: {exc}"
            ) from exc
    else:
        logger.info("使用已有 PDB: %s", raw_pdb)

    # =====================
    # Step 2. 清理 PDB
    # =====================
    clean_path = resolve_path(config["structures"]["cleaned_pdb_file"])
    kept_atoms = clean_pdb(raw_pdb, clean_path)
    if kept_atoms == 0:
        raise RuntimeError(f"清理后的 PDB 不含蛋白 ATOM 记录: {clean_path}")
    logger.info("输出 clean PDB: %s", clean_path)
    logger.info("保留蛋白原子数: %s", kept_atoms)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
