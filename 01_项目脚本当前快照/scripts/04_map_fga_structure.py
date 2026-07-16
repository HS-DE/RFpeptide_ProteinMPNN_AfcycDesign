from __future__ import annotations

import argparse

from common import append_run_header, load_config, read_fasta_sequence, setup_logger, write_csv
from pdb_utils import parse_residues, residue_sequence, smith_waterman, write_visible_pdb


def main() -> int:
    parser = argparse.ArgumentParser(description="识别 3GHG 中的 FGA 链并建立 PDB-UniProt 映射。")
    parser.add_argument("--config", default="config/project.yaml")
    args = parser.parse_args()

    logger = setup_logger("04_map_fga_structure")
    append_run_header(logger, "04_map_fga_structure.py")
    config = load_config(args.config)

    # =====================
    # Step 1. 读取序列和结构
    # =====================
    fga_sequence = read_fasta_sequence("data/input/FGA_full_length_1_866.fasta")
    clean_pdb = config["structures"]["cleaned_pdb_file"]
    chains = parse_residues(clean_pdb)
    if not chains:
        raise RuntimeError(f"clean PDB 中没有可解析蛋白链: {clean_pdb}")
    logger.info("解析到 PDB 链: %s", ", ".join(sorted(chains)))

    # =====================
    # Step 2. 对每条 PDB 链做局部序列比对
    # =====================
    selected = {}
    for chain_id, residues in chains.items():
        seq = residue_sequence(residues)
        score, identity, pairs = smith_waterman(seq, fga_sequence)
        aligned = len(pairs)
        logger.info("链 %s: aligned=%s identity=%.3f score=%s", chain_id, aligned, identity, score)
        if aligned >= 80 and identity >= 0.80:
            selected[chain_id] = {"residues": residues, "identity": identity, "pairs": pairs}
    if not selected:
        raise RuntimeError("无法可靠识别 3GHG 中的 FGA 链；请检查 PDB 或补充链注释。")

    # =====================
    # Step 3. 输出 PDB residue 到 UniProt residue 的映射
    # =====================
    rows = []
    for chain_id, info in selected.items():
        residues = info["residues"]
        confidence = "high" if info["identity"] >= 0.95 else "medium"
        for pdb_idx, uni_idx in info["pairs"]:
            res = residues[pdb_idx]
            rows.append(
                {
                    "pdb_id": config["structures"]["primary_pdb"],
                    "chain_id": chain_id,
                    "pdb_residue_number": res["pdb_residue_number"],
                    "pdb_residue_name": res["pdb_residue_name"],
                    "uniprot_id": config["project"]["target_uniprot"],
                    "uniprot_residue_number": uni_idx + 1,
                    "uniprot_residue_name": fga_sequence[uni_idx],
                    "mapping_confidence": confidence,
                    "is_visible": True,
                }
            )
    fields = [
        "pdb_id",
        "chain_id",
        "pdb_residue_number",
        "pdb_residue_name",
        "uniprot_id",
        "uniprot_residue_number",
        "uniprot_residue_name",
        "mapping_confidence",
        "is_visible",
    ]
    write_csv("data/annotations/FGA_structure_mapping.csv", rows, fields)
    logger.info("输出映射行数: %s", len(rows))

    # =====================
    # Step 4. 输出可见 FGA 链 PDB
    # =====================
    kept_atoms = write_visible_pdb(clean_pdb, "data/structures/prepared/FGA_visible_regions.pdb", list(selected))
    logger.info("输出 FGA_visible_regions.pdb，原子数: %s", kept_atoms)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
