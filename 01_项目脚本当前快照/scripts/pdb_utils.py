from __future__ import annotations

import math
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, Iterable, List, Sequence, Tuple

from common import resolve_path


AA3_TO_1 = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
    "SEC": "U",
    "PYL": "O",
    "MSE": "M",
}


def is_protein_residue(resname: str) -> bool:
    return resname.strip().upper() in AA3_TO_1


def clean_pdb(input_pdb: str | Path, output_pdb: str | Path) -> int:
    src = resolve_path(input_pdb)
    out = resolve_path(output_pdb)
    out.parent.mkdir(parents=True, exist_ok=True)
    if not src.exists():
        raise FileNotFoundError(f"缺少 PDB 文件: {src}")
    kept = 0
    saw_end = False
    with src.open("r", encoding="utf-8", errors="replace") as inp, out.open("w", encoding="utf-8", newline="\n") as handle:
        for line in inp:
            if line.startswith("END"):
                saw_end = True
            if line.startswith(("ATOM  ", "TER", "END")):
                if line.startswith("ATOM  "):
                    resname = line[17:20].strip()
                    if not is_protein_residue(resname):
                        continue
                    kept += 1
                handle.write(line.rstrip() + "\n")
        if kept and not saw_end:
            handle.write("END\n")
    return kept


def parse_residues(pdb_path: str | Path) -> Dict[str, List[Dict[str, Any]]]:
    path = resolve_path(pdb_path)
    chains: Dict[str, Dict[Tuple[str, str], Dict[str, Any]]] = defaultdict(dict)
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith("ATOM  "):
                continue
            atom = line[12:16].strip()
            resname = line[17:20].strip().upper()
            chain_id = (line[21].strip() or "_")
            resseq = line[22:26].strip()
            icode = line[26].strip()
            if not is_protein_residue(resname):
                continue
            key = (resseq, icode)
            residue = chains[chain_id].setdefault(
                key,
                {
                    "chain_id": chain_id,
                    "pdb_residue_number": resseq + icode,
                    "pdb_resseq": resseq,
                    "pdb_icode": icode,
                    "pdb_residue_name": resname,
                    "aa": AA3_TO_1.get(resname, "X"),
                    "atoms": {},
                },
            )
            try:
                residue["atoms"][atom] = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            except ValueError:
                continue
    return {chain: list(res_map.values()) for chain, res_map in chains.items()}


def residue_sequence(residues: Sequence[Mapping[str, Any]]) -> str:
    return "".join(str(res.get("aa", "X")) for res in residues)


def smith_waterman(query: str, target: str) -> Tuple[int, float, List[Tuple[int, int]]]:
    # query 是 PDB 可见序列，target 是 UniProt full-length；返回匹配 residue 下标对。
    m, n = len(query), len(target)
    prev = [0] * (n + 1)
    matrix: List[List[int]] = [[0] * (n + 1) for _ in range(m + 1)]
    pointer: List[List[int]] = [[0] * (n + 1) for _ in range(m + 1)]
    best_score = 0
    best_pos = (0, 0)
    for i in range(1, m + 1):
        current = [0] * (n + 1)
        for j in range(1, n + 1):
            match = matrix[i - 1][j - 1] + (3 if query[i - 1] == target[j - 1] else -2)
            delete = matrix[i - 1][j] - 2
            insert = current[j - 1] - 2
            score = max(0, match, delete, insert)
            current[j] = score
            matrix[i][j] = score
            if score == 0:
                pointer[i][j] = 0
            elif score == match:
                pointer[i][j] = 1
            elif score == delete:
                pointer[i][j] = 2
            else:
                pointer[i][j] = 3
            if score > best_score:
                best_score = score
                best_pos = (i, j)
        prev = current

    i, j = best_pos
    pairs: List[Tuple[int, int]] = []
    matches = 0
    aligned = 0
    while i > 0 and j > 0 and matrix[i][j] > 0:
        move = pointer[i][j]
        if move == 1:
            pairs.append((i - 1, j - 1))
            aligned += 1
            if query[i - 1] == target[j - 1]:
                matches += 1
            i -= 1
            j -= 1
        elif move == 2:
            i -= 1
        elif move == 3:
            j -= 1
        else:
            break
    pairs.reverse()
    identity = matches / aligned if aligned else 0.0
    return best_score, identity, pairs


def ca_coord(residue: Mapping[str, Any]) -> Tuple[float, float, float] | None:
    atoms = residue.get("atoms", {})
    if "CA" in atoms:
        return atoms["CA"]
    if atoms:
        return next(iter(atoms.values()))
    return None


def distance(a: Tuple[float, float, float], b: Tuple[float, float, float]) -> float:
    return math.sqrt(sum((a[i] - b[i]) ** 2 for i in range(3)))


def centroid(coords: Iterable[Tuple[float, float, float]]) -> Tuple[float, float, float]:
    points = list(coords)
    if not points:
        return (0.0, 0.0, 0.0)
    return tuple(sum(p[i] for p in points) / len(points) for i in range(3))  # type: ignore[return-value]


def write_visible_pdb(input_pdb: str | Path, output_pdb: str | Path, chain_ids: Sequence[str]) -> int:
    src = resolve_path(input_pdb)
    out = resolve_path(output_pdb)
    out.parent.mkdir(parents=True, exist_ok=True)
    chain_set = set(chain_ids)
    kept = 0
    with src.open("r", encoding="utf-8", errors="replace") as inp, out.open("w", encoding="utf-8", newline="\n") as handle:
        for line in inp:
            if line.startswith("ATOM  ") and (line[21].strip() or "_") in chain_set:
                handle.write(line.rstrip() + "\n")
                kept += 1
        handle.write("END\n")
    return kept


def write_patch_pdb(input_pdb: str | Path, output_pdb: str | Path, chain_residues: Dict[str, set[str]]) -> int:
    src = resolve_path(input_pdb)
    out = resolve_path(output_pdb)
    out.parent.mkdir(parents=True, exist_ok=True)
    kept = 0
    with src.open("r", encoding="utf-8", errors="replace") as inp, out.open("w", encoding="utf-8", newline="\n") as handle:
        for line in inp:
            if not line.startswith("ATOM  "):
                continue
            chain_id = line[21].strip() or "_"
            resid = line[22:26].strip() + line[26].strip()
            if resid in chain_residues.get(chain_id, set()):
                handle.write(line.rstrip() + "\n")
                kept += 1
        handle.write("END\n")
    return kept
