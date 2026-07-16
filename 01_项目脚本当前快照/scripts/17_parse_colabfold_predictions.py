from __future__ import annotations

import argparse
import json
import math
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Tuple

from common import append_run_header, load_config, project_root, read_csv, read_fasta_sequence, resolve_path, setup_logger, write_csv
from ranking import normalized_complex_score


SEED_FIELDS = [
    "colabfold_job_id",
    "peptide_id",
    "core_sequence",
    "patch_id",
    "seed",
    "iptm",
    "ptm",
    "complex_plddt",
    "interface_pae",
    "peptide_plddt",
    "interface_contacts",
    "patch_contact_count",
    "patch_consistency_flag",
    "cys_sg_distance",
    "cys_cys_geometry",
    "scores_json",
    "structure_file",
    "notes",
]

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
}


def _to_float(value: Any, default: float = float("nan")) -> float:
    try:
        if value in ("", None):
            return default
        return float(value)
    except (TypeError, ValueError):
        return default


def _mean(values: Iterable[float], default: float = float("nan")) -> float:
    vals = [v for v in values if not math.isnan(v)]
    if not vals:
        return default
    return sum(vals) / len(vals)


def _first_file(root: Path, patterns: List[str]) -> Path | None:
    if not root.exists():
        return None
    for pattern in patterns:
        matches = sorted(root.glob(pattern))
        if matches:
            return matches[0]
    for pattern in patterns:
        matches = sorted(root.rglob(pattern))
        if matches:
            return matches[0]
    return None


def _load_json(path: Path | None) -> Dict[str, Any]:
    if path is None or not path.exists():
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {}


def _scale_plddt(value: float) -> float:
    if math.isnan(value):
        return value
    if value <= 1.0:
        return value * 100.0
    return value


def _matrix_from_json(conf: Mapping[str, Any], keys: Iterable[str]) -> Any:
    for key in keys:
        if key in conf:
            return conf[key]
    return None


def _interface_pae(conf: Mapping[str, Any], target_len: int, peptide_len: int) -> float:
    matrix = _matrix_from_json(conf, ["pae", "predicted_aligned_error"])
    if matrix is None:
        return float("nan")
    try:
        total = target_len + peptide_len
        block = [row[target_len:total] for row in matrix[:target_len]]
        vals = [float(v) for row in block for v in row]
        return _mean(vals)
    except Exception:
        return float("nan")


def _peptide_plddt(conf: Mapping[str, Any], target_len: int, peptide_len: int, fallback: float) -> float:
    arr = _matrix_from_json(conf, ["plddt", "predicted_lddt"])
    if arr is None:
        return _scale_plddt(fallback)
    try:
        total = target_len + peptide_len
        return _scale_plddt(_mean(float(v) for v in arr[target_len:total]))
    except Exception:
        return _scale_plddt(fallback)


def _parse_hotspot_positions(value: str) -> List[int]:
    out: List[int] = []
    for part in str(value).replace(";", ",").split(","):
        part = part.strip()
        if not part:
            continue
        try:
            out.append(int(part))
        except ValueError:
            digits = "".join(ch for ch in part if ch.isdigit())
            if digits:
                out.append(int(digits))
    return out


def _parse_pdb_atoms(path: Path | None) -> Dict[str, Dict[int, List[Tuple[str, Tuple[float, float, float]]]]]:
    chains: Dict[str, Dict[int, List[Tuple[str, Tuple[float, float, float]]]]] = defaultdict(lambda: defaultdict(list))
    if path is None or not path.exists() or path.suffix.lower() != ".pdb":
        return chains
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            atom = line[12:16].strip()
            if atom.startswith("H"):
                continue
            resname = line[17:20].strip()
            if resname not in AA3_TO_1:
                continue
            chain = line[21].strip() or "A"
            try:
                resseq = int(line[22:26])
                xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            except ValueError:
                continue
            chains[chain][resseq].append((atom, xyz))
    return chains


def _dist(a: Tuple[float, float, float], b: Tuple[float, float, float]) -> float:
    return math.sqrt(sum((a[i] - b[i]) ** 2 for i in range(3)))


def _contact_stats(
    path: Path | None,
    peptide_len: int,
    hotspot_positions: List[int],
    cutoff: float,
) -> Tuple[int, int, float, str, str]:
    chains = _parse_pdb_atoms(path)
    target = chains.get("A", {})
    peptide = chains.get("B", {})
    if not target or not peptide:
        return 0, 0, float("nan"), "fail", "missing"

    contact_pairs = set()
    patch_contact_count = 0
    hotspot_set = set(hotspot_positions)
    for tres, tatoms in target.items():
        for pres, patoms in peptide.items():
            hit = False
            for _, txyz in tatoms:
                for _, pxyz in patoms:
                    if _dist(txyz, pxyz) <= cutoff:
                        hit = True
                        break
                if hit:
                    break
            if hit:
                contact_pairs.add((tres, pres))
                if tres in hotspot_set:
                    patch_contact_count += 1

    patch_flag = "pass" if patch_contact_count > 0 else "fail"
    sg_first = None
    sg_last = None
    for atom_name, xyz in peptide.get(1, []):
        if atom_name == "SG":
            sg_first = xyz
            break
    for atom_name, xyz in peptide.get(peptide_len, []):
        if atom_name == "SG":
            sg_last = xyz
            break
    if sg_first is None or sg_last is None:
        cys_dist = float("nan")
        cys_flag = "missing"
    else:
        cys_dist = _dist(sg_first, sg_last)
        cys_flag = "pass" if 1.8 <= cys_dist <= 2.4 else "fail"
    return len(contact_pairs), patch_contact_count, cys_dist, patch_flag, cys_flag


def _score_pass(row: Mapping[str, Any], thresholds: Mapping[str, Any]) -> bool:
    return (
        _to_float(row.get("best_interface_pae"), 999.0) <= float(thresholds["max_interface_pae"])
        and _to_float(row.get("mean_peptide_plddt"), 0.0) >= float(thresholds["min_peptide_plddt"])
        and _to_float(row.get("interface_contacts"), 0.0) >= float(thresholds["min_interface_contacts"])
        and str(row.get("patch_consistency_flag", "")).lower() == "pass"
        and str(row.get("cys_cys_geometry", "")).lower() == "pass"
    )


def main() -> int:
    parser = argparse.ArgumentParser(description="Parse ColabFold/AlphaFold-Multimer outputs into FGA score tables.")
    parser.add_argument("--config", default="config/project.yaml")
    parser.add_argument("--jobs-csv", default="results/colabfold_predictions_top30_seed1/colabfold_jobs.csv")
    parser.add_argument("--output-root", default="results/colabfold_predictions_top30_seed1")
    parser.add_argument("--target-fasta", default="data/input/FGA_chain_36_866.fasta")
    parser.add_argument("--target-length", type=int, default=None)
    parser.add_argument("--contact-cutoff", type=float, default=5.0)
    args = parser.parse_args()

    logger = setup_logger("17_parse_colabfold_predictions")
    append_run_header(logger, "17_parse_colabfold_predictions.py")
    config = load_config(args.config)
    thresholds = config["scoring_thresholds"]
    target_len = args.target_length or len(read_fasta_sequence(args.target_fasta))
    output_root = resolve_path(args.output_root)

    seed_rows: List[Dict[str, Any]] = []
    for job in read_csv(args.jobs_csv):
        out_dir = resolve_path(job.get("output_dir", ""))
        scores_path = _first_file(out_dir, ["*_scores_rank_001*.json", "*scores*.json"])
        structure_path = _first_file(out_dir, ["*_unrelaxed_rank_001*.pdb", "*rank_001*.pdb", "*.pdb"])
        conf = _load_json(scores_path)
        peptide_len = int(job.get("peptide_length") or len(job.get("core_sequence", "")))
        iptm = _to_float(conf.get("iptm"), float("nan"))
        ptm = _to_float(conf.get("ptm"), float("nan"))
        complex_plddt = _scale_plddt(_mean(float(v) for v in conf.get("plddt", []) if str(v) != "")) if conf.get("plddt") else float("nan")
        interface_pae = _interface_pae(conf, target_len, peptide_len)
        peptide_plddt = _peptide_plddt(conf, target_len, peptide_len, complex_plddt)
        contacts, patch_contacts, cys_dist, patch_flag, cys_flag = _contact_stats(
            structure_path,
            peptide_len,
            _parse_hotspot_positions(job.get("hotspot_target_positions", "")),
            args.contact_cutoff,
        )
        notes = []
        if not conf:
            notes.append("missing_scores_json")
        if math.isnan(interface_pae):
            notes.append("missing_or_unusable_pae")
        if structure_path is None:
            notes.append("missing_structure")
        if not notes:
            notes.append("parsed_real_colabfold_output")
        seed_rows.append(
            {
                "colabfold_job_id": job.get("colabfold_job_id", ""),
                "peptide_id": job.get("peptide_id", ""),
                "core_sequence": job.get("core_sequence", ""),
                "patch_id": job.get("patch_id", ""),
                "seed": job.get("seed", ""),
                "iptm": "" if math.isnan(iptm) else round(iptm, 4),
                "ptm": "" if math.isnan(ptm) else round(ptm, 4),
                "complex_plddt": "" if math.isnan(complex_plddt) else round(complex_plddt, 2),
                "interface_pae": "" if math.isnan(interface_pae) else round(interface_pae, 3),
                "peptide_plddt": "" if math.isnan(peptide_plddt) else round(peptide_plddt, 2),
                "interface_contacts": contacts,
                "patch_contact_count": patch_contacts,
                "patch_consistency_flag": patch_flag,
                "cys_sg_distance": "" if math.isnan(cys_dist) else round(cys_dist, 3),
                "cys_cys_geometry": cys_flag,
                "scores_json": "" if scores_path is None else str(scores_path.relative_to(project_root())),
                "structure_file": "" if structure_path is None else str(structure_path.relative_to(project_root())),
                "notes": ";".join(notes),
            }
        )

    write_csv(output_root / "FGA_colabfold_seed_scores.csv", seed_rows, SEED_FIELDS)

    grouped: Dict[str, List[Dict[str, Any]]] = defaultdict(list)
    for row in seed_rows:
        if row.get("scores_json"):
            grouped[str(row.get("peptide_id", ""))].append(row)

    summary_rows: List[Dict[str, Any]] = []
    for peptide_id, rows in grouped.items():
        best = max(rows, key=lambda r: _to_float(r.get("iptm"), -1.0))
        summary = {
            "peptide_id": peptide_id,
            "core_sequence": best.get("core_sequence", ""),
            "patch_id": best.get("patch_id", ""),
            "n_seeds": len(rows),
            "best_seed": best.get("seed", ""),
            "mean_iptm": round(_mean(_to_float(r.get("iptm"), float("nan")) for r in rows), 4),
            "best_iptm": best.get("iptm", ""),
            "mean_interface_pae": round(_mean((_to_float(r.get("interface_pae"), float("nan")) for r in rows), default=999.0), 3),
            "best_interface_pae": best.get("interface_pae", ""),
            "mean_peptide_plddt": round(_mean(_to_float(r.get("peptide_plddt"), float("nan")) for r in rows), 2),
            "interface_contacts": best.get("interface_contacts", ""),
            "pose_consistency_rmsd": "not_computed",
            "patch_consistency_flag": "pass" if any(str(r.get("patch_consistency_flag", "")).lower() == "pass" for r in rows) else "fail",
            "cys_cys_geometry": "pass" if any(str(r.get("cys_cys_geometry", "")).lower() == "pass" for r in rows) else "fail",
            "notes": "aggregated_real_colabfold_outputs",
        }
        summary["complex_score"] = round(normalized_complex_score(summary), 4)
        summary["complex_score_pass"] = _score_pass(summary, thresholds)
        summary_rows.append(summary)

    write_csv(output_root / "FGA_colabfold_complex_prediction_summary.csv", summary_rows, SUMMARY_FIELDS)
    logger.info("Parsed ColabFold seed rows: %s", len(seed_rows))
    logger.info("Wrote ColabFold peptide summaries: %s", len(summary_rows))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
