from __future__ import annotations

import argparse
import json
import math
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Sequence, Tuple

from common import append_run_header, load_config, project_root, read_csv, resolve_path, setup_logger, write_csv
from ranking import normalized_complex_score


SEED_FIELDS = [
    "boltz_job_id",
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
    "confidence_json",
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


def _to_float(value: Any, default: float = 0.0) -> float:
    try:
        if value in ("", None):
            return default
        return float(value)
    except (TypeError, ValueError):
        return default


def _mean(values: Iterable[float], default: float = 0.0) -> float:
    vals = [v for v in values if not math.isnan(v)]
    if not vals:
        return default
    return sum(vals) / len(vals)


def _first_file(root: Path, patterns: Sequence[str]) -> Path | None:
    for pattern in patterns:
        try:
            hits = sorted(root.rglob(pattern))
        except OSError:
            continue
        if hits:
            return hits[0]
    return None


def _load_json(path: Path | None) -> Dict[str, Any]:
    if path is None or not path.exists():
        return {}
    try:
        parsed = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}
    return parsed if isinstance(parsed, dict) else {}


def _scale_plddt(value: float) -> float:
    if value <= 1.0:
        return value * 100.0
    return value


def _npz_array(path: Path | None, ndim: int) -> Any:
    if path is None or not path.exists():
        return None
    try:
        import numpy as np  # type: ignore

        data = np.load(path)
        for key in data.files:
            arr = data[key]
            if arr.ndim == ndim:
                return arr
    except Exception:
        return None
    return None


def _interface_pae(path: Path | None, target_len: int, peptide_len: int) -> float:
    arr = _npz_array(path, 2)
    if arr is None:
        return float("nan")
    total = target_len + peptide_len
    if arr.shape[0] < total or arr.shape[1] < total:
        return float("nan")
    block = arr[:target_len, target_len:total]
    return float(block.mean())


def _peptide_plddt(path: Path | None, target_len: int, peptide_len: int, fallback: float) -> float:
    arr = _npz_array(path, 1)
    total = target_len + peptide_len
    if arr is not None and arr.shape[0] >= total:
        return _scale_plddt(float(arr[target_len:total].mean()))
    return _scale_plddt(fallback)


def _parse_hotspot_positions(value: str) -> List[int]:
    out: List[int] = []
    for part in str(value).split(","):
        part = part.strip()
        if not part:
            continue
        try:
            out.append(int(part))
        except ValueError:
            continue
    return sorted(set(out))


def _parse_pdb_atoms(path: Path | None) -> Dict[str, Dict[int, List[Tuple[str, Tuple[float, float, float]]]]]:
    atoms: Dict[str, Dict[int, List[Tuple[str, Tuple[float, float, float]]]]] = defaultdict(lambda: defaultdict(list))
    if path is None or not path.exists() or path.suffix.lower() != ".pdb":
        return atoms
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            atom_name = line[12:16].strip()
            if atom_name.startswith("H"):
                continue
            chain_id = line[21].strip()
            if not chain_id:
                continue
            try:
                res_id = int(line[22:26])
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                continue
            atoms[chain_id][res_id].append((atom_name, (x, y, z)))
    return atoms


def _dist(a: Tuple[float, float, float], b: Tuple[float, float, float]) -> float:
    return math.sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2)


def _contact_stats(
    structure: Path | None,
    target_chain: str,
    peptide_chain: str,
    peptide_len: int,
    hotspot_positions: Sequence[int],
    cutoff: float,
) -> Tuple[int, int, float, str, str]:
    atoms = _parse_pdb_atoms(structure)
    target = atoms.get(target_chain, {})
    peptide = atoms.get(peptide_chain, {})
    if not target or not peptide:
        return 0, 0, float("nan"), "missing_structure_or_chain", "missing"

    contact_pairs = set()
    target_contacts = set()
    for t_res, t_atoms in target.items():
        for p_res, p_atoms in peptide.items():
            hit = False
            for _, t_xyz in t_atoms:
                for _, p_xyz in p_atoms:
                    if _dist(t_xyz, p_xyz) <= cutoff:
                        hit = True
                        break
                if hit:
                    break
            if hit:
                contact_pairs.add((t_res, p_res))
                target_contacts.add(t_res)

    hotspot_set = set(hotspot_positions)
    patch_contact_count = len(target_contacts & hotspot_set) if hotspot_set else 0
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
    parser = argparse.ArgumentParser(description="Parse Boltz prediction outputs into FGA complex prediction score tables.")
    parser.add_argument("--config", default="config/project.yaml")
    parser.add_argument("--jobs-csv", default="results/boltz_predictions/boltz_jobs.csv")
    parser.add_argument("--output-root", default="results/boltz_predictions")
    parser.add_argument("--target-length", type=int, default=None)
    parser.add_argument("--target-fasta", default="data/input/FGA_chain_36_866.fasta")
    parser.add_argument("--contact-cutoff", type=float, default=5.0)
    args = parser.parse_args()

    logger = setup_logger("15_parse_boltz_predictions")
    append_run_header(logger, "15_parse_boltz_predictions.py")
    config = load_config(args.config)
    thresholds = config["scoring_thresholds"]

    if args.target_length is None:
        from common import read_fasta_sequence

        target_len = len(read_fasta_sequence(args.target_fasta))
    else:
        target_len = args.target_length

    jobs = read_csv(args.jobs_csv)
    output_root = resolve_path(args.output_root)
    seed_rows: List[Dict[str, Any]] = []

    for job in jobs:
        job_id = job.get("boltz_job_id", "")
        out_dir = resolve_path(job.get("output_dir", ""))
        confidence_path = _first_file(out_dir, ["confidence_*.json", "*confidence*.json"])
        structure_path = _first_file(out_dir, ["*_model_0.pdb", "*.pdb", "*_model_0.cif", "*.cif"])
        pae_path = _first_file(out_dir, ["pae_*_model_0.npz", "*pae*.npz"])
        plddt_path = _first_file(out_dir, ["plddt_*_model_0.npz", "*plddt*.npz"])
        conf = _load_json(confidence_path)

        core = job.get("core_sequence", "")
        peptide_len = len(core)
        iptm = _to_float(conf.get("iptm"), float("nan"))
        ptm = _to_float(conf.get("ptm"), float("nan"))
        complex_plddt = _scale_plddt(_to_float(conf.get("complex_plddt"), float("nan")))
        interface_pae = _interface_pae(pae_path, target_len, peptide_len)
        peptide_plddt = _peptide_plddt(plddt_path, target_len, peptide_len, complex_plddt)
        hotspot_positions = _parse_hotspot_positions(job.get("hotspot_target_positions", ""))
        contacts, patch_contacts, cys_dist, patch_flag, cys_flag = _contact_stats(
            structure_path,
            job.get("target_chain_id", "A"),
            job.get("peptide_chain_id", "B"),
            peptide_len,
            hotspot_positions,
            args.contact_cutoff,
        )

        notes = []
        if not conf:
            notes.append("missing_confidence_json")
        if math.isnan(interface_pae):
            notes.append("missing_or_unusable_pae")
        if structure_path is None:
            notes.append("missing_structure")
        if structure_path is not None and structure_path.suffix.lower() != ".pdb":
            notes.append("structure_contact_checks_require_pdb")

        seed_rows.append(
            {
                "boltz_job_id": job_id,
                "peptide_id": job.get("peptide_id", ""),
                "core_sequence": core,
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
                "confidence_json": "" if confidence_path is None else str(confidence_path.relative_to(project_root())),
                "structure_file": "" if structure_path is None else str(structure_path.relative_to(project_root())),
                "notes": ";".join(notes) if notes else "parsed_real_boltz_output",
            }
        )

    write_csv(output_root / "FGA_boltz_seed_scores.csv", seed_rows, SEED_FIELDS)

    grouped: Dict[str, List[Dict[str, Any]]] = defaultdict(list)
    for row in seed_rows:
        if row.get("confidence_json"):
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
            "notes": "aggregated_real_boltz_outputs",
        }
        summary["complex_score"] = round(normalized_complex_score(summary), 4)
        summary["complex_score_pass"] = _score_pass(summary, thresholds)
        summary_rows.append(summary)

    write_csv(output_root / "FGA_boltz_complex_prediction_summary.csv", summary_rows, SUMMARY_FIELDS)
    logger.info("Parsed Boltz seed rows: %s", len(seed_rows))
    logger.info("Wrote Boltz peptide summaries: %s", len(summary_rows))
    if not summary_rows:
        logger.warning("No complete Boltz confidence outputs were found; summary table is empty.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
