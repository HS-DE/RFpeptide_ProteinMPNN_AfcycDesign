from __future__ import annotations

import argparse
import math
import os
from collections import defaultdict
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np

from common import append_run_header, read_csv, resolve_path, rows_to_markdown, setup_logger, write_csv, write_markdown
from pdb_utils import parse_residues, residue_sequence


MODEL_FIELDS = [
    "stage5_control_job_id",
    "control_group",
    "target_msa_mode",
    "msa_rows_input",
    "use_mlm",
    "mlm_replace_fraction",
    "seed",
    "model_name",
    "prediction_pdb",
    "target_sequence_match",
    "target_global_CA_RMSD_A",
    "Site_2_local_CA_RMSD_A",
    "hotspot_local_CA_RMSD_A",
    "local_rmsd_alignment_basis",
    "target_mean_pLDDT_fraction",
    "Site_2_mean_pLDDT_fraction",
    "hotspot_mean_pLDDT_fraction",
    "target_mean_pLDDT_100",
    "Site_2_mean_pLDDT_100",
    "hotspot_mean_pLDDT_100",
    "ptm",
    "target_recovery_pass",
    "target_recovery_failure_reasons",
]

GROUP_FIELDS = [
    "control_group",
    "target_msa_mode",
    "use_mlm",
    "mlm_replace_fraction",
    "seeds_planned",
    "seeds_completed",
    "models_completed",
    "models_passing_target_recovery",
    "median_target_global_CA_RMSD_A",
    "median_Site_2_local_CA_RMSD_A",
    "median_hotspot_local_CA_RMSD_A",
    "median_target_mean_pLDDT_100",
    "median_Site_2_mean_pLDDT_100",
    "median_hotspot_mean_pLDDT_100",
    "status",
    "notes",
]


def _resolve_mixed_path(value: str | Path) -> Path:
    text = str(value).strip().replace("\\", "/")
    if os.name == "nt" and text.startswith("/mnt/") and len(text) > 7 and text[6] == "/":
        return Path(f"{text[5].upper()}:/{text[7:]}")
    if os.name != "nt" and len(text) >= 3 and text[1] == ":" and text[2] == "/":
        return Path(f"/mnt/{text[0].lower()}{text[2:]}")
    path = Path(text)
    return path if path.is_absolute() else resolve_path(path)


def _float(value: Any, default: float = float("nan")) -> float:
    try:
        parsed = float(str(value).strip())
    except (TypeError, ValueError):
        return default
    return parsed if math.isfinite(parsed) else default


def _rounded(value: float, digits: int = 3) -> float | str:
    return round(value, digits) if math.isfinite(value) else ""


def _median(rows: Sequence[Mapping[str, Any]], field: str) -> float | str:
    values = [_float(row.get(field, "")) for row in rows]
    finite = [value for value in values if math.isfinite(value)]
    return round(float(np.median(finite)), 3) if finite else ""


def _ca_coords(residues: Sequence[Mapping[str, Any]]) -> np.ndarray:
    points = []
    for residue in residues:
        atoms = residue.get("atoms", {})
        if "CA" not in atoms:
            raise RuntimeError("A target residue lacks a CA atom")
        points.append(atoms["CA"])
    return np.asarray(points, dtype=float)


def _kabsch(mobile: np.ndarray, reference: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    if mobile.shape != reference.shape or mobile.shape[0] < 3:
        raise RuntimeError(f"Cannot align target CA arrays with shapes {mobile.shape} and {reference.shape}")
    mobile_center = mobile.mean(axis=0)
    reference_center = reference.mean(axis=0)
    covariance = (mobile - mobile_center).T @ (reference - reference_center)
    u_matrix, _, vt_matrix = np.linalg.svd(covariance)
    rotation = vt_matrix.T @ u_matrix.T
    if np.linalg.det(rotation) < 0:
        vt_matrix[-1, :] *= -1
        rotation = vt_matrix.T @ u_matrix.T
    translation = reference_center - mobile_center @ rotation.T
    aligned = mobile @ rotation.T + translation
    rmsd = float(np.sqrt(np.mean(np.sum((aligned - reference) ** 2, axis=1))))
    return rotation, translation, rmsd


def _subset_rmsd(
    mobile: np.ndarray,
    reference: np.ndarray,
    indices_1based: Sequence[int],
    rotation: np.ndarray,
    translation: np.ndarray,
) -> float:
    indices = np.asarray([int(index) - 1 for index in indices_1based], dtype=int)
    aligned = mobile[indices] @ rotation.T + translation
    delta = aligned - reference[indices]
    return float(np.sqrt(np.mean(np.sum(delta * delta, axis=1))))


def _mapping_indices(mapping_csv: Path, expected_length: int) -> tuple[list[int], list[int]]:
    rows = read_csv(mapping_csv)
    if len(rows) != expected_length:
        raise RuntimeError(f"Stage 0 mapping has {len(rows)} rows, expected {expected_length}")
    site = [
        int(row["rfpeptides_residue_number"])
        for row in rows
        if str(row.get("is_target_site_residue", "")).strip().lower() == "true"
    ]
    hotspots = [
        int(row["rfpeptides_residue_number"])
        for row in rows
        if str(row.get("is_selected_hotspot", "")).strip().lower() == "true"
    ]
    if not site or not hotspots:
        raise RuntimeError("Stage 0 mapping lacks Site_2 or hotspot indices")
    return site, hotspots


def _model_row(
    *,
    job: Mapping[str, str],
    metric: Mapping[str, str],
    prediction_pdb: Path,
    reference_residues: Sequence[Mapping[str, Any]],
    reference_sequence: str,
    site_indices: Sequence[int],
    hotspot_indices: Sequence[int],
    args: argparse.Namespace,
) -> dict[str, Any]:
    prediction_chains = parse_residues(prediction_pdb)
    if "A" not in prediction_chains:
        raise RuntimeError(f"Predicted target chain A is missing: {prediction_pdb}")
    predicted_residues = prediction_chains["A"]
    sequence_match = residue_sequence(predicted_residues) == reference_sequence
    predicted_ca = _ca_coords(predicted_residues)
    reference_ca = _ca_coords(reference_residues)
    rotation, translation, global_rmsd = _kabsch(predicted_ca, reference_ca)
    site_rmsd = _subset_rmsd(predicted_ca, reference_ca, site_indices, rotation, translation)
    hotspot_rmsd = _subset_rmsd(predicted_ca, reference_ca, hotspot_indices, rotation, translation)

    target_plddt_fraction = _float(metric.get("plddt_mean_fraction", ""))
    site_plddt_fraction = _float(metric.get("plddt_site2_mean_fraction", ""))
    hotspot_plddt_fraction = _float(metric.get("plddt_hotspot_mean_fraction", ""))
    target_plddt_100 = _float(metric.get("plddt_mean_100", ""))
    site_plddt_100 = _float(metric.get("plddt_site2_mean_100", ""))
    hotspot_plddt_100 = _float(metric.get("plddt_hotspot_mean_100", ""))
    failure_reasons: list[str] = []
    if not sequence_match:
        failure_reasons.append("target_sequence_mismatch")
    if global_rmsd > args.max_target_global_rmsd:
        failure_reasons.append("high_target_global_CA_RMSD")
    if site_rmsd > args.max_site2_local_rmsd:
        failure_reasons.append("high_Site_2_local_CA_RMSD")
    if hotspot_rmsd > args.max_hotspot_local_rmsd:
        failure_reasons.append("high_hotspot_local_CA_RMSD")
    if not math.isfinite(target_plddt_100) or target_plddt_100 < args.min_target_plddt:
        failure_reasons.append("low_target_mean_pLDDT")
    if not math.isfinite(site_plddt_100) or site_plddt_100 < args.min_site2_plddt:
        failure_reasons.append("low_Site_2_mean_pLDDT")

    return {
        "stage5_control_job_id": job["stage5_control_job_id"],
        "control_group": job["control_group"],
        "target_msa_mode": job["target_msa_mode"],
        "msa_rows_input": metric.get("msa_rows_input", job.get("msa_rows_expected", "")),
        "use_mlm": job["use_mlm"],
        "mlm_replace_fraction": job["mlm_replace_fraction"],
        "seed": job["seed"],
        "model_name": metric.get("model_name", prediction_pdb.stem),
        "prediction_pdb": prediction_pdb,
        "target_sequence_match": str(sequence_match).lower(),
        "target_global_CA_RMSD_A": _rounded(global_rmsd),
        "Site_2_local_CA_RMSD_A": _rounded(site_rmsd),
        "hotspot_local_CA_RMSD_A": _rounded(hotspot_rmsd),
        "local_rmsd_alignment_basis": "after_global_target_CA_alignment",
        "target_mean_pLDDT_fraction": _rounded(target_plddt_fraction, 5),
        "Site_2_mean_pLDDT_fraction": _rounded(site_plddt_fraction, 5),
        "hotspot_mean_pLDDT_fraction": _rounded(hotspot_plddt_fraction, 5),
        "target_mean_pLDDT_100": _rounded(target_plddt_100),
        "Site_2_mean_pLDDT_100": _rounded(site_plddt_100),
        "hotspot_mean_pLDDT_100": _rounded(hotspot_plddt_100),
        "ptm": metric.get("ptm", ""),
        "target_recovery_pass": str(not failure_reasons).lower(),
        "target_recovery_failure_reasons": ";".join(failure_reasons),
    }


def _group_rows(
    controls: Sequence[Mapping[str, str]],
    model_rows: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    grouped: dict[str, list[Mapping[str, Any]]] = defaultdict(list)
    for row in model_rows:
        grouped[str(row["control_group"])].append(row)
    summaries: list[dict[str, Any]] = []
    for control in controls:
        rows = grouped.get(str(control["control_group"]), [])
        seeds = {str(row.get("seed", "")) for row in rows}
        pass_count = sum(row.get("target_recovery_pass") == "true" for row in rows)
        if str(control.get("status", "")).startswith("blocked"):
            status = str(control["status"])
        elif not rows:
            status = "prepared_not_run"
        elif pass_count:
            status = "has_target_recovery_models"
        else:
            status = "target_not_recovered_under_thresholds"
        summaries.append(
            {
                "control_group": control["control_group"],
                "target_msa_mode": control["target_msa_mode"],
                "use_mlm": control["use_mlm"],
                "mlm_replace_fraction": control["mlm_replace_fraction"],
                "seeds_planned": control["seeds_planned"],
                "seeds_completed": len(seeds),
                "models_completed": len(rows),
                "models_passing_target_recovery": pass_count,
                "median_target_global_CA_RMSD_A": _median(rows, "target_global_CA_RMSD_A"),
                "median_Site_2_local_CA_RMSD_A": _median(rows, "Site_2_local_CA_RMSD_A"),
                "median_hotspot_local_CA_RMSD_A": _median(rows, "hotspot_local_CA_RMSD_A"),
                "median_target_mean_pLDDT_100": _median(rows, "target_mean_pLDDT_100"),
                "median_Site_2_mean_pLDDT_100": _median(rows, "Site_2_mean_pLDDT_100"),
                "median_hotspot_mean_pLDDT_100": _median(rows, "hotspot_mean_pLDDT_100"),
                "status": status,
                "notes": "Target globally aligned by CA; Site_2/hotspot RMSD measured in that common frame.",
            }
        )
    return summaries


def _report(
    model_rows: Sequence[Mapping[str, Any]],
    groups: Sequence[Mapping[str, Any]],
    args: argparse.Namespace,
) -> str:
    group_columns = [
        "control_group",
        "target_msa_mode",
        "use_mlm",
        "seeds_completed",
        "models_completed",
        "models_passing_target_recovery",
        "median_target_global_CA_RMSD_A",
        "median_Site_2_local_CA_RMSD_A",
        "median_hotspot_local_CA_RMSD_A",
        "median_target_mean_pLDDT_100",
        "median_Site_2_mean_pLDDT_100",
        "status",
    ]
    model_columns = [
        "control_group",
        "seed",
        "model_name",
        "target_global_CA_RMSD_A",
        "Site_2_local_CA_RMSD_A",
        "hotspot_local_CA_RMSD_A",
        "target_mean_pLDDT_100",
        "Site_2_mean_pLDDT_100",
        "hotspot_mean_pLDDT_100",
        "target_recovery_pass",
        "target_recovery_failure_reasons",
    ]
    return f"""# Stage 5 Target-Only Control Results

These are sequence-based target recovery controls. No peptide, target template,
initial guess, or cyclic offset is used. pLDDT is reported on the 0-100 scale.

Local Site_2 and hotspot RMSDs are measured after a single global target CA
alignment, so they report local distortion in a common target reference frame.

## Descriptive Recovery Thresholds

```text
max_target_global_CA_RMSD_A: {args.max_target_global_rmsd:g}
max_Site_2_local_CA_RMSD_A: {args.max_site2_local_rmsd:g}
max_hotspot_local_CA_RMSD_A: {args.max_hotspot_local_rmsd:g}
min_target_mean_pLDDT_100: {args.min_target_plddt:g}
min_Site_2_mean_pLDDT_100: {args.min_site2_plddt:g}
```

These thresholds are control diagnostics, not peptide pass/fail criteria.

## Group Comparison

{rows_to_markdown(groups, group_columns, "No control groups were found.")}

## Model Results

{rows_to_markdown(model_rows, model_columns, "No target-control predictions have been run.")}
"""


def main() -> int:
    parser = argparse.ArgumentParser(description="Collect Stage 5 target-only recovery controls.")
    parser.add_argument(
        "--control-root",
        default="results/rfpeptides_article_route_clean_20260623_stage5_target_controls_v1/07_structure_validation_target_controls",
    )
    parser.add_argument("--stage0-root", default="results/rfpeptides_article_route_clean_20260615_fpocket")
    parser.add_argument("--max-target-global-rmsd", type=float, default=3.0)
    parser.add_argument("--max-site2-local-rmsd", type=float, default=2.0)
    parser.add_argument("--max-hotspot-local-rmsd", type=float, default=2.0)
    parser.add_argument("--min-target-plddt", type=float, default=70.0)
    parser.add_argument("--min-site2-plddt", type=float, default=70.0)
    args = parser.parse_args()

    logger = setup_logger("29_collect_stage5_target_controls")
    append_run_header(logger, "29_collect_stage5_target_controls.py")
    control_root = _resolve_mixed_path(args.control_root)
    stage0_root = _resolve_mixed_path(args.stage0_root)
    controls = read_csv(control_root / "FGA_rfpeptides_stage5_target_control_manifest.csv")
    jobs = read_csv(control_root / "inputs" / "FGA_rfpeptides_stage5_target_control_jobs.csv")
    if not controls or not jobs:
        raise RuntimeError(f"Missing target-control manifest or jobs CSV: {control_root}")

    reference_pdb = stage0_root / "00_target_inputs" / "RFpep_Site_2_target.pdb"
    reference_chains = parse_residues(reference_pdb)
    if "A" not in reference_chains:
        raise RuntimeError("Stage 0 reference target chain A is missing")
    reference_residues = reference_chains["A"]
    reference_sequence = residue_sequence(reference_residues)
    site_indices, hotspot_indices = _mapping_indices(
        stage0_root / "00_target_inputs" / "RFpep_Site_2_crop_renumbering_mapping.csv",
        len(reference_residues),
    )

    model_rows: list[dict[str, Any]] = []
    for job in jobs:
        prediction_dir = _resolve_mixed_path(job["prediction_output_dir"])
        metrics = read_csv(prediction_dir / "model_metrics.csv")
        for metric in metrics:
            prediction_pdb = _resolve_mixed_path(metric.get("prediction_pdb", ""))
            if not prediction_pdb.is_file():
                logger.warning("Skipping missing target-control PDB: %s", prediction_pdb)
                continue
            try:
                model_rows.append(
                    _model_row(
                        job=job,
                        metric=metric,
                        prediction_pdb=prediction_pdb,
                        reference_residues=reference_residues,
                        reference_sequence=reference_sequence,
                        site_indices=site_indices,
                        hotspot_indices=hotspot_indices,
                        args=args,
                    )
                )
            except Exception as exc:
                logger.error("Failed to parse %s: %s: %s", prediction_pdb, exc.__class__.__name__, exc)

    groups = _group_rows(controls, model_rows)
    write_csv(control_root / "FGA_rfpeptides_stage5_target_control_model_results.csv", model_rows, MODEL_FIELDS)
    write_csv(control_root / "FGA_rfpeptides_stage5_target_control_group_summary.csv", groups, GROUP_FIELDS)
    write_markdown(
        control_root / "FGA_rfpeptides_stage5_target_control_results.md",
        _report(model_rows, groups, args),
    )
    logger.info("Target-control models parsed: %s", len(model_rows))
    logger.info("Target-control model recovery passes: %s", sum(row["target_recovery_pass"] == "true" for row in model_rows))
    logger.info("Output directory: %s", control_root)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
