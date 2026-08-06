from __future__ import annotations

import argparse
import json
import math
import os
from collections import defaultdict
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np

from common import (
    ROUTE_PROVENANCE_FIELDS,
    add_route_provenance,
    assert_active_route_path,
    append_run_header,
    load_route_manifest,
    read_csv,
    resolve_path,
    route_provenance_fields,
    rows_to_markdown,
    setup_logger,
    sha256_file,
    validate_route_project_config,
    validate_row_route_provenance,
    write_csv,
    write_markdown,
)
from pdb_utils import parse_residues, residue_sequence


MODEL_FIELDS = [
    "stage5_context_job_id",
    "context_id",
    "template_mode",
    "seed",
    "model_name",
    "prediction_pdb",
    "context_sequence_match",
    "fga_global_CA_RMSD_A",
    "Site_2_local_CA_RMSD_A",
    "hotspot_local_CA_RMSD_A",
    "context_global_CA_RMSD_A",
    "context_after_fga_alignment_CA_RMSD_A",
    "partner_after_fga_alignment_CA_RMSD_A",
    "local_rmsd_alignment_basis",
    "fga_mean_pLDDT_100",
    "Site_2_mean_pLDDT_100",
    "hotspot_mean_pLDDT_100",
    "partner_mean_pLDDT_100",
    "ptm",
    "iptm",
    "context_interchain_pae_mean_A",
    "fga_target_recovery_pass",
    "fga_target_recovery_failure_reasons",
    "native_context_assembly_status",
] + ROUTE_PROVENANCE_FIELDS

SUMMARY_FIELDS = [
    "context_id",
    "context_description",
    "template_mode",
    "template_sequence_masked",
    "chain_lengths",
    "models_planned",
    "models_completed",
    "models_passing_fga_target_recovery",
    "median_fga_global_CA_RMSD_A",
    "best_fga_global_CA_RMSD_A",
    "median_Site_2_local_CA_RMSD_A",
    "best_Site_2_local_CA_RMSD_A",
    "median_hotspot_local_CA_RMSD_A",
    "best_hotspot_local_CA_RMSD_A",
    "median_context_after_fga_alignment_CA_RMSD_A",
    "median_fga_mean_pLDDT_100",
    "median_Site_2_mean_pLDDT_100",
    "median_hotspot_mean_pLDDT_100",
    "status",
    "interpretation",
] + ROUTE_PROVENANCE_FIELDS


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


def _finite_values(rows: Sequence[Mapping[str, Any]], field: str) -> list[float]:
    return [value for value in (_float(row.get(field, "")) for row in rows) if math.isfinite(value)]


def _median(rows: Sequence[Mapping[str, Any]], field: str) -> float | str:
    values = _finite_values(rows, field)
    return round(float(np.median(values)), 3) if values else ""


def _best(rows: Sequence[Mapping[str, Any]], field: str) -> float | str:
    values = _finite_values(rows, field)
    return round(min(values), 3) if values else ""


def _ca_coords(residues: Sequence[Mapping[str, Any]], label: str) -> np.ndarray:
    points = []
    for residue in residues:
        atoms = residue.get("atoms", {})
        if "CA" not in atoms:
            raise RuntimeError(f"{label} residue lacks CA")
        points.append(atoms["CA"])
    return np.asarray(points, dtype=float)


def _kabsch(mobile: np.ndarray, reference: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    if mobile.shape != reference.shape or mobile.shape[0] < 3:
        raise RuntimeError(f"Cannot align CA arrays with shapes {mobile.shape} and {reference.shape}")
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


def _rmsd_after_transform(
    mobile: np.ndarray,
    reference: np.ndarray,
    rotation: np.ndarray,
    translation: np.ndarray,
) -> float:
    aligned = mobile @ rotation.T + translation
    delta = aligned - reference
    return float(np.sqrt(np.mean(np.sum(delta * delta, axis=1))))


def _subset_rmsd(
    mobile: np.ndarray,
    reference: np.ndarray,
    indices_1based: Sequence[int],
    rotation: np.ndarray,
    translation: np.ndarray,
) -> float:
    indices = np.asarray([int(index) - 1 for index in indices_1based], dtype=int)
    return _rmsd_after_transform(mobile[indices], reference[indices], rotation, translation)


def _mapping_indices(mapping_csv: Path) -> tuple[list[int], list[int]]:
    rows = read_csv(mapping_csv)
    site = [
        int(row["context_residue_number"])
        for row in rows
        if row.get("context_chain") == "A" and str(row.get("is_target_site_residue", "")).lower() == "true"
    ]
    hotspots = [
        int(row["context_residue_number"])
        for row in rows
        if row.get("context_chain") == "A" and str(row.get("is_selected_hotspot", "")).lower() == "true"
    ]
    if not site or len(hotspots) != 4:
        raise RuntimeError(f"Invalid Site_2/hotspot mapping in {mapping_csv}")
    return site, hotspots


def _validate_metadata(metadata_path: Path, job: Mapping[str, str]) -> None:
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    expected = {
        "stage5_context_job_id": job["stage5_context_job_id"],
        "context_id": job["context_id"],
        "protocol_hash": job["protocol_hash"],
        "context_pdb_sha256": job["context_pdb_sha256"],
    }
    for field, value in expected.items():
        if str(metadata.get(field, "")) != str(value):
            raise RuntimeError(f"Run metadata mismatch for {field}: {metadata_path}")
    if metadata.get("template_input_verified") is not True:
        raise RuntimeError(f"Run metadata lacks verified template input: {metadata_path}")
    if metadata.get("peptide_loaded_by_prediction_runner") is not False:
        raise RuntimeError(f"Target-context control unexpectedly loaded a peptide: {metadata_path}")


def _model_row(
    *,
    job: Mapping[str, str],
    metric: Mapping[str, str],
    prediction_pdb: Path,
    reference_chains: Mapping[str, Sequence[Mapping[str, Any]]],
    site_indices: Sequence[int],
    hotspot_indices: Sequence[int],
    args: argparse.Namespace,
) -> dict[str, Any]:
    predicted_chains = parse_residues(prediction_pdb)
    reference_order = list(reference_chains)
    if list(predicted_chains) != reference_order:
        raise RuntimeError(
            f"Prediction chain order {list(predicted_chains)} != reference chain order {reference_order}"
        )
    sequence_match = all(
        residue_sequence(predicted_chains[chain]) == residue_sequence(reference_chains[chain])
        for chain in reference_order
    )
    reference_arrays = {
        chain: _ca_coords(reference_chains[chain], f"reference chain {chain}") for chain in reference_order
    }
    predicted_arrays = {
        chain: _ca_coords(predicted_chains[chain], f"prediction chain {chain}") for chain in reference_order
    }
    for chain in reference_order:
        if predicted_arrays[chain].shape != reference_arrays[chain].shape:
            raise RuntimeError(f"Prediction/reference length mismatch for chain {chain}")

    fga_rotation, fga_translation, fga_rmsd = _kabsch(predicted_arrays["A"], reference_arrays["A"])
    site_rmsd = _subset_rmsd(
        predicted_arrays["A"], reference_arrays["A"], site_indices, fga_rotation, fga_translation
    )
    hotspot_rmsd = _subset_rmsd(
        predicted_arrays["A"], reference_arrays["A"], hotspot_indices, fga_rotation, fga_translation
    )
    all_predicted = np.concatenate([predicted_arrays[chain] for chain in reference_order], axis=0)
    all_reference = np.concatenate([reference_arrays[chain] for chain in reference_order], axis=0)
    _, _, context_global_rmsd = _kabsch(all_predicted, all_reference)
    context_after_fga = _rmsd_after_transform(all_predicted, all_reference, fga_rotation, fga_translation)
    partner_after_fga: float | str = ""
    if len(reference_order) > 1:
        partner_predicted = np.concatenate([predicted_arrays[chain] for chain in reference_order[1:]], axis=0)
        partner_reference = np.concatenate([reference_arrays[chain] for chain in reference_order[1:]], axis=0)
        partner_after_fga = _rmsd_after_transform(
            partner_predicted, partner_reference, fga_rotation, fga_translation
        )

    fga_plddt = _float(metric.get("plddt_fga_mean_100", ""))
    site_plddt = _float(metric.get("plddt_site2_mean_100", ""))
    hotspot_plddt = _float(metric.get("plddt_hotspot_mean_100", ""))
    failures: list[str] = []
    if not sequence_match:
        failures.append("context_sequence_mismatch")
    if fga_rmsd > args.max_fga_global_rmsd:
        failures.append("high_fga_global_CA_RMSD")
    if site_rmsd > args.max_site2_local_rmsd:
        failures.append("high_Site_2_local_CA_RMSD")
    if hotspot_rmsd > args.max_hotspot_local_rmsd:
        failures.append("high_hotspot_local_CA_RMSD")
    if not math.isfinite(fga_plddt) or fga_plddt < args.min_fga_plddt:
        failures.append("low_fga_mean_pLDDT")
    if not math.isfinite(site_plddt) or site_plddt < args.min_site2_plddt:
        failures.append("low_Site_2_mean_pLDDT")
    if isinstance(partner_after_fga, float):
        assembly_status = (
            "native_context_assembly_pass"
            if partner_after_fga <= args.max_partner_after_fga_rmsd
            else "native_context_assembly_not_recovered"
        )
    else:
        assembly_status = "not_applicable_single_chain"

    return {
        "stage5_context_job_id": job["stage5_context_job_id"],
        "context_id": job["context_id"],
        "template_mode": job["template_mode"],
        "seed": job["seed"],
        "model_name": metric.get("model_name", prediction_pdb.stem),
        "prediction_pdb": prediction_pdb,
        "context_sequence_match": str(sequence_match).lower(),
        "fga_global_CA_RMSD_A": _rounded(fga_rmsd),
        "Site_2_local_CA_RMSD_A": _rounded(site_rmsd),
        "hotspot_local_CA_RMSD_A": _rounded(hotspot_rmsd),
        "context_global_CA_RMSD_A": _rounded(context_global_rmsd),
        "context_after_fga_alignment_CA_RMSD_A": _rounded(context_after_fga),
        "partner_after_fga_alignment_CA_RMSD_A": (
            _rounded(partner_after_fga) if isinstance(partner_after_fga, float) else ""
        ),
        "local_rmsd_alignment_basis": "after_global_FGA_chain_A_CA_alignment",
        "fga_mean_pLDDT_100": _rounded(fga_plddt),
        "Site_2_mean_pLDDT_100": _rounded(site_plddt),
        "hotspot_mean_pLDDT_100": _rounded(hotspot_plddt),
        "partner_mean_pLDDT_100": metric.get("plddt_partner_mean_100", ""),
        "ptm": metric.get("ptm", ""),
        "iptm": metric.get("iptm", ""),
        "context_interchain_pae_mean_A": metric.get("context_interchain_pae_mean_A", ""),
        "fga_target_recovery_pass": str(not failures).lower(),
        "fga_target_recovery_failure_reasons": ";".join(failures),
        "native_context_assembly_status": assembly_status,
    }


def _summaries(
    controls: Sequence[Mapping[str, str]],
    models: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    grouped: dict[str, list[Mapping[str, Any]]] = defaultdict(list)
    for row in models:
        grouped[str(row["context_id"])].append(row)
    output: list[dict[str, Any]] = []
    for control in controls:
        rows = grouped.get(str(control["context_id"]), [])
        passed = sum(row.get("fga_target_recovery_pass") == "true" for row in rows)
        if not rows:
            status = "prepared_not_run"
            interpretation = "No predictions were available."
        elif passed:
            status = "target_recovery_observed"
            interpretation = "At least one model recovered FGA Site_2 under the diagnostic thresholds."
        else:
            status = "target_not_recovered_under_thresholds"
            interpretation = "This context did not recover FGA Site_2 under the diagnostic thresholds."
        output.append(
            {
                "context_id": control["context_id"],
                "context_description": control["context_description"],
                "template_mode": control["template_mode"],
                "template_sequence_masked": control["template_sequence_masked"],
                "chain_lengths": control["chain_lengths"],
                "models_planned": control["model_predictions_planned"],
                "models_completed": len(rows),
                "models_passing_fga_target_recovery": passed,
                "median_fga_global_CA_RMSD_A": _median(rows, "fga_global_CA_RMSD_A"),
                "best_fga_global_CA_RMSD_A": _best(rows, "fga_global_CA_RMSD_A"),
                "median_Site_2_local_CA_RMSD_A": _median(rows, "Site_2_local_CA_RMSD_A"),
                "best_Site_2_local_CA_RMSD_A": _best(rows, "Site_2_local_CA_RMSD_A"),
                "median_hotspot_local_CA_RMSD_A": _median(rows, "hotspot_local_CA_RMSD_A"),
                "best_hotspot_local_CA_RMSD_A": _best(rows, "hotspot_local_CA_RMSD_A"),
                "median_context_after_fga_alignment_CA_RMSD_A": _median(
                    rows, "context_after_fga_alignment_CA_RMSD_A"
                ),
                "median_fga_mean_pLDDT_100": _median(rows, "fga_mean_pLDDT_100"),
                "median_Site_2_mean_pLDDT_100": _median(rows, "Site_2_mean_pLDDT_100"),
                "median_hotspot_mean_pLDDT_100": _median(rows, "hotspot_mean_pLDDT_100"),
                "status": status,
                "interpretation": interpretation,
            }
        )
    return output


def _report(
    models: Sequence[Mapping[str, Any]],
    summaries: Sequence[Mapping[str, Any]],
    args: argparse.Namespace,
) -> str:
    summary_columns = [
        "context_id", "template_mode", "chain_lengths", "models_completed",
        "models_passing_fga_target_recovery", "median_fga_global_CA_RMSD_A",
        "median_Site_2_local_CA_RMSD_A", "median_hotspot_local_CA_RMSD_A",
        "median_fga_mean_pLDDT_100", "status",
    ]
    model_columns = [
        "context_id", "model_name", "fga_global_CA_RMSD_A", "Site_2_local_CA_RMSD_A",
        "hotspot_local_CA_RMSD_A", "context_after_fga_alignment_CA_RMSD_A",
        "fga_mean_pLDDT_100", "Site_2_mean_pLDDT_100", "fga_target_recovery_pass",
        "fga_target_recovery_failure_reasons",
    ]
    return f"""# Stage 5 Expanded Target Context Control Results

These are target-only diagnostic controls. No peptide, peptide template, cyclic
offset, or initial guess was used. Therefore these results diagnose whether the
target representation can be recovered; they do not validate or reject peptide
candidates.

Site_2 and hotspot RMSDs are measured after one global FGA chain-A CA alignment.
The thresholds below are internal diagnostics, not experimental standards.

```text
max_fga_global_CA_RMSD_A: {args.max_fga_global_rmsd:g}
max_Site_2_local_CA_RMSD_A: {args.max_site2_local_rmsd:g}
max_hotspot_local_CA_RMSD_A: {args.max_hotspot_local_rmsd:g}
min_fga_mean_pLDDT_100: {args.min_fga_plddt:g}
min_Site_2_mean_pLDDT_100: {args.min_site2_plddt:g}
max_partner_after_FGA_alignment_CA_RMSD_A: {args.max_partner_after_fga_rmsd:g}
```

## Context Comparison

{rows_to_markdown(summaries, summary_columns, "No target contexts were found.")}

## Model Results

{rows_to_markdown(models, model_columns, "No target-context predictions have been run.")}
"""


def main() -> int:
    parser = argparse.ArgumentParser(description="Collect expanded-target AfCycDesign recovery controls.")
    parser.add_argument("--control-root", required=True)
    parser.add_argument("--project-config", required=True)
    parser.add_argument("--max-fga-global-rmsd", type=float, default=3.0)
    parser.add_argument("--max-site2-local-rmsd", type=float, default=2.0)
    parser.add_argument("--max-hotspot-local-rmsd", type=float, default=2.0)
    parser.add_argument("--min-fga-plddt", type=float, default=70.0)
    parser.add_argument("--min-site2-plddt", type=float, default=70.0)
    parser.add_argument("--max-partner-after-fga-rmsd", type=float, default=4.0)
    args = parser.parse_args()

    logger = setup_logger("33_collect_stage5_target_context_controls")
    append_run_header(logger, "33_collect_stage5_target_context_controls.py")
    control_root = assert_active_route_path(
        _resolve_mixed_path(args.control_root), "Stage 33 target-context control root"
    )
    route_manifest_path, route_manifest, route_manifest_sha256 = load_route_manifest(control_root.parent)
    validate_route_project_config(args.project_config, route_manifest)
    provenance = route_provenance_fields(route_manifest_path, route_manifest, route_manifest_sha256)
    manifest_csv = assert_active_route_path(
        control_root / "FGA_rfpeptides_stage5_target_context_control_manifest.csv",
        "Stage 33 context manifest CSV",
    )
    jobs_csv = assert_active_route_path(
        control_root / "inputs" / "FGA_rfpeptides_stage5_target_context_control_jobs.csv",
        "Stage 33 context jobs CSV",
    )
    controls = read_csv(manifest_csv)
    jobs = read_csv(jobs_csv)
    if not controls or not jobs:
        raise RuntimeError("Target-context manifest or job table is empty")
    for row in controls:
        validate_row_route_provenance(row, provenance, "Stage 33 context manifest row")
    for row in jobs:
        validate_row_route_provenance(row, provenance, "Stage 33 context job row")

    models: list[dict[str, Any]] = []
    for job in jobs:
        reference_pdb = assert_active_route_path(job["context_pdb"], "Stage 33 context reference PDB")
        if sha256_file(reference_pdb) != str(job["context_pdb_sha256"]):
            raise RuntimeError(f"Context PDB changed after task preparation: {reference_pdb}")
        reference_chains = parse_residues(reference_pdb)
        mapping_csv = reference_pdb.parent.parent / "mappings" / f"{job['context_id']}_mapping.csv"
        assert_active_route_path(mapping_csv, "Stage 33 context mapping CSV")
        site_indices, hotspot_indices = _mapping_indices(mapping_csv)
        prediction_dir = _resolve_mixed_path(job["prediction_output_dir"])
        metrics_csv = prediction_dir / "model_metrics.csv"
        metadata_json = prediction_dir / "run_metadata.json"
        if not metrics_csv.is_file() and not metadata_json.is_file():
            continue
        if not metrics_csv.is_file() or not metadata_json.is_file():
            raise RuntimeError(f"Incomplete target-context output directory: {prediction_dir}")
        assert_active_route_path(metrics_csv, "Stage 33 model metrics CSV")
        assert_active_route_path(metadata_json, "Stage 33 run metadata JSON")
        _validate_metadata(metadata_json, job)
        metrics = read_csv(metrics_csv)
        if len(metrics) != int(job["models_per_seed"]):
            raise RuntimeError(f"Expected {job['models_per_seed']} model rows in {metrics_csv}, observed {len(metrics)}")
        for metric in metrics:
            prediction_pdb = assert_active_route_path(
                _resolve_mixed_path(metric["prediction_pdb"]), "Stage 33 prediction PDB"
            )
            models.append(
                _model_row(
                    job=job,
                    metric=metric,
                    prediction_pdb=prediction_pdb,
                    reference_chains=reference_chains,
                    site_indices=site_indices,
                    hotspot_indices=hotspot_indices,
                    args=args,
                )
            )

    add_route_provenance(models, provenance)
    summaries = _summaries(controls, models)
    add_route_provenance(summaries, provenance)
    write_csv(
        control_root / "FGA_rfpeptides_stage5_target_context_control_model_results.csv",
        models,
        MODEL_FIELDS,
    )
    write_csv(
        control_root / "FGA_rfpeptides_stage5_target_context_control_summary.csv",
        summaries,
        SUMMARY_FIELDS,
    )
    write_markdown(
        control_root / "FGA_rfpeptides_stage5_target_context_control_results.md",
        _report(models, summaries, args),
    )
    logger.info("Target-context model predictions parsed: %s", len(models))
    logger.info(
        "Target-context FGA recovery passes: %s",
        sum(row["fga_target_recovery_pass"] == "true" for row in models),
    )
    logger.info("Output directory: %s", control_root)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
