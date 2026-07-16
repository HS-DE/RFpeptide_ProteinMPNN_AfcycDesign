from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
from typing import Any, Mapping


LEGAL_AA = set("ACDEFGHIKLMNPQRSTVWY")
CACHE_IDENTITY_FIELDS = (
    "stage5_control_job_id",
    "control_group",
    "target_sequence",
    "target_sequence_hash",
    "target_a3m_sha1",
    "target_msa_mode",
    "seed",
    "requested_recycles",
    "forward_passes",
    "colabdesign_commit",
    "protocol_hash",
    "template_mode",
    "use_initial_guess",
    "use_mlm",
    "mlm_replace_fraction",
)


def _sha1_text(value: str) -> str:
    return hashlib.sha1(value.encode("utf-8")).hexdigest()


def _sha1_file(path: Path) -> str:
    digest = hashlib.sha1()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _validate_sequence(value: Any, label: str) -> str:
    sequence = str(value).strip().upper()
    if not sequence:
        raise RuntimeError(f"{label} is empty")
    invalid = sorted(set(sequence) - LEGAL_AA)
    if invalid:
        raise RuntimeError(f"{label} contains unsupported amino acids: {','.join(invalid)}")
    return sequence


def _read_a3m_records(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header = ""
    chunks: list[str] = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header:
                    records.append((header, "".join(chunks)))
                header = line[1:].strip() or f"row_{len(records) + 1}"
                chunks = []
            elif not line.startswith("#"):
                if not header:
                    raise RuntimeError(f"A3M sequence encountered before a header: {path}")
                chunks.append(line)
    if header:
        records.append((header, "".join(chunks)))
    if not records:
        raise RuntimeError(f"Missing A3M records: {path}")
    return records


def _a3m_query_match_sequence(a3m_sequence: str) -> str:
    return "".join(ch for ch in a3m_sequence if not ch.islower() and ch != "-").upper()


def _read_job_spec(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        spec = json.load(handle)
    required = [
        "stage5_control_job_id",
        "control_group",
        "target_sequence",
        "target_sequence_length",
        "target_sequence_hash",
        "target_a3m",
        "target_a3m_sha1",
        "target_msa_mode",
        "msa_rows_expected",
        "seed",
        "models_per_seed",
        "requested_recycles",
        "forward_passes",
        "site2_target_indices_1based",
        "hotspot_target_indices_1based",
        "use_mlm",
        "mlm_replace_fraction",
        "template_mode",
        "use_initial_guess",
        "validation_test_type",
        "protocol_hash",
        "colabdesign_commit",
        "prediction_output_dir",
    ]
    missing = [field for field in required if field not in spec]
    if missing:
        raise RuntimeError(f"Control job spec is missing required fields: {','.join(missing)}")
    if spec["template_mode"] != "none" or bool(spec["use_initial_guess"]):
        raise RuntimeError("Target-only controls require template=none and use_initial_guess=false")
    if spec["validation_test_type"] != "sequence_based_target_recovery_control":
        raise RuntimeError("Unexpected validation_test_type for target-only control")
    if str(spec["target_msa_mode"]) not in {"single_sequence", "unpaired_homolog_msa"}:
        raise RuntimeError("target_msa_mode must be single_sequence or unpaired_homolog_msa")
    if int(spec["models_per_seed"]) != 5:
        raise RuntimeError("Target-only controls require all five AlphaFold model parameter sets")

    target_sequence = _validate_sequence(spec["target_sequence"], "target_sequence")
    if int(spec["target_sequence_length"]) != len(target_sequence):
        raise RuntimeError("target_sequence_length does not match target_sequence")
    if str(spec["target_sequence_hash"]) != _sha1_text(target_sequence)[:8]:
        raise RuntimeError("target_sequence_hash does not match target_sequence")
    requested_recycles = int(spec["requested_recycles"])
    if int(spec["forward_passes"]) != requested_recycles + 1:
        raise RuntimeError("forward_passes must equal requested_recycles + 1")
    for field in ("site2_target_indices_1based", "hotspot_target_indices_1based"):
        indices = [int(value) for value in spec[field]]
        if not indices or min(indices) < 1 or max(indices) > len(target_sequence):
            raise RuntimeError(f"{field} contains invalid target indices")
        spec[field] = indices
    spec["target_sequence"] = target_sequence
    return spec


def _validated_plddt_fraction(value: Any) -> Any:
    import numpy as np

    plddt = np.asarray(value, dtype=float)
    if plddt.size == 0 or not np.isfinite(plddt).all():
        raise RuntimeError("model.aux['plddt'] is empty or contains non-finite values")
    if float(plddt.min()) < -1e-6 or float(plddt.max()) > 1.000001:
        raise RuntimeError(
            "Expected model.aux['plddt'] on the 0-1 fraction scale; "
            f"observed range {float(plddt.min())} to {float(plddt.max())}"
        )
    return plddt


def _write_metrics(path: Path, rows: list[Mapping[str, Any]]) -> None:
    fields = [
        "stage5_control_job_id",
        "control_group",
        "protocol_hash",
        "target_msa_mode",
        "msa_rows_input",
        "use_mlm",
        "mlm_replace_fraction",
        "seed",
        "model_name",
        "requested_recycles",
        "forward_passes",
        "plddt_mean_fraction",
        "plddt_site2_mean_fraction",
        "plddt_hotspot_mean_fraction",
        "plddt_mean_100",
        "plddt_site2_mean_100",
        "plddt_hotspot_mean_100",
        "ptm",
        "ranking_confidence",
        "prediction_pdb",
        "prediction_npz",
        "template_mode",
        "use_initial_guess",
        "validation_test_type",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _verify_colabdesign_source(spec: Mapping[str, Any]) -> tuple[Path, str]:
    source_value = os.environ.get("COLABDESIGN_GAMMA_SOURCE", "").strip()
    if not source_value:
        raise RuntimeError("COLABDESIGN_GAMMA_SOURCE is required")
    source_dir = Path(os.path.expandvars(os.path.expanduser(source_value))).resolve()
    marker = source_dir / ".stage5_colabdesign_commit"
    if not marker.is_file():
        raise RuntimeError(f"ColabDesign commit marker is missing: {marker}")
    loaded_commit = marker.read_text(encoding="utf-8").strip()
    if loaded_commit != str(spec["colabdesign_commit"]):
        raise RuntimeError(f"ColabDesign source commit is {loaded_commit}, expected {spec['colabdesign_commit']}")
    return source_dir, loaded_commit


def _completed_output_is_valid(
    output_dir: Path,
    spec: Mapping[str, Any],
    source_dir: Path,
    loaded_commit: str,
) -> bool:
    metadata_path = output_dir / "run_metadata.json"
    metrics_path = output_dir / "model_metrics.csv"
    if not metadata_path.is_file() or not metrics_path.is_file():
        return False
    try:
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        with metrics_path.open("r", encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle))
    except (OSError, ValueError, csv.Error):
        return False
    for field in CACHE_IDENTITY_FIELDS:
        if metadata.get(field) != spec.get(field):
            return False
    if metadata.get("loaded_colabdesign_source") != str(source_dir):
        return False
    if metadata.get("loaded_colabdesign_commit") != loaded_commit:
        return False
    if metadata.get("input_identity_verified") is not True:
        return False
    if int(metadata.get("prediction_count", 0)) != int(spec["models_per_seed"]):
        return False
    if len(rows) != int(spec["models_per_seed"]):
        return False
    if len({row.get("model_name", "") for row in rows}) != int(spec["models_per_seed"]):
        return False
    for row in rows:
        if row.get("stage5_control_job_id") != spec["stage5_control_job_id"]:
            return False
        if row.get("protocol_hash") != spec["protocol_hash"]:
            return False
        pdb_name = Path(str(row.get("prediction_pdb", ""))).name
        npz_name = Path(str(row.get("prediction_npz", ""))).name
        if not pdb_name or not npz_name:
            return False
        if not (output_dir / pdb_name).is_file() or not (output_dir / npz_name).is_file():
            return False
    return True


def _run(spec: Mapping[str, Any], af_params: Path) -> None:
    os.environ.setdefault("XLA_PYTHON_CLIENT_PREALLOCATE", "false")
    source_dir, loaded_commit = _verify_colabdesign_source(spec)
    target_sequence = str(spec["target_sequence"])
    target_a3m = Path(str(spec["target_a3m"])).resolve()
    if not target_a3m.is_file():
        raise RuntimeError(f"Missing target A3M: {target_a3m}")
    actual_a3m_sha1 = _sha1_file(target_a3m)
    if actual_a3m_sha1 != str(spec["target_a3m_sha1"]):
        raise RuntimeError("Target A3M content changed after job preparation; regenerate control jobs")
    records = _read_a3m_records(target_a3m)
    if _a3m_query_match_sequence(records[0][1]) != target_sequence:
        raise RuntimeError("Target A3M query does not match the Stage 0 target sequence")
    if int(spec["msa_rows_expected"]) != len(records):
        raise RuntimeError(
            f"Target A3M row count changed: observed {len(records)}, expected {spec['msa_rows_expected']}"
        )

    output_dir = Path(str(spec["prediction_output_dir"]))
    output_dir.mkdir(parents=True, exist_ok=True)
    if _completed_output_is_valid(output_dir, spec, source_dir, loaded_commit):
        print(f"[Stage5 target control] validated cache; skipping {spec['stage5_control_job_id']}", flush=True)
        return

    import numpy as np
    import colabdesign
    from colabdesign import clear_mem, mk_af_model
    from colabdesign.af.contrib import predict

    loaded_module_path = Path(colabdesign.__file__).resolve()
    if source_dir != loaded_module_path and source_dir not in loaded_module_path.parents:
        raise RuntimeError(f"Imported colabdesign from {loaded_module_path}, expected source {source_dir}")

    msa, deletion_matrix = predict.parse_a3m(str(target_a3m))
    msa_rows_input = int(msa.shape[0])
    if msa_rows_input != len(records):
        raise RuntimeError(f"Parsed A3M rows={msa_rows_input}, file records={len(records)}")

    use_mlm = bool(spec["use_mlm"])
    replace_fraction = float(spec["mlm_replace_fraction"])
    requested_recycles = int(spec["requested_recycles"])
    forward_passes = int(spec["forward_passes"])
    seed = int(spec["seed"])
    site_indices = np.asarray(spec["site2_target_indices_1based"], dtype=int) - 1
    hotspot_indices = np.asarray(spec["hotspot_target_indices_1based"], dtype=int) - 1

    print(
        f"[Stage5 target control] starting {spec['stage5_control_job_id']}: "
        f"msa_mode={spec['target_msa_mode']}, msa_rows={msa_rows_input}, use_mlm={use_mlm}, "
        f"requested_recycles={requested_recycles}, forward_passes={forward_passes}",
        flush=True,
    )
    clear_mem()
    model = mk_af_model(
        data_dir=str(af_params),
        use_mlm=use_mlm,
        num_msa=512,
        num_extra_msa=1024,
        num_templates=1,
        use_cluster_profile=True,
        model_type="alphafold2_multimer_v3",
        use_templates=False,
        use_batch_as_template=False,
        use_dgram=True,
        protocol="hallucination",
        best_metric="multi",
        optimize_seq=False,
        debug=True,
        clear_prev=False,
        query_bias=True,
    )
    model.prep_inputs([len(target_sequence)], copies=1, seed=seed)
    model.set_msa(msa, deletion_matrix)
    if use_mlm:
        model.set_opt("mlm", replace_fraction=replace_fraction)

    model_names = list(model._model_names)
    if len(model_names) != 5:
        raise RuntimeError(f"Expected five AlphaFold model parameter sets, found {model_names}")
    model.set_seed(seed)
    metrics: list[dict[str, Any]] = []
    for model_name in model_names:
        print(f"[Stage5 target control] running {model_name}", flush=True)
        model._inputs.pop("prev", None)
        for _ in range(forward_passes):
            model.predict(dropout=False, models=[model_name], verbose=False)
            model._inputs["prev"] = model.aux["prev"]

        safe_model = str(model_name).replace("/", "_")
        pdb_path = output_dir / f"{safe_model}_seed{seed:02d}.pdb"
        npz_path = output_dir / f"{safe_model}_seed{seed:02d}.npz"
        model.save_current_pdb(str(pdb_path))
        plddt = _validated_plddt_fraction(model.aux["plddt"])
        log = dict(model.aux.get("log", {}))
        np.savez_compressed(
            npz_path,
            plddt_fraction=plddt.astype(np.float16),
            pae=np.asarray(model.aux["pae"], dtype=np.float16),
            atom_positions=np.asarray(model.aux["atom_positions"], dtype=np.float32),
            target_length=np.array(len(target_sequence)),
            site2_target_indices_1based=np.asarray(spec["site2_target_indices_1based"], dtype=np.int16),
            hotspot_target_indices_1based=np.asarray(spec["hotspot_target_indices_1based"], dtype=np.int16),
            seed=np.array(seed),
            model_name=np.array(str(model_name)),
        )
        metrics.append(
            {
                "stage5_control_job_id": spec["stage5_control_job_id"],
                "control_group": spec["control_group"],
                "protocol_hash": spec["protocol_hash"],
                "target_msa_mode": spec["target_msa_mode"],
                "msa_rows_input": msa_rows_input,
                "use_mlm": str(use_mlm).lower(),
                "mlm_replace_fraction": replace_fraction,
                "seed": seed,
                "model_name": model_name,
                "requested_recycles": requested_recycles,
                "forward_passes": forward_passes,
                "plddt_mean_fraction": round(float(plddt.mean()), 5),
                "plddt_site2_mean_fraction": round(float(plddt[site_indices].mean()), 5),
                "plddt_hotspot_mean_fraction": round(float(plddt[hotspot_indices].mean()), 5),
                "plddt_mean_100": round(float(plddt.mean()) * 100.0, 3),
                "plddt_site2_mean_100": round(float(plddt[site_indices].mean()) * 100.0, 3),
                "plddt_hotspot_mean_100": round(float(plddt[hotspot_indices].mean()) * 100.0, 3),
                "ptm": log.get("ptm", ""),
                "ranking_confidence": log.get("multi", ""),
                "prediction_pdb": pdb_path,
                "prediction_npz": npz_path,
                "template_mode": "none",
                "use_initial_guess": "false",
                "validation_test_type": "sequence_based_target_recovery_control",
            }
        )
    _write_metrics(output_dir / "model_metrics.csv", metrics)
    metadata = {
        **dict(spec),
        "af_params": str(af_params),
        "loaded_colabdesign_source": str(source_dir),
        "loaded_colabdesign_module": str(loaded_module_path),
        "loaded_colabdesign_commit": loaded_commit,
        "input_identity_verified": True,
        "msa_rows_input": msa_rows_input,
        "models_evaluated": model_names,
        "prediction_count": len(metrics),
        "plddt_source_scale": "fraction_0_1",
        "plddt_report_scales": ["fraction_0_1", "0_100"],
        "target_template_used": False,
        "initial_guess_used": False,
        "peptide_included": False,
        "cyclic_offset_applied": False,
    }
    with (output_dir / "run_metadata.json").open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(metadata, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(f"[Stage5 target control] completed {spec['stage5_control_job_id']}: models={len(metrics)}", flush=True)


def main() -> int:
    parser = argparse.ArgumentParser(description="Run one Stage 5 target-only recovery control job.")
    parser.add_argument("--job-spec", required=True)
    parser.add_argument("--af-params", required=True)
    args = parser.parse_args()
    spec = _read_job_spec(Path(args.job_spec))
    _run(spec, Path(args.af_params))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
