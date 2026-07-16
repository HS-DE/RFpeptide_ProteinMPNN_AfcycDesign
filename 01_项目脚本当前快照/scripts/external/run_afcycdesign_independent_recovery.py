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
    "stage5_job_id",
    "target_sequence",
    "peptide_sequence",
    "target_sequence_hash",
    "peptide_sequence_hash",
    "protocol_hash",
    "seed",
    "requested_recycles",
    "forward_passes",
    "colabdesign_commit",
    "cyclic_chain_index",
    "template_mode",
    "use_initial_guess",
    "target_msa_mode",
    "peptide_msa_mode",
    "use_mlm",
    "mlm_replace_fraction",
)


def _sequence_hash(sequence: str) -> str:
    return hashlib.sha1(sequence.encode("utf-8")).hexdigest()[:8]


def _validate_sequence(sequence: Any, label: str) -> str:
    normalized = str(sequence).strip().upper()
    if not normalized:
        raise RuntimeError(f"{label} is empty")
    invalid = sorted(set(normalized) - LEGAL_AA)
    if invalid:
        raise RuntimeError(f"{label} contains unsupported amino acids: {','.join(invalid)}")
    return normalized


def _read_job_spec(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        spec = json.load(handle)
    required = [
        "stage5_job_id",
        "target_sequence",
        "target_sequence_length",
        "target_sequence_hash",
        "peptide_sequence",
        "peptide_sequence_length",
        "peptide_sequence_hash",
        "protocol_hash",
        "cyclic_chain_index",
        "seed",
        "models_per_seed",
        "requested_recycles",
        "forward_passes",
        "colabdesign_commit",
        "site2_target_indices_1based",
        "hotspot_target_indices_1based",
        "target_msa_mode",
        "peptide_msa_mode",
        "use_mlm",
        "mlm_replace_fraction",
        "template_mode",
        "use_initial_guess",
        "validation_test_type",
        "prediction_output_dir",
    ]
    missing = [key for key in required if key not in spec]
    if missing:
        raise RuntimeError(f"Job spec is missing required fields: {','.join(missing)}")
    if spec["template_mode"] != "none" or bool(spec["use_initial_guess"]):
        raise RuntimeError("This runner only permits template-free, no-initial-guess independent recovery.")
    if spec["validation_test_type"] != "independent_recovery":
        raise RuntimeError("This runner only permits validation_test_type=independent_recovery.")
    if int(spec["cyclic_chain_index"]) != 1:
        raise RuntimeError("Expected target chain index 0 and cyclic peptide chain index 1.")
    if int(spec["models_per_seed"]) != 5:
        raise RuntimeError("Independent recovery requires all five AlphaFold model parameter sets.")
    target_sequence = _validate_sequence(spec["target_sequence"], "target_sequence")
    peptide_sequence = _validate_sequence(spec["peptide_sequence"], "peptide_sequence")
    if int(spec["target_sequence_length"]) != len(target_sequence):
        raise RuntimeError("target_sequence_length does not match target_sequence")
    if int(spec["peptide_sequence_length"]) != len(peptide_sequence):
        raise RuntimeError("peptide_sequence_length does not match peptide_sequence")
    if str(spec["target_sequence_hash"]) != _sequence_hash(target_sequence):
        raise RuntimeError("target_sequence_hash does not match target_sequence")
    if str(spec["peptide_sequence_hash"]) != _sequence_hash(peptide_sequence):
        raise RuntimeError("peptide_sequence_hash does not match peptide_sequence")
    if str(spec["target_msa_mode"]) != "single_sequence" or str(spec["peptide_msa_mode"]) != "single_sequence":
        raise RuntimeError("This Stage 5A runner requires single-sequence target and peptide inputs")
    requested_recycles = int(spec["requested_recycles"])
    if int(spec["forward_passes"]) != requested_recycles + 1:
        raise RuntimeError("forward_passes must equal requested_recycles + 1")
    for field in ("site2_target_indices_1based", "hotspot_target_indices_1based"):
        values = [int(value) for value in spec[field]]
        if not values or min(values) < 1 or max(values) > len(target_sequence):
            raise RuntimeError(f"{field} contains invalid target residue indices")
        spec[field] = values
    spec["target_sequence"] = target_sequence
    spec["peptide_sequence"] = peptide_sequence
    return spec


def _bidirectional_interface_pae_values(
    pae: Any,
    target_length: int,
    target_indices_1based: list[int] | None = None,
) -> Any:
    import numpy as np

    array = np.asarray(pae, dtype=float)
    if array.ndim != 2 or array.shape[0] != array.shape[1] or array.shape[0] <= target_length:
        raise RuntimeError(f"Unexpected PAE shape for target_length={target_length}: {array.shape}")
    if target_indices_1based is None:
        target_indices = np.arange(target_length, dtype=int)
    else:
        target_indices = np.asarray([index - 1 for index in target_indices_1based], dtype=int)
    peptide_indices = np.arange(target_length, array.shape[0], dtype=int)
    left = array[np.ix_(target_indices, peptide_indices)]
    right = array[np.ix_(peptide_indices, target_indices)]
    return np.concatenate([left.reshape(-1), right.reshape(-1)])


def _pae_mean(values: Any) -> float:
    import numpy as np

    array = np.asarray(values, dtype=float)
    return float(array.mean()) if array.size else float("nan")


def _pae_median(values: Any) -> float:
    import numpy as np

    array = np.asarray(values, dtype=float)
    return float(np.median(array)) if array.size else float("nan")


def _validated_plddt_fraction(plddt: Any) -> Any:
    import numpy as np

    values = np.asarray(plddt, dtype=float)
    if values.size == 0 or not np.isfinite(values).all():
        raise RuntimeError("model.aux['plddt'] is empty or contains non-finite values")
    if float(values.min()) < -1e-6 or float(values.max()) > 1.000001:
        raise RuntimeError(
            "Expected model.aux['plddt'] on the ColabDesign 0-1 fraction scale; "
            f"observed range {float(values.min())} to {float(values.max())}"
        )
    return values


def _write_metrics(path: Path, rows: list[Mapping[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "stage5_job_id",
        "stage5_candidate_id",
        "protocol_hash",
        "peptide_sequence_hash",
        "seed",
        "model_name",
        "requested_recycles",
        "forward_passes",
        "plddt_mean_fraction",
        "plddt_peptide_mean_fraction",
        "plddt_mean_100",
        "plddt_peptide_mean_100",
        "ptm",
        "iptm",
        "ranking_confidence",
        "interface_pae_global_mean_A",
        "interface_pae_site2_mean_A",
        "interface_pae_site2_median_A",
        "interface_pae_hotspot_mean_A",
        "target_msa_mode",
        "peptide_msa_mode",
        "msa_rows_input",
        "use_mlm",
        "mlm_replace_fraction",
        "prediction_pdb",
        "prediction_npz",
        "validation_test_type",
        "template_mode",
        "use_initial_guess",
        "cyclic_topology_encoding",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _verify_colabdesign_source_marker(spec: Mapping[str, Any]) -> tuple[Path, str]:
    source_value = os.environ.get("COLABDESIGN_GAMMA_SOURCE", "").strip()
    if not source_value:
        raise RuntimeError("COLABDESIGN_GAMMA_SOURCE is required for runtime commit verification")
    source_dir = Path(os.path.expandvars(os.path.expanduser(source_value))).resolve()
    marker = source_dir / ".stage5_colabdesign_commit"
    if not marker.is_file():
        raise RuntimeError(f"ColabDesign commit marker is missing: {marker}")
    loaded_commit = marker.read_text(encoding="utf-8").strip()
    expected_commit = str(spec["colabdesign_commit"]).strip()
    if loaded_commit != expected_commit:
        raise RuntimeError(f"ColabDesign source commit is {loaded_commit}, expected {expected_commit}")
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
        with metadata_path.open("r", encoding="utf-8") as handle:
            metadata = json.load(handle)
        with metrics_path.open("r", encoding="utf-8", newline="") as handle:
            metrics = list(csv.DictReader(handle))
    except (OSError, ValueError, csv.Error):
        return False
    for field in CACHE_IDENTITY_FIELDS:
        if field not in metadata or metadata[field] != spec[field]:
            return False
    if metadata.get("loaded_colabdesign_commit") != loaded_commit:
        return False
    if metadata.get("loaded_colabdesign_source") != str(source_dir):
        return False
    if metadata.get("input_identity_verified") is not True:
        return False
    if int(metadata.get("prediction_count", 0)) != int(spec["models_per_seed"]):
        return False
    if len(metrics) != int(spec["models_per_seed"]) or len({row.get("model_name", "") for row in metrics}) != int(
        spec["models_per_seed"]
    ):
        return False
    for row in metrics:
        if row.get("stage5_job_id") != spec["stage5_job_id"]:
            return False
        if row.get("protocol_hash") != spec["protocol_hash"]:
            return False
        if row.get("peptide_sequence_hash") != spec["peptide_sequence_hash"]:
            return False
        pdb_name = Path(str(row.get("prediction_pdb", ""))).name
        npz_name = Path(str(row.get("prediction_npz", ""))).name
        if not pdb_name or not npz_name:
            return False
        if not (output_dir / pdb_name).is_file() or not (output_dir / npz_name).is_file():
            return False
    return True


def _audit_imports() -> None:
    from colabdesign import mk_af_model  # noqa: F401
    from colabdesign.af.contrib import predict  # noqa: F401
    from colabdesign.af.contrib.cyclic import add_cyclic_offset  # noqa: F401

    print("PASS: AfCycDesign gamma prediction imports are available.")


def _run_prediction(spec: Mapping[str, Any], af_params: Path) -> None:
    os.environ.setdefault("XLA_PYTHON_CLIENT_PREALLOCATE", "false")

    target_sequence = str(spec["target_sequence"]).strip().upper()
    peptide_sequence = str(spec["peptide_sequence"]).strip().upper()
    lengths = [len(target_sequence), len(peptide_sequence)]
    seed = int(spec["seed"])
    requested_recycles = int(spec["requested_recycles"])
    forward_passes = int(spec["forward_passes"])
    use_mlm = bool(spec["use_mlm"])
    mlm_replace_fraction = float(spec["mlm_replace_fraction"])
    source_dir, loaded_commit = _verify_colabdesign_source_marker(spec)
    output_dir = Path(str(spec["prediction_output_dir"]))
    output_dir.mkdir(parents=True, exist_ok=True)
    if _completed_output_is_valid(output_dir, spec, source_dir, loaded_commit):
        print(
            f"[Stage5 AfCycDesign] skipping completed job {spec['stage5_job_id']}: "
            f"validated 5 PDB/NPZ model pairs in {output_dir}",
            flush=True,
        )
        return

    print(
        f"[Stage5 AfCycDesign] starting {spec['stage5_job_id']}: "
        f"target_length={lengths[0]}, peptide_length={lengths[1]}, "
        f"seed={seed}, requested_recycles={requested_recycles}, forward_passes={forward_passes}, "
        f"use_mlm={use_mlm}, mlm_replace_fraction={mlm_replace_fraction}",
        flush=True,
    )
    print(f"[Stage5 AfCycDesign] output: {output_dir}", flush=True)

    import numpy as np
    import colabdesign
    from colabdesign import clear_mem, mk_af_model
    from colabdesign.af.contrib import predict
    from colabdesign.af.contrib.cyclic import add_cyclic_offset

    loaded_module_path = Path(colabdesign.__file__).resolve()
    if source_dir != loaded_module_path and source_dir not in loaded_module_path.parents:
        raise RuntimeError(
            f"Imported colabdesign from {loaded_module_path}, outside expected pinned source {source_dir}"
        )

    a3m_path = output_dir / "single_sequence.a3m"
    with a3m_path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(f">{spec['stage5_job_id']}\n{target_sequence}{peptide_sequence}\n")
    msa, deletion_matrix = predict.parse_a3m(str(a3m_path))
    msa_rows_input = int(msa.shape[0])
    if msa_rows_input != 1:
        raise RuntimeError(f"Stage 5A single-sequence protocol expected one MSA row, observed {msa_rows_input}")

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
    model.prep_inputs(lengths, copies=1, seed=seed)
    model.set_msa(msa, deletion_matrix)
    add_cyclic_offset(model, [int(spec["cyclic_chain_index"])])
    if use_mlm:
        model.set_opt("mlm", replace_fraction=mlm_replace_fraction)

    model_names = list(model._model_names)
    if len(model_names) != 5:
        raise RuntimeError(f"Expected five AlphaFold model parameter sets, found {len(model_names)}: {model_names}")

    metrics: list[dict[str, Any]] = []
    model.set_seed(seed)
    for model_name in model_names:
        print(f"[Stage5 AfCycDesign] running {model_name}", flush=True)
        model._inputs.pop("prev", None)
        for _ in range(forward_passes):
            model.predict(dropout=False, models=[model_name], verbose=False)
            model._inputs["prev"] = model.aux["prev"]

        safe_model = str(model_name).replace("/", "_")
        pdb_path = output_dir / f"{safe_model}_seed{seed:02d}.pdb"
        npz_path = output_dir / f"{safe_model}_seed{seed:02d}.npz"
        model.save_current_pdb(str(pdb_path))

        plddt = _validated_plddt_fraction(model.aux["plddt"])
        pae = np.asarray(model.aux["pae"], dtype=float)
        global_pae_values = _bidirectional_interface_pae_values(pae, lengths[0])
        site2_pae_values = _bidirectional_interface_pae_values(
            pae,
            lengths[0],
            list(spec["site2_target_indices_1based"]),
        )
        hotspot_pae_values = _bidirectional_interface_pae_values(
            pae,
            lengths[0],
            list(spec["hotspot_target_indices_1based"]),
        )
        log = dict(model.aux.get("log", {}))
        iptm = log.get("i_ptm", log.get("iptm", ""))
        ptm = log.get("ptm", "")
        ranking = log.get("multi", "")
        np.savez_compressed(
            npz_path,
            plddt_fraction=plddt.astype(np.float16),
            pae=pae.astype(np.float16),
            atom_positions=np.asarray(model.aux["atom_positions"], dtype=np.float32),
            target_length=np.array(lengths[0]),
            peptide_length=np.array(lengths[1]),
            site2_target_indices_1based=np.asarray(spec["site2_target_indices_1based"], dtype=np.int16),
            hotspot_target_indices_1based=np.asarray(spec["hotspot_target_indices_1based"], dtype=np.int16),
            seed=np.array(seed),
            model_name=np.array(str(model_name)),
        )
        metrics.append(
            {
                "stage5_job_id": spec["stage5_job_id"],
                "stage5_candidate_id": spec["stage5_candidate_id"],
                "protocol_hash": spec["protocol_hash"],
                "peptide_sequence_hash": spec["peptide_sequence_hash"],
                "seed": seed,
                "model_name": model_name,
                "requested_recycles": requested_recycles,
                "forward_passes": forward_passes,
                "plddt_mean_fraction": round(float(plddt.mean()), 5),
                "plddt_peptide_mean_fraction": round(float(plddt[lengths[0] :].mean()), 5),
                "plddt_mean_100": round(float(plddt.mean()) * 100.0, 3),
                "plddt_peptide_mean_100": round(float(plddt[lengths[0] :].mean()) * 100.0, 3),
                "ptm": ptm,
                "iptm": iptm,
                "ranking_confidence": ranking,
                "interface_pae_global_mean_A": round(_pae_mean(global_pae_values), 5),
                "interface_pae_site2_mean_A": round(_pae_mean(site2_pae_values), 5),
                "interface_pae_site2_median_A": round(_pae_median(site2_pae_values), 5),
                "interface_pae_hotspot_mean_A": round(_pae_mean(hotspot_pae_values), 5),
                "target_msa_mode": spec["target_msa_mode"],
                "peptide_msa_mode": spec["peptide_msa_mode"],
                "msa_rows_input": msa_rows_input,
                "use_mlm": str(use_mlm).lower(),
                "mlm_replace_fraction": mlm_replace_fraction,
                "prediction_pdb": pdb_path,
                "prediction_npz": npz_path,
                "validation_test_type": "independent_recovery",
                "template_mode": "none",
                "use_initial_guess": "false",
                "cyclic_topology_encoding": "peptide_chain_relative_position_cyclic_offset",
            }
        )
        print(
            f"[Stage5 AfCycDesign] finished {model_name}: "
            f"peptide_pLDDT_100={metrics[-1]['plddt_peptide_mean_100']}, "
            f"ipTM={metrics[-1]['iptm']}, "
            f"global_interface_PAE_A={metrics[-1]['interface_pae_global_mean_A']}, "
            f"hotspot_interface_PAE_A={metrics[-1]['interface_pae_hotspot_mean_A']}",
            flush=True,
        )

    _write_metrics(output_dir / "model_metrics.csv", metrics)
    with (output_dir / "run_metadata.json").open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(
            {
                **dict(spec),
                "af_params": str(af_params),
                "loaded_colabdesign_source": str(source_dir),
                "loaded_colabdesign_module": str(loaded_module_path),
                "loaded_colabdesign_commit": loaded_commit,
                "input_identity_verified": True,
                "models_evaluated": model_names,
                "prediction_count": len(metrics),
                "msa_rows_input": msa_rows_input,
                "plddt_source_scale": "fraction_0_1",
                "plddt_report_scales": ["fraction_0_1", "0_100"],
                "pae_units": "angstrom",
                "design_pose_loaded_by_prediction_runner": False,
                "explicit_terminal_cn_bond_record_used": False,
                "cyclic_encoding": "relative_position_cyclic_offset",
            },
            handle,
            indent=2,
            sort_keys=True,
        )
        handle.write("\n")
    print(
        f"[Stage5 AfCycDesign] completed {spec['stage5_job_id']}: "
        f"predictions={len(metrics)}",
        flush=True,
    )


def main() -> int:
    parser = argparse.ArgumentParser(description="Run one gated AfCycDesign independent-recovery seed job.")
    parser.add_argument("--job-spec")
    parser.add_argument("--af-params")
    parser.add_argument("--audit-imports", action="store_true")
    args = parser.parse_args()

    if args.audit_imports:
        _audit_imports()
        return 0
    if not args.job_spec or not args.af_params:
        parser.error("--job-spec and --af-params are required unless --audit-imports is used")
    if os.environ.get("RUN_STAGE5_PREDICTIONS") != "YES":
        raise RuntimeError("Prediction is review-gated. Set RUN_STAGE5_PREDICTIONS=YES only after manual approval.")

    spec = _read_job_spec(Path(args.job_spec))
    af_params = Path(os.path.expandvars(os.path.expanduser(args.af_params)))
    if not af_params.is_dir():
        raise RuntimeError(f"AlphaFold parameter directory does not exist: {af_params}")
    _run_prediction(spec, af_params)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
