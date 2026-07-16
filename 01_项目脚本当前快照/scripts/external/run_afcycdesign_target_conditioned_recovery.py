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
    "stage5B_job_id",
    "stage5B_candidate_id",
    "target_sequence",
    "target_sequence_length",
    "target_sequence_hash",
    "peptide_sequence",
    "peptide_sequence_length",
    "peptide_sequence_hash",
    "target_template_sha1",
    "seed",
    "protocol_hash",
    "protocol_version",
    "requested_recycles",
    "forward_passes",
    "colabdesign_commit",
    "cyclic_chain_index",
    "cyclic_topology_encoding",
    "template_mode",
    "use_templates",
    "use_batch_as_template",
    "template_sequence_masked",
    "template_sidechains_masked",
    "template_interchain_features_masked",
    "use_initial_guess",
    "target_msa_mode",
    "peptide_msa_mode",
    "use_mlm",
    "mlm_replace_fraction",
    "use_dropout",
    "validation_test_type",
)

METRIC_FIELDS = [
    "stage5B_job_id",
    "stage5B_candidate_id",
    "protocol_hash",
    "peptide_sequence_hash",
    "target_template_sha1",
    "seed",
    "model_name",
    "requested_recycles",
    "forward_passes",
    "target_template_path",
    "target_template_residues_covered",
    "target_template_coverage_fraction",
    "peptide_template_residues_covered",
    "peptide_template_atom_mask_sum",
    "template_input_verified",
    "plddt_mean_fraction",
    "plddt_target_mean_fraction",
    "plddt_peptide_mean_fraction",
    "plddt_site2_mean_fraction",
    "plddt_hotspot_mean_fraction",
    "plddt_mean_100",
    "plddt_target_mean_100",
    "plddt_peptide_mean_100",
    "plddt_site2_mean_100",
    "plddt_hotspot_mean_100",
    "ptm",
    "iptm",
    "ranking_confidence",
    "ranking_confidence_source",
    "interface_pae_global_mean_A",
    "interface_pae_site2_mean_A",
    "interface_pae_site2_median_A",
    "interface_pae_hotspot_mean_A",
    "target_msa_mode",
    "peptide_msa_mode",
    "msa_rows_input",
    "use_mlm",
    "use_dropout",
    "prediction_pdb",
    "prediction_npz",
    "validation_test_type",
    "template_mode",
    "use_initial_guess",
    "cyclic_topology_encoding",
]


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


def _pdb_chain_sequences(path: Path) -> dict[str, str]:
    aa3_to_aa1 = {
        "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
        "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "MSE": "M",
        "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    }
    seen: dict[str, set[tuple[str, str]]] = {}
    sequences: dict[str, list[str]] = {}
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith("ATOM  "):
                continue
            chain = line[21].strip() or "_"
            key = (line[22:26].strip(), line[26].strip())
            seen.setdefault(chain, set())
            sequences.setdefault(chain, [])
            if key in seen[chain]:
                continue
            seen[chain].add(key)
            resname = line[17:20].strip().upper()
            if resname not in aa3_to_aa1:
                raise RuntimeError(f"Unsupported residue {resname} in {path}")
            sequences[chain].append(aa3_to_aa1[resname])
    return {chain: "".join(sequence) for chain, sequence in sequences.items()}


def _read_job_spec(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        spec = json.load(handle)
    required = set(CACHE_IDENTITY_FIELDS) | {
        "target_template_path",
        "target_template_chain",
        "target_template_expected_coverage",
        "peptide_template_expected_coverage",
        "reference_design_pdb_for_posthoc",
        "reference_design_loaded_by_prediction_runner",
        "site2_target_indices_1based",
        "hotspot_target_indices_1based",
        "models_per_seed",
        "model_names_expected",
        "prediction_output_dir",
        "template_padding_mode",
    }
    missing = sorted(required - set(spec))
    if missing:
        raise RuntimeError(f"Stage 5B job spec is missing required fields: {','.join(missing)}")

    target = _validate_sequence(spec["target_sequence"], "target_sequence")
    peptide = _validate_sequence(spec["peptide_sequence"], "peptide_sequence")
    if len(target) != 86 or int(spec["target_sequence_length"]) != len(target):
        raise RuntimeError("Stage 5B requires an 86-aa target sequence")
    if int(spec["peptide_sequence_length"]) != len(peptide):
        raise RuntimeError("peptide_sequence_length does not match peptide_sequence")
    if str(spec["target_sequence_hash"]) != _sha1_text(target)[:8]:
        raise RuntimeError("target_sequence_hash does not match target_sequence")
    if str(spec["peptide_sequence_hash"]) != _sha1_text(peptide)[:8]:
        raise RuntimeError("peptide_sequence_hash does not match peptide_sequence")
    if str(spec["validation_test_type"]) != "target_structure_conditioned_recovery":
        raise RuntimeError("Stage 5B validation_test_type must be target_structure_conditioned_recovery")
    if str(spec["template_mode"]) != "target_only":
        raise RuntimeError("Stage 5B template_mode must be target_only")
    if not bool(spec["use_templates"]) or bool(spec["use_batch_as_template"]):
        raise RuntimeError("Stage 5B requires use_templates=true and use_batch_as_template=false")
    if not all(bool(spec[field]) for field in (
        "template_sequence_masked", "template_sidechains_masked", "template_interchain_features_masked"
    )):
        raise RuntimeError("All Stage 5B target-template masking controls must be enabled")
    if bool(spec["use_initial_guess"]):
        raise RuntimeError("Stage 5B forbids an initial guess")
    if bool(spec["reference_design_loaded_by_prediction_runner"]):
        raise RuntimeError("Stage 4 reference design must not be loaded by the prediction runner")
    if int(spec["cyclic_chain_index"]) != 1:
        raise RuntimeError("Cyclic offset must be applied only to the second peptide chain")
    if str(spec["target_msa_mode"]) != "single_sequence" or str(spec["peptide_msa_mode"]) != "single_sequence":
        raise RuntimeError("Stage 5B v1 requires single-sequence target and peptide inputs")
    if bool(spec["use_mlm"]) or float(spec["mlm_replace_fraction"]) != 0.0:
        raise RuntimeError("Stage 5B v1 requires use_mlm=false")
    if not bool(spec["use_dropout"]):
        raise RuntimeError("Stage 5B v1 uses dropout=true so the seed ensemble is stochastic")
    if int(spec["requested_recycles"]) != 6 or int(spec["forward_passes"]) != 7:
        raise RuntimeError("Stage 5B v1 requires requested_recycles=6 and forward_passes=7")
    if int(spec["models_per_seed"]) != 5:
        raise RuntimeError("Stage 5B requires all five AlphaFold multimer parameter sets")
    expected_models = [f"model_{index}_multimer_v3" for index in range(1, 6)]
    if list(spec["model_names_expected"]) != expected_models:
        raise RuntimeError("model_names_expected does not list all five multimer_v3 models in order")
    if int(spec["target_template_expected_coverage"]) != len(target) or int(spec["peptide_template_expected_coverage"]) != 0:
        raise RuntimeError("Unexpected target/peptide template coverage requirements")
    for field in ("site2_target_indices_1based", "hotspot_target_indices_1based"):
        values = [int(value) for value in spec[field]]
        if not values or min(values) < 1 or max(values) > len(target):
            raise RuntimeError(f"{field} contains invalid target indices")
        spec[field] = values

    target_template = Path(str(spec["target_template_path"]))
    if not target_template.is_file():
        raise RuntimeError(f"Target template does not exist: {target_template}")
    if _sha1_file(target_template) != str(spec["target_template_sha1"]):
        raise RuntimeError("target_template_sha1 does not match the target template file")
    template_sequences = _pdb_chain_sequences(target_template)
    if set(template_sequences) != {"A"}:
        raise RuntimeError(f"Target template must contain only chain A; observed {sorted(template_sequences)}")
    if template_sequences["A"] != target:
        raise RuntimeError("Target template sequence does not match target_sequence")

    spec["target_sequence"] = target
    spec["peptide_sequence"] = peptide
    return spec


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
        raise RuntimeError(f"Loaded ColabDesign commit {loaded_commit} does not match the job spec")
    return source_dir, loaded_commit


def _build_padded_target_template(spec: Mapping[str, Any]) -> tuple[dict[str, Any], dict[str, Any]]:
    import numpy as np
    from colabdesign.af.alphafold.common import residue_constants
    from colabdesign.af.prep import prep_pdb

    target_length = int(spec["target_sequence_length"])
    peptide_length = int(spec["peptide_sequence_length"])
    total_length = target_length + peptide_length
    target_path = Path(str(spec["target_template_path"]))
    target_batch = prep_pdb(str(target_path), chain=str(spec["target_template_chain"]), ignore_missing=True)["batch"]
    if int(target_batch["aatype"].shape[0]) != target_length:
        raise RuntimeError(f"Target template parser returned {target_batch['aatype'].shape[0]} residues, expected {target_length}")

    aatype = np.full((total_length,), 21, dtype=np.int32)
    atom_positions = np.zeros((total_length,) + tuple(target_batch["all_atom_positions"].shape[1:]), dtype=np.float32)
    atom_mask = np.zeros((total_length,) + tuple(target_batch["all_atom_mask"].shape[1:]), dtype=np.float32)
    aatype[:target_length] = np.asarray(target_batch["aatype"], dtype=np.int32)
    atom_positions[:target_length] = np.asarray(target_batch["all_atom_positions"], dtype=np.float32)
    atom_mask[:target_length] = np.asarray(target_batch["all_atom_mask"], dtype=np.float32)

    ca_index = residue_constants.atom_order["CA"]
    target_ca_coverage = int(np.count_nonzero(atom_mask[:target_length, ca_index] > 0.5))
    peptide_atom_mask_sum = float(atom_mask[target_length:].sum())
    peptide_position_nonzero_count = int(np.count_nonzero(atom_positions[target_length:]))
    peptide_unknown_count = int(np.count_nonzero(aatype[target_length:] == 21))
    if target_ca_coverage != target_length:
        raise RuntimeError(f"Target template CA coverage is {target_ca_coverage}/{target_length}")
    if peptide_atom_mask_sum != 0.0 or peptide_position_nonzero_count != 0:
        raise RuntimeError("Peptide template coordinates or atom masks are non-zero")
    if peptide_unknown_count != peptide_length:
        raise RuntimeError("Every peptide template aatype must be unknown")
    if not np.isfinite(atom_positions).all() or not np.isfinite(atom_mask).all():
        raise RuntimeError("Template tensors contain non-finite values")

    batch = {
        "aatype": aatype,
        "all_atom_positions": atom_positions,
        "all_atom_mask": atom_mask,
    }
    audit = {
        "template_total_length": total_length,
        "target_template_residues_covered": target_ca_coverage,
        "target_template_coverage_fraction": target_ca_coverage / target_length,
        "peptide_template_residues_covered": 0,
        "peptide_template_atom_mask_sum": peptide_atom_mask_sum,
        "peptide_template_position_nonzero_count": peptide_position_nonzero_count,
        "peptide_template_unknown_aatype_count": peptide_unknown_count,
        "template_input_verified": True,
    }
    return batch, audit


def _verify_cyclic_offset(model: Any, target_length: int, peptide_length: int) -> None:
    import numpy as np

    if list(model._lengths) != [target_length, peptide_length]:
        raise RuntimeError(f"Unexpected model chain lengths: {model._lengths}")
    offset = np.asarray(model._inputs.get("offset"))
    linear = np.asarray(model._inputs["residue_index"])[:, None] - np.asarray(model._inputs["residue_index"])[None, :]
    if offset.shape != linear.shape:
        raise RuntimeError("Cyclic offset matrix shape is invalid")
    if not np.array_equal(offset[:target_length, :target_length], linear[:target_length, :target_length]):
        raise RuntimeError("Cyclic offset incorrectly changed the target-chain block")
    if not np.array_equal(offset[:target_length, target_length:], linear[:target_length, target_length:]):
        raise RuntimeError("Cyclic offset incorrectly changed target-to-peptide offsets")
    if not np.array_equal(offset[target_length:, :target_length], linear[target_length:, :target_length]):
        raise RuntimeError("Cyclic offset incorrectly changed peptide-to-target offsets")
    peptide_block = offset[target_length:, target_length:]
    if peptide_length > 2 and np.array_equal(peptide_block, linear[target_length:, target_length:]):
        raise RuntimeError("Peptide-chain block was not cyclically encoded")
    if peptide_length > 1 and abs(int(peptide_block[0, -1])) != 1:
        raise RuntimeError("Peptide terminal positions are not adjacent in the cyclic offset matrix")


def _verify_cyclic_offset_static(spec: Mapping[str, Any]) -> None:
    import numpy as np
    from colabdesign.af.contrib.cyclic import add_cyclic_offset

    target_length = int(spec["target_sequence_length"])
    peptide_length = int(spec["peptide_sequence_length"])

    class DummyModel:
        pass

    model = DummyModel()
    model._lengths = [target_length, peptide_length]
    model._inputs = {"residue_index": np.concatenate([np.arange(target_length), np.arange(peptide_length)])}
    add_cyclic_offset(model, [int(spec["cyclic_chain_index"])])
    _verify_cyclic_offset(model, target_length, peptide_length)


def _bidirectional_pae_values(pae: Any, target_length: int, target_indices_1based: list[int] | None = None) -> Any:
    import numpy as np

    array = np.asarray(pae, dtype=float)
    if array.ndim != 2 or array.shape[0] != array.shape[1] or array.shape[0] <= target_length:
        raise RuntimeError(f"Unexpected PAE shape: {array.shape}")
    target_indices = (
        np.arange(target_length, dtype=int)
        if target_indices_1based is None
        else np.asarray([index - 1 for index in target_indices_1based], dtype=int)
    )
    peptide_indices = np.arange(target_length, array.shape[0], dtype=int)
    return np.concatenate(
        [array[np.ix_(target_indices, peptide_indices)].reshape(-1), array[np.ix_(peptide_indices, target_indices)].reshape(-1)]
    )


def _validated_plddt_fraction(value: Any) -> Any:
    import numpy as np

    array = np.asarray(value, dtype=float)
    if not array.size or not np.isfinite(array).all() or float(array.min()) < -1e-6 or float(array.max()) > 1.000001:
        raise RuntimeError("Expected finite model.aux['plddt'] values on the 0-1 fraction scale")
    return array


def _scalar(value: Any) -> float | str:
    try:
        return round(float(value), 6)
    except (TypeError, ValueError):
        return ""


def _write_metrics(path: Path, rows: list[Mapping[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=METRIC_FIELDS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


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
            metrics = list(csv.DictReader(handle))
    except (OSError, ValueError, csv.Error):
        return False
    for field in CACHE_IDENTITY_FIELDS:
        if metadata.get(field) != spec.get(field):
            return False
    if metadata.get("loaded_colabdesign_commit") != loaded_commit:
        return False
    if metadata.get("loaded_colabdesign_source") != str(source_dir):
        return False
    if metadata.get("template_input_verified") is not True:
        return False
    if int(metadata.get("target_template_residues_covered", 0)) != int(spec["target_sequence_length"]):
        return False
    if int(metadata.get("peptide_template_residues_covered", -1)) != 0:
        return False
    expected_models = set(spec["model_names_expected"])
    if len(metrics) != int(spec["models_per_seed"]) or {row.get("model_name", "") for row in metrics} != expected_models:
        return False
    for row in metrics:
        if row.get("stage5B_job_id") != spec["stage5B_job_id"] or row.get("protocol_hash") != spec["protocol_hash"]:
            return False
        if row.get("target_template_sha1") != spec["target_template_sha1"]:
            return False
        pdb_path = Path(str(row.get("prediction_pdb", "")))
        npz_path = Path(str(row.get("prediction_npz", "")))
        if not pdb_path.is_file() or not npz_path.is_file():
            return False
    return True


def _audit_imports() -> None:
    from colabdesign import mk_af_model  # noqa: F401
    from colabdesign.af.contrib import predict  # noqa: F401
    from colabdesign.af.contrib.cyclic import add_cyclic_offset  # noqa: F401
    from colabdesign.af.prep import prep_pdb  # noqa: F401

    print("PASS: Stage 5B gamma prediction, target-template, and cyclic-offset imports are available.")


def _preflight(spec: Mapping[str, Any], af_params: Path) -> None:
    source_dir, loaded_commit = _verify_colabdesign_source(spec)
    if not af_params.is_dir():
        raise RuntimeError(f"AlphaFold parameter directory does not exist: {af_params}")
    import colabdesign

    loaded_module = Path(colabdesign.__file__).resolve()
    if source_dir != loaded_module and source_dir not in loaded_module.parents:
        raise RuntimeError(f"Imported colabdesign from {loaded_module}, expected pinned source {source_dir}")
    _, audit = _build_padded_target_template(spec)
    _verify_cyclic_offset_static(spec)
    print(
        f"PASS: {spec['stage5B_candidate_id']} target-only padded template "
        f"target_coverage={audit['target_template_residues_covered']}/86, peptide_coverage=0, "
        f"peptide_length={spec['peptide_sequence_length']}, commit={loaded_commit}",
        flush=True,
    )


def _run_prediction(spec: Mapping[str, Any], af_params: Path) -> None:
    os.environ.setdefault("XLA_PYTHON_CLIENT_PREALLOCATE", "false")
    import numpy as np
    import colabdesign
    from colabdesign import clear_mem, mk_af_model
    from colabdesign.af.contrib import predict
    from colabdesign.af.contrib.cyclic import add_cyclic_offset

    source_dir, loaded_commit = _verify_colabdesign_source(spec)
    loaded_module = Path(colabdesign.__file__).resolve()
    if source_dir != loaded_module and source_dir not in loaded_module.parents:
        raise RuntimeError(f"Imported colabdesign from {loaded_module}, expected pinned source {source_dir}")

    output_dir = Path(str(spec["prediction_output_dir"]))
    output_dir.mkdir(parents=True, exist_ok=True)
    if _completed_output_is_valid(output_dir, spec, source_dir, loaded_commit):
        print(f"[Stage5B] skipping cache-valid completed job {spec['stage5B_job_id']}", flush=True)
        return

    target = str(spec["target_sequence"])
    peptide = str(spec["peptide_sequence"])
    lengths = [len(target), len(peptide)]
    seed = int(spec["seed"])
    template_batch, template_audit = _build_padded_target_template(spec)
    _verify_cyclic_offset_static(spec)
    print(
        f"[Stage5B] starting {spec['stage5B_job_id']}: target={len(target)}, peptide={len(peptide)}, "
        f"seed={seed}, recycles={spec['requested_recycles']}, dropout=true, MLM=false",
        flush=True,
    )

    a3m_path = output_dir / "single_sequence.a3m"
    with a3m_path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(f">{spec['stage5B_job_id']}\n{target}{peptide}\n")
    msa, deletion_matrix = predict.parse_a3m(str(a3m_path))
    if int(msa.shape[0]) != 1:
        raise RuntimeError(f"Expected one MSA row, observed {msa.shape[0]}")

    clear_mem()
    model = mk_af_model(
        data_dir=str(af_params),
        use_mlm=False,
        num_msa=512,
        num_extra_msa=1024,
        num_templates=1,
        use_cluster_profile=True,
        model_type="alphafold2_multimer_v3",
        use_templates=True,
        use_batch_as_template=False,
        use_dgram=False,
        protocol="hallucination",
        best_metric="multi",
        optimize_seq=False,
        debug=True,
        clear_prev=True,
        query_bias=True,
    )
    model.prep_inputs(lengths, copies=1, seed=seed)
    model.set_msa(msa, deletion_matrix)
    model.set_template(batch=template_batch, n=0)
    model.set_opt("template", rm_seq=True, rm_sc=True, rm_ic=True)
    add_cyclic_offset(model, [int(spec["cyclic_chain_index"])])
    _verify_cyclic_offset(model, lengths[0], lengths[1])

    template_mask = np.asarray(model._inputs["template_all_atom_mask"])[0]
    if int(np.count_nonzero(template_mask[lengths[0] :])) != 0:
        raise RuntimeError("Actual model input contains peptide template atom masks")
    if bool(model._args.get("use_initial_guess")):
        raise RuntimeError("Actual model unexpectedly enabled an initial guess")
    if bool(model._args.get("use_batch_as_template")):
        raise RuntimeError("Actual model unexpectedly enabled use_batch_as_template")

    model_names = list(model._model_names)
    if model_names != list(spec["model_names_expected"]):
        raise RuntimeError(f"Loaded model parameter sets do not match the job spec: {model_names}")

    target_indices = np.arange(lengths[0], dtype=int)
    site_indices = np.asarray([int(index) - 1 for index in spec["site2_target_indices_1based"]], dtype=int)
    hotspot_indices = np.asarray([int(index) - 1 for index in spec["hotspot_target_indices_1based"]], dtype=int)
    metrics: list[dict[str, Any]] = []
    for model_name in model_names:
        print(f"[Stage5B] running seed={seed:02d} {model_name}", flush=True)
        model.predict(
            num_recycles=int(spec["requested_recycles"]),
            dropout=True,
            models=[model_name],
            seed=seed,
            verbose=False,
        )
        safe_model = str(model_name).replace("/", "_")
        pdb_path = output_dir / f"{safe_model}_seed{seed:02d}.pdb"
        npz_path = output_dir / f"{safe_model}_seed{seed:02d}.npz"
        model.save_current_pdb(str(pdb_path))

        plddt = _validated_plddt_fraction(model.aux["plddt"])
        pae = np.asarray(model.aux["pae"], dtype=float)
        global_pae = _bidirectional_pae_values(pae, lengths[0])
        site_pae = _bidirectional_pae_values(pae, lengths[0], list(spec["site2_target_indices_1based"]))
        hotspot_pae = _bidirectional_pae_values(pae, lengths[0], list(spec["hotspot_target_indices_1based"]))
        log = dict(model.aux.get("log", {}))
        ptm_value = _scalar(log.get("ptm", ""))
        iptm_value = _scalar(log.get("i_ptm", log.get("iptm", "")))
        ranking_value = _scalar(log.get("multi", ""))
        ranking_source = "model_aux_log_multi"
        if ranking_value == "" and isinstance(ptm_value, float) and isinstance(iptm_value, float):
            ranking_value = round((0.8 * iptm_value) + (0.2 * ptm_value), 6)
            ranking_source = "derived_0.8_iPTM_plus_0.2_pTM"
        atom_positions = np.asarray(model.aux["atom_positions"], dtype=np.float32)
        np.savez_compressed(
            npz_path,
            plddt_fraction=plddt.astype(np.float16),
            pae=pae.astype(np.float16),
            atom_positions=atom_positions,
            target_length=np.array(lengths[0]),
            peptide_length=np.array(lengths[1]),
            site2_target_indices_1based=np.asarray(spec["site2_target_indices_1based"], dtype=np.int16),
            hotspot_target_indices_1based=np.asarray(spec["hotspot_target_indices_1based"], dtype=np.int16),
            seed=np.array(seed),
            model_name=np.array(str(model_name)),
        )
        metrics.append(
            {
                "stage5B_job_id": spec["stage5B_job_id"],
                "stage5B_candidate_id": spec["stage5B_candidate_id"],
                "protocol_hash": spec["protocol_hash"],
                "peptide_sequence_hash": spec["peptide_sequence_hash"],
                "target_template_sha1": spec["target_template_sha1"],
                "seed": seed,
                "model_name": model_name,
                "requested_recycles": spec["requested_recycles"],
                "forward_passes": spec["forward_passes"],
                "target_template_path": spec["target_template_path"],
                "target_template_residues_covered": template_audit["target_template_residues_covered"],
                "target_template_coverage_fraction": template_audit["target_template_coverage_fraction"],
                "peptide_template_residues_covered": 0,
                "peptide_template_atom_mask_sum": 0.0,
                "template_input_verified": "true",
                "plddt_mean_fraction": round(float(plddt.mean()), 5),
                "plddt_target_mean_fraction": round(float(plddt[target_indices].mean()), 5),
                "plddt_peptide_mean_fraction": round(float(plddt[lengths[0] :].mean()), 5),
                "plddt_site2_mean_fraction": round(float(plddt[site_indices].mean()), 5),
                "plddt_hotspot_mean_fraction": round(float(plddt[hotspot_indices].mean()), 5),
                "plddt_mean_100": round(float(plddt.mean()) * 100.0, 3),
                "plddt_target_mean_100": round(float(plddt[target_indices].mean()) * 100.0, 3),
                "plddt_peptide_mean_100": round(float(plddt[lengths[0] :].mean()) * 100.0, 3),
                "plddt_site2_mean_100": round(float(plddt[site_indices].mean()) * 100.0, 3),
                "plddt_hotspot_mean_100": round(float(plddt[hotspot_indices].mean()) * 100.0, 3),
                "ptm": ptm_value,
                "iptm": iptm_value,
                "ranking_confidence": ranking_value,
                "ranking_confidence_source": ranking_source,
                "interface_pae_global_mean_A": round(float(global_pae.mean()), 5),
                "interface_pae_site2_mean_A": round(float(site_pae.mean()), 5),
                "interface_pae_site2_median_A": round(float(np.median(site_pae)), 5),
                "interface_pae_hotspot_mean_A": round(float(hotspot_pae.mean()), 5),
                "target_msa_mode": spec["target_msa_mode"],
                "peptide_msa_mode": spec["peptide_msa_mode"],
                "msa_rows_input": int(msa.shape[0]),
                "use_mlm": "false",
                "use_dropout": "true",
                "prediction_pdb": pdb_path,
                "prediction_npz": npz_path,
                "validation_test_type": spec["validation_test_type"],
                "template_mode": spec["template_mode"],
                "use_initial_guess": "false",
                "cyclic_topology_encoding": spec["cyclic_topology_encoding"],
            }
        )

    _write_metrics(output_dir / "model_metrics.csv", metrics)
    metadata = {
        **dict(spec),
        "loaded_colabdesign_source": str(source_dir),
        "loaded_colabdesign_commit": loaded_commit,
        "template_input_verified": True,
        "target_template_residues_covered": template_audit["target_template_residues_covered"],
        "target_template_coverage_fraction": template_audit["target_template_coverage_fraction"],
        "peptide_template_residues_covered": 0,
        "peptide_template_atom_mask_sum": 0.0,
        "prediction_count": len(metrics),
        "model_names_completed": model_names,
        "plddt_source_scale": "fraction_0_1",
        "plddt_report_scale": "0_100",
        "pae_units": "angstrom",
        "target_template_loaded_by_prediction_runner": True,
        "peptide_template_loaded_by_prediction_runner": False,
        "reference_design_loaded_by_prediction_runner": False,
        "explicit_terminal_cn_bond_record_used": False,
    }
    with (output_dir / "run_metadata.json").open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(metadata, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(f"[Stage5B] completed {spec['stage5B_job_id']}: models={len(metrics)}", flush=True)


def main() -> int:
    parser = argparse.ArgumentParser(description="Run one Stage 5B target-structure-conditioned AfCycDesign seed job.")
    parser.add_argument("--job-spec")
    parser.add_argument("--af-params")
    parser.add_argument("--audit-imports", action="store_true")
    parser.add_argument("--preflight-only", action="store_true")
    args = parser.parse_args()

    if args.audit_imports:
        _audit_imports()
        return 0
    if not args.job_spec or not args.af_params:
        parser.error("--job-spec and --af-params are required")
    spec = _read_job_spec(Path(args.job_spec))
    af_params = Path(os.path.expandvars(os.path.expanduser(args.af_params)))
    if args.preflight_only:
        _preflight(spec, af_params)
        return 0
    if os.environ.get("RUN_STAGE5B_PREDICTIONS") != "YES":
        raise RuntimeError("Stage 5B prediction is review-gated. Set RUN_STAGE5B_PREDICTIONS=YES after review.")
    if not af_params.is_dir():
        raise RuntimeError(f"AlphaFold parameter directory does not exist: {af_params}")
    _run_prediction(spec, af_params)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
