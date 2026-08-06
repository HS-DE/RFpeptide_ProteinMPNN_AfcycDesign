from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
from typing import Any, Mapping

import run_afcycdesign_target_context_control as context_runner


PROTOCOL_VERSION = "stage5B_v2_C1_C3_full_target_template_top5_v1"
VALIDATION_TEST_TYPE = "target_context_structure_conditioned_recovery"
LEGAL_AA = set("ACDEFGHIKLMNPQRSTVWY")
STAGE4_IDENTITY_FIELDS = (
    "global_backbone_id",
    "backbone_family_id",
    "source_batch_label",
    "source_stage4_design_id",
    "source_stage4_scored_pdb_sha256",
    "stage4_run_id",
    "stage4_protocol_identity_sha256",
    "stage4_run_group_id",
    "stage4_scores_csv",
    "stage4_scores_csv_sha256",
    "stage4_top_candidates_csv",
    "stage4_top_candidates_csv_sha256",
    "stage4_validation_selection_rank",
)

CACHE_IDENTITY_FIELDS = (
    "stage5B_v2_job_id",
    "stage5B_v2_candidate_context_id",
    "stage5B_v2_candidate_id",
    "protocol_hash",
    "protocol_version",
    "stage5_campaign_id",
    "validation_test_type",
    "context_id",
    "context_pdb_sha256",
    "context_mapping_csv_sha256",
    "context_chain_ids",
    "context_chain_lengths",
    "context_chain_sequences",
    "stage0_crop_context_indices_1based",
    "site2_fga_indices_1based",
    "hotspot_fga_indices_1based",
    "fga_chain_index",
    "peptide_chain_id",
    "peptide_sequence",
    "peptide_sequence_length",
    "peptide_sequence_hash",
    "cyclic_chain_index",
    "cyclic_topology_encoding",
    "template_mode",
    "use_templates",
    "use_batch_as_template",
    "template_sequence_masked",
    "template_sidechains_masked",
    "template_interchain_features_masked",
    "template_padding_mode",
    "target_template_expected_coverage",
    "peptide_template_expected_coverage",
    "use_initial_guess",
    "reference_design_pdb_sha256",
    "reference_target_chain",
    "reference_peptide_chain",
    "reference_design_loaded_by_prediction_runner",
    "target_msa_mode",
    "peptide_msa_mode",
    "msa_rows_expected",
    "use_mlm",
    "use_dropout",
    "seed",
    "requested_recycles",
    "forward_passes",
    "models_per_seed",
    "model_names_expected",
    "colabdesign_commit",
    "route_run_id",
    "route_manifest_sha256",
    *STAGE4_IDENTITY_FIELDS,
)

METRIC_FIELDS = [
    "stage5B_v2_job_id",
    "stage5B_v2_candidate_context_id",
    "stage5B_v2_candidate_id",
    "context_id",
    "protocol_hash",
    "peptide_sequence_hash",
    "seed",
    "model_name",
    "requested_recycles",
    "forward_passes",
    "context_template_residues_covered",
    "context_template_coverage_fraction",
    "peptide_template_residues_covered",
    "peptide_template_atom_mask_sum",
    "template_input_verified",
    "template_mode",
    "template_sequence_masked",
    "template_sidechains_masked",
    "template_interchain_features_masked",
    "target_msa_mode",
    "peptide_msa_mode",
    "msa_rows_input",
    "use_mlm",
    "use_dropout",
    "plddt_mean_fraction",
    "plddt_context_mean_fraction",
    "plddt_fga_mean_fraction",
    "plddt_site2_mean_fraction",
    "plddt_hotspot_mean_fraction",
    "plddt_partner_mean_fraction",
    "plddt_peptide_mean_fraction",
    "plddt_mean_100",
    "plddt_context_mean_100",
    "plddt_fga_mean_100",
    "plddt_site2_mean_100",
    "plddt_hotspot_mean_100",
    "plddt_partner_mean_100",
    "plddt_peptide_mean_100",
    "ptm",
    "iptm",
    "ranking_confidence",
    "ranking_confidence_source",
    "interface_pae_context_mean_A",
    "interface_pae_fga_mean_A",
    "interface_pae_site2_mean_A",
    "interface_pae_site2_median_A",
    "interface_pae_hotspot_mean_A",
    "context_interchain_pae_mean_A",
    "prediction_pdb",
    "prediction_npz",
    "validation_test_type",
    "cyclic_topology_encoding",
]


def _sha1_text(value: str) -> str:
    return hashlib.sha1(value.encode("utf-8")).hexdigest()


def _validate_sequence(value: Any, label: str) -> str:
    sequence = str(value).strip().upper()
    if not sequence:
        raise RuntimeError(f"{label} is empty")
    invalid = sorted(set(sequence) - LEGAL_AA)
    if invalid:
        raise RuntimeError(f"{label} contains unsupported amino acids: {','.join(invalid)}")
    return sequence


def _as_int_list(value: Any, label: str) -> list[int]:
    try:
        values = [int(item) for item in value]
    except (TypeError, ValueError) as exc:
        raise RuntimeError(f"{label} is not an integer list") from exc
    if not values:
        raise RuntimeError(f"{label} is empty")
    return values


def _read_job_spec(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        spec = json.load(handle)
    required = set(CACHE_IDENTITY_FIELDS) | {
        "context_pdb",
        "context_mapping_csv",
        "reference_design_pdb_for_posthoc",
        "prediction_output_dir",
        "route_manifest_path",
    }
    missing = sorted(required - set(spec))
    if missing:
        raise RuntimeError(f"Stage 5B-v2 job spec is missing fields: {','.join(missing)}")
    if str(spec["protocol_version"]) != PROTOCOL_VERSION:
        raise RuntimeError("Unexpected Stage 5B-v2 protocol version")
    if str(spec["validation_test_type"]) != VALIDATION_TEST_TYPE:
        raise RuntimeError("Unexpected Stage 5B-v2 validation_test_type")
    if str(spec["context_id"]) not in {
        "C1_crop86_full_template",
        "C3_native_GHI301_full_template",
    }:
        raise RuntimeError("Only C1 and C3 are valid Stage 5B-v2 contexts")
    if str(spec["template_mode"]) != "target_context_full":
        raise RuntimeError("Stage 5B-v2 requires template_mode=target_context_full")
    if not bool(spec["use_templates"]) or bool(spec["use_batch_as_template"]):
        raise RuntimeError("Stage 5B-v2 requires an explicit target template only")
    for field in (
        "template_sequence_masked",
        "template_sidechains_masked",
        "template_interchain_features_masked",
    ):
        if bool(spec[field]):
            raise RuntimeError(f"Stage 5B-v2 full target template forbids {field}=true")
    if bool(spec["use_initial_guess"]):
        raise RuntimeError("Stage 5B-v2 forbids an initial guess")
    if bool(spec["reference_design_loaded_by_prediction_runner"]):
        raise RuntimeError("Stage 4 reference coordinates are posthoc-only")
    if int(spec["peptide_template_expected_coverage"]) != 0:
        raise RuntimeError("Peptide template coverage must be exactly zero")
    if str(spec["target_msa_mode"]) != "single_sequence" or str(spec["peptide_msa_mode"]) != "single_sequence":
        raise RuntimeError("Stage 5B-v2 uses single-sequence target and peptide inputs")
    if int(spec["msa_rows_expected"]) != 1 or bool(spec["use_mlm"]) or bool(spec["use_dropout"]):
        raise RuntimeError("Stage 5B-v2 requires one MSA row, MLM=false, and dropout=false")
    if int(spec["requested_recycles"]) != 6 or int(spec["forward_passes"]) != 7:
        raise RuntimeError("Stage 5B-v2 requires six requested recycles and seven forward passes")
    if int(spec["models_per_seed"]) != 5 or list(spec["model_names_expected"]) != [
        f"model_{index}_multimer_v3" for index in range(1, 6)
    ]:
        raise RuntimeError("Stage 5B-v2 requires all five multimer-v3 model parameter sets")

    chain_ids = [str(value) for value in spec["context_chain_ids"]]
    chain_lengths = _as_int_list(spec["context_chain_lengths"], "context_chain_lengths")
    chain_sequences = [
        _validate_sequence(value, f"context_chain_sequences[{index}]")
        for index, value in enumerate(spec["context_chain_sequences"])
    ]
    if len(chain_ids) != len(chain_lengths) or len(chain_ids) != len(chain_sequences):
        raise RuntimeError("Context chain IDs, lengths, and sequences have inconsistent sizes")
    if chain_ids != [chr(ord("A") + index) for index in range(len(chain_ids))]:
        raise RuntimeError("Context chains must be consecutive starting at A")
    if [len(sequence) for sequence in chain_sequences] != chain_lengths:
        raise RuntimeError("Context chain sequence lengths do not match context_chain_lengths")
    if int(spec["fga_chain_index"]) != 0:
        raise RuntimeError("FGA must be context chain index 0")
    expected_peptide_chain = chr(ord("A") + len(chain_ids))
    if str(spec["peptide_chain_id"]) != expected_peptide_chain:
        raise RuntimeError("Peptide chain ID does not follow the context chains")
    if int(spec["cyclic_chain_index"]) != len(chain_ids):
        raise RuntimeError("Cyclic offset index must identify only the peptide chain")
    if str(spec["cyclic_topology_encoding"]) != "peptide_chain_relative_position_cyclic_offset":
        raise RuntimeError("Unexpected cyclic topology encoding")

    peptide = _validate_sequence(spec["peptide_sequence"], "peptide_sequence")
    if len(peptide) != int(spec["peptide_sequence_length"]):
        raise RuntimeError("peptide_sequence_length does not match peptide_sequence")
    if _sha1_text(peptide)[:8] != str(spec["peptide_sequence_hash"]):
        raise RuntimeError("peptide_sequence_hash mismatch")
    crop = _as_int_list(spec["stage0_crop_context_indices_1based"], "stage0_crop_context_indices_1based")
    site = _as_int_list(spec["site2_fga_indices_1based"], "site2_fga_indices_1based")
    hotspots = _as_int_list(spec["hotspot_fga_indices_1based"], "hotspot_fga_indices_1based")
    if len(crop) != 86 or len(hotspots) != 4:
        raise RuntimeError("Stage 5B-v2 requires the 86-aa crop and four hotspots")
    if min(crop + site + hotspots) < 1 or max(crop + site + hotspots) > chain_lengths[0]:
        raise RuntimeError("Crop/Site_2/hotspot indices lie outside FGA chain A")

    context_pdb = Path(str(spec["context_pdb"]))
    mapping_csv = Path(str(spec["context_mapping_csv"]))
    if not context_pdb.is_file() or not mapping_csv.is_file():
        raise RuntimeError("Context PDB or mapping CSV is missing")
    if context_runner._sha256_file(context_pdb) != str(spec["context_pdb_sha256"]):
        raise RuntimeError("context_pdb_sha256 mismatch")
    if context_runner._sha256_file(mapping_csv) != str(spec["context_mapping_csv_sha256"]):
        raise RuntimeError("context_mapping_csv_sha256 mismatch")
    observed_ids, observed_sequences, observed_ca = context_runner._pdb_chain_sequences(context_pdb)
    if observed_ids != chain_ids or observed_sequences != chain_sequences or observed_ca != chain_lengths:
        raise RuntimeError("Context PDB chain identity, sequence, or CA coverage mismatch")
    context_runner._validate_context_mapping(
        mapping_csv,
        {
            "chain_ids": chain_ids,
            "chain_lengths": chain_lengths,
            "chain_sequences": chain_sequences,
            "site2_fga_indices_1based": site,
            "hotspot_fga_indices_1based": hotspots,
        },
    )
    if int(spec["target_template_expected_coverage"]) != sum(chain_lengths):
        raise RuntimeError("target_template_expected_coverage mismatch")

    route_manifest = Path(str(spec["route_manifest_path"]))
    if not route_manifest.is_file():
        raise RuntimeError(f"Route manifest does not exist: {route_manifest}")
    if context_runner._sha256_file(route_manifest) != str(spec["route_manifest_sha256"]):
        raise RuntimeError("route_manifest_sha256 mismatch")
    payload = json.loads(route_manifest.read_text(encoding="utf-8"))
    if str(payload.get("run_id")) != str(spec["route_run_id"]):
        raise RuntimeError("route_run_id does not match route manifest")

    # The reference path is identity-bound for collection, but this runner
    # deliberately does not open it or load any of its coordinates.
    if not str(spec["reference_design_pdb_for_posthoc"]).strip() or not str(spec["reference_design_pdb_sha256"]).strip():
        raise RuntimeError("Posthoc reference identity is incomplete")
    if (
        not str(spec["reference_target_chain"]).strip()
        or not str(spec["reference_peptide_chain"]).strip()
        or str(spec["reference_target_chain"]) == str(spec["reference_peptide_chain"])
    ):
        raise RuntimeError("Posthoc reference chain roles are incomplete or ambiguous")
    spec["context_chain_ids"] = chain_ids
    spec["context_chain_lengths"] = chain_lengths
    spec["context_chain_sequences"] = chain_sequences
    spec["stage0_crop_context_indices_1based"] = crop
    spec["site2_fga_indices_1based"] = site
    spec["hotspot_fga_indices_1based"] = hotspots
    spec["peptide_sequence"] = peptide
    return spec


def _build_padded_template(spec: Mapping[str, Any]) -> tuple[dict[str, Any], dict[str, Any]]:
    import numpy as np
    from colabdesign.af.alphafold.common import residue_constants

    context_batch, context_audit = context_runner._build_template_batch(
        {
            "context_pdb": spec["context_pdb"],
            "chain_ids": spec["context_chain_ids"],
            "chain_lengths": spec["context_chain_lengths"],
        }
    )
    context_length = sum(spec["context_chain_lengths"])
    peptide_length = int(spec["peptide_sequence_length"])
    total_length = context_length + peptide_length
    aatype = np.full((total_length,), 21, dtype=np.int32)
    positions = np.zeros((total_length,) + tuple(context_batch["all_atom_positions"].shape[1:]), dtype=np.float32)
    mask = np.zeros((total_length,) + tuple(context_batch["all_atom_mask"].shape[1:]), dtype=np.float32)
    aatype[:context_length] = np.asarray(context_batch["aatype"], dtype=np.int32)
    positions[:context_length] = np.asarray(context_batch["all_atom_positions"], dtype=np.float32)
    mask[:context_length] = np.asarray(context_batch["all_atom_mask"], dtype=np.float32)
    ca_index = residue_constants.atom_order["CA"]
    context_coverage = int(np.count_nonzero(mask[:context_length, ca_index] > 0.5))
    peptide_mask_sum = float(mask[context_length:].sum())
    peptide_position_count = int(np.count_nonzero(positions[context_length:]))
    peptide_unknown_count = int(np.count_nonzero(aatype[context_length:] == 21))
    if context_coverage != context_length:
        raise RuntimeError(f"Context template CA coverage is {context_coverage}/{context_length}")
    if peptide_mask_sum != 0.0 or peptide_position_count != 0 or peptide_unknown_count != peptide_length:
        raise RuntimeError("Peptide template coordinates, masks, or residue identities leaked into the template")
    if not np.isfinite(positions).all() or not np.isfinite(mask).all():
        raise RuntimeError("Padded template contains non-finite values")
    return {
        "aatype": aatype,
        "all_atom_positions": positions,
        "all_atom_mask": mask,
    }, {
        **context_audit,
        "template_total_length": total_length,
        "context_template_residues_covered": context_coverage,
        "context_template_coverage_fraction": context_coverage / context_length,
        "peptide_template_residues_covered": 0,
        "peptide_template_atom_mask_sum": peptide_mask_sum,
        "peptide_template_position_nonzero_count": peptide_position_count,
        "peptide_template_unknown_aatype_count": peptide_unknown_count,
        "template_input_verified": True,
    }


def _verify_cyclic_offset(model: Any, lengths: list[int], cyclic_chain_index: int) -> None:
    import numpy as np

    if list(model._lengths) != lengths:
        raise RuntimeError(f"Unexpected model chain lengths: {model._lengths}")
    offset = np.asarray(model._inputs.get("offset"))
    residue_index = np.asarray(model._inputs["residue_index"])
    linear = residue_index[:, None] - residue_index[None, :]
    if offset.shape != linear.shape:
        raise RuntimeError("Cyclic offset matrix shape is invalid")
    starts = np.cumsum([0] + lengths)
    start, end = int(starts[cyclic_chain_index]), int(starts[cyclic_chain_index + 1])
    outside = np.ones(offset.shape, dtype=bool)
    outside[start:end, start:end] = False
    if not np.array_equal(offset[outside], linear[outside]):
        raise RuntimeError("Cyclic offset modified a target or target-peptide block")
    peptide_block = offset[start:end, start:end]
    if lengths[cyclic_chain_index] > 2 and np.array_equal(peptide_block, linear[start:end, start:end]):
        raise RuntimeError("Peptide block was not cyclically encoded")
    if lengths[cyclic_chain_index] > 1 and abs(int(peptide_block[0, -1])) != 1:
        raise RuntimeError("Peptide termini are not adjacent in the cyclic offset matrix")


def _verify_cyclic_offset_static(spec: Mapping[str, Any]) -> None:
    import numpy as np
    from colabdesign.af.contrib.cyclic import add_cyclic_offset

    lengths = list(spec["context_chain_lengths"]) + [int(spec["peptide_sequence_length"])]

    class DummyModel:
        pass

    model = DummyModel()
    model._lengths = lengths
    model._inputs = {"residue_index": np.concatenate([np.arange(length) for length in lengths])}
    add_cyclic_offset(model, [int(spec["cyclic_chain_index"])])
    _verify_cyclic_offset(model, lengths, int(spec["cyclic_chain_index"]))


def _bidirectional_pae_values(pae: Any, left_indices: Any, right_indices: Any) -> Any:
    import numpy as np

    array = np.asarray(pae, dtype=float)
    left = np.asarray(left_indices, dtype=int)
    right = np.asarray(right_indices, dtype=int)
    if array.ndim != 2 or array.shape[0] != array.shape[1] or not len(left) or not len(right):
        raise RuntimeError(f"Cannot compute interface PAE from shape {array.shape}")
    return np.concatenate(
        [array[np.ix_(left, right)].reshape(-1), array[np.ix_(right, left)].reshape(-1)]
    )


def _context_interchain_pae_mean(pae: Any, lengths: list[int]) -> float | str:
    import numpy as np

    if len(lengths) < 2:
        return ""
    starts = np.cumsum([0] + lengths)
    values = []
    for left in range(len(lengths)):
        for right in range(left + 1, len(lengths)):
            left_idx = np.arange(starts[left], starts[left + 1], dtype=int)
            right_idx = np.arange(starts[right], starts[right + 1], dtype=int)
            values.append(_bidirectional_pae_values(pae, left_idx, right_idx))
    return round(float(np.concatenate(values).mean()), 5)


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
    if any(metadata.get(field) != spec.get(field) for field in CACHE_IDENTITY_FIELDS):
        return False
    if metadata.get("loaded_colabdesign_commit") != loaded_commit or metadata.get("loaded_colabdesign_source") != str(source_dir):
        return False
    if metadata.get("template_input_verified") is not True:
        return False
    if int(metadata.get("context_template_residues_covered", 0)) != sum(spec["context_chain_lengths"]):
        return False
    if int(metadata.get("peptide_template_residues_covered", -1)) != 0:
        return False
    if metadata.get("reference_design_loaded_by_prediction_runner") is not False:
        return False
    expected_models = set(spec["model_names_expected"])
    if len(metrics) != int(spec["models_per_seed"]) or {row.get("model_name") for row in metrics} != expected_models:
        return False
    for row in metrics:
        if row.get("stage5B_v2_job_id") != spec["stage5B_v2_job_id"]:
            return False
        if row.get("protocol_hash") != spec["protocol_hash"]:
            return False
        if not Path(str(row.get("prediction_pdb", ""))).is_file() or not Path(str(row.get("prediction_npz", ""))).is_file():
            return False
    return True


def _preflight(spec: Mapping[str, Any], af_params: Path) -> None:
    source_dir, loaded_commit = context_runner._verify_colabdesign_source(spec)
    if not af_params.is_dir():
        raise RuntimeError(f"AlphaFold parameter directory does not exist: {af_params}")
    import colabdesign

    loaded_module = Path(colabdesign.__file__).resolve()
    if source_dir != loaded_module and source_dir not in loaded_module.parents:
        raise RuntimeError(f"Imported colabdesign from {loaded_module}, expected {source_dir}")
    _, audit = _build_padded_template(spec)
    _verify_cyclic_offset_static(spec)
    print(
        f"PASS: {spec['stage5B_v2_candidate_context_id']} context_coverage="
        f"{audit['context_template_residues_covered']}/{sum(spec['context_chain_lengths'])}, "
        f"peptide_coverage=0, cyclic_chain_index={spec['cyclic_chain_index']}, commit={loaded_commit}",
        flush=True,
    )


def _run_prediction(spec: Mapping[str, Any], af_params: Path) -> None:
    os.environ.setdefault("XLA_PYTHON_CLIENT_PREALLOCATE", "false")
    import numpy as np
    import colabdesign
    from colabdesign import clear_mem, mk_af_model
    from colabdesign.af.alphafold.common import residue_constants
    from colabdesign.af.contrib import predict
    from colabdesign.af.contrib.cyclic import add_cyclic_offset

    source_dir, loaded_commit = context_runner._verify_colabdesign_source(spec)
    loaded_module = Path(colabdesign.__file__).resolve()
    if source_dir != loaded_module and source_dir not in loaded_module.parents:
        raise RuntimeError(f"Imported colabdesign from {loaded_module}, expected pinned source {source_dir}")
    output_dir = Path(str(spec["prediction_output_dir"]))
    output_dir.mkdir(parents=True, exist_ok=True)
    if _completed_output_is_valid(output_dir, spec, source_dir, loaded_commit):
        print(f"[Stage5B-v2] skipping cache-valid job {spec['stage5B_v2_job_id']}", flush=True)
        return

    context_lengths = list(spec["context_chain_lengths"])
    context_sequences = list(spec["context_chain_sequences"])
    peptide = str(spec["peptide_sequence"])
    lengths = context_lengths + [len(peptide)]
    sequence = "".join(context_sequences) + peptide
    context_length = sum(context_lengths)
    peptide_indices = np.arange(context_length, context_length + len(peptide), dtype=int)
    fga_indices = np.arange(context_lengths[0], dtype=int)
    site_indices = np.asarray([value - 1 for value in spec["site2_fga_indices_1based"]], dtype=int)
    hotspot_indices = np.asarray([value - 1 for value in spec["hotspot_fga_indices_1based"]], dtype=int)
    partner_indices = np.arange(context_lengths[0], context_length, dtype=int)
    seed = int(spec["seed"])
    template_batch, template_audit = _build_padded_template(spec)
    _verify_cyclic_offset_static(spec)
    print(
        f"[Stage5B-v2] starting {spec['stage5B_v2_job_id']}: context={context_lengths}, "
        f"peptide={len(peptide)}, seed={seed}, recycles={spec['requested_recycles']}",
        flush=True,
    )

    a3m_path = output_dir / "single_sequence.a3m"
    with a3m_path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(f">{spec['stage5B_v2_job_id']}\n{sequence}\n")
    msa, deletion_matrix = predict.parse_a3m(str(a3m_path))
    if tuple(msa.shape) != (1, len(sequence)):
        raise RuntimeError(f"Unexpected single-sequence MSA shape: {msa.shape}")

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
    model.set_opt("template", rm_seq=False, rm_sc=False, rm_ic=False)
    add_cyclic_offset(model, [int(spec["cyclic_chain_index"])])
    _verify_cyclic_offset(model, lengths, int(spec["cyclic_chain_index"]))
    if bool(model._args.get("use_initial_guess")) or bool(model._args.get("use_batch_as_template")):
        raise RuntimeError("Actual model enabled a forbidden initial guess or batch template")
    template_mask = np.asarray(model._inputs["template_all_atom_mask"])[0]
    ca_index = residue_constants.atom_order["CA"]
    if int(np.count_nonzero(template_mask[:context_length, ca_index] > 0.5)) != context_length:
        raise RuntimeError("Actual model input lacks full target-context CA template coverage")
    if int(np.count_nonzero(template_mask[context_length:])) != 0:
        raise RuntimeError("Actual model input contains peptide template atom masks")
    model_names = list(model._model_names)
    if model_names != list(spec["model_names_expected"]):
        raise RuntimeError(f"Loaded model parameter sets do not match the spec: {model_names}")

    metrics: list[dict[str, Any]] = []
    for model_name in model_names:
        print(f"[Stage5B-v2] running seed={seed:02d} {model_name}", flush=True)
        model.predict(
            num_recycles=int(spec["requested_recycles"]),
            dropout=False,
            models=[model_name],
            seed=seed,
            verbose=False,
        )
        safe_model = str(model_name).replace("/", "_")
        pdb_path = output_dir / f"{safe_model}_seed{seed:02d}.pdb"
        npz_path = output_dir / f"{safe_model}_seed{seed:02d}.npz"
        model.save_current_pdb(str(pdb_path))
        plddt = context_runner._validated_plddt_fraction(model.aux["plddt"])
        pae = np.asarray(model.aux["pae"], dtype=float)
        atom_positions = np.asarray(model.aux["atom_positions"], dtype=np.float32)
        if len(plddt) != len(sequence) or pae.shape != (len(sequence), len(sequence)):
            raise RuntimeError("Model output dimensions do not match the context-plus-peptide input")
        context_indices = np.arange(context_length, dtype=int)
        context_pae = _bidirectional_pae_values(pae, context_indices, peptide_indices)
        fga_pae = _bidirectional_pae_values(pae, fga_indices, peptide_indices)
        site_pae = _bidirectional_pae_values(pae, site_indices, peptide_indices)
        hotspot_pae = _bidirectional_pae_values(pae, hotspot_indices, peptide_indices)
        log = dict(model.aux.get("log", {}))
        ptm = context_runner._scalar(log.get("ptm", ""))
        iptm = context_runner._scalar(log.get("i_ptm", log.get("iptm", "")))
        ranking = context_runner._scalar(log.get("multi", ""))
        ranking_source = "model_aux_log_multi"
        if ranking == "" and isinstance(ptm, float) and isinstance(iptm, float):
            ranking = round((0.8 * iptm) + (0.2 * ptm), 6)
            ranking_source = "derived_0.8_iPTM_plus_0.2_pTM"
        partner_fraction: float | str = ""
        if len(partner_indices):
            partner_fraction = round(float(plddt[partner_indices].mean()), 5)
        np.savez_compressed(
            npz_path,
            plddt_fraction=plddt.astype(np.float16),
            pae=pae.astype(np.float16),
            atom_positions=atom_positions,
            context_chain_lengths=np.asarray(context_lengths, dtype=np.int16),
            peptide_length=np.asarray(len(peptide), dtype=np.int16),
            stage0_crop_context_indices_1based=np.asarray(spec["stage0_crop_context_indices_1based"], dtype=np.int16),
            site2_fga_indices_1based=np.asarray(spec["site2_fga_indices_1based"], dtype=np.int16),
            hotspot_fga_indices_1based=np.asarray(spec["hotspot_fga_indices_1based"], dtype=np.int16),
            seed=np.asarray(seed),
            model_name=np.asarray(str(model_name)),
        )
        metrics.append(
            {
                "stage5B_v2_job_id": spec["stage5B_v2_job_id"],
                "stage5B_v2_candidate_context_id": spec["stage5B_v2_candidate_context_id"],
                "stage5B_v2_candidate_id": spec["stage5B_v2_candidate_id"],
                "context_id": spec["context_id"],
                "protocol_hash": spec["protocol_hash"],
                "peptide_sequence_hash": spec["peptide_sequence_hash"],
                "seed": seed,
                "model_name": model_name,
                "requested_recycles": spec["requested_recycles"],
                "forward_passes": spec["forward_passes"],
                "context_template_residues_covered": template_audit["context_template_residues_covered"],
                "context_template_coverage_fraction": template_audit["context_template_coverage_fraction"],
                "peptide_template_residues_covered": 0,
                "peptide_template_atom_mask_sum": 0.0,
                "template_input_verified": "true",
                "template_mode": spec["template_mode"],
                "template_sequence_masked": "false",
                "template_sidechains_masked": "false",
                "template_interchain_features_masked": "false",
                "target_msa_mode": spec["target_msa_mode"],
                "peptide_msa_mode": spec["peptide_msa_mode"],
                "msa_rows_input": int(msa.shape[0]),
                "use_mlm": "false",
                "use_dropout": "false",
                "plddt_mean_fraction": round(float(plddt.mean()), 5),
                "plddt_context_mean_fraction": round(float(plddt[context_indices].mean()), 5),
                "plddt_fga_mean_fraction": round(float(plddt[fga_indices].mean()), 5),
                "plddt_site2_mean_fraction": round(float(plddt[site_indices].mean()), 5),
                "plddt_hotspot_mean_fraction": round(float(plddt[hotspot_indices].mean()), 5),
                "plddt_partner_mean_fraction": partner_fraction,
                "plddt_peptide_mean_fraction": round(float(plddt[peptide_indices].mean()), 5),
                "plddt_mean_100": round(float(plddt.mean()) * 100.0, 3),
                "plddt_context_mean_100": round(float(plddt[context_indices].mean()) * 100.0, 3),
                "plddt_fga_mean_100": round(float(plddt[fga_indices].mean()) * 100.0, 3),
                "plddt_site2_mean_100": round(float(plddt[site_indices].mean()) * 100.0, 3),
                "plddt_hotspot_mean_100": round(float(plddt[hotspot_indices].mean()) * 100.0, 3),
                "plddt_partner_mean_100": round(float(partner_fraction) * 100.0, 3) if isinstance(partner_fraction, float) else "",
                "plddt_peptide_mean_100": round(float(plddt[peptide_indices].mean()) * 100.0, 3),
                "ptm": ptm,
                "iptm": iptm,
                "ranking_confidence": ranking,
                "ranking_confidence_source": ranking_source,
                "interface_pae_context_mean_A": round(float(context_pae.mean()), 5),
                "interface_pae_fga_mean_A": round(float(fga_pae.mean()), 5),
                "interface_pae_site2_mean_A": round(float(site_pae.mean()), 5),
                "interface_pae_site2_median_A": round(float(np.median(site_pae)), 5),
                "interface_pae_hotspot_mean_A": round(float(hotspot_pae.mean()), 5),
                "context_interchain_pae_mean_A": _context_interchain_pae_mean(pae, context_lengths),
                "prediction_pdb": pdb_path,
                "prediction_npz": npz_path,
                "validation_test_type": spec["validation_test_type"],
                "cyclic_topology_encoding": spec["cyclic_topology_encoding"],
            }
        )

    _write_metrics(output_dir / "model_metrics.csv", metrics)
    metadata = {
        **dict(spec),
        "loaded_colabdesign_source": str(source_dir),
        "loaded_colabdesign_commit": loaded_commit,
        **template_audit,
        "prediction_count": len(metrics),
        "model_names_completed": model_names,
        "plddt_source_scale": "fraction_0_1",
        "plddt_report_scale": "0_100",
        "pae_units": "angstrom",
        "target_context_template_loaded_by_prediction_runner": True,
        "peptide_template_loaded_by_prediction_runner": False,
        "reference_design_loaded_by_prediction_runner": False,
        "initial_guess_loaded_by_prediction_runner": False,
        "cyclic_offset_applied": True,
        "cyclic_offset_chain_index": int(spec["cyclic_chain_index"]),
        "explicit_terminal_cn_bond_record_used": False,
    }
    with (output_dir / "run_metadata.json").open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(metadata, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(f"[Stage5B-v2] completed {spec['stage5B_v2_job_id']}: models={len(metrics)}", flush=True)


def main() -> int:
    parser = argparse.ArgumentParser(description="Run one Stage 5B-v2 C1/C3 target-context-conditioned prediction job.")
    parser.add_argument("--job-spec", required=True)
    parser.add_argument("--af-params", required=True)
    parser.add_argument("--preflight-only", action="store_true")
    args = parser.parse_args()
    spec = _read_job_spec(Path(args.job_spec))
    af_params = Path(os.path.expandvars(os.path.expanduser(args.af_params)))
    if args.preflight_only:
        _preflight(spec, af_params)
        return 0
    if os.environ.get("RUN_STAGE5B_V2_PREDICTIONS") != "YES":
        raise RuntimeError("Stage 5B-v2 is review-gated. Set RUN_STAGE5B_V2_PREDICTIONS=YES after review.")
    if not af_params.is_dir():
        raise RuntimeError(f"AlphaFold parameter directory does not exist: {af_params}")
    _run_prediction(spec, af_params)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
