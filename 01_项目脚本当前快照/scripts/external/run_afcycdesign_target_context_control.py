from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
from typing import Any, Mapping


LEGAL_AA = set("ACDEFGHIKLMNPQRSTVWY")
AA3_TO_AA1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "MSE": "M", "PHE": "F",
    "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y",
    "VAL": "V",
}

CACHE_IDENTITY_FIELDS = (
    "stage5_context_job_id",
    "context_id",
    "protocol_hash",
    "protocol_version",
    "context_pdb_sha256",
    "context_mapping_csv_sha256",
    "source_cleaned_pdb_sha256",
    "chain_ids",
    "chain_lengths",
    "chain_sequences",
    "fga_chain_index",
    "site2_fga_indices_1based",
    "hotspot_fga_indices_1based",
    "template_mode",
    "use_templates",
    "use_batch_as_template",
    "template_sequence_masked",
    "template_sidechains_masked",
    "template_interchain_features_masked",
    "target_msa_mode",
    "msa_rows_expected",
    "peptide_included",
    "use_initial_guess",
    "use_dropout",
    "cyclic_topology_encoding",
    "validation_test_type",
    "seed",
    "requested_recycles",
    "forward_passes",
    "models_per_seed",
    "model_names_expected",
    "colabdesign_commit",
    "route_run_id",
    "route_manifest_sha256",
)

METRIC_FIELDS = [
    "stage5_context_job_id",
    "context_id",
    "protocol_hash",
    "seed",
    "model_name",
    "requested_recycles",
    "forward_passes",
    "template_mode",
    "template_sequence_masked",
    "template_sidechains_masked",
    "template_interchain_features_masked",
    "context_template_residues_covered",
    "context_template_coverage_fraction",
    "template_input_verified",
    "target_msa_mode",
    "msa_rows_input",
    "peptide_included",
    "use_initial_guess",
    "use_dropout",
    "plddt_mean_fraction",
    "plddt_fga_mean_fraction",
    "plddt_site2_mean_fraction",
    "plddt_hotspot_mean_fraction",
    "plddt_partner_mean_fraction",
    "plddt_mean_100",
    "plddt_fga_mean_100",
    "plddt_site2_mean_100",
    "plddt_hotspot_mean_100",
    "plddt_partner_mean_100",
    "ptm",
    "iptm",
    "ranking_confidence",
    "ranking_confidence_source",
    "context_interchain_pae_mean_A",
    "prediction_pdb",
    "prediction_npz",
    "validation_test_type",
]


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
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


def _pdb_chain_sequences(path: Path) -> tuple[list[str], list[str], list[int]]:
    order: list[str] = []
    residues: dict[str, list[str]] = {}
    seen: dict[str, set[tuple[str, str]]] = {}
    ca_seen: dict[str, set[tuple[str, str]]] = {}
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith("ATOM  "):
                continue
            chain = line[21].strip() or "_"
            if chain not in residues:
                order.append(chain)
                residues[chain] = []
                seen[chain] = set()
                ca_seen[chain] = set()
            key = (line[22:26].strip(), line[26].strip())
            if line[12:16].strip() == "CA":
                ca_seen[chain].add(key)
            if key in seen[chain]:
                continue
            seen[chain].add(key)
            resname = line[17:20].strip().upper()
            if resname not in AA3_TO_AA1:
                raise RuntimeError(f"Unsupported residue {resname} in {path}")
            residues[chain].append(AA3_TO_AA1[resname])
    sequences = ["".join(residues[chain]) for chain in order]
    ca_counts = [len(ca_seen[chain]) for chain in order]
    if ca_counts != [len(sequence) for sequence in sequences]:
        raise RuntimeError(f"Context PDB lacks complete CA coverage: chains={order}, CA={ca_counts}")
    return order, sequences, ca_counts


def _validate_context_mapping(path: Path, spec: Mapping[str, Any]) -> None:
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != sum(spec["chain_lengths"]):
        raise RuntimeError("Context mapping row count does not match total context length")
    expected_sequence_by_chain = dict(zip(spec["chain_ids"], spec["chain_sequences"]))
    observed_sequence_by_chain: dict[str, list[str]] = {chain: [] for chain in spec["chain_ids"]}
    observed_numbers: dict[str, list[int]] = {chain: [] for chain in spec["chain_ids"]}
    site_indices: list[int] = []
    hotspot_indices: list[int] = []
    for row in rows:
        chain = str(row.get("context_chain", ""))
        if chain not in observed_sequence_by_chain:
            raise RuntimeError(f"Unexpected chain {chain!r} in context mapping")
        observed_sequence_by_chain[chain].append(str(row.get("context_residue_aa", "")))
        observed_numbers[chain].append(int(row["context_residue_number"]))
        if chain == "A" and str(row.get("is_target_site_residue", "")).lower() == "true":
            site_indices.append(int(row["context_residue_number"]))
        if chain == "A" and str(row.get("is_selected_hotspot", "")).lower() == "true":
            hotspot_indices.append(int(row["context_residue_number"]))
    for chain, expected_sequence in expected_sequence_by_chain.items():
        if "".join(observed_sequence_by_chain[chain]) != expected_sequence:
            raise RuntimeError(f"Context mapping sequence mismatch for chain {chain}")
        if observed_numbers[chain] != list(range(1, len(expected_sequence) + 1)):
            raise RuntimeError(f"Context mapping numbering is not contiguous for chain {chain}")
    if site_indices != list(spec["site2_fga_indices_1based"]):
        raise RuntimeError("Context mapping Site_2 indices do not match the job spec")
    if hotspot_indices != list(spec["hotspot_fga_indices_1based"]):
        raise RuntimeError("Context mapping hotspot indices do not match the job spec")


def _read_job_spec(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        spec = json.load(handle)
    required = set(CACHE_IDENTITY_FIELDS) | {
        "context_pdb",
        "context_mapping_csv",
        "source_cleaned_pdb",
        "prediction_output_dir",
        "route_manifest_path",
    }
    missing = sorted(required - set(spec))
    if missing:
        raise RuntimeError(f"Target-context job spec is missing fields: {','.join(missing)}")
    if str(spec["protocol_version"]) != "stage5_target_context_controls_v1":
        raise RuntimeError("Unexpected target-context protocol version")
    if str(spec["validation_test_type"]) != "target_context_template_recovery_control":
        raise RuntimeError("Unexpected target-context validation_test_type")
    if not bool(spec["use_templates"]) or bool(spec["use_batch_as_template"]):
        raise RuntimeError("Target-context controls require explicit templates without use_batch_as_template")
    if bool(spec["peptide_included"]):
        raise RuntimeError("Target-context controls must not include a peptide")
    if bool(spec["use_initial_guess"]):
        raise RuntimeError("Target-context controls forbid an initial guess")
    if bool(spec["use_dropout"]):
        raise RuntimeError("Target-context v1 uses dropout=false for deterministic context comparison")
    if str(spec["target_msa_mode"]) != "single_sequence" or int(spec["msa_rows_expected"]) != 1:
        raise RuntimeError("Target-context v1 requires a one-row single-sequence MSA")
    if str(spec["cyclic_topology_encoding"]) != "not_applicable_no_peptide":
        raise RuntimeError("Cyclic topology must be marked not applicable when peptide is absent")
    if int(spec["requested_recycles"]) < 1:
        raise RuntimeError("requested_recycles must be positive")
    if int(spec["forward_passes"]) != int(spec["requested_recycles"]) + 1:
        raise RuntimeError("forward_passes must equal requested_recycles + 1")
    if int(spec["models_per_seed"]) != 5:
        raise RuntimeError("All five AlphaFold multimer-v3 parameter sets are required")
    expected_models = [f"model_{index}_multimer_v3" for index in range(1, 6)]
    if list(spec["model_names_expected"]) != expected_models:
        raise RuntimeError("model_names_expected is not the complete ordered multimer-v3 set")

    chain_ids = [str(value) for value in spec["chain_ids"]]
    chain_lengths = [int(value) for value in spec["chain_lengths"]]
    chain_sequences = [_validate_sequence(value, f"chain_sequences[{index}]") for index, value in enumerate(spec["chain_sequences"])]
    if not chain_ids or len(chain_ids) != len(chain_lengths) or len(chain_ids) != len(chain_sequences):
        raise RuntimeError("chain_ids, chain_lengths, and chain_sequences have inconsistent lengths")
    if len(set(chain_ids)) != len(chain_ids):
        raise RuntimeError("Context chain IDs are not unique")
    if [len(sequence) for sequence in chain_sequences] != chain_lengths:
        raise RuntimeError("Context chain sequence lengths do not match chain_lengths")
    if int(spec["fga_chain_index"]) != 0 or chain_ids[0] != "A":
        raise RuntimeError("The first context chain must be the FGA chain A")
    for field in ("site2_fga_indices_1based", "hotspot_fga_indices_1based"):
        values = [int(value) for value in spec[field]]
        if not values or min(values) < 1 or max(values) > chain_lengths[0]:
            raise RuntimeError(f"{field} contains invalid FGA indices")
        spec[field] = values
    if len(spec["hotspot_fga_indices_1based"]) != 4:
        raise RuntimeError("Exactly four Site_2 hotspots are required")

    context_pdb = Path(str(spec["context_pdb"]))
    mapping_csv = Path(str(spec["context_mapping_csv"]))
    source_pdb = Path(str(spec["source_cleaned_pdb"]))
    for label, source in (
        ("context_pdb", context_pdb),
        ("context_mapping_csv", mapping_csv),
        ("source_cleaned_pdb", source_pdb),
    ):
        if not source.is_file():
            raise RuntimeError(f"{label} does not exist: {source}")
    if _sha256_file(context_pdb) != str(spec["context_pdb_sha256"]):
        raise RuntimeError("context_pdb_sha256 mismatch")
    if _sha256_file(mapping_csv) != str(spec["context_mapping_csv_sha256"]):
        raise RuntimeError("context_mapping_csv_sha256 mismatch")
    if _sha256_file(source_pdb) != str(spec["source_cleaned_pdb_sha256"]):
        raise RuntimeError("source_cleaned_pdb_sha256 mismatch")
    observed_ids, observed_sequences, observed_ca = _pdb_chain_sequences(context_pdb)
    if observed_ids != chain_ids:
        raise RuntimeError(f"Context PDB chain order {observed_ids} != spec {chain_ids}")
    if observed_sequences != chain_sequences:
        raise RuntimeError("Context PDB sequences do not match the job spec")
    if observed_ca != chain_lengths:
        raise RuntimeError("Context PDB CA coverage does not match chain lengths")

    route_manifest = Path(str(spec["route_manifest_path"]))
    if not route_manifest.is_file():
        raise RuntimeError(f"Route manifest does not exist: {route_manifest}")
    if _sha256_file(route_manifest) != str(spec["route_manifest_sha256"]):
        raise RuntimeError("route_manifest_sha256 mismatch")
    manifest_payload = json.loads(route_manifest.read_text(encoding="utf-8"))
    if str(manifest_payload.get("run_id", "")) != str(spec["route_run_id"]):
        raise RuntimeError("route_run_id does not match the route manifest")

    spec["chain_ids"] = chain_ids
    spec["chain_lengths"] = chain_lengths
    spec["chain_sequences"] = chain_sequences
    _validate_context_mapping(mapping_csv, spec)
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


def _build_template_batch(spec: Mapping[str, Any]) -> tuple[dict[str, Any], dict[str, Any]]:
    import numpy as np
    from colabdesign.af.alphafold.common import residue_constants
    from colabdesign.af.prep import prep_pdb

    context_pdb = Path(str(spec["context_pdb"]))
    chain_selector = ",".join(spec["chain_ids"])
    template = prep_pdb(str(context_pdb), chain=chain_selector, ignore_missing=True)["batch"]
    expected_length = sum(int(value) for value in spec["chain_lengths"])
    observed_length = int(template["aatype"].shape[0])
    if observed_length != expected_length:
        raise RuntimeError(f"Template parser returned {observed_length} residues, expected {expected_length}")
    atom_mask = np.asarray(template["all_atom_mask"], dtype=float)
    atom_positions = np.asarray(template["all_atom_positions"], dtype=float)
    ca_index = residue_constants.atom_order["CA"]
    ca_coverage = int(np.count_nonzero(atom_mask[:, ca_index] > 0.5))
    if ca_coverage != expected_length:
        raise RuntimeError(f"Context template CA coverage is {ca_coverage}/{expected_length}")
    if not np.isfinite(atom_mask).all() or not np.isfinite(atom_positions).all():
        raise RuntimeError("Context template tensors contain non-finite values")
    return template, {
        "context_template_residues_covered": ca_coverage,
        "context_template_coverage_fraction": ca_coverage / expected_length,
        "template_input_verified": True,
    }


def _validated_plddt_fraction(value: Any) -> Any:
    import numpy as np

    array = np.asarray(value, dtype=float)
    if not array.size or not np.isfinite(array).all() or float(array.min()) < -1e-6 or float(array.max()) > 1.000001:
        raise RuntimeError("Expected finite model.aux['plddt'] values on the 0-1 scale")
    return array


def _scalar(value: Any) -> float | str:
    try:
        return round(float(value), 6)
    except (TypeError, ValueError):
        return ""


def _interchain_pae_mean(pae: Any, lengths: list[int]) -> float | str:
    import numpy as np

    if len(lengths) < 2:
        return ""
    values = []
    starts = np.cumsum([0] + lengths)
    for left in range(len(lengths)):
        for right in range(left + 1, len(lengths)):
            a = np.asarray(pae[starts[left] : starts[left + 1], starts[right] : starts[right + 1]], dtype=float).ravel()
            b = np.asarray(pae[starts[right] : starts[right + 1], starts[left] : starts[left + 1]], dtype=float).ravel()
            values.extend([a, b])
    combined = np.concatenate(values)
    return round(float(combined.mean()), 5)


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
    if metadata.get("loaded_colabdesign_commit") != loaded_commit:
        return False
    if metadata.get("loaded_colabdesign_source") != str(source_dir):
        return False
    if metadata.get("template_input_verified") is not True:
        return False
    if int(metadata.get("context_template_residues_covered", 0)) != sum(spec["chain_lengths"]):
        return False
    expected_models = set(spec["model_names_expected"])
    if len(metrics) != int(spec["models_per_seed"]) or {row.get("model_name", "") for row in metrics} != expected_models:
        return False
    for row in metrics:
        if row.get("stage5_context_job_id") != spec["stage5_context_job_id"]:
            return False
        if row.get("protocol_hash") != spec["protocol_hash"]:
            return False
        if not Path(str(row.get("prediction_pdb", ""))).is_file():
            return False
        if not Path(str(row.get("prediction_npz", ""))).is_file():
            return False
    return True


def _preflight(spec: Mapping[str, Any], af_params: Path) -> None:
    source_dir, loaded_commit = _verify_colabdesign_source(spec)
    if not af_params.is_dir():
        raise RuntimeError(f"AlphaFold parameter directory does not exist: {af_params}")
    import colabdesign

    loaded_module = Path(colabdesign.__file__).resolve()
    if source_dir != loaded_module and source_dir not in loaded_module.parents:
        raise RuntimeError(f"Imported colabdesign from {loaded_module}, expected {source_dir}")
    _, audit = _build_template_batch(spec)
    print(
        f"PASS: {spec['context_id']} chains={spec['chain_ids']} lengths={spec['chain_lengths']} "
        f"template_coverage={audit['context_template_residues_covered']}/{sum(spec['chain_lengths'])} "
        f"mode={spec['template_mode']} commit={loaded_commit}",
        flush=True,
    )


def _run_prediction(spec: Mapping[str, Any], af_params: Path) -> None:
    os.environ.setdefault("XLA_PYTHON_CLIENT_PREALLOCATE", "false")
    import numpy as np
    import colabdesign
    from colabdesign import clear_mem, mk_af_model
    from colabdesign.af.alphafold.common import residue_constants
    from colabdesign.af.contrib import predict

    source_dir, loaded_commit = _verify_colabdesign_source(spec)
    loaded_module = Path(colabdesign.__file__).resolve()
    if source_dir != loaded_module and source_dir not in loaded_module.parents:
        raise RuntimeError(f"Imported colabdesign from {loaded_module}, expected {source_dir}")
    output_dir = Path(str(spec["prediction_output_dir"]))
    output_dir.mkdir(parents=True, exist_ok=True)
    if _completed_output_is_valid(output_dir, spec, source_dir, loaded_commit):
        print(f"[Stage5 context] skipping cache-valid job {spec['stage5_context_job_id']}", flush=True)
        return

    lengths = list(spec["chain_lengths"])
    sequences = list(spec["chain_sequences"])
    sequence = "".join(sequences)
    seed = int(spec["seed"])
    template_batch, template_audit = _build_template_batch(spec)
    print(
        f"[Stage5 context] starting {spec['stage5_context_job_id']}: lengths={lengths}, "
        f"seed={seed}, recycles={spec['requested_recycles']}, mode={spec['template_mode']}",
        flush=True,
    )

    a3m_path = output_dir / "single_sequence.a3m"
    with a3m_path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(f">{spec['stage5_context_job_id']}\n{sequence}\n")
    msa, deletion_matrix = predict.parse_a3m(str(a3m_path))
    if int(msa.shape[0]) != 1 or int(msa.shape[1]) != len(sequence):
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
    model.set_opt(
        "template",
        rm_seq=bool(spec["template_sequence_masked"]),
        rm_sc=bool(spec["template_sidechains_masked"]),
        rm_ic=bool(spec["template_interchain_features_masked"]),
    )
    if list(model._lengths) != lengths:
        raise RuntimeError(f"Actual model chain lengths {model._lengths} != spec {lengths}")
    if bool(model._args.get("use_initial_guess")):
        raise RuntimeError("Actual model unexpectedly enabled an initial guess")
    if bool(model._args.get("use_batch_as_template")):
        raise RuntimeError("Actual model unexpectedly enabled use_batch_as_template")
    template_mask = np.asarray(model._inputs["template_all_atom_mask"])[0]
    ca_index = residue_constants.atom_order["CA"]
    if int(np.count_nonzero(template_mask[:, ca_index] > 0.5)) != sum(lengths):
        raise RuntimeError("Actual model template input lacks complete CA coverage")
    model_names = list(model._model_names)
    if model_names != list(spec["model_names_expected"]):
        raise RuntimeError(f"Loaded model sets do not match the job spec: {model_names}")

    fga_length = lengths[0]
    fga_indices = np.arange(fga_length, dtype=int)
    site_indices = np.asarray([int(value) - 1 for value in spec["site2_fga_indices_1based"]], dtype=int)
    hotspot_indices = np.asarray([int(value) - 1 for value in spec["hotspot_fga_indices_1based"]], dtype=int)
    partner_indices = np.arange(fga_length, sum(lengths), dtype=int)
    metrics: list[dict[str, Any]] = []
    for model_name in model_names:
        print(f"[Stage5 context] running seed={seed:02d} {model_name}", flush=True)
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
        plddt = _validated_plddt_fraction(model.aux["plddt"])
        pae = np.asarray(model.aux["pae"], dtype=float)
        atom_positions = np.asarray(model.aux["atom_positions"], dtype=np.float32)
        if len(plddt) != sum(lengths) or pae.shape != (sum(lengths), sum(lengths)):
            raise RuntimeError("Model output dimensions do not match context length")
        log = dict(model.aux.get("log", {}))
        ptm = _scalar(log.get("ptm", ""))
        iptm = _scalar(log.get("i_ptm", log.get("iptm", "")))
        ranking = _scalar(log.get("multi", ""))
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
            chain_lengths=np.asarray(lengths, dtype=np.int16),
            site2_fga_indices_1based=np.asarray(spec["site2_fga_indices_1based"], dtype=np.int16),
            hotspot_fga_indices_1based=np.asarray(spec["hotspot_fga_indices_1based"], dtype=np.int16),
            seed=np.array(seed),
            model_name=np.array(str(model_name)),
        )
        metrics.append(
            {
                "stage5_context_job_id": spec["stage5_context_job_id"],
                "context_id": spec["context_id"],
                "protocol_hash": spec["protocol_hash"],
                "seed": seed,
                "model_name": model_name,
                "requested_recycles": spec["requested_recycles"],
                "forward_passes": spec["forward_passes"],
                "template_mode": spec["template_mode"],
                "template_sequence_masked": str(spec["template_sequence_masked"]).lower(),
                "template_sidechains_masked": str(spec["template_sidechains_masked"]).lower(),
                "template_interchain_features_masked": str(spec["template_interchain_features_masked"]).lower(),
                "context_template_residues_covered": template_audit["context_template_residues_covered"],
                "context_template_coverage_fraction": template_audit["context_template_coverage_fraction"],
                "template_input_verified": "true",
                "target_msa_mode": spec["target_msa_mode"],
                "msa_rows_input": int(msa.shape[0]),
                "peptide_included": "false",
                "use_initial_guess": "false",
                "use_dropout": "false",
                "plddt_mean_fraction": round(float(plddt.mean()), 5),
                "plddt_fga_mean_fraction": round(float(plddt[fga_indices].mean()), 5),
                "plddt_site2_mean_fraction": round(float(plddt[site_indices].mean()), 5),
                "plddt_hotspot_mean_fraction": round(float(plddt[hotspot_indices].mean()), 5),
                "plddt_partner_mean_fraction": partner_fraction,
                "plddt_mean_100": round(float(plddt.mean()) * 100.0, 3),
                "plddt_fga_mean_100": round(float(plddt[fga_indices].mean()) * 100.0, 3),
                "plddt_site2_mean_100": round(float(plddt[site_indices].mean()) * 100.0, 3),
                "plddt_hotspot_mean_100": round(float(plddt[hotspot_indices].mean()) * 100.0, 3),
                "plddt_partner_mean_100": (
                    round(float(partner_fraction) * 100.0, 3) if isinstance(partner_fraction, float) else ""
                ),
                "ptm": ptm,
                "iptm": iptm,
                "ranking_confidence": ranking,
                "ranking_confidence_source": ranking_source,
                "context_interchain_pae_mean_A": _interchain_pae_mean(pae, lengths),
                "prediction_pdb": pdb_path,
                "prediction_npz": npz_path,
                "validation_test_type": spec["validation_test_type"],
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
        "context_template_loaded_by_prediction_runner": True,
        "peptide_loaded_by_prediction_runner": False,
        "cyclic_offset_applied": False,
    }
    with (output_dir / "run_metadata.json").open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(metadata, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(f"[Stage5 context] completed {spec['stage5_context_job_id']}: models={len(metrics)}", flush=True)


def main() -> int:
    parser = argparse.ArgumentParser(description="Run one expanded-target AfCycDesign recovery control.")
    parser.add_argument("--job-spec", required=True)
    parser.add_argument("--af-params", required=True)
    parser.add_argument("--preflight-only", action="store_true")
    args = parser.parse_args()
    spec = _read_job_spec(Path(args.job_spec))
    af_params = Path(os.path.expandvars(os.path.expanduser(args.af_params)))
    if args.preflight_only:
        _preflight(spec, af_params)
        return 0
    if os.environ.get("RUN_STAGE5_CONTEXT_CONTROLS") != "YES":
        raise RuntimeError(
            "Target-context prediction is review-gated. Set RUN_STAGE5_CONTEXT_CONTROLS=YES after review."
        )
    if not af_params.is_dir():
        raise RuntimeError(f"AlphaFold parameter directory does not exist: {af_params}")
    _run_prediction(spec, af_params)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
