from __future__ import annotations

import csv
import copy
import hashlib
import importlib.util
import json
import logging
import os
import shutil
import sys
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional


DEFAULT_CONFIG: Dict[str, Any] = {
    "project": {
        "name": "fga_cyclic_peptide_design",
        "target_gene": "FGA",
        "target_uniprot": "P02671",
        "target_description": "Human fibrinogen alpha chain; design should focus on FGA regions exposed in native fibrinogen.",
    },
    "input": {
        "excel_file": "data/input/高丰度蛋白信息.xlsx",
        "gene_column": "Gene",
        "uniprot_column": "UniprotID",
        "sequence_column": "Sequence",
        "abundance_column": "estimated_ng_per_ml",
    },
    "target_regions": {
        "full_length": {
            "start": 1,
            "end": 866,
            "use_for_design": False,
            "priority": "record",
            "note": "Full-length precursor; keep for record only.",
        },
        "extracellular": {
            "start": 20,
            "end": 866,
            "use_for_design": True,
            "priority": "secondary",
            "note": "Signal peptide removed.",
        },
        "main_chain": {
            "start": 36,
            "end": 866,
            "use_for_design": True,
            "preferred": True,
            "priority": "preferred",
            "note": "Avoids fibrinopeptide A as primary target.",
        },
    },
    "structures": {
        "primary_pdb": "3GHG",
        "primary_pdb_file": "data/structures/raw/3GHG.pdb",
        "cleaned_pdb_file": "data/structures/prepared/fibrinogen_3GHG_clean.pdb",
        "alphafold_pdb_file": "data/structures/raw/AF-P02671-F1-model_v4.pdb",
        "prefer_native_complex": True,
    },
    "peptide_design": {
        "scheme": "A",
        "cyclization": "Cys-Cys disulfide",
        "final_format_prefix": "Biotin-PEG4-GSG-",
        "final_format_suffix": "-NH2",
        "core_length_min": 10,
        "core_length_max": 18,
        "preferred_core_lengths": [12, 14, 16],
        "terminal_residue": "C",
        "forbid_internal_cys": True,
    },
    "generation": {
        "total_raw_designs_target": 5000,
        "patches": {
            "Patch_A": {"description": "Stable visible exposed FGA surface in 3GHG", "n_designs": 2000, "priority": "high"},
            "Patch_B": {"description": "FGA 36-200 visible exposed surface", "n_designs": 2000, "priority": "medium"},
            "Patch_C": {"description": "FGA C-terminal / alphaC-related exploratory region", "n_designs": 1000, "priority": "exploratory"},
        },
        "length_distribution": {10: 0.10, 12: 0.30, 14: 0.30, 16: 0.20, 18: 0.10},
    },
    "sequence_filters": {
        "max_hydrophobic_run": 4,
        "net_charge_min": -3,
        "net_charge_max": 3,
        "max_w_count": 1,
        "max_m_count": 1,
        "forbid_low_complexity": True,
        "forbid_poly_basic": True,
        "forbid_poly_acidic": True,
    },
    "complex_prediction": {
        "run_prediction": True,
        "prediction_engines": ["colabfold", "boltz2_optional"],
        "seeds_per_candidate": 5,
        "require_patch_consistency": True,
    },
    "scoring_thresholds": {
        "max_interface_pae": 10.0,
        "min_peptide_plddt": 70.0,
        "min_interface_contacts": 8,
        "min_iptm_soft": 0.50,
        "min_iptm_preferred": 0.65,
        "require_cys_geometry_pass": True,
    },
    "negative_screen": {
        "enabled": True,
        "targets": ["ALB", "APOA1", "TF", "A2M", "C3", "IGG_FC"],
        "purpose": "Remove obviously sticky non-specific peptides.",
    },
    "ranking": {"top_n_candidates": 50, "top_n_synthesis_priority": 10},
    "report": {"language": "zh-CN", "include_warnings": True, "do_not_claim_experimental_validation": True},
}


REQUIRED_DIRS = [
    "config",
    "data/input",
    "data/structures/raw",
    "data/structures/prepared",
    "data/annotations",
    "scripts",
    "notebooks",
    "tests",
    "results/raw_designs",
    "results/raw_designs/colabdesign_jobs",
    "results/raw_designs/rfdiffusion_jobs_optional",
    "results/complex_predictions",
    "results/filtered",
    "results/final",
    "logs",
    "env",
]


ACTIVE_ROUTE_DENY_MARKERS = (
    "_archived_invalid_",
    "03_旧错误路线生成脚本",
    "04_旧n20假preflight样例_禁止运行",
    "旧错误路线",
    "禁止运行",
)

ROUTE_NAME = "rfpeptides_head_to_tail"
ROUTE_PROTOCOL_VERSION = "rfpeptides_head_to_tail_v1"
HOTSPOT_MAPPING_VERSION = "site2_stage0_trb_tensor_pdb_complex_global_v1"
ROUTE_CYCLIZATION = "head_to_tail_amide"
ROUTE_PROVENANCE_FIELDS = [
    "route_run_id",
    "route_batch_id",
    "route_manifest_path",
    "route_manifest_sha256",
    "route_protocol_version",
    "hotspot_mapping_version",
    "route_cyclization",
]
SOURCE_ROUTE_PROVENANCE_FIELDS = [
    "source_run_id",
    "source_batch_id",
    "source_route_manifest",
    "source_route_manifest_sha256",
]


def project_root() -> Path:
    return Path(__file__).resolve().parents[1]


def resolve_path(path: str | Path, root: Optional[Path] = None) -> Path:
    p = Path(path)
    if p.is_absolute():
        return p
    return (root or project_root()) / p


def _mixed_path_variants(path: str | Path) -> list[str]:
    text = str(path).strip().replace("\\", "/")
    variants = [text]
    if len(text) >= 3 and text[1] == ":" and text[2] == "/":
        variants.append(f"/mnt/{text[0].lower()}{text[2:]}")
    elif text.startswith("/mnt/") and len(text) > 7 and text[6] == "/":
        variants.append(f"{text[5].upper()}:{text[6:]}")
    return variants


def _check_active_route_markers(path: str | Path, label: str) -> None:
    for variant in _mixed_path_variants(path):
        folded = variant.casefold()
        for marker in ACTIVE_ROUTE_DENY_MARKERS:
            if marker.casefold() in folded:
                raise RuntimeError(
                    f"Active-route input {label} is blocked: path={path}; "
                    f"normalized={variant}; matched_marker={marker}"
                )


def assert_active_route_path(
    path: str | Path,
    label: str,
    *,
    must_exist: bool = True,
) -> Path:
    """Resolve a production path and reject archived or explicitly forbidden locations."""

    _check_active_route_markers(path, label)
    text = str(path).strip().replace("\\", "/")
    if os.name == "nt" and text.startswith("/mnt/") and len(text) > 7 and text[6] == "/":
        text = f"{text[5].upper()}:{text[6:]}"
    elif os.name != "nt" and len(text) >= 3 and text[1] == ":" and text[2] == "/":
        text = f"/mnt/{text[0].lower()}{text[2:]}"
    candidate = Path(text)
    if not candidate.is_absolute():
        candidate = resolve_path(candidate)
    if must_exist and not candidate.exists():
        raise FileNotFoundError(f"Missing active-route input {label}: {candidate}")

    if candidate.exists():
        resolved = candidate.resolve(strict=True)
    else:
        parent = candidate.parent
        missing_parts = [candidate.name]
        while not parent.exists() and parent != parent.parent:
            missing_parts.append(parent.name)
            parent = parent.parent
        resolved_parent = parent.resolve(strict=True) if parent.exists() else parent.resolve()
        resolved = resolved_parent.joinpath(*reversed(missing_parts))
    _check_active_route_markers(resolved, label)
    return resolved


def sha256_file(path: str | Path) -> str:
    source = assert_active_route_path(path, "SHA-256 source file")
    if not source.is_file():
        raise RuntimeError(f"SHA-256 source is not a file: {source}")
    digest = hashlib.sha256()
    with source.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def canonical_json_sha256(value: Any) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _route_identity_payload(manifest: Mapping[str, Any]) -> dict[str, Any]:
    payload = copy.deepcopy({key: value for key, value in manifest.items() if key not in {"created_at", "run_id"}})
    payload.pop("project_config", None)
    for site in payload.get("stage0_sites", []):
        if isinstance(site, dict):
            site.pop("target_pdb", None)
            site.pop("mapping_csv", None)
    for source in payload.get("source_route_manifests", []):
        if isinstance(source, dict):
            source.pop("manifest_path", None)
            source.pop("path", None)
    return payload


def finalize_route_manifest(payload: Mapping[str, Any]) -> dict[str, Any]:
    manifest = dict(payload)
    manifest.setdefault("schema_version", 1)
    manifest.setdefault("route_name", ROUTE_NAME)
    manifest.setdefault("route_protocol_version", ROUTE_PROTOCOL_VERSION)
    manifest.setdefault("hotspot_mapping_version", HOTSPOT_MAPPING_VERSION)
    manifest.setdefault("production_route", True)
    manifest.setdefault("cyclization", ROUTE_CYCLIZATION)
    identity_sha256 = canonical_json_sha256(_route_identity_payload(manifest))
    manifest["run_id"] = f"rfp_{identity_sha256[:16]}"
    manifest.setdefault("created_at", datetime.now().astimezone().isoformat(timespec="seconds"))
    return manifest


def write_route_manifest(run_root: str | Path, payload: Mapping[str, Any]) -> tuple[Path, dict[str, Any], str]:
    root = assert_active_route_path(run_root, "route manifest output root", must_exist=False)
    root.mkdir(parents=True, exist_ok=True)
    manifest = finalize_route_manifest(payload)
    path = root / "route_manifest.json"
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True, ensure_ascii=False)
        handle.write("\n")
    return path, manifest, sha256_file(path)


def validate_route_manifest(manifest_path: str | Path, manifest: Mapping[str, Any]) -> None:
    required = {
        "schema_version",
        "route_name",
        "route_protocol_version",
        "hotspot_mapping_version",
        "production_route",
        "cyclization",
        "run_id",
        "batch_id",
        "site_labels",
        "protocol_peptide_length_min",
        "protocol_peptide_length_max",
        "run_peptide_length_min",
        "run_peptide_length_max",
        "num_designs_requested",
        "project_config",
        "project_config_sha256",
        "effective_project_config_sha256",
        "stage0_sites",
    }
    missing = sorted(required - set(manifest))
    if missing:
        raise RuntimeError(f"Route manifest is missing required fields {missing}: {manifest_path}")
    if manifest.get("schema_version") != 1:
        raise RuntimeError(f"Unsupported route manifest schema: {manifest.get('schema_version')}")
    expected_fixed = {
        "route_name": ROUTE_NAME,
        "route_protocol_version": ROUTE_PROTOCOL_VERSION,
        "hotspot_mapping_version": HOTSPOT_MAPPING_VERSION,
        "production_route": True,
        "cyclization": ROUTE_CYCLIZATION,
    }
    for key, expected in expected_fixed.items():
        if manifest.get(key) != expected:
            raise RuntimeError(f"Route manifest {key}={manifest.get(key)!r}, expected {expected!r}")
    expected_run_id = f"rfp_{canonical_json_sha256(_route_identity_payload(manifest))[:16]}"
    if manifest.get("run_id") != expected_run_id:
        raise RuntimeError(f"Route manifest run_id mismatch: {manifest.get('run_id')} != {expected_run_id}")

    site_labels = manifest.get("site_labels")
    stage0_sites = manifest.get("stage0_sites")
    if not isinstance(site_labels, list) or len(site_labels) != 1:
        raise RuntimeError("Current production route requires exactly one site label")
    if not isinstance(stage0_sites, list) or len(stage0_sites) != 1:
        raise RuntimeError("Current production route requires exactly one Stage 0 site provenance record")
    if stage0_sites[0].get("site_label") != site_labels[0]:
        raise RuntimeError("Route manifest Stage 0 site label does not match site_labels")

    protocol_min = int(manifest["protocol_peptide_length_min"])
    protocol_max = int(manifest["protocol_peptide_length_max"])
    run_min = int(manifest["run_peptide_length_min"])
    run_max = int(manifest["run_peptide_length_max"])
    if (protocol_min, protocol_max) != (12, 24):
        raise RuntimeError(f"Active route protocol length range must be 12-24, got {protocol_min}-{protocol_max}")
    if run_min < protocol_min or run_max > protocol_max or run_max < run_min:
        raise RuntimeError(f"Run peptide length range {run_min}-{run_max} is outside protocol range")
    if int(manifest["num_designs_requested"]) < 1:
        raise RuntimeError("Route manifest num_designs_requested must be positive")

    config_path = assert_active_route_path(manifest["project_config"], "route manifest project config")
    if sha256_file(config_path) != str(manifest["project_config_sha256"]):
        raise RuntimeError("Route manifest project config SHA-256 mismatch")
    if not str(manifest["effective_project_config_sha256"]).strip():
        raise RuntimeError("Route manifest effective project config SHA-256 is empty")

    site = stage0_sites[0]
    for key in [
        "target_pdb",
        "target_pdb_sha256",
        "mapping_csv",
        "mapping_csv_sha256",
        "normalized_hotspots",
        "normalized_hotspots_sha256",
    ]:
        if key not in site:
            raise RuntimeError(f"Route manifest Stage 0 provenance is missing {key}")
    target_pdb = assert_active_route_path(site["target_pdb"], "route manifest Stage 0 target PDB")
    mapping_csv = assert_active_route_path(site["mapping_csv"], "route manifest Stage 0 mapping CSV")
    if sha256_file(target_pdb) != str(site["target_pdb_sha256"]):
        raise RuntimeError("Route manifest Stage 0 target PDB SHA-256 mismatch")
    if sha256_file(mapping_csv) != str(site["mapping_csv_sha256"]):
        raise RuntimeError("Route manifest Stage 0 mapping CSV SHA-256 mismatch")
    hotspots = site["normalized_hotspots"]
    if not isinstance(hotspots, list) or not hotspots:
        raise RuntimeError("Route manifest normalized hotspot list is empty or invalid")
    if canonical_json_sha256(hotspots) != str(site["normalized_hotspots_sha256"]):
        raise RuntimeError("Route manifest normalized hotspot SHA-256 mismatch")


def load_route_manifest(run_root: str | Path) -> tuple[Path, dict[str, Any], str]:
    root = assert_active_route_path(run_root, "route manifest run root")
    path = assert_active_route_path(root / "route_manifest.json", "route manifest JSON")
    try:
        manifest = json.loads(path.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, UnicodeDecodeError) as exc:
        raise RuntimeError(f"Invalid route manifest JSON: {path}") from exc
    if not isinstance(manifest, dict):
        raise RuntimeError(f"Route manifest root must be an object: {path}")
    validate_route_manifest(path, manifest)
    return path, manifest, sha256_file(path)


def route_provenance_fields(
    manifest_path: str | Path,
    manifest: Mapping[str, Any],
    manifest_sha256: str,
) -> dict[str, str]:
    return {
        "route_run_id": str(manifest["run_id"]),
        "route_batch_id": str(manifest["batch_id"]),
        "route_manifest_path": str(manifest_path),
        "route_manifest_sha256": str(manifest_sha256),
        "route_protocol_version": str(manifest["route_protocol_version"]),
        "hotspot_mapping_version": str(manifest["hotspot_mapping_version"]),
        "route_cyclization": str(manifest["cyclization"]),
    }


def validate_row_route_provenance(row: Mapping[str, Any], expected: Mapping[str, str], label: str) -> None:
    for key in ROUTE_PROVENANCE_FIELDS:
        observed = str(row.get(key, "")).strip()
        if key == "route_manifest_path":
            observed_path = assert_active_route_path(observed, f"{label} route manifest")
            observed_digest = sha256_file(observed_path)
            row_digest = str(row.get("route_manifest_sha256", "")).strip()
            if observed_digest != row_digest:
                raise RuntimeError(
                    f"{label} route manifest file SHA-256 mismatch: observed={observed_digest!r}, "
                    f"row={row_digest!r}"
                )
            continue
        if observed != str(expected[key]):
            raise RuntimeError(
                f"{label} route provenance mismatch for {key}: observed={observed!r}, expected={expected[key]!r}"
            )


def validate_source_route_provenance(row: Mapping[str, Any], label: str) -> None:
    missing = [field for field in SOURCE_ROUTE_PROVENANCE_FIELDS if not str(row.get(field, "")).strip()]
    if missing:
        raise RuntimeError(f"{label} is missing source route provenance fields: {missing}")
    manifest_path = assert_active_route_path(row["source_route_manifest"], f"{label} source route manifest")
    actual_sha256 = sha256_file(manifest_path)
    recorded_sha256 = str(row["source_route_manifest_sha256"]).strip()
    if actual_sha256 != recorded_sha256:
        raise RuntimeError(
            f"{label} source route manifest SHA-256 mismatch: actual={actual_sha256}, recorded={recorded_sha256}"
        )
    _, manifest, loaded_sha256 = load_route_manifest(manifest_path.parent)
    if loaded_sha256 != recorded_sha256:
        raise RuntimeError(f"{label} source route manifest changed during validation")
    if str(manifest["run_id"]) != str(row["source_run_id"]):
        raise RuntimeError(f"{label} source run_id does not match its source route manifest")
    if str(manifest["batch_id"]) != str(row["source_batch_id"]):
        raise RuntimeError(f"{label} source batch_id does not match its source route manifest")


def add_route_provenance(rows: Iterable[dict[str, Any]], provenance: Mapping[str, str]) -> None:
    for row in rows:
        row.update(provenance)


def deep_update(base: Dict[str, Any], override: Mapping[str, Any]) -> Dict[str, Any]:
    out = dict(base)
    for key, value in override.items():
        if isinstance(value, Mapping) and isinstance(out.get(key), Mapping):
            out[key] = deep_update(dict(out[key]), value)
        else:
            out[key] = value
    return out


def load_config(config_path: str | Path) -> Dict[str, Any]:
    path = resolve_path(config_path)
    if not path.exists():
        raise FileNotFoundError(f"缺少配置文件: {path}")
    try:
        import yaml  # type: ignore

        with path.open("r", encoding="utf-8") as handle:
            loaded = yaml.safe_load(handle) or {}
        return deep_update(DEFAULT_CONFIG, loaded)
    except ModuleNotFoundError:
        # 当前基础系统可能还没有 pyyaml；完整环境中会按 env/environment.yml 安装。
        # 为了允许 prepare 流程继续生成目录和基础输出，这里回退到内置默认配置。
        return DEFAULT_CONFIG


def load_active_route_config(config_path: str | Path) -> tuple[Dict[str, Any], str, str]:
    path = assert_active_route_path(config_path, "active route project config")
    try:
        import yaml  # type: ignore
    except ModuleNotFoundError as exc:
        raise RuntimeError("PyYAML is required for active-route production configuration") from exc
    try:
        with path.open("r", encoding="utf-8") as handle:
            loaded = yaml.safe_load(handle)
    except Exception as exc:
        raise RuntimeError(f"Could not parse active-route project config: {path}") from exc
    if not isinstance(loaded, dict):
        raise RuntimeError("Active-route project config must be a YAML mapping")

    expected = {
        ("status",): "active",
        ("production_use",): True,
        ("route", "name"): ROUTE_NAME,
        ("route", "route_protocol_version"): ROUTE_PROTOCOL_VERSION,
        ("route", "hotspot_mapping_version"): HOTSPOT_MAPPING_VERSION,
        ("peptide_design", "cyclization"): ROUTE_CYCLIZATION,
        ("peptide_design", "length_min"): 12,
        ("peptide_design", "length_max"): 24,
        ("peptide_design", "required_terminal_residue"): None,
        ("peptide_design", "require_cys_pair"): False,
    }
    for keys, required_value in expected.items():
        value: Any = loaded
        for key in keys:
            if not isinstance(value, Mapping) or key not in value:
                raise RuntimeError(f"Active-route project config is missing {'.'.join(keys)}")
            value = value[key]
        if value != required_value:
            raise RuntimeError(
                f"Active-route project config {'.'.join(keys)}={value!r}, expected {required_value!r}"
            )

    shared_paths = [
        ("project", "name"),
        ("project", "target_gene"),
        ("project", "target_uniprot"),
        ("input", "excel_file"),
        ("target_regions", "main_chain", "start"),
        ("target_regions", "main_chain", "end"),
        ("structures", "cleaned_pdb_file"),
        ("pipeline", "stages"),
    ]
    for keys in shared_paths:
        value = loaded
        for key in keys:
            if not isinstance(value, Mapping) or key not in value:
                raise RuntimeError(f"Active-route project config is missing shared field {'.'.join(keys)}")
            value = value[key]
        if value is None or value == "":
            raise RuntimeError(f"Active-route shared field {'.'.join(keys)} is empty")

    peptide = loaded["peptide_design"]
    if peptide.get("terminal_residue") not in {None, ""}:
        raise RuntimeError("Active head-to-tail route must not require a terminal residue")
    if "disulfide" in str(peptide.get("cyclization", "")).casefold():
        raise RuntimeError("Active head-to-tail route must not use disulfide cyclization")
    scoring = loaded.get("scoring_thresholds", {})
    if isinstance(scoring, Mapping) and scoring.get("require_cys_geometry_pass") not in {None, False}:
        raise RuntimeError("Active head-to-tail route must disable Cys-specific geometry gating")
    omit_aas = str(peptide.get("omit_aas", ""))
    if omit_aas != "CX":
        raise RuntimeError("Active route currently requires peptide_design.omit_aas='CX'")

    raw_sha256 = sha256_file(path)
    effective_sha256 = canonical_json_sha256(loaded)
    return dict(loaded), raw_sha256, effective_sha256


def validate_route_project_config(config_path: str | Path, manifest: Mapping[str, Any]) -> Dict[str, Any]:
    config, raw_sha256, effective_sha256 = load_active_route_config(config_path)
    if raw_sha256 != str(manifest.get("project_config_sha256", "")):
        raise RuntimeError("Active project config file SHA-256 does not match route manifest")
    if effective_sha256 != str(manifest.get("effective_project_config_sha256", "")):
        raise RuntimeError("Active effective project config SHA-256 does not match route manifest")
    manifest_config = assert_active_route_path(manifest.get("project_config", ""), "manifest project config")
    supplied_config = assert_active_route_path(config_path, "supplied active project config")
    if manifest_config.resolve() != supplied_config.resolve():
        raise RuntimeError("Supplied active project config path differs from route manifest")
    return config


def ensure_project_dirs(root: Optional[Path] = None) -> None:
    base = root or project_root()
    for rel in REQUIRED_DIRS:
        (base / rel).mkdir(parents=True, exist_ok=True)


def setup_logger(name: str) -> logging.Logger:
    ensure_project_dirs()
    logger = logging.getLogger(name)
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    formatter = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s")
    log_path = project_root() / "logs" / f"{name}.log"
    file_handler = logging.FileHandler(log_path, encoding="utf-8")
    file_handler.setFormatter(formatter)
    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)
    return logger


def import_available(module_name: str) -> bool:
    return importlib.util.find_spec(module_name) is not None


def executable_available(name: str) -> bool:
    return shutil.which(name) is not None


def clean_sequence(seq: str) -> str:
    return "".join(ch for ch in str(seq).upper() if ch.isalpha())


def write_fasta(path: str | Path, header: str, sequence: str, width: int = 70) -> None:
    out = resolve_path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    sequence = clean_sequence(sequence)
    with out.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(f">{header}\n")
        for i in range(0, len(sequence), width):
            handle.write(sequence[i : i + width] + "\n")


def read_fasta_sequence(path: str | Path) -> str:
    fasta = resolve_path(path)
    if not fasta.exists():
        raise FileNotFoundError(f"缺少 FASTA 文件: {fasta}")
    return clean_sequence("".join(line.strip() for line in fasta.read_text(encoding="utf-8").splitlines() if not line.startswith(">")))


def read_fasta_header(path: str | Path) -> str:
    fasta = resolve_path(path)
    if not fasta.exists():
        raise FileNotFoundError(f"缺少 FASTA 文件: {fasta}")
    for line in fasta.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            return line[1:].strip()
    return ""


def write_csv(path: str | Path, rows: Iterable[Mapping[str, Any]], fieldnames: List[str]) -> None:
    out = resolve_path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})


def read_csv(path: str | Path) -> List[Dict[str, str]]:
    src = resolve_path(path)
    if not src.exists():
        return []
    with src.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_markdown(path: str | Path, text: str) -> None:
    out = resolve_path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(text.rstrip() + "\n", encoding="utf-8")


def append_run_header(logger: logging.Logger, script_name: str) -> None:
    logger.info("脚本: %s", script_name)
    logger.info("运行时间: %s", datetime.now().isoformat(timespec="seconds"))
    logger.info("项目目录: %s", project_root())


def bool_text(value: Any) -> str:
    return "true" if bool(value) else "false"


def rows_to_markdown(rows: List[Mapping[str, Any]], columns: List[str], empty_text: str) -> str:
    if not rows:
        return empty_text
    header = "| " + " | ".join(columns) + " |"
    sep = "| " + " | ".join(["---"] * len(columns)) + " |"
    lines = [header, sep]
    for row in rows:
        lines.append("| " + " | ".join(str(row.get(col, "")) for col in columns) + " |")
    return "\n".join(lines)
