from __future__ import annotations

import csv
import importlib.util
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
    candidate = Path(path)
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
