from __future__ import annotations

import argparse
import importlib.util
import math
import os
import shutil
from collections import Counter
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

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
    validate_route_project_config,
    validate_row_route_provenance,
    write_csv,
    write_markdown,
    write_route_manifest,
)


SCRIPT_DIR = Path(__file__).resolve().parent
STAGE3C_PATH = SCRIPT_DIR / "23_collect_proteinmpnn_sequences.py"
HYDROPHOBIC_AAS = set("AVILMFWY")
AROMATIC_AAS = set("FWY")
BASIC_AAS = set("KR")
ACIDIC_AAS = set("DE")


STAGE4_FIELDS = [
    "stage4_design_id",
    "source_stage3d1_design_id",
    "source_sequence_design_id",
    "backbone_id",
    "site_label",
    "site_id",
    "input_pdb",
    "scored_pdb",
    "peptide_sequence",
    "peptide_length",
    "net_charge_approx",
    "hydrophobic_fraction",
    "aromatic_count",
    "pro_count",
    "gly_count",
    "sequence_liability_notes",
    "sequence_liability_count",
    "score_mode",
    "pyrosetta_score_status",
    "complex_rosetta_total_score",
    "separated_rosetta_total_score",
    "ddg_proxy_no_repack",
    "ddg_proxy_status",
    "CMS",
    "SAP",
    "target_contact_status",
    "target_site_recovery_status",
    "hotspot_recovery_status",
    "macrocycle_geometry_status",
    "clash_status",
    "num_target_contacts",
    "num_target_site_contacts",
    "num_hotspot_contacts",
    "peptide_target_min_distance",
    "peptide_site_min_distance",
    "peptide_hotspot_min_distance",
    "closest_target_residue",
    "closest_site_residue",
    "closest_hotspot_residue",
    "macrocycle_terminal_cn_distance",
    "peptide_radius_of_gyration",
    "detached_or_collapsed_flag",
    "stage4_energy_rank",
    "stage4_priority_rank",
    "stage4_priority_class",
    "stage4_priority_notes",
    "pass_stage4_qc",
    "stage4_failure_reasons",
    "notes",
] + ROUTE_PROVENANCE_FIELDS

TOP_VALIDATION_FIELDS = [
    "stage4_priority_rank",
    "stage4_priority_class",
    "stage4_design_id",
    "peptide_sequence",
    "scored_pdb",
    "ddg_proxy_no_repack",
    "site_contacts",
    "hotspot_contacts",
    "macrocycle_terminal_cn_distance",
    "clash_status",
    "sequence_liability_notes",
    "reason_selected",
] + ROUTE_PROVENANCE_FIELDS


def _load_stage3c_module() -> Any:
    spec = importlib.util.spec_from_file_location("stage3c_qc", STAGE3C_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not import Stage 3C QC helper: {STAGE3C_PATH}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _split_csv(value: str) -> list[str]:
    return [item.strip() for item in str(value or "").split(",") if item.strip()]


def _safe_token(value: str) -> str:
    keep = []
    for ch in str(value):
        if ch.isalnum() or ch in {"_", "-", "."}:
            keep.append(ch)
        else:
            keep.append("_")
    return "".join(keep).strip("_") or "item"


def _resolve_mixed_path(value: str | Path) -> Path:
    text = str(value).strip().replace("\\", "/")
    if not text:
        return resolve_path(text)
    if os.name == "nt" and text.startswith("/mnt/") and len(text) > 7 and text[6] == "/":
        return Path(f"{text[5].upper()}:/{text[7:]}")
    if os.name != "nt" and len(text) >= 3 and text[1] == ":" and text[2] == "/":
        return Path(f"/mnt/{text[0].lower()}{text[2:]}")
    path = Path(text)
    if path.is_absolute():
        return path
    return resolve_path(path)


def _read_required_csv(path: Path) -> list[dict[str, str]]:
    rows = read_csv(path)
    if not rows:
        raise RuntimeError(f"Missing or empty CSV: {path}")
    return rows


def _lookup_rows(rows: Iterable[Mapping[str, str]], keys: Sequence[str]) -> dict[str, dict[str, str]]:
    lookup: dict[str, dict[str, str]] = {}
    for row in rows:
        item = dict(row)
        for key in keys:
            value = item.get(key, "")
            if value:
                lookup[value] = item
    return lookup


def _parse_float(value: Any) -> float | None:
    try:
        if str(value).strip() == "":
            return None
        out = float(str(value).strip())
    except ValueError:
        return None
    if math.isnan(out) or math.isinf(out):
        return None
    return out


def _parse_int(value: Any, default: int = 0) -> int:
    try:
        return int(str(value).strip())
    except ValueError:
        return default


def _sequence_properties(sequence: str) -> dict[str, Any]:
    seq = "".join(aa for aa in str(sequence).upper() if aa.isalpha())
    length = len(seq)
    if not seq:
        return {
            "peptide_length": 0,
            "net_charge_approx": "",
            "hydrophobic_fraction": "",
            "aromatic_count": 0,
            "pro_count": 0,
            "gly_count": 0,
            "sequence_liability_notes": "missing_sequence",
            "sequence_liability_count": 1,
        }

    hydrophobic_count = sum(1 for aa in seq if aa in HYDROPHOBIC_AAS)
    aromatic_count = sum(1 for aa in seq if aa in AROMATIC_AAS)
    pro_count = seq.count("P")
    gly_count = seq.count("G")
    net_charge = sum(1.0 for aa in seq if aa in BASIC_AAS)
    net_charge += 0.1 * seq.count("H")
    net_charge -= sum(1.0 for aa in seq if aa in ACIDIC_AAS)
    hydrophobic_fraction = hydrophobic_count / length

    notes = []
    if hydrophobic_fraction >= 0.60:
        notes.append("warning_high_hydrophobic_fraction")
    if pro_count >= max(4, math.ceil(length * 0.30)):
        notes.append("warning_high_pro_fraction")
    if gly_count >= max(4, math.ceil(length * 0.30)):
        notes.append("warning_high_gly_fraction")
    most_common_count = max(Counter(seq).values())
    if most_common_count >= max(5, math.ceil(length * 0.45)):
        notes.append("warning_low_complexity_single_residue_dominant")
    if any(seq[idx] == seq[idx + 1] == seq[idx + 2] == seq[idx + 3] for idx in range(max(0, length - 3))):
        notes.append("warning_low_complexity_four_residue_run")

    return {
        "peptide_length": length,
        "net_charge_approx": round(net_charge, 3),
        "hydrophobic_fraction": round(hydrophobic_fraction, 3),
        "aromatic_count": aromatic_count,
        "pro_count": pro_count,
        "gly_count": gly_count,
        "sequence_liability_notes": ";".join(notes) if notes else "none",
        "sequence_liability_count": len(notes),
    }


def _load_pyrosetta() -> Any:
    try:
        import pyrosetta  # type: ignore
    except ModuleNotFoundError as exc:
        raise RuntimeError(
            "PyRosetta is required for Stage 4 score-only interface scoring. "
            "Run in the proteinmpnn_binder_design environment."
        ) from exc
    return pyrosetta


def _init_pyrosetta(pyrosetta: Any) -> None:
    pyrosetta.init("-beta_nov16 -mute all -use_terminal_residues true")


def _xyz_tuple(xyz: Any) -> tuple[float, float, float]:
    try:
        return float(xyz.x), float(xyz.y), float(xyz.z)
    except TypeError:
        return float(xyz[0]), float(xyz[1]), float(xyz[2])


def _pose_chain(pose: Any, pose_idx: int) -> str:
    pdb_info = pose.pdb_info()
    chain = pdb_info.chain(pose_idx) if pdb_info is not None else ""
    return str(chain).strip() or "_"


def _chain_positions(pose: Any, chain_id: str) -> list[int]:
    return [idx for idx in range(1, pose.total_residue() + 1) if _pose_chain(pose, idx) == chain_id]


def _translate_positions(pyrosetta: Any, pose: Any, positions: Sequence[int], separation_distance: float) -> None:
    atom_id_cls = pyrosetta.rosetta.core.id.AtomID
    xyz_vector_cls = pyrosetta.rosetta.numeric.xyzVector_double_t
    for res_idx in positions:
        residue = pose.residue(res_idx)
        for atom_idx in range(1, residue.natoms() + 1):
            atom_id = atom_id_cls(atom_idx, res_idx)
            x, y, z = _xyz_tuple(pose.xyz(atom_id))
            pose.set_xyz(atom_id, xyz_vector_cls(x + separation_distance, y, z))


def _score_complex_and_separated(
    *,
    pyrosetta: Any,
    input_pdb: Path,
    peptide_chain: str,
    separation_distance: float,
) -> dict[str, Any]:
    pose = pyrosetta.pose_from_pdb(str(input_pdb))
    peptide_positions = _chain_positions(pose, peptide_chain)
    if not peptide_positions:
        raise RuntimeError(f"Peptide chain {peptide_chain} was not found in {input_pdb}")

    scorefxn = pyrosetta.get_fa_scorefxn()
    complex_score = float(scorefxn(pose))
    separated_pose = pose.clone()
    _translate_positions(pyrosetta, separated_pose, peptide_positions, separation_distance)
    separated_score = float(scorefxn(separated_pose))
    return {
        "complex_rosetta_total_score": round(complex_score, 3),
        "separated_rosetta_total_score": round(separated_score, 3),
        "ddg_proxy_no_repack": round(complex_score - separated_score, 3),
    }


def _qc_scored_pdb(
    *,
    qchelper: Any,
    scored_pdb: Path,
    backbone_row: Mapping[str, str],
    site_numbers: set[str],
    hotspot_numbers: set[str],
    args: argparse.Namespace,
) -> dict[str, Any]:
    return qchelper._qc_row_for_relaxed_pdb(
        relaxed_pdb=scored_pdb,
        backbone_row=backbone_row,
        site_numbers=site_numbers,
        hotspot_numbers=hotspot_numbers,
        contact_cutoff=args.contact_cutoff,
        site_near_distance=args.site_near_distance,
        hotspot_near_distance=args.hotspot_near_distance,
        severe_clash_distance=args.severe_clash_distance,
        min_target_contacts=args.min_target_contacts,
        min_site_contacts=args.min_site_contacts,
        min_hotspot_contacts=args.min_hotspot_contacts,
        macrocycle_pass_distance=args.macrocycle_pass_distance,
        macrocycle_warn_distance=args.macrocycle_warn_distance,
        forbidden_aas=args.forbidden_aas,
    )


def _detached_or_collapsed_flag(qc: Mapping[str, Any], min_peptide_rog: float) -> str:
    flags = []
    if qc.get("target_contact_status") not in {"target_contact_pass", "target_contact_low_count"}:
        flags.append(str(qc.get("target_contact_status", "detached_from_target_crop")))
    if qc.get("target_site_recovery_status") != "site_contact_pass":
        flags.append(str(qc.get("target_site_recovery_status", "site_not_recovered")))
    if qc.get("hotspot_recovery_status") != "hotspot_contact_pass":
        flags.append(str(qc.get("hotspot_recovery_status", "hotspot_not_recovered")))
    rog = _parse_float(qc.get("peptide_radius_of_gyration", ""))
    if rog is not None and rog < min_peptide_rog:
        flags.append("peptide_collapsed_low_radius_of_gyration")
    return ";".join(flags) if flags else "pass_basic_pose_geometry"


def _ddg_proxy_status(ddg_proxy: float | None) -> str:
    if ddg_proxy is None:
        return "not_available"
    if ddg_proxy <= 0.0:
        return "favorable_or_neutral_interface_proxy"
    return "unfavorable_interface_proxy"


def _stage4_failure_reasons(
    *,
    qc: Mapping[str, Any],
    score_status: str,
    detached_flag: str,
    ddg_status: str,
    require_favorable_interface: bool,
) -> list[str]:
    reasons: list[str] = []
    if score_status != "success":
        reasons.append(score_status)
    if qc.get("sequence_status") != "pass_sequence":
        reasons.append(str(qc.get("sequence_status", "sequence_not_pass")))
    if qc.get("target_site_recovery_status") != "site_contact_pass":
        reasons.append(str(qc.get("target_site_recovery_status", "site_not_pass")))
    if qc.get("hotspot_recovery_status") != "hotspot_contact_pass":
        reasons.append(str(qc.get("hotspot_recovery_status", "hotspot_not_pass")))
    if qc.get("macrocycle_geometry_status") != "pass_head_to_tail_macrocycle":
        reasons.append(str(qc.get("macrocycle_geometry_status", "macrocycle_not_pass")))
    if qc.get("clash_status") != "pass_no_severe_clash":
        reasons.append(str(qc.get("clash_status", "clash_not_pass")))
    if detached_flag != "pass_basic_pose_geometry":
        reasons.append(detached_flag)
    if require_favorable_interface and ddg_status != "favorable_or_neutral_interface_proxy":
        reasons.append(ddg_status)
    return reasons


def _assign_energy_ranks(rows: list[dict[str, Any]]) -> None:
    ranked = sorted(
        [row for row in rows if _parse_float(row.get("ddg_proxy_no_repack", "")) is not None],
        key=lambda row: (
            float(row.get("ddg_proxy_no_repack", 999999.0)),
            float(row.get("complex_rosetta_total_score", 999999.0)),
            str(row.get("stage4_design_id", "")),
        ),
    )
    for rank, row in enumerate(ranked, start=1):
        row["stage4_energy_rank"] = rank


def _priority_hard_qc_pass(row: Mapping[str, Any]) -> bool:
    return (
        row.get("macrocycle_geometry_status") == "pass_head_to_tail_macrocycle"
        and row.get("target_site_recovery_status") == "site_contact_pass"
        and row.get("hotspot_recovery_status") == "hotspot_contact_pass"
        and row.get("clash_status") == "pass_no_severe_clash"
        and row.get("detached_or_collapsed_flag") == "pass_basic_pose_geometry"
    )


def _priority_sort_key(row: Mapping[str, Any]) -> tuple[Any, ...]:
    hard_qc_rank = 0 if _priority_hard_qc_pass(row) else 1
    ddg = _parse_float(row.get("ddg_proxy_no_repack", ""))
    return (
        hard_qc_rank,
        ddg if ddg is not None else 999999.0,
        -_parse_int(row.get("num_hotspot_contacts", 0)),
        -_parse_int(row.get("num_target_site_contacts", 0)),
        _parse_int(row.get("sequence_liability_count", 0)),
        str(row.get("stage4_design_id", "")),
    )


def _priority_notes(row: Mapping[str, Any]) -> str:
    notes = []
    if _priority_hard_qc_pass(row):
        notes.append("hard geometry/contact QC passed")
    else:
        notes.append("hard geometry/contact QC not fully passed")
    if row.get("ddg_proxy_status") == "favorable_or_neutral_interface_proxy":
        notes.append("no-repack ddG proxy favorable_or_neutral")
    elif row.get("ddg_proxy_status") == "unfavorable_interface_proxy":
        notes.append("no-repack ddG proxy weak_or_not_clearly_favorable")
    else:
        notes.append("no-repack ddG proxy not_available")
    seq_notes = str(row.get("sequence_liability_notes", ""))
    if seq_notes and seq_notes != "none":
        notes.append(f"sequence warnings: {seq_notes}")
    else:
        notes.append("no obvious sequence liability warning")
    return "; ".join(notes)


def _assign_priority(rows: list[dict[str, Any]], top_validation_count: int) -> None:
    ranked = sorted(rows, key=_priority_sort_key)
    for rank, row in enumerate(ranked, start=1):
        row["stage4_priority_rank"] = rank
        if not _priority_hard_qc_pass(row):
            row["stage4_priority_class"] = "priority_3_low_priority"
        elif rank <= min(3, top_validation_count):
            row["stage4_priority_class"] = "priority_1_validation_ready"
        elif rank <= top_validation_count:
            row["stage4_priority_class"] = "priority_2_backup"
        else:
            row["stage4_priority_class"] = "priority_3_low_priority"
        row["stage4_priority_notes"] = _priority_notes(row)


def _top_validation_rows(rows: list[dict[str, Any]], top_validation_count: int) -> list[dict[str, Any]]:
    selected = [
        row
        for row in sorted(rows, key=lambda item: _parse_int(item.get("stage4_priority_rank", 999999), 999999))
        if _priority_hard_qc_pass(row)
    ][:top_validation_count]
    top_rows = []
    for row in selected:
        ddg = _parse_float(row.get("ddg_proxy_no_repack", ""))
        if ddg is None:
            energy_note = "interface-energy proxy is not available"
        elif ddg <= 0:
            energy_note = "interface-energy proxy is favorable or neutral"
        else:
            energy_note = "interface-energy proxy is weak or not clearly favorable"
        top_rows.append(
            {
                "stage4_priority_rank": row.get("stage4_priority_rank", ""),
                "stage4_priority_class": row.get("stage4_priority_class", ""),
                "stage4_design_id": row.get("stage4_design_id", ""),
                "peptide_sequence": row.get("peptide_sequence", ""),
                "scored_pdb": row.get("scored_pdb", ""),
                "ddg_proxy_no_repack": row.get("ddg_proxy_no_repack", ""),
                "site_contacts": row.get("num_target_site_contacts", ""),
                "hotspot_contacts": row.get("num_hotspot_contacts", ""),
                "macrocycle_terminal_cn_distance": row.get("macrocycle_terminal_cn_distance", ""),
                "clash_status": row.get("clash_status", ""),
                "sequence_liability_notes": row.get("sequence_liability_notes", ""),
                "reason_selected": (
                    "geometry/contact/macrocycle QC passed; no severe clash; "
                    f"{energy_note}; selected for downstream structure validation."
                ),
            }
        )
    return top_rows


def _status_lines(rows: list[Mapping[str, Any]], field: str) -> str:
    counts = Counter(str(row.get(field, "")) for row in rows)
    if not counts:
        return "- none: 0"
    return "\n".join(f"- {key or 'blank'}: {counts[key]}" for key in sorted(counts))


def _summary_markdown(
    *,
    rows: list[Mapping[str, Any]],
    top_validation_rows: list[Mapping[str, Any]],
    output_dir: Path,
    args: argparse.Namespace,
) -> str:
    pass_rows = [row for row in rows if row.get("pass_stage4_qc") == "true"]
    top_rows = sorted(
        rows,
        key=lambda row: (
            row.get("pass_stage4_qc") != "true",
            float(row.get("stage4_priority_rank") or 999999.0),
            float(row.get("peptide_hotspot_min_distance") or 999.0),
            str(row.get("stage4_design_id", "")),
        ),
    )[: args.top_report]
    ddg_values = [_parse_float(row.get("ddg_proxy_no_repack", "")) for row in rows]
    ddg_values = [value for value in ddg_values if value is not None]
    all_ddg_positive = bool(ddg_values) and all(value > 0.0 for value in ddg_values)
    ddg_interpretation = (
        "All parsed ddg_proxy_no_repack values are positive. Geometry/contact/macrocycle QC passed, "
        "but the interface-energy proxy is weak or not clearly favorable; selected rows are only inputs "
        "for downstream structure validation."
        if all_ddg_positive
        else "ddg_proxy_no_repack is used only as a lightweight ranking proxy for downstream validation triage."
    )
    columns = [
        "stage4_design_id",
        "peptide_sequence",
        "ddg_proxy_no_repack",
        "ddg_proxy_status",
        "num_target_site_contacts",
        "num_hotspot_contacts",
        "peptide_hotspot_min_distance",
        "macrocycle_terminal_cn_distance",
        "clash_status",
        "detached_or_collapsed_flag",
        "stage4_energy_rank",
        "stage4_priority_rank",
        "stage4_priority_class",
        "sequence_liability_notes",
        "pass_stage4_qc",
        "stage4_failure_reasons",
    ]
    return f"""# FGA RFpeptides Stage 4 Rosetta Score-Only Interface Scoring

Status: Stage 3D-1 repack-only structures were scored without ordinary
FastRelax or backbone minimization.

Important rule: Stage 4 ranks and filters structures for downstream validation.
It does not make final peptide candidates. Ordinary unconstrained FastRelax is
not used here because it opened/damaged the RFpep_Site_2 macrocycle control.

`ddg_proxy_no_repack` is a no-repack separated-state proxy. It is not a formal
experimental binding energy and is not a final affinity judgment.

Method:

```text
input_stage3d1_pass_csv: {args.stage3d1_pass_csv}
score_mode: score_only_no_repack_no_minimization
separation_distance_A_for_ddg_proxy: {args.separation_distance:g}
require_favorable_interface_proxy: {args.require_favorable_interface}
top_validation_candidates: {args.top_validation_candidates}
CMS: not_available
SAP: not_available
```

Output directory:

```text
{output_dir}
```

## Counts

```text
input_rows_scored: {len(rows)}
pass_stage4_qc: {len(pass_rows)}
```

## QC Status Counts

Target-site recovery:

{_status_lines(rows, "target_site_recovery_status")}

Hotspot recovery:

{_status_lines(rows, "hotspot_recovery_status")}

Macrocycle geometry:

{_status_lines(rows, "macrocycle_geometry_status")}

Clash:

{_status_lines(rows, "clash_status")}

ddG proxy:

{_status_lines(rows, "ddg_proxy_status")}

Sequence liability notes:

{_status_lines(rows, "sequence_liability_notes")}

## ddG Proxy Interpretation

{ddg_interpretation}

## Recommended Top Validation Candidates

{rows_to_markdown(top_validation_rows, TOP_VALIDATION_FIELDS, "No top validation candidates were selected.")}

## Ranked Stage 4 Rows

{rows_to_markdown(top_rows, columns, "No Stage 4 rows were generated.")}
"""


def main() -> int:
    parser = argparse.ArgumentParser(description="Stage 4 Rosetta score-only interface scoring for Stage 3D-1 RFpeptides outputs.")
    parser.add_argument("--stage0-root", required=True)
    parser.add_argument("--stage3-root", required=True)
    parser.add_argument("--project-config", required=True)
    parser.add_argument("--output-root", default="", help="Defaults to --stage3-root.")
    parser.add_argument("--selected-backbones", required=True)
    parser.add_argument("--stage2-pass-csv", default="")
    parser.add_argument("--stage3d1-pass-csv", default="")
    parser.add_argument("--output-pdb-dir", default="")
    parser.add_argument("--separation-distance", type=float, default=1000.0)
    parser.add_argument("--require-favorable-interface", action="store_true")
    parser.add_argument("--min-peptide-rog", type=float, default=1.0)
    parser.add_argument("--contact-cutoff", type=float, default=5.0)
    parser.add_argument("--site-near-distance", type=float, default=6.0)
    parser.add_argument("--hotspot-near-distance", type=float, default=8.0)
    parser.add_argument("--severe-clash-distance", type=float, default=1.2)
    parser.add_argument("--min-target-contacts", type=int, default=1)
    parser.add_argument("--min-site-contacts", type=int, default=1)
    parser.add_argument("--min-hotspot-contacts", type=int, default=1)
    parser.add_argument("--macrocycle-pass-distance", type=float, default=2.0)
    parser.add_argument("--macrocycle-warn-distance", type=float, default=3.0)
    parser.add_argument("--forbidden-aas", default="CX")
    parser.add_argument("--top-report", type=int, default=50)
    parser.add_argument("--top-pymol", type=int, default=20)
    parser.add_argument("--top-validation-candidates", type=int, default=5)
    parser.add_argument(
        "--pymol-path-style",
        choices=["windows", "wsl", "native"],
        default="windows",
        help="Path style for generated PyMOL review scripts. Use windows for Windows PyMOL.",
    )
    args = parser.parse_args()

    logger = setup_logger("25_stage4_rosetta_interface_scoring")
    append_run_header(logger, "25_stage4_rosetta_interface_scoring.py")

    if args.separation_distance <= 0:
        raise RuntimeError("--separation-distance must be > 0")
    if args.contact_cutoff <= 0:
        raise RuntimeError("--contact-cutoff must be > 0")
    if args.site_near_distance < args.contact_cutoff:
        raise RuntimeError("--site-near-distance must be >= --contact-cutoff")
    if args.hotspot_near_distance < args.contact_cutoff:
        raise RuntimeError("--hotspot-near-distance must be >= --contact-cutoff")
    if args.macrocycle_warn_distance < args.macrocycle_pass_distance:
        raise RuntimeError("--macrocycle-warn-distance must be >= --macrocycle-pass-distance")
    if args.top_validation_candidates <= 0:
        raise RuntimeError("--top-validation-candidates must be > 0")

    qchelper = _load_stage3c_module()
    stage0_root = _resolve_mixed_path(args.stage0_root)
    stage3_root = _resolve_mixed_path(args.stage3_root)
    output_root = _resolve_mixed_path(args.output_root) if args.output_root else stage3_root
    output_dir = output_root / "06_rosetta_scoring"
    output_pdb_dir = _resolve_mixed_path(args.output_pdb_dir) if args.output_pdb_dir else output_dir / "scored_pdbs"
    stage2_pass_csv = (
        _resolve_mixed_path(args.stage2_pass_csv)
        if args.stage2_pass_csv
        else stage3_root / "03_backbone_qc" / "FGA_rfpeptides_backbones_qc_pass.csv"
    )
    stage3d1_pass_csv = (
        _resolve_mixed_path(args.stage3d1_pass_csv)
        if args.stage3d1_pass_csv
        else stage3_root / "05_proteinmpnn_sequences" / "FGA_rfpeptides_stage3D1_sidechain_repack_qc_pass.csv"
    )
    assert_active_route_path(stage0_root, "Stage 25 Stage 0 root")
    assert_active_route_path(stage3_root, "Stage 25 Stage 3 root")
    assert_active_route_path(output_root, "Stage 25 output root", must_exist=False)
    assert_active_route_path(stage2_pass_csv, "Stage 25 Stage 2 pass CSV")
    assert_active_route_path(stage3d1_pass_csv, "Stage 25 Stage 3D-1 pass CSV")
    route_manifest_path, route_manifest, route_manifest_sha256 = load_route_manifest(stage3_root)
    validate_route_project_config(args.project_config, route_manifest)
    source_route_provenance = route_provenance_fields(route_manifest_path, route_manifest, route_manifest_sha256)
    if output_root.resolve() != stage3_root.resolve():
        route_manifest_path, route_manifest, route_manifest_sha256 = write_route_manifest(output_root, route_manifest)
    route_provenance = route_provenance_fields(route_manifest_path, route_manifest, route_manifest_sha256)
    args.stage3d1_pass_csv = str(stage3d1_pass_csv)

    selected_backbones = set(_split_csv(args.selected_backbones))
    if not selected_backbones:
        raise RuntimeError("--selected-backbones must not be empty")

    stage2_pass_lookup = _lookup_rows(_read_required_csv(stage2_pass_csv), ["design_id"])
    stage3d1_rows = [
        row
        for row in _read_required_csv(stage3d1_pass_csv)
        if (not selected_backbones or row.get("backbone_id", "") in selected_backbones)
        and str(row.get("pass_stage3d1_qc", "")).strip().lower() == "true"
    ]
    if not stage3d1_rows:
        raise RuntimeError("No Stage 3D-1 pass rows were selected for Stage 4 scoring.")

    pyrosetta = _load_pyrosetta()
    _init_pyrosetta(pyrosetta)
    output_pdb_dir.mkdir(parents=True, exist_ok=True)

    result_rows: list[dict[str, Any]] = []
    pymol_rows: list[dict[str, Any]] = []
    last_site_numbers: set[str] = set()
    last_hotspot_numbers: set[str] = set()
    last_target_chain = "A"
    last_peptide_chain = "B"

    for input_row in stage3d1_rows:
        backbone_id = str(input_row.get("backbone_id", "")).strip()
        backbone_row = stage2_pass_lookup.get(backbone_id)
        if backbone_row is None:
            raise RuntimeError(f"Missing Stage 2 pass row for backbone: {backbone_id}")
        validate_row_route_provenance(backbone_row, source_route_provenance, f"Stage 25 Stage 2 row {backbone_id}")
        validate_row_route_provenance(input_row, source_route_provenance, f"Stage 25 Stage 3D-1 row {backbone_id}")

        site_label = str(backbone_row.get("site_label", ""))
        site_mapping_csv = stage0_root / "00_target_inputs" / f"{_safe_token(site_label)}_crop_renumbering_mapping.csv"
        site_numbers, hotspot_numbers = qchelper._load_site_mapping(site_mapping_csv)
        last_site_numbers = site_numbers
        last_hotspot_numbers = hotspot_numbers
        target_chain = str(backbone_row.get("target_chain", "")).strip() or "A"
        peptide_chain = str(backbone_row.get("peptide_chain", "")).strip() or "B"
        last_target_chain = target_chain
        last_peptide_chain = peptide_chain

        input_pdb = _resolve_mixed_path(str(input_row.get("repacked_pdb", "")))
        assert_active_route_path(input_pdb, "Stage 25 Stage 3D-1 repacked PDB")
        if not input_pdb.exists():
            raise RuntimeError(f"Missing Stage 3D-1 repacked PDB: {input_pdb}")
        source_stage3d1_id = str(input_row.get("stage3d1_design_id", input_pdb.stem)).strip() or input_pdb.stem
        stage4_id = f"{source_stage3d1_id}_stage4_score"
        scored_pdb = output_pdb_dir / f"{_safe_token(stage4_id)}.pdb"
        shutil.copy2(input_pdb, scored_pdb)

        score_metrics: dict[str, Any]
        score_status = "success"
        notes = "Score-only; no ordinary FastRelax, no repack, no backbone minimization."
        try:
            score_metrics = _score_complex_and_separated(
                pyrosetta=pyrosetta,
                input_pdb=scored_pdb,
                peptide_chain=peptide_chain,
                separation_distance=args.separation_distance,
            )
        except Exception as exc:
            score_status = f"failed:{exc.__class__.__name__}"
            score_metrics = {
                "complex_rosetta_total_score": "",
                "separated_rosetta_total_score": "",
                "ddg_proxy_no_repack": "",
            }
            notes = str(exc)

        qc = _qc_scored_pdb(
            qchelper=qchelper,
            scored_pdb=scored_pdb,
            backbone_row=backbone_row,
            site_numbers=site_numbers,
            hotspot_numbers=hotspot_numbers,
            args=args,
        )
        pymol_rows.append(qc)
        ddg_proxy = _parse_float(score_metrics.get("ddg_proxy_no_repack", ""))
        ddg_status = _ddg_proxy_status(ddg_proxy)
        detached_flag = _detached_or_collapsed_flag(qc, args.min_peptide_rog)
        failure_reasons = _stage4_failure_reasons(
            qc=qc,
            score_status=score_status,
            detached_flag=detached_flag,
            ddg_status=ddg_status,
            require_favorable_interface=args.require_favorable_interface,
        )
        peptide_sequence = str(input_row.get("peptide_sequence", qc.get("peptide_sequence", ""))).strip()
        sequence_metrics = _sequence_properties(peptide_sequence)

        result_rows.append(
            {
                "stage4_design_id": stage4_id,
                "source_stage3d1_design_id": source_stage3d1_id,
                "source_sequence_design_id": input_row.get("source_sequence_design_id", ""),
                "backbone_id": backbone_id,
                "site_label": backbone_row.get("site_label", ""),
                "site_id": backbone_row.get("site_id", ""),
                "input_pdb": input_pdb,
                "scored_pdb": scored_pdb,
                "peptide_sequence": peptide_sequence,
                **sequence_metrics,
                "score_mode": "score_only_no_repack_no_minimization",
                "pyrosetta_score_status": score_status,
                **score_metrics,
                "ddg_proxy_status": ddg_status,
                "CMS": "not_available",
                "SAP": "not_available",
                "target_contact_status": qc.get("target_contact_status", ""),
                "target_site_recovery_status": qc.get("target_site_recovery_status", ""),
                "hotspot_recovery_status": qc.get("hotspot_recovery_status", ""),
                "macrocycle_geometry_status": qc.get("macrocycle_geometry_status", ""),
                "clash_status": qc.get("clash_status", ""),
                "num_target_contacts": qc.get("num_target_contacts", ""),
                "num_target_site_contacts": qc.get("num_target_site_contacts", ""),
                "num_hotspot_contacts": qc.get("num_hotspot_contacts", ""),
                "peptide_target_min_distance": qc.get("peptide_target_min_distance", ""),
                "peptide_site_min_distance": qc.get("peptide_site_min_distance", ""),
                "peptide_hotspot_min_distance": qc.get("peptide_hotspot_min_distance", ""),
                "closest_target_residue": qc.get("closest_target_residue", ""),
                "closest_site_residue": qc.get("closest_site_residue", ""),
                "closest_hotspot_residue": qc.get("closest_hotspot_residue", ""),
                "macrocycle_terminal_cn_distance": qc.get("macrocycle_terminal_cn_distance", ""),
                "peptide_radius_of_gyration": qc.get("peptide_radius_of_gyration", ""),
                "detached_or_collapsed_flag": detached_flag,
                "stage4_energy_rank": "",
                "pass_stage4_qc": "true" if not failure_reasons else "false",
                "stage4_failure_reasons": ";".join(failure_reasons),
                "notes": notes,
            }
        )

    add_route_provenance(result_rows, route_provenance)
    _assign_energy_ranks(result_rows)
    _assign_priority(result_rows, args.top_validation_candidates)
    top_validation_rows = _top_validation_rows(result_rows, args.top_validation_candidates)
    add_route_provenance(top_validation_rows, route_provenance)
    pass_rows = [row for row in result_rows if row.get("pass_stage4_qc") == "true"]
    write_csv(output_dir / "FGA_rfpeptides_stage4_rosetta_interface_scores.csv", result_rows, STAGE4_FIELDS)
    write_csv(output_dir / "FGA_rfpeptides_stage4_rosetta_interface_scores_pass.csv", pass_rows, STAGE4_FIELDS)
    write_csv(
        output_dir / "FGA_rfpeptides_stage4_top_validation_candidates.csv",
        top_validation_rows,
        TOP_VALIDATION_FIELDS,
    )
    write_markdown(
        output_dir / "FGA_rfpeptides_stage4_rosetta_interface_scores.md",
        _summary_markdown(rows=result_rows, top_validation_rows=top_validation_rows, output_dir=output_dir, args=args),
    )
    qchelper._write_pymol_review(
        rows=pymol_rows,
        output_path=output_dir / "RFpep_Site_2_stage4_rosetta_score_only_review.pml",
        target_chain=last_target_chain,
        peptide_chain=last_peptide_chain,
        site_numbers=last_site_numbers,
        hotspot_numbers=last_hotspot_numbers,
        top_n=args.top_pymol,
        pymol_path_style=args.pymol_path_style,
    )

    logger.info("Stage 4 input rows scored: %s", len(result_rows))
    logger.info("Stage 4 passed QC: %s", len(pass_rows))
    logger.info("Stage 4 top validation candidates: %s", len(top_validation_rows))
    logger.info("Output directory: %s", output_dir)
    if not pass_rows:
        logger.warning("No Stage 4 outputs passed score-only QC.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
