from __future__ import annotations

import argparse
import importlib.util
import math
import os
import sys
from collections import Counter
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

from common import assert_active_route_path, append_run_header, read_csv, resolve_path, rows_to_markdown, setup_logger, write_csv, write_markdown
from pdb_utils import parse_residues, residue_sequence


SCRIPT_DIR = Path(__file__).resolve().parent
STAGE3C_PATH = SCRIPT_DIR / "23_collect_proteinmpnn_sequences.py"


def _load_stage3c_module() -> Any:
    spec = importlib.util.spec_from_file_location("stage3c_qc", STAGE3C_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not import Stage 3C QC helper: {STAGE3C_PATH}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


STAGE3D1_FIELDS = [
    "stage3d1_design_id",
    "source_sequence_design_id",
    "backbone_id",
    "site_label",
    "site_id",
    "input_pdb",
    "repacked_pdb",
    "peptide_sequence",
    "sequence_duplicate_status",
    "repack_status",
    "repack_mode",
    "input_pose_total_score",
    "repacked_pose_total_score",
    "repack_score_delta",
    "repack_residue_count",
    "peptide_repack_residue_count",
    "target_repack_residue_count",
    "target_repack_residue_labels",
    "pre_target_site_recovery_status",
    "pre_hotspot_recovery_status",
    "pre_macrocycle_geometry_status",
    "pre_clash_status",
    "pre_peptide_site_min_distance",
    "pre_peptide_hotspot_min_distance",
    "pre_macrocycle_terminal_cn_distance",
    "post_target_site_recovery_status",
    "post_hotspot_recovery_status",
    "post_macrocycle_geometry_status",
    "post_clash_status",
    "post_num_target_site_contacts",
    "post_num_hotspot_contacts",
    "post_peptide_site_min_distance",
    "post_peptide_hotspot_min_distance",
    "post_macrocycle_terminal_cn_distance",
    "post_rosetta_pose_total",
    "pass_stage3d1_qc",
    "stage3d1_failure_reasons",
    "notes",
]


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


def _residue_sort_key(value: str) -> tuple[int, str]:
    digits = "".join(ch for ch in str(value) if ch.isdigit() or ch == "-")
    suffix = "".join(ch for ch in str(value) if not (ch.isdigit() or ch == "-"))
    try:
        return int(digits), suffix
    except ValueError:
        return 0, suffix


def _selected_stage3c_rows(
    rows: list[dict[str, str]],
    selected_backbones: set[str],
    include_duplicate_sequences: bool,
    require_contact_pass: bool,
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    selected: list[dict[str, str]] = []
    skipped: list[dict[str, str]] = []
    seen_sequences: set[str] = set()
    for row in rows:
        if selected_backbones and row.get("backbone_id", "") not in selected_backbones:
            continue
        sequence = row.get("peptide_sequence", "").strip()
        if not sequence:
            skipped.append({**row, "skip_reason": "missing_peptide_sequence"})
            continue
        if not include_duplicate_sequences and sequence in seen_sequences:
            skipped.append({**row, "skip_reason": "duplicate_peptide_sequence"})
            continue
        if row.get("sequence_status", "") != "pass_sequence":
            skipped.append({**row, "skip_reason": "sequence_status_not_pass"})
            continue
        if require_contact_pass:
            required = {
                "target_site_recovery_status": "site_contact_pass",
                "hotspot_recovery_status": "hotspot_contact_pass",
                "macrocycle_geometry_status": "pass_head_to_tail_macrocycle",
            }
            failed = [field for field, expected in required.items() if row.get(field, "") != expected]
            if failed:
                skipped.append({**row, "skip_reason": "failed_pre_repack_filter:" + ",".join(failed)})
                continue
        selected.append(row)
        seen_sequences.add(sequence)
    return selected, skipped


def _load_pyrosetta() -> Any:
    try:
        import pyrosetta  # type: ignore
    except ModuleNotFoundError as exc:
        raise RuntimeError(
            "PyRosetta is required for Stage 3D-1 side-chain repack. "
            "Run in the proteinmpnn_binder_design environment."
        ) from exc
    return pyrosetta


def _init_pyrosetta(pyrosetta: Any) -> None:
    pyrosetta.init(
        "-beta_nov16 -mute all -use_terminal_residues true "
        "-ex1 -ex2aro -packing:use_input_sc"
    )


def _xyz_tuple(xyz: Any) -> tuple[float, float, float]:
    try:
        return float(xyz.x), float(xyz.y), float(xyz.z)
    except TypeError:
        return float(xyz[0]), float(xyz[1]), float(xyz[2])


def _residue_atom_coords(residue: Any, heavy_only: bool = True) -> list[tuple[float, float, float]]:
    coords: list[tuple[float, float, float]] = []
    for atom_idx in range(1, residue.natoms() + 1):
        if heavy_only and residue.atom_is_hydrogen(atom_idx):
            continue
        coords.append(_xyz_tuple(residue.xyz(atom_idx)))
    return coords


def _sq_distance(left: tuple[float, float, float], right: tuple[float, float, float]) -> float:
    return sum((left[idx] - right[idx]) ** 2 for idx in range(3))


def _pose_chain(pose: Any, pose_idx: int) -> str:
    pdb_info = pose.pdb_info()
    chain = pdb_info.chain(pose_idx) if pdb_info is not None else ""
    return str(chain).strip() or "_"


def _pose_residue_label(pose: Any, pose_idx: int) -> str:
    pdb_info = pose.pdb_info()
    chain = _pose_chain(pose, pose_idx)
    number = pdb_info.number(pose_idx) if pdb_info is not None else pose_idx
    icode = pdb_info.icode(pose_idx).strip() if pdb_info is not None else ""
    return f"{chain}{number}{icode}"


def _pose_number(pose: Any, pose_idx: int) -> str:
    pdb_info = pose.pdb_info()
    if pdb_info is None:
        return str(pose_idx)
    icode = pdb_info.icode(pose_idx).strip()
    return f"{pdb_info.number(pose_idx)}{icode}"


def _mapped_site_hotspots(
    *,
    qchelper: Any,
    stage0_root: Path,
    backbone_row: Mapping[str, str],
    input_pdb: Path,
    target_chain: str,
) -> tuple[set[str], set[str], set[str]]:
    site_label = str(backbone_row.get("site_label", "")).strip()
    site_mapping_csv = stage0_root / "00_target_inputs" / f"{_safe_token(site_label)}_crop_renumbering_mapping.csv"
    site_numbers, hotspot_numbers = qchelper._load_site_mapping(site_mapping_csv)
    input_chains = parse_residues(input_pdb)
    input_target_residues = list(input_chains.get(target_chain, []))
    source_backbone_pdb = _resolve_mixed_path(str(backbone_row.get("rf_pdb", "")))
    assert_active_route_path(source_backbone_pdb, "Stage 24 source backbone PDB")
    target_number_map, _, _ = qchelper._target_number_mapping(
        source_backbone_pdb=source_backbone_pdb,
        target_chain=target_chain,
        relaxed_target_residues=input_target_residues,
    )
    mapped_site_numbers = qchelper._apply_number_mapping(site_numbers, target_number_map)
    mapped_hotspot_numbers = qchelper._apply_number_mapping(hotspot_numbers, target_number_map)
    return site_numbers, mapped_site_numbers, mapped_hotspot_numbers


def _repack_positions(
    *,
    pose: Any,
    peptide_chain: str,
    target_chain: str,
    target_repack_radius: float,
    mapped_hotspot_numbers: set[str],
    include_hotspots: bool,
) -> tuple[set[int], set[int], set[int]]:
    peptide_positions = {idx for idx in range(1, pose.total_residue() + 1) if _pose_chain(pose, idx) == peptide_chain}
    target_positions = {idx for idx in range(1, pose.total_residue() + 1) if _pose_chain(pose, idx) == target_chain}
    peptide_coords: list[tuple[float, float, float]] = []
    for idx in sorted(peptide_positions):
        peptide_coords.extend(_residue_atom_coords(pose.residue(idx), heavy_only=True))

    radius_sq = target_repack_radius * target_repack_radius
    interface_target_positions: set[int] = set()
    for idx in sorted(target_positions):
        target_coords = _residue_atom_coords(pose.residue(idx), heavy_only=True)
        if any(_sq_distance(left, right) <= radius_sq for left in target_coords for right in peptide_coords):
            interface_target_positions.add(idx)

    if include_hotspots:
        for idx in sorted(target_positions):
            if _pose_number(pose, idx) in mapped_hotspot_numbers:
                interface_target_positions.add(idx)

    return peptide_positions | interface_target_positions, peptide_positions, interface_target_positions


def _run_repack(
    *,
    pyrosetta: Any,
    input_pdb: Path,
    output_pdb: Path,
    peptide_chain: str,
    target_chain: str,
    target_repack_radius: float,
    mapped_hotspot_numbers: set[str],
    include_hotspots: bool,
) -> dict[str, Any]:
    pose = pyrosetta.pose_from_pdb(str(input_pdb))
    scorefxn = pyrosetta.get_fa_scorefxn()
    before_score = float(scorefxn(pose))
    repack_positions, peptide_positions, target_positions = _repack_positions(
        pose=pose,
        peptide_chain=peptide_chain,
        target_chain=target_chain,
        target_repack_radius=target_repack_radius,
        mapped_hotspot_numbers=mapped_hotspot_numbers,
        include_hotspots=include_hotspots,
    )
    if not peptide_positions:
        raise RuntimeError(f"Peptide chain {peptide_chain} was not found in {input_pdb}")
    if not target_positions:
        raise RuntimeError(f"No target interface residues selected for {input_pdb}")

    task = pyrosetta.standard_packer_task(pose)
    task.restrict_to_repacking()
    task.or_include_current(True)
    for idx in range(1, pose.total_residue() + 1):
        if idx not in repack_positions:
            task.nonconst_residue_task(idx).prevent_repacking()

    mover = pyrosetta.rosetta.protocols.minimization_packing.PackRotamersMover(scorefxn, task)
    mover.apply(pose)
    after_score = float(scorefxn(pose))

    output_pdb.parent.mkdir(parents=True, exist_ok=True)
    pose.dump_pdb(str(output_pdb))
    return {
        "input_pose_total_score": round(before_score, 3),
        "repacked_pose_total_score": round(after_score, 3),
        "repack_score_delta": round(after_score - before_score, 3),
        "repack_residue_count": len(repack_positions),
        "peptide_repack_residue_count": len(peptide_positions),
        "target_repack_residue_count": len(target_positions),
        "target_repack_residue_labels": ",".join(_pose_residue_label(pose, idx) for idx in sorted(target_positions)),
    }


def _qc_repacked_pdb(
    *,
    qchelper: Any,
    repacked_pdb: Path,
    backbone_row: Mapping[str, str],
    site_numbers: set[str],
    hotspot_numbers: set[str],
    args: argparse.Namespace,
) -> dict[str, Any]:
    return qchelper._qc_row_for_relaxed_pdb(
        relaxed_pdb=repacked_pdb,
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


def _stage3d1_failure_reasons(post_qc: Mapping[str, Any], input_sequence: str) -> list[str]:
    reasons: list[str] = []
    if str(post_qc.get("peptide_sequence", "")) != input_sequence:
        reasons.append("peptide_sequence_changed")
    if post_qc.get("sequence_status") != "pass_sequence":
        reasons.append(str(post_qc.get("sequence_status", "sequence_not_pass")))
    if post_qc.get("target_site_recovery_status") != "site_contact_pass":
        reasons.append(str(post_qc.get("target_site_recovery_status", "site_not_pass")))
    if post_qc.get("hotspot_recovery_status") != "hotspot_contact_pass":
        reasons.append(str(post_qc.get("hotspot_recovery_status", "hotspot_not_pass")))
    if post_qc.get("macrocycle_geometry_status") != "pass_head_to_tail_macrocycle":
        reasons.append(str(post_qc.get("macrocycle_geometry_status", "macrocycle_not_pass")))
    if post_qc.get("clash_status") != "pass_no_severe_clash":
        reasons.append(str(post_qc.get("clash_status", "clash_not_pass")))
    return reasons


def _summary_markdown(
    *,
    rows: list[Mapping[str, Any]],
    skipped_rows: list[Mapping[str, Any]],
    output_dir: Path,
    args: argparse.Namespace,
) -> str:
    pass_rows = [row for row in rows if row.get("pass_stage3d1_qc") == "true"]
    columns = [
        "stage3d1_design_id",
        "peptide_sequence",
        "post_target_site_recovery_status",
        "post_hotspot_recovery_status",
        "post_macrocycle_geometry_status",
        "post_clash_status",
        "post_peptide_hotspot_min_distance",
        "post_macrocycle_terminal_cn_distance",
        "pass_stage3d1_qc",
        "stage3d1_failure_reasons",
    ]
    skip_counts = Counter(str(row.get("skip_reason", "")) for row in skipped_rows)
    skip_lines = "\n".join(f"- {key or 'blank'}: {skip_counts[key]}" for key in sorted(skip_counts)) or "- none: 0"
    return f"""# FGA RFpeptides Stage 3D-1 Side-Chain Repack QC

Status: side-chain repack-only cleanup completed for ProteinMPNN-only Stage 3
outputs. This is not sequence redesign and not FastRelax.

Method:

```text
fixed sequence: true
backbone minimization: false
FastRelax: false
peptide chain repack: all residues
target chain repack: target residues within {args.target_repack_radius:g} A of peptide, plus mapped hotspots
```

Output directory:

```text
{output_dir}
```

Important rule: Stage 3D-1 passes only if repacking removes severe clashes
while preserving RFpep_Site_2 contact, hotspot contact, and head-to-tail
macrocycle geometry. Passing this step does not make a final peptide candidate.

## Counts

```text
input_rows_repacked: {len(rows)}
pass_stage3d1_qc: {len(pass_rows)}
skipped_rows: {len(skipped_rows)}
```

Skipped input rows:

{skip_lines}

## Repacked Rows

{rows_to_markdown(rows, columns, "No Stage 3D-1 rows were generated.")}
"""


def main() -> int:
    parser = argparse.ArgumentParser(description="Stage 3D-1 side-chain repack-only cleanup for ProteinMPNN-only outputs.")
    parser.add_argument("--stage0-root", required=True)
    parser.add_argument("--stage3-root", required=True)
    parser.add_argument("--output-root", default="", help="Defaults to --stage3-root.")
    parser.add_argument("--selected-backbones", required=True)
    parser.add_argument("--stage2-pass-csv", default="")
    parser.add_argument("--stage3c-qc-csv", default="")
    parser.add_argument("--output-pdb-dir", default="")
    parser.add_argument("--target-repack-radius", type=float, default=8.0)
    parser.add_argument("--include-duplicate-sequences", action="store_true")
    parser.add_argument("--no-require-stage3c-contact-pass", action="store_true")
    parser.add_argument("--no-include-hotspots", action="store_true")
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
    parser.add_argument(
        "--pymol-path-style",
        choices=["windows", "wsl", "native"],
        default="windows",
        help="Path style for generated PyMOL review scripts. Use windows for Windows PyMOL.",
    )
    args = parser.parse_args()

    logger = setup_logger("24_stage3d1_sidechain_repack")
    append_run_header(logger, "24_stage3d1_sidechain_repack.py")

    if args.target_repack_radius <= 0:
        raise RuntimeError("--target-repack-radius must be > 0")
    if args.contact_cutoff <= 0:
        raise RuntimeError("--contact-cutoff must be > 0")
    if args.site_near_distance < args.contact_cutoff:
        raise RuntimeError("--site-near-distance must be >= --contact-cutoff")
    if args.hotspot_near_distance < args.contact_cutoff:
        raise RuntimeError("--hotspot-near-distance must be >= --contact-cutoff")
    if args.macrocycle_warn_distance < args.macrocycle_pass_distance:
        raise RuntimeError("--macrocycle-warn-distance must be >= --macrocycle-pass-distance")

    qchelper = _load_stage3c_module()
    stage0_root = _resolve_mixed_path(args.stage0_root)
    stage3_root = _resolve_mixed_path(args.stage3_root)
    output_root = _resolve_mixed_path(args.output_root) if args.output_root else stage3_root
    output_dir = output_root / "05_proteinmpnn_sequences"
    output_pdb_dir = _resolve_mixed_path(args.output_pdb_dir) if args.output_pdb_dir else output_dir / "stage3d1_repack_only_pdbs"
    stage2_pass_csv = (
        _resolve_mixed_path(args.stage2_pass_csv)
        if args.stage2_pass_csv
        else stage3_root / "03_backbone_qc" / "FGA_rfpeptides_backbones_qc_pass.csv"
    )
    stage3c_qc_csv = (
        _resolve_mixed_path(args.stage3c_qc_csv)
        if args.stage3c_qc_csv
        else output_dir / "FGA_rfpeptides_stage3_proteinmpnn_only_sequences_qc.csv"
    )
    assert_active_route_path(stage0_root, "Stage 24 Stage 0 root")
    assert_active_route_path(stage3_root, "Stage 24 Stage 3 root")
    assert_active_route_path(output_root, "Stage 24 output root", must_exist=False)
    assert_active_route_path(stage2_pass_csv, "Stage 24 Stage 2 pass CSV")
    assert_active_route_path(stage3c_qc_csv, "Stage 24 Stage 3C QC CSV")

    selected_backbones = set(_split_csv(args.selected_backbones))
    if not selected_backbones:
        raise RuntimeError("--selected-backbones must not be empty")

    stage2_pass_lookup = _lookup_rows(_read_required_csv(stage2_pass_csv), ["design_id"])
    stage3c_rows = _read_required_csv(stage3c_qc_csv)
    selected_rows, skipped_rows = _selected_stage3c_rows(
        stage3c_rows,
        selected_backbones,
        include_duplicate_sequences=args.include_duplicate_sequences,
        require_contact_pass=not args.no_require_stage3c_contact_pass,
    )
    if not selected_rows:
        raise RuntimeError("No Stage 3C rows were selected for Stage 3D-1 repack.")

    pyrosetta = _load_pyrosetta()
    _init_pyrosetta(pyrosetta)

    output_pdb_dir.mkdir(parents=True, exist_ok=True)
    result_rows: list[dict[str, Any]] = []
    post_qc_rows: list[dict[str, Any]] = []
    last_site_numbers: set[str] = set()
    last_hotspot_numbers: set[str] = set()
    last_target_chain = "A"
    last_peptide_chain = "B"
    for idx, input_row in enumerate(selected_rows, start=1):
        backbone_id = str(input_row.get("backbone_id", "")).strip()
        backbone_row = stage2_pass_lookup.get(backbone_id)
        if backbone_row is None:
            raise RuntimeError(f"Missing Stage 2 pass row for backbone: {backbone_id}")
        peptide_chain = str(input_row.get("peptide_chain") or backbone_row.get("peptide_chain") or "B").strip()
        target_chain = str(input_row.get("target_chain") or backbone_row.get("target_chain") or "A").strip()
        last_target_chain = target_chain
        last_peptide_chain = peptide_chain

        input_pdb = _resolve_mixed_path(str(input_row.get("relaxed_pdb", "")))
        assert_active_route_path(input_pdb, "Stage 24 Stage 3C input PDB")
        if not input_pdb.exists():
            raise RuntimeError(f"Missing Stage 3C input PDB: {input_pdb}")
        input_sequence = str(input_row.get("peptide_sequence", "")).strip()
        source_id = str(input_row.get("sequence_design_id", input_pdb.stem)).strip() or input_pdb.stem
        stage3d1_id = f"{source_id}_stage3d1_repack"
        repacked_pdb = output_pdb_dir / f"{_safe_token(stage3d1_id)}.pdb"

        site_numbers, mapped_site_numbers, mapped_hotspot_numbers = _mapped_site_hotspots(
            qchelper=qchelper,
            stage0_root=stage0_root,
            backbone_row=backbone_row,
            input_pdb=input_pdb,
            target_chain=target_chain,
        )
        site_label = str(backbone_row.get("site_label", ""))
        site_mapping_csv = stage0_root / "00_target_inputs" / f"{_safe_token(site_label)}_crop_renumbering_mapping.csv"
        original_site_numbers, original_hotspot_numbers = qchelper._load_site_mapping(site_mapping_csv)
        last_site_numbers = mapped_site_numbers
        last_hotspot_numbers = mapped_hotspot_numbers

        try:
            repack_metrics = _run_repack(
                pyrosetta=pyrosetta,
                input_pdb=input_pdb,
                output_pdb=repacked_pdb,
                peptide_chain=peptide_chain,
                target_chain=target_chain,
                target_repack_radius=args.target_repack_radius,
                mapped_hotspot_numbers=mapped_hotspot_numbers,
                include_hotspots=not args.no_include_hotspots,
            )
            repack_status = "success"
            notes = ""
        except Exception as exc:
            repack_metrics = {
                "input_pose_total_score": "",
                "repacked_pose_total_score": "",
                "repack_score_delta": "",
                "repack_residue_count": "",
                "peptide_repack_residue_count": "",
                "target_repack_residue_count": "",
                "target_repack_residue_labels": "",
            }
            repack_status = f"failed:{exc.__class__.__name__}"
            notes = str(exc)

        if repack_status == "success":
            post_qc = _qc_repacked_pdb(
                qchelper=qchelper,
                repacked_pdb=repacked_pdb,
                backbone_row=backbone_row,
                site_numbers=original_site_numbers,
                hotspot_numbers=original_hotspot_numbers,
                args=args,
            )
            post_qc_rows.append(post_qc)
            failure_reasons = _stage3d1_failure_reasons(post_qc, input_sequence)
        else:
            post_qc = {}
            failure_reasons = [repack_status]

        result_rows.append(
            {
                "stage3d1_design_id": stage3d1_id,
                "source_sequence_design_id": source_id,
                "backbone_id": backbone_id,
                "site_label": backbone_row.get("site_label", ""),
                "site_id": backbone_row.get("site_id", ""),
                "input_pdb": input_pdb,
                "repacked_pdb": repacked_pdb if repack_status == "success" else "",
                "peptide_sequence": input_sequence,
                "sequence_duplicate_status": "included_duplicate" if args.include_duplicate_sequences else "unique_sequence",
                "repack_status": repack_status,
                "repack_mode": "side_chain_repack_only_no_backbone_minimization",
                **repack_metrics,
                "pre_target_site_recovery_status": input_row.get("target_site_recovery_status", ""),
                "pre_hotspot_recovery_status": input_row.get("hotspot_recovery_status", ""),
                "pre_macrocycle_geometry_status": input_row.get("macrocycle_geometry_status", ""),
                "pre_clash_status": input_row.get("clash_status", ""),
                "pre_peptide_site_min_distance": input_row.get("peptide_site_min_distance", ""),
                "pre_peptide_hotspot_min_distance": input_row.get("peptide_hotspot_min_distance", ""),
                "pre_macrocycle_terminal_cn_distance": input_row.get("macrocycle_terminal_cn_distance", ""),
                "post_target_site_recovery_status": post_qc.get("target_site_recovery_status", ""),
                "post_hotspot_recovery_status": post_qc.get("hotspot_recovery_status", ""),
                "post_macrocycle_geometry_status": post_qc.get("macrocycle_geometry_status", ""),
                "post_clash_status": post_qc.get("clash_status", ""),
                "post_num_target_site_contacts": post_qc.get("num_target_site_contacts", ""),
                "post_num_hotspot_contacts": post_qc.get("num_hotspot_contacts", ""),
                "post_peptide_site_min_distance": post_qc.get("peptide_site_min_distance", ""),
                "post_peptide_hotspot_min_distance": post_qc.get("peptide_hotspot_min_distance", ""),
                "post_macrocycle_terminal_cn_distance": post_qc.get("macrocycle_terminal_cn_distance", ""),
                "post_rosetta_pose_total": post_qc.get("rosetta_pose_total", ""),
                "pass_stage3d1_qc": "true" if not failure_reasons else "false",
                "stage3d1_failure_reasons": ";".join(failure_reasons),
                "notes": notes,
            }
        )

    pass_rows = [row for row in result_rows if row.get("pass_stage3d1_qc") == "true"]
    write_csv(output_dir / "FGA_rfpeptides_stage3D1_sidechain_repack_qc.csv", result_rows, STAGE3D1_FIELDS)
    write_csv(output_dir / "FGA_rfpeptides_stage3D1_sidechain_repack_qc_pass.csv", pass_rows, STAGE3D1_FIELDS)
    write_csv(output_dir / "FGA_rfpeptides_stage3D1_sidechain_repack_skipped_inputs.csv", skipped_rows, list(skipped_rows[0].keys()) if skipped_rows else ["skip_reason"])
    write_markdown(
        output_dir / "FGA_rfpeptides_stage3D1_sidechain_repack_qc.md",
        _summary_markdown(rows=result_rows, skipped_rows=skipped_rows, output_dir=output_dir, args=args),
    )
    if post_qc_rows:
        qchelper._write_pymol_review(
            rows=post_qc_rows,
            output_path=output_dir / "RFpep_Site_2_stage3D1_sidechain_repack_review.pml",
            target_chain=last_target_chain,
            peptide_chain=last_peptide_chain,
            site_numbers=last_site_numbers,
            hotspot_numbers=last_hotspot_numbers,
            top_n=20,
            pymol_path_style=args.pymol_path_style,
        )

    logger.info("Stage 3D-1 input rows repacked: %s", len(result_rows))
    logger.info("Stage 3D-1 passed QC: %s", len(pass_rows))
    logger.info("Skipped Stage 3C rows: %s", len(skipped_rows))
    logger.info("Output directory: %s", output_dir)
    if not pass_rows:
        logger.warning("No Stage 3D-1 outputs passed QC.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
