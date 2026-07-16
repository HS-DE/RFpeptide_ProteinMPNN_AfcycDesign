from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Sequence

from common import append_run_header, clean_sequence, executable_available, load_config, project_root, read_csv, read_fasta_sequence, resolve_path, setup_logger, write_csv


JOB_FIELDS = [
    "boltz_job_id",
    "peptide_id",
    "core_sequence",
    "patch_id",
    "seed",
    "input_yaml",
    "output_dir",
    "command",
    "status",
    "msa_mode",
    "target_chain_id",
    "peptide_chain_id",
    "template_pdb",
    "template_chain_id",
    "hotspot_uniprot_residues",
    "hotspot_target_positions",
    "disulfide_bond",
    "notes",
]

SKIP_FIELDS = [
    "raw_id",
    "patch_id",
    "core_sequence",
    "reason",
]


def _truthy(value: Any) -> bool:
    return str(value).lower() in {"true", "1", "yes", "pass"}


def _safe_id(value: str) -> str:
    safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value).strip())
    return safe.strip("_") or "unnamed"


def _json_notes(row: Mapping[str, str]) -> Dict[str, Any]:
    raw = row.get("notes", "")
    if not raw:
        return {}
    try:
        parsed = json.loads(raw)
    except json.JSONDecodeError:
        return {}
    return parsed if isinstance(parsed, dict) else {}


def _comma_ints(value: str) -> List[int]:
    vals: List[int] = []
    for part in str(value).split(","):
        part = part.strip()
        if not part:
            continue
        try:
            vals.append(int(part))
        except ValueError:
            continue
    return vals


def _load_patch_hotspots() -> Dict[str, List[int]]:
    rows = read_csv("data/annotations/FGA_epitope_candidates.csv")
    out: Dict[str, List[int]] = {}
    for row in rows:
        patch_id = row.get("patch_id", "")
        residues = _comma_ints(row.get("representative_hotspot_uniprot_residues", ""))
        if patch_id and residues:
            out[patch_id] = sorted(set(residues))
    return out


def _target_positions_from_uniprot(uniprot_residues: Sequence[int], target_region_start: int) -> List[int]:
    return [res - target_region_start + 1 for res in uniprot_residues if res >= target_region_start]


def _pdb_has_seqres(path: str) -> bool:
    if not path:
        return False
    pdb_path = resolve_path(path)
    if not pdb_path.exists() or pdb_path.suffix.lower() != ".pdb":
        return True
    with pdb_path.open("r", encoding="utf-8", errors="ignore") as handle:
        return any(line.startswith("SEQRES") for line in handle)


def _boltz_template_chain_id(template_pdb: str, chain_id: str) -> str:
    if not template_pdb or not chain_id:
        return chain_id
    path = resolve_path(template_pdb)
    if path.suffix.lower() != ".pdb":
        return chain_id
    if any(ch.isdigit() for ch in chain_id):
        return chain_id
    # Boltz converts PDB chains to mmCIF subchains with a numeric suffix
    # during template parsing; the first subchain of PDB chain J becomes J1.
    return f"{chain_id}1"


def _candidate_skip_reason(row: Mapping[str, str], require_generation_guards: bool) -> str:
    core = clean_sequence(row.get("core_sequence", ""))
    if not _truthy(row.get("sequence_filter_pass", "")):
        return "sequence_filter_pass=false"
    if len(core) < 2 or not core.startswith("C") or not core.endswith("C"):
        return "not_terminal_cys_cys"
    if core[1:-1].count("C"):
        return "internal_cys_present"
    if not require_generation_guards:
        return ""

    notes = _json_notes(row)
    if notes.get("final_sequence_changed") is not True:
        return "final_sequence_changed_not_true"
    if notes.get("terminal_cys_mutation_lock") is not True:
        return "terminal_cys_mutation_lock_not_true"
    return ""


def _yaml_text(
    *,
    target_sequence: str,
    peptide_sequence: str,
    target_chain_id: str,
    peptide_chain_id: str,
    msa_mode: str,
    template_pdb: str,
    template_chain_id: str,
    hotspot_target_positions: Sequence[int],
    pocket_max_distance: float,
    pocket_force: bool,
) -> str:
    lines = [
        "version: 1",
        "sequences:",
        "  - protein:",
        f"      id: {target_chain_id}",
        f"      sequence: {target_sequence}",
    ]
    if msa_mode == "empty":
        lines.append("      msa: empty")

    lines.extend(
        [
            "  - protein:",
            f"      id: {peptide_chain_id}",
            f"      sequence: {peptide_sequence}",
        ]
    )
    if msa_mode == "empty":
        lines.append("      msa: empty")

    lines.extend(
        [
            "constraints:",
            "  - bond:",
            f"      atom1: [{peptide_chain_id}, 1, SG]",
            f"      atom2: [{peptide_chain_id}, {len(peptide_sequence)}, SG]",
        ]
    )

    if hotspot_target_positions:
        lines.extend(
            [
                "  - pocket:",
                f"      binder: {peptide_chain_id}",
                "      contacts:",
            ]
        )
        for pos in hotspot_target_positions:
            lines.append(f"        - [{target_chain_id}, {pos}]")
        lines.extend(
            [
                f"      max_distance: {pocket_max_distance:g}",
                f"      force: {'true' if pocket_force else 'false'}",
            ]
        )

    if template_pdb and template_chain_id:
        lines.extend(
            [
                "templates:",
                f"  - pdb: {template_pdb}",
                f"    chain_id: {target_chain_id}",
                f"    template_id: {template_chain_id}",
            ]
        )

    return "\n".join(lines) + "\n"


def _build_command(input_yaml: str, output_dir: str, seed: int, msa_mode: str, use_potentials: bool, output_format: str) -> str:
    parts = [
        "boltz",
        "predict",
        input_yaml,
        "--out_dir",
        output_dir,
        "--model",
        "boltz2",
        "--seed",
        str(seed),
        "--output_format",
        output_format,
        "--write_full_pae",
    ]
    if msa_mode == "server":
        parts.append("--use_msa_server")
    if use_potentials:
        parts.append("--use_potentials")
    return " ".join(parts)


def _iter_selected_candidates(rows: Iterable[Mapping[str, str]], max_candidates: int | None) -> Iterable[Mapping[str, str]]:
    count = 0
    for row in rows:
        yield row
        count += 1
        if max_candidates is not None and count >= max_candidates:
            break


def main() -> int:
    parser = argparse.ArgumentParser(description="Prepare Boltz YAML inputs and job table for hard-filtered FGA cyclic peptide candidates.")
    parser.add_argument("--config", default="config/project.yaml")
    parser.add_argument("--candidate-csv", default="results/filtered/FGA_hard_filtered_candidates.csv")
    parser.add_argument("--target-fasta", default="data/input/FGA_chain_36_866.fasta")
    parser.add_argument("--output-root", default="results/boltz_predictions")
    parser.add_argument("--seeds", type=int, default=None)
    parser.add_argument("--max-candidates", type=int, default=None)
    parser.add_argument("--msa-mode", choices=["server", "empty"], default="server")
    parser.add_argument("--target-chain-id", default="A")
    parser.add_argument("--peptide-chain-id", default="B")
    parser.add_argument("--target-region-start", type=int, default=36)
    parser.add_argument("--template-pdb", default="")
    parser.add_argument("--default-template-chain-id", default="J")
    parser.add_argument("--no-template", action="store_true")
    parser.add_argument("--no-pocket-constraints", action="store_true")
    parser.add_argument("--pocket-max-distance", type=float, default=8.0)
    parser.add_argument("--pocket-force", action="store_true")
    parser.add_argument("--no-potentials", action="store_true")
    parser.add_argument("--output-format", choices=["pdb", "mmcif"], default="pdb")
    parser.add_argument("--allow-missing-generation-guards", action="store_true")
    args = parser.parse_args()

    logger = setup_logger("14_prepare_boltz_prediction_jobs")
    append_run_header(logger, "14_prepare_boltz_prediction_jobs.py")
    config = load_config(args.config)
    seeds = args.seeds if args.seeds is not None else int(config["complex_prediction"]["seeds_per_candidate"])
    if seeds < 1:
        raise ValueError("--seeds must be >= 1")

    output_root = resolve_path(args.output_root)
    yaml_dir = output_root / "inputs" / "yaml"
    output_dir_root = output_root / "outputs"
    yaml_dir.mkdir(parents=True, exist_ok=True)
    output_dir_root.mkdir(parents=True, exist_ok=True)

    target_sequence = read_fasta_sequence(args.target_fasta)
    candidates = read_csv(args.candidate_csv)
    patch_hotspots = _load_patch_hotspots()
    boltz_available = executable_available("boltz")
    template_skip_reason = ""
    template_pdb = "" if args.no_template else str(args.template_pdb).replace("\\", "/").strip()
    if template_pdb and not _pdb_has_seqres(template_pdb):
        logger.warning(
            "Template PDB %s has no SEQRES records; skipping template to avoid Boltz template parsing failure.",
            template_pdb,
        )
        template_skip_reason = "template_pdb_missing_seqres"
        template_pdb = ""
    elif args.no_template:
        template_skip_reason = "no_template_requested"
    elif not template_pdb:
        template_skip_reason = "no_template_configured"

    rows: List[Dict[str, Any]] = []
    skipped: List[Dict[str, str]] = []
    selected = _iter_selected_candidates(candidates, args.max_candidates)
    accepted_idx = 0
    for cand in selected:
        reason = _candidate_skip_reason(cand, require_generation_guards=not args.allow_missing_generation_guards)
        if reason:
            skipped.append(
                {
                    "raw_id": cand.get("raw_id", ""),
                    "patch_id": cand.get("patch_id", ""),
                    "core_sequence": cand.get("core_sequence", ""),
                    "reason": reason,
                }
            )
            continue

        peptide_id = _safe_id(cand.get("raw_id", ""))
        core = clean_sequence(cand.get("core_sequence", ""))
        patch_id = cand.get("patch_id", "")
        notes = _json_notes(cand)
        source_template_chain_id = str(notes.get("target_chain") or args.default_template_chain_id)
        template_chain_id = _boltz_template_chain_id(template_pdb, source_template_chain_id)
        hotspot_uniprot = patch_hotspots.get(patch_id, [])
        hotspot_target_positions = [] if args.no_pocket_constraints else _target_positions_from_uniprot(hotspot_uniprot, args.target_region_start)
        accepted_idx += 1
        job_prefix = _safe_id(f"{patch_id}_boltz_{accepted_idx:05d}")

        for seed in range(1, seeds + 1):
            job_id = f"{job_prefix}_seed{seed}"
            yaml_rel = f"{args.output_root}/inputs/yaml/{job_id}.yaml"
            out_rel = f"{args.output_root}/outputs/{job_id}"
            yaml_path = resolve_path(yaml_rel)
            yaml_path.write_text(
                _yaml_text(
                    target_sequence=target_sequence,
                    peptide_sequence=core,
                    target_chain_id=args.target_chain_id,
                    peptide_chain_id=args.peptide_chain_id,
                    msa_mode=args.msa_mode,
                    template_pdb=template_pdb,
                    template_chain_id=template_chain_id,
                    hotspot_target_positions=hotspot_target_positions,
                    pocket_max_distance=args.pocket_max_distance,
                    pocket_force=args.pocket_force,
                ),
                encoding="utf-8",
                newline="\n",
            )
            rows.append(
                {
                    "boltz_job_id": job_id,
                    "peptide_id": peptide_id,
                    "core_sequence": core,
                    "patch_id": patch_id,
                    "seed": seed,
                    "input_yaml": yaml_rel.replace("\\", "/"),
                    "output_dir": out_rel.replace("\\", "/"),
                    "command": _build_command(
                        yaml_rel.replace("\\", "/"),
                        out_rel.replace("\\", "/"),
                        seed,
                        args.msa_mode,
                        use_potentials=not args.no_potentials,
                        output_format=args.output_format,
                    ),
                    "status": "ready" if boltz_available else "pending_external_model",
                    "msa_mode": args.msa_mode,
                    "target_chain_id": args.target_chain_id,
                    "peptide_chain_id": args.peptide_chain_id,
                    "template_pdb": template_pdb,
                    "template_chain_id": template_chain_id if template_pdb else "",
                    "hotspot_uniprot_residues": ",".join(str(x) for x in hotspot_uniprot),
                    "hotspot_target_positions": ",".join(str(x) for x in hotspot_target_positions),
                    "disulfide_bond": f"{args.peptide_chain_id}:1:SG-{args.peptide_chain_id}:{len(core)}:SG",
                    "notes": json.dumps(
                        {
                            "source_candidate_csv": args.candidate_csv,
                            "source_candidate_raw_id": cand.get("raw_id", ""),
                            "source_template_chain_id": source_template_chain_id,
                            "template_skip_reason": template_skip_reason,
                            "uses_template": bool(template_pdb),
                            "uses_pocket_constraints": not args.no_pocket_constraints,
                            "uses_disulfide_bond_constraint": True,
                        },
                        sort_keys=True,
                    ),
                }
            )

    write_csv(output_root / "boltz_jobs.csv", rows, JOB_FIELDS)
    write_csv(output_root / "boltz_job_skipped_candidates.csv", skipped, SKIP_FIELDS)
    logger.info("Wrote Boltz jobs: %s", len(rows))
    logger.info("Skipped candidates: %s", len(skipped))
    if not boltz_available:
        logger.warning("boltz executable not detected in this shell; jobs are written as pending external model tasks.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
