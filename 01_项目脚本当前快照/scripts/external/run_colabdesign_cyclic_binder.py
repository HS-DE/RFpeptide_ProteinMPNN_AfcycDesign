#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import types
from pathlib import Path

import numpy as np


def add_cyclic_offset(model, bug_fix: bool = True) -> None:
    """Apply the cyclic peptide complex offset used by the official notebook."""

    def cyclic_offset(length: int) -> np.ndarray:
        i = np.arange(length)
        ij = np.stack([i, i + length], -1)
        offset = i[:, None] - i[None, :]
        c_offset = np.abs(ij[:, None, :, None] - ij[None, :, None, :]).min((2, 3))
        if bug_fix:
            mask = c_offset < np.abs(offset)
            c_offset[mask] = -c_offset[mask]
        return c_offset * np.sign(offset)

    idx = model._inputs["residue_index"]
    offset = np.array(idx[:, None] - idx[None, :])
    if model.protocol != "binder":
        raise RuntimeError("Cyclic peptide offset is only supported for binder protocol.")
    offset[model._target_len :, model._target_len :] = cyclic_offset(model._binder_len)
    model._inputs["offset"] = offset


def require_af_params(data_dir: Path, num_models: int, use_multimer: bool) -> None:
    model_suffix = "multimer_v3" if use_multimer else "ptm"
    missing = []
    for idx in range(1, num_models + 1):
        name = f"model_{idx}_{model_suffix}"
        candidates = [
            data_dir / "params" / f"params_{name}.npz",
            data_dir / f"params_{name}.npz",
            data_dir / "params" / f"{name}.npz",
            data_dir / f"{name}.npz",
        ]
        if not any(path.is_file() for path in candidates):
            missing.append(name)
    if missing:
        raise FileNotFoundError(
            "Missing AlphaFold parameter files for ColabDesign: "
            + ", ".join(missing)
            + f". Expected under {data_dir}/params or {data_dir}."
        )


def cys_cys_scheme_pass(seq: str) -> bool:
    return len(seq) >= 2 and seq[0] == "C" and seq[-1] == "C" and "C" not in seq[1:-1]


def salted_seed(seed: int, salt: str) -> int:
    if not salt:
        return seed
    digest = hashlib.blake2b(f"{salt}:{seed}".encode("utf-8"), digest_size=8).digest()
    return int.from_bytes(digest, "little") % (2**32)


def initial_sequence(length: int, seed: int, mode: str, alphabet: str, salt: str = "") -> str:
    if mode == "fixed":
        return "C" + ("X" * (length - 2)) + "C"

    residues = [aa for aa in alphabet.upper() if aa.isalpha() and aa != "C"]
    if not residues:
        raise ValueError("--init_alphabet must contain at least one non-Cys amino acid.")
    rng = np.random.default_rng(salted_seed(seed, salt))
    middle = "".join(rng.choice(residues, size=length - 2))
    return f"C{middle}C"


def cys_cys_position_bias(length: int, cys_index: int, alphabet_size: int = 20) -> np.ndarray:
    bias = np.zeros((length, alphabet_size), dtype=np.float32)
    bias[1:-1, cys_index] -= 1e6
    bias[0, :] -= 1e6
    bias[-1, :] -= 1e6
    bias[0, cys_index] = 1e6
    bias[-1, cys_index] = 1e6
    return bias


def lock_terminal_cys_mutations(model, cys_index: int) -> None:
    """Force semigreedy mutation proposals to keep binder termini as Cys.

    The position bias used at restart steers logits and mutation sampling, but
    ColabDesign's semigreedy mutator can still choose a terminal position. This
    wrapper repairs every proposed mutant before prediction/scoring, so accepted
    semigreedy moves are evaluated as Cys-Cys cyclic peptides.
    """

    original_mutate = model._mutate

    def locked_mutate(self, seq, *args, **kwargs):
        mutant = np.array(original_mutate(seq, *args, **kwargs), copy=True)
        if mutant.shape[-1] < 2:
            raise ValueError("Cys-Cys terminal lock requires a binder length of at least 2.")
        mutant[..., 0] = cys_index
        mutant[..., -1] = cys_index
        return mutant

    model._mutate = types.MethodType(locked_mutate, model)


def best_log(model) -> dict:
    aux = model._tmp.get("best", {}).get("aux", None)
    if aux is None:
        aux = getattr(model, "aux", {})
    log = aux.get("log", {}) if isinstance(aux, dict) else {}
    out = {}
    for key, value in log.items():
        try:
            out[key] = float(value)
        except Exception:
            out[key] = str(value)
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description="Run ColabDesign cyclic binder hallucination for one FGA job.")
    parser.add_argument("--target_pdb", required=True)
    parser.add_argument("--target_chain", default="A")
    parser.add_argument("--hotspot", default="")
    parser.add_argument("--peptide_length", type=int, required=True)
    parser.add_argument("--num_designs", type=int, required=True)
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--data_dir", required=True, help="Directory containing AlphaFold params/ files.")
    parser.add_argument("--job_id", default="colabdesign_job")
    parser.add_argument("--patch_id", default="")
    parser.add_argument("--target_region", default="FGA_chain_36_866")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--num_models", type=int, default=1)
    parser.add_argument("--num_recycles", type=int, default=0)
    parser.add_argument("--use_multimer", action="store_true")
    parser.add_argument("--pssm_iters", type=int, default=120)
    parser.add_argument("--greedy_iters", type=int, default=32)
    parser.add_argument("--learning_rate", type=float, default=0.1)
    parser.add_argument("--dropout", action="store_true")
    parser.add_argument("--init_mode", choices=["fixed", "random"], default="random")
    parser.add_argument("--init_alphabet", default="ADEFGHIKLMNPQRSTVWY")
    parser.add_argument(
        "--init_seed_salt",
        default="",
        help="Optional salt mixed into random initialization; defaults to job/patch/target metadata.",
    )
    parser.add_argument("--bugfix_cyclic_offset", action="store_true", default=True)
    args = parser.parse_args()

    if args.peptide_length < 2:
        raise ValueError("--peptide_length must be at least 2 for Cys-Cys cyclic peptides.")

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    data_dir = Path(args.data_dir).expanduser().resolve()
    require_af_params(data_dir, args.num_models, args.use_multimer)

    # Import after fast parameter checks so missing params fail before JAX starts.
    from colabdesign import clear_mem, mk_afdesign_model
    from colabdesign.af.alphafold.common import residue_constants

    clear_mem()
    if args.use_multimer:
        model_names = [f"model_{idx}_multimer_v3" for idx in range(1, args.num_models + 1)]
    else:
        model_names = [f"model_{idx}_ptm" for idx in range(1, args.num_models + 1)]

    model = mk_afdesign_model(
        protocol="binder",
        use_multimer=args.use_multimer,
        num_recycles=args.num_recycles,
        recycle_mode="sample",
        model_names=model_names,
        data_dir=str(data_dir),
    )
    model.prep_inputs(
        pdb_filename=args.target_pdb,
        chain=args.target_chain,
        binder_len=args.peptide_length,
        hotspot=args.hotspot or None,
        ignore_missing=False,
    )
    add_cyclic_offset(model, bug_fix=args.bugfix_cyclic_offset)

    fields = [
        "raw_id",
        "method",
        "job_id",
        "patch_id",
        "target_region",
        "core_sequence",
        "core_length",
        "raw_score",
        "source_file",
        "notes",
    ]
    rows = []
    changed_sequences = 0
    flags = {
        "num_recycles": args.num_recycles,
        "models": model._model_names[: args.num_models],
        "dropout": args.dropout,
    }
    init_seed_salt = args.init_seed_salt or f"{args.job_id}:{args.patch_id}:{args.target_chain}:{args.hotspot}"
    alphabet_size = int(getattr(model, "_args", {}).get("alphabet_size", 20))
    cys_bias = cys_cys_position_bias(args.peptide_length, residue_constants.restype_order["C"], alphabet_size)
    lock_terminal_cys_mutations(model, residue_constants.restype_order["C"])

    for design_idx in range(1, args.num_designs + 1):
        seed = args.seed + design_idx - 1
        init_seq = initial_sequence(args.peptide_length, seed, args.init_mode, args.init_alphabet, init_seed_salt)
        model.restart(seed=seed, seq=init_seq, bias=cys_bias)
        model.set_optimizer(optimizer="sgd", learning_rate=args.learning_rate, norm_seq_grad=True)
        model.design_pssm_semigreedy(args.pssm_iters, args.greedy_iters, **flags)

        seqs = model.get_seqs()
        seq = seqs[0] if seqs else ""
        pdb_path = out_dir / f"{args.job_id}_design_{design_idx:04d}.pdb"
        model.save_pdb(str(pdb_path))
        log = best_log(model)
        accepted = cys_cys_scheme_pass(seq)
        sequence_changed = seq != init_seq
        if sequence_changed:
            changed_sequences += 1
        status = "accepted_cys_cys_scheme_A" if accepted else "rejected_non_scheme_A"
        notes = {
            "status": status,
            "seed": seed,
            "init_mode": args.init_mode,
            "init_sequence": init_seq,
            "init_seed_salt": init_seed_salt,
            "final_sequence_changed": sequence_changed,
            "terminal_cys_mutation_lock": True,
            "target_chain": args.target_chain,
            "hotspot": args.hotspot,
            "colabdesign_log": log,
        }
        if accepted:
            rows.append(
                {
                    "raw_id": f"{args.job_id}_raw_{len(rows) + 1:04d}",
                    "method": "ColabDesign-cyclic-binder",
                    "job_id": args.job_id,
                    "patch_id": args.patch_id,
                    "target_region": args.target_region,
                    "core_sequence": seq,
                    "core_length": len(seq),
                    "raw_score": log.get("loss", ""),
                    "source_file": str(pdb_path.name),
                    "notes": json.dumps(notes, sort_keys=True),
                }
            )
        else:
            reject_path = out_dir / "rejected_non_scheme_A.csv"
            with reject_path.open("a", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=["job_id", "design_index", "sequence", "notes"])
                if handle.tell() == 0:
                    writer.writeheader()
                writer.writerow(
                    {
                        "job_id": args.job_id,
                        "design_index": design_idx,
                        "sequence": seq,
                        "notes": json.dumps(notes, sort_keys=True),
                    }
                )

        csv_path = out_dir / "candidates.csv"
        with csv_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)

    summary = {
        "job_id": args.job_id,
        "patch_id": args.patch_id,
        "num_designs_requested": args.num_designs,
        "num_scheme_A_candidates": len(rows),
        "num_changed_sequences": changed_sequences,
        "terminal_cys_mutation_lock": True,
        "output_csv": str(out_dir / "candidates.csv"),
    }
    (out_dir / "run_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
