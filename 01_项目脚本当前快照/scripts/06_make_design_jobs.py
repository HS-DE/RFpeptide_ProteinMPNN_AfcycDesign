from __future__ import annotations

import argparse
from pathlib import Path

from common import append_run_header, executable_available, load_config, project_root, read_csv, resolve_path, setup_logger, write_csv


def _first_chain(chain_list: str) -> str:
    return next((c.strip() for c in str(chain_list).split(",") if c.strip()), "A")


def _prefix_hotspots(chain_id: str, pdb_residue_numbers: str) -> str:
    labels = []
    seen = set()
    for value in str(pdb_residue_numbers).split(","):
        value = value.strip()
        if not value:
            continue
        label = value if value[0].isalpha() else f"{chain_id}{value}"
        if label not in seen:
            seen.add(label)
            labels.append(label)
    return ",".join(labels)


def _job_script(
    *,
    cleaned_pdb_file: str,
    job_id: str,
    patch_id: str,
    target_chain: str,
    hotspot: str,
    length: int,
    n_designs: int,
    out_dir: str,
) -> str:
    return "\n".join(
        [
            "#!/usr/bin/env bash",
            "set -euo pipefail",
            "",
            "# Production ColabDesign job generated from chain-qualified FGA patch residues.",
            "# Override NUM_DESIGNS/START_SEED/JOB_ID/OUT_DIR/LCF_PY/AF_PARAMS/PROJECT_DIR as needed.",
            'PROJECT_DIR="${PROJECT_DIR:-/mnt/c/SH/fga_cyclic_peptide_design}"',
            'LCF_PY="${LCF_PY:-$HOME/fga_model_envs/colabdesign-py310/.pixi/envs/default/bin/python}"',
            'AF_PARAMS="${AF_PARAMS:-$HOME/fga_model_envs/af_params}"',
            f'BASE_JOB_ID="{job_id}"',
            'JOB_ID="${JOB_ID:-$BASE_JOB_ID}"',
            'OUT_DIR="${OUT_DIR:-results/raw_designs/colabdesign_outputs/${JOB_ID}}"',
            f'NUM_DESIGNS="${{NUM_DESIGNS:-{n_designs}}}"',
            'START_SEED="${START_SEED:-0}"',
            'PSSM_ITERS="${PSSM_ITERS:-80}"',
            'GREEDY_ITERS="${GREEDY_ITERS:-32}"',
            'INIT_MODE="${INIT_MODE:-random}"',
            'INIT_SEED_SALT="${INIT_SEED_SALT:-$JOB_ID}"',
            'cd "$PROJECT_DIR"',
            "",
            'export CUDA_VISIBLE_DEVICES="${CUDA_VISIBLE_DEVICES:-0}"',
            'export CUDA_HOME="${CUDA_HOME:-/usr/local/cuda}"',
            'export LD_LIBRARY_PATH="${CUDA_HOME}/lib64:${CUDA_HOME}/targets/x86_64-linux/lib:/usr/lib/x86_64-linux-gnu:/usr/lib/wsl/lib:${LD_LIBRARY_PATH:-}"',
            'export XLA_PYTHON_CLIENT_PREALLOCATE="${XLA_PYTHON_CLIENT_PREALLOCATE:-false}"',
            'export XLA_PYTHON_CLIENT_MEM_FRACTION="${XLA_PYTHON_CLIENT_MEM_FRACTION:-0.30}"',
            'export XLA_PYTHON_CLIENT_ALLOCATOR="${XLA_PYTHON_CLIENT_ALLOCATOR:-platform}"',
            'export XLA_FLAGS="${XLA_FLAGS:---xla_gpu_enable_triton_gemm=false --xla_gpu_autotune_level=0}"',
            'export JAX_PLATFORMS="${JAX_PLATFORMS:-cuda}"',
            'export PYTHONUNBUFFERED="${PYTHONUNBUFFERED:-1}"',
            'export PYTHONFAULTHANDLER="${PYTHONFAULTHANDLER:-1}"',
            'export TF_FORCE_UNIFIED_MEMORY="${TF_FORCE_UNIFIED_MEMORY:-0}"',
            "",
            '"$LCF_PY" -X faulthandler scripts/external/run_colabdesign_cyclic_binder.py \\',
            f'  --target_pdb "{cleaned_pdb_file}" \\',
            f'  --target_chain "{target_chain}" \\',
            f'  --hotspot "{hotspot}" \\',
            f"  --peptide_length {length} \\",
            '  --num_designs "$NUM_DESIGNS" \\',
            '  --job_id "$JOB_ID" \\',
            f'  --patch_id "{patch_id}" \\',
            '  --seed "$START_SEED" \\',
            '  --data_dir "$AF_PARAMS" \\',
            "  --num_models 1 \\",
            "  --num_recycles 0 \\",
            '  --pssm_iters "$PSSM_ITERS" \\',
            '  --greedy_iters "$GREEDY_ITERS" \\',
            '  --init_mode "$INIT_MODE" \\',
            '  --init_seed_salt "$INIT_SEED_SALT" \\',
            '  --output_dir "$OUT_DIR"',
            "",
        ]
    )


def main() -> int:
    parser = argparse.ArgumentParser(description="Generate FGA cyclic peptide design jobs.")
    parser.add_argument("--config", default="config/project.yaml")
    args = parser.parse_args()

    logger = setup_logger("06_make_design_jobs")
    append_run_header(logger, "06_make_design_jobs.py")
    config = load_config(args.config)

    patches = read_csv("data/annotations/FGA_epitope_candidates.csv")
    if not patches:
        raise RuntimeError("Missing FGA_epitope_candidates.csv; run 05_select_surface_patches.py first.")

    root = project_root()
    colab_dir = root / "results/raw_designs/colabdesign_jobs"
    rf_dir = root / "results/raw_designs/rfdiffusion_jobs_optional"
    colab_dir.mkdir(parents=True, exist_ok=True)
    rf_dir.mkdir(parents=True, exist_ok=True)

    target_pdb = resolve_path(config["structures"]["cleaned_pdb_file"])
    cleaned_pdb_file = config["structures"]["cleaned_pdb_file"]
    length_distribution = {int(k): v for k, v in config["generation"]["length_distribution"].items()}
    colab_available = executable_available("colabdesign_cyclic_binder") or executable_available("afdesign")
    rf_available = executable_available("run_inference.py")

    rows = []
    for patch in patches:
        patch_cfg = config["generation"]["patches"].get(patch["patch_id"], {})
        n_total = int(patch_cfg.get("n_designs", 100))
        target_chain = patch.get("representative_chain") or _first_chain(patch.get("chain_id", "A"))
        hotspot = patch.get("representative_hotspot_residues") or _prefix_hotspots(target_chain, patch.get("pdb_residue_numbers", ""))
        hotspot_all_chains = patch.get("hotspot_pdb_residues") or _prefix_hotspots(target_chain, patch.get("pdb_residue_numbers", ""))
        hotspot_uniprot = patch.get("uniprot_residue_numbers", "")

        for length, fraction in length_distribution.items():
            n_designs = max(1, int(round(n_total * float(fraction))))
            job_id = f"{patch['patch_id']}_L{length}_colabdesign"
            out_dir = f"results/raw_designs/colabdesign_outputs/{job_id}"
            command = (
                "$LCF_PY scripts/external/run_colabdesign_cyclic_binder.py "
                f"--target_pdb {cleaned_pdb_file} "
                f"--target_chain {target_chain} "
                f"--hotspot {hotspot} "
                f"--peptide_length {length} "
                f"--num_designs {n_designs} "
                f"--job_id {job_id} "
                f"--patch_id {patch['patch_id']} "
                f"--data_dir $AF_PARAMS "
                f"--output_dir {out_dir}"
            )
            status = "ready" if colab_available else "pending_external_model"
            script_path = colab_dir / f"{job_id}.sh"
            script_path.write_text(
                _job_script(
                    cleaned_pdb_file=cleaned_pdb_file,
                    job_id=job_id,
                    patch_id=patch["patch_id"],
                    target_chain=target_chain,
                    hotspot=hotspot,
                    length=length,
                    n_designs=n_designs,
                    out_dir=out_dir,
                ),
                encoding="utf-8",
                newline="\n",
            )
            rows.append(
                {
                    "job_id": job_id,
                    "target_pdb": str(target_pdb),
                    "patch_id": patch["patch_id"],
                    "target_chain": target_chain,
                    "hotspot_residues": hotspot,
                    "hotspot_all_chains": hotspot_all_chains,
                    "hotspot_uniprot_residues": hotspot_uniprot,
                    "peptide_length": length,
                    "n_designs": n_designs,
                    "method": "ColabDesign-cyclic-binder",
                    "output_dir": out_dir,
                    "command_or_notebook": command,
                    "status": status,
                }
            )

        rf_job_id = f"{patch['patch_id']}_rfdiffusion_optional"
        rf_yaml = rf_dir / f"{rf_job_id}.yaml"
        rf_yaml.write_text(
            "\n".join(
                [
                    "# RFdiffusion/RFpeptides optional macrocycle design template.",
                    "# Fill in model-specific parameters in an RFdiffusion/RFpeptides environment.",
                    f"target_pdb: {target_pdb.as_posix()}",
                    f"patch_id: {patch['patch_id']}",
                    f"target_chain: {target_chain}",
                    f"hotspot_pdb_residues: \"{hotspot}\"",
                    f"hotspot_all_chains: \"{hotspot_all_chains}\"",
                    f"hotspot_uniprot_residues: \"{hotspot_uniprot}\"",
                    "cyclization: Cys-Cys disulfide",
                    "status: pending_external_model",
                    "",
                ]
            ),
            encoding="utf-8",
            newline="\n",
        )
        rows.append(
            {
                "job_id": rf_job_id,
                "target_pdb": str(target_pdb),
                "patch_id": patch["patch_id"],
                "target_chain": target_chain,
                "hotspot_residues": hotspot,
                "hotspot_all_chains": hotspot_all_chains,
                "hotspot_uniprot_residues": hotspot_uniprot,
                "peptide_length": "10-18",
                "n_designs": patch_cfg.get("n_designs", 100),
                "method": "RFdiffusion/RFpeptides optional",
                "output_dir": f"results/raw_designs/rfdiffusion_outputs/{rf_job_id}",
                "command_or_notebook": str(rf_yaml.relative_to(root)),
                "status": "ready" if rf_available else "pending_external_model",
            }
        )

    fields = [
        "job_id",
        "target_pdb",
        "patch_id",
        "target_chain",
        "hotspot_residues",
        "hotspot_all_chains",
        "hotspot_uniprot_residues",
        "peptide_length",
        "n_designs",
        "method",
        "output_dir",
        "command_or_notebook",
        "status",
    ]
    write_csv("results/raw_designs/design_jobs.csv", rows, fields)
    logger.info("Wrote design_jobs.csv, jobs=%s", len(rows))
    if not colab_available:
        logger.warning("ColabDesign executable not detected in this shell; jobs were generated for external/GPU execution.")
    if not rf_available:
        logger.warning("RFdiffusion/RFpeptides executable not detected; optional job templates only.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
