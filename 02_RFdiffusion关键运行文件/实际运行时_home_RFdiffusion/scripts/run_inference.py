#!/usr/bin/env python
"""
Inference script.

To run with base.yaml as the config,

> python run_inference.py

To specify a different config,

> python run_inference.py --config-name symmetry

where symmetry can be the filename of any other config (without .yaml extension)
See https://hydra.cc/docs/advanced/hydra-command-line-flags/ for more options.

"""

import re
import os, time, pickle
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import torch
from omegaconf import OmegaConf
import hydra
import logging
from rfdiffusion.util import writepdb_multi, writepdb
import rfdiffusion.util as rfd_util
from rfdiffusion.inference import utils as iu
from hydra.core.hydra_config import HydraConfig
import numpy as np
import random
import glob


def make_deterministic(seed=0):
    torch.manual_seed(seed)
    np.random.seed(seed)
    random.seed(seed)


def sha256_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sha256_text(value):
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def canonical_hotspots(values):
    parsed = []
    for value in values or []:
        text = str(value).strip()
        match = re.fullmatch(r"([A-Za-z])(\-?\d+)", text)
        if match is None:
            raise RuntimeError(f"Invalid hotspot residue label: {value}")
        parsed.append((match.group(1), int(match.group(2))))
    if len(set(parsed)) != len(parsed):
        raise RuntimeError(f"Duplicate hotspot residue labels are not allowed: {values}")
    return ",".join(f"{chain}{number}" for chain, number in sorted(parsed))


def required_environment(name):
    value = os.environ.get(name, "").strip()
    if not value:
        raise RuntimeError(f"Required Stage 1 provenance environment variable is missing: {name}")
    return value


@hydra.main(version_base=None, config_path="../config/inference", config_name="base")
def main(conf: HydraConfig) -> None:
    log = logging.getLogger(__name__)
    if conf.inference.deterministic:
        make_deterministic()

    # Check for available GPU and print result of check
    if torch.cuda.is_available():
        device_name = torch.cuda.get_device_name(torch.cuda.current_device())
        log.info(f"Found GPU with device_name {device_name}. Will run RFdiffusion on {device_name}")
    else:
        log.info("////////////////////////////////////////////////")
        log.info("///// NO GPU DETECTED! Falling back to CPU /////")
        log.info("////////////////////////////////////////////////")

    # Initialize sampler and target/contig.
    sampler = iu.sampler_selector(conf)

    provenance_required = os.environ.get("RFPEPTIDES_REQUIRE_PROVENANCE_CLOSURE", "") == "1"
    input_pdb_path = Path(sampler.inf_conf.input_pdb).resolve()
    input_pdb_sha256 = sha256_file(input_pdb_path)
    configured_hotspots_normalized = canonical_hotspots(sampler.ppi_conf.hotspot_res or [])
    configured_hotspots_sha256 = sha256_text(configured_hotspots_normalized)
    stage0_mapping_csv = os.environ.get("STAGE0_MAPPING_CSV", "").strip()
    stage0_mapping_csv_sha256 = ""
    if stage0_mapping_csv:
        stage0_mapping_csv = str(Path(stage0_mapping_csv).resolve())
        stage0_mapping_csv_sha256 = sha256_file(stage0_mapping_csv)

    if provenance_required:
        expected_target_sha256 = required_environment("STAGE0_TARGET_PDB_SHA256")
        expected_mapping_sha256 = required_environment("STAGE0_MAPPING_CSV_SHA256")
        expected_hotspots_sha256 = required_environment("STAGE0_HOTSPOTS_SHA256")
        expected_hotspots_normalized = required_environment("STAGE0_HOTSPOTS_NORMALIZED")
        required_environment("STAGE0_MAPPING_CSV")
        if input_pdb_sha256 != expected_target_sha256:
            raise RuntimeError(
                "Stage 0 target PDB hash changed before inference: "
                f"observed={input_pdb_sha256} expected={expected_target_sha256}"
            )
        if stage0_mapping_csv_sha256 != expected_mapping_sha256:
            raise RuntimeError(
                "Stage 0 mapping CSV hash changed before inference: "
                f"observed={stage0_mapping_csv_sha256} expected={expected_mapping_sha256}"
            )
        if configured_hotspots_normalized != expected_hotspots_normalized:
            raise RuntimeError(
                "Configured hotspots differ from the locked normalized Stage 0 list: "
                f"configured={configured_hotspots_normalized} expected={expected_hotspots_normalized}"
            )
        if configured_hotspots_sha256 != expected_hotspots_sha256:
            raise RuntimeError(
                "Configured hotspot hash differs from the locked Stage 0 hash: "
                f"observed={configured_hotspots_sha256} expected={expected_hotspots_sha256}"
            )

    # Loop over number of designs to sample.
    design_startnum = sampler.inf_conf.design_startnum
    if sampler.inf_conf.design_startnum == -1:
        existing = glob.glob(sampler.inf_conf.output_prefix + "*.pdb")
        indices = [-1]
        for e in existing:
            print(e)
            m = re.match(".*_(\d+)\.pdb$", e)
            print(m)
            if not m:
                continue
            m = m.groups()[0]
            indices.append(int(m))
        design_startnum = max(indices) + 1

    for i_des in range(design_startnum, design_startnum + sampler.inf_conf.num_designs):
        if conf.inference.deterministic:
            make_deterministic(i_des)

        start_time = time.time()
        out_prefix = f"{sampler.inf_conf.output_prefix}_{i_des}"
        log.info(f"Making design {out_prefix}")
        if sampler.inf_conf.cautious and os.path.exists(out_prefix + ".pdb"):
            log.info(
                f"(cautious mode) Skipping this design because {out_prefix}.pdb already exists."
            )
            continue

        x_init, seq_init = sampler.sample_init()
        model_runner_module = sys.modules[sampler.__class__.__module__]
        runtime_root = Path(iu.__file__).resolve().parents[2]
        runtime_paths = {
            "inference_utils": Path(iu.__file__).resolve(),
            "model_runners": Path(model_runner_module.__file__).resolve(),
            "run_inference": Path(__file__).resolve(),
            "util": Path(rfd_util.__file__).resolve(),
        }
        runtime_audit = {
            "runtime_root": str(runtime_root),
            "runtime_git_commit": subprocess.check_output(
                ["git", "-C", str(runtime_root), "rev-parse", "HEAD"], text=True
            ).strip(),
            "inference_entrypoint": str(runtime_paths["run_inference"]),
            "rfdiffusion_util": str(runtime_paths["util"]),
            "rfdiffusion_inference_utils": str(runtime_paths["inference_utils"]),
            "rfdiffusion_model_runners": str(runtime_paths["model_runners"]),
            "runtime_source_sha256": {
                key: sha256_file(path) for key, path in runtime_paths.items()
            },
            "binderlen": int(sampler.binderlen),
            "hotspot_0idx": [int(value) for value in (sampler.hotspot_0idx or [])],
            "contigmap_derived_hotspot_0idx": [],
            "model_hotspot_tensor_0idx": [],
            "chain_idx": [str(value) for value in sampler.chain_idx],
            "idx_pdb": [int(value) for value in sampler.idx_pdb],
            "cyclic_residue_indices": [
                int(value)
                for value in torch.where(sampler.cyclic_reses)[0].detach().cpu().tolist()
            ],
            "input_pdb": str(input_pdb_path),
            "input_pdb_sha256": input_pdb_sha256,
            "stage0_target_pdb_sha256": input_pdb_sha256,
            "stage0_mapping_csv": stage0_mapping_csv,
            "stage0_mapping_csv_sha256": stage0_mapping_csv_sha256,
            "stage0_hotspots_normalized": configured_hotspots_normalized,
            "stage0_hotspots_sha256": configured_hotspots_sha256,
            "provenance_closure_required": provenance_required,
            "requested_hotspots": [str(value) for value in (sampler.ppi_conf.hotspot_res or [])],
            "contigs": [str(value) for value in (sampler.contig_conf.contigs or [])],
            "cyclic": bool(sampler.inf_conf.cyclic),
            "cyc_chains": str(sampler.inf_conf.cyc_chains),
        }
        denoised_xyz_stack = []
        px0_xyz_stack = []
        seq_stack = []
        plddt_stack = []

        x_t = torch.clone(x_init)
        seq_t = torch.clone(seq_init)
        # Loop over number of reverse diffusion time steps.
        for t in range(int(sampler.t_step_input), sampler.inf_conf.final_step - 1, -1):
            px0, x_t, seq_t, plddt = sampler.sample_step(
                t=t, x_t=x_t, seq_init=seq_t, final_step=sampler.inf_conf.final_step
            )
            px0_xyz_stack.append(px0)
            denoised_xyz_stack.append(x_t)
            seq_stack.append(seq_t)
            plddt_stack.append(plddt[0])  # remove singleton leading dimension

        runtime_audit["contigmap_derived_hotspot_0idx"] = [
            int(value) for value in sampler.contigmap_derived_hotspot_0idx
        ]
        runtime_audit["model_hotspot_tensor_0idx"] = [
            int(value) for value in sampler.model_hotspot_tensor_0idx
        ]
        expected_runtime_hotspots = sorted(runtime_audit["hotspot_0idx"])
        if runtime_audit["contigmap_derived_hotspot_0idx"] != expected_runtime_hotspots:
            raise RuntimeError(
                "Final ContigMap hotspot audit disagrees with helper mapping: "
                f"contigmap={runtime_audit['contigmap_derived_hotspot_0idx']} "
                f"helper={expected_runtime_hotspots}"
            )
        if runtime_audit["model_hotspot_tensor_0idx"] != expected_runtime_hotspots:
            raise RuntimeError(
                "Final model hotspot tensor audit disagrees with helper mapping: "
                f"tensor={runtime_audit['model_hotspot_tensor_0idx']} "
                f"helper={expected_runtime_hotspots}"
            )
        runtime_audit["runtime_audit_finalized_after_model_hotspot_tensor"] = True

        # Flip order for better visualization in pymol
        denoised_xyz_stack = torch.stack(denoised_xyz_stack)
        denoised_xyz_stack = torch.flip(
            denoised_xyz_stack,
            [
                0,
            ],
        )
        px0_xyz_stack = torch.stack(px0_xyz_stack)
        px0_xyz_stack = torch.flip(
            px0_xyz_stack,
            [
                0,
            ],
        )

        # For logging -- don't flip
        plddt_stack = torch.stack(plddt_stack)

        # Save outputs
        os.makedirs(os.path.dirname(out_prefix), exist_ok=True)
        final_seq = seq_stack[-1]

        # Output glycines, except for motif region
        final_seq = torch.where(
            torch.argmax(seq_init, dim=-1) == 21, 7, torch.argmax(seq_init, dim=-1)
        )  # 7 is glycine

        bfacts = torch.ones_like(final_seq.squeeze())
        # make bfact=0 for diffused coordinates
        bfacts[torch.where(torch.argmax(seq_init, dim=-1) == 21, True, False)] = 0
        # pX0 last step
        out = f"{out_prefix}.pdb"

        # Now don't output sidechains
        writepdb(
            out,
            denoised_xyz_stack[0, :, :4],
            final_seq,
            sampler.binderlen,
            chain_idx=sampler.chain_idx,
            bfacts=bfacts,
            idx_pdb=sampler.idx_pdb
        )

        # run metadata
        trb = dict(
            config=OmegaConf.to_container(sampler._conf, resolve=True),
            plddt=plddt_stack.cpu().numpy(),
            device=torch.cuda.get_device_name(torch.cuda.current_device())
            if torch.cuda.is_available()
            else "CPU",
            time=time.time() - start_time,
            runtime_audit=runtime_audit,
        )
        if hasattr(sampler, "contig_map"):
            for key, value in sampler.contig_map.get_mappings().items():
                trb[key] = value
        audit_json_path = f"{out_prefix}.runtime_audit.json"
        with open(audit_json_path, "w", encoding="utf-8") as audit_handle:
            json.dump(runtime_audit, audit_handle, indent=2, sort_keys=True)
            audit_handle.write("\n")
        with open(f"{out_prefix}.trb", "wb") as f_out:
            pickle.dump(trb, f_out)
        with open(audit_json_path, "r", encoding="utf-8") as audit_handle:
            json_audit_roundtrip = json.load(audit_handle)
        with open(f"{out_prefix}.trb", "rb") as f_in:
            trb_audit_roundtrip = pickle.load(f_in).get("runtime_audit")
        if json_audit_roundtrip != trb_audit_roundtrip:
            raise RuntimeError(
                "Final runtime audit JSON and TRB runtime_audit differ; output provenance is invalid."
            )
        log.info(
            "Final runtime audit: binderlen=%s helper=%s contigmap=%s tensor=%s chain_counts=%s",
            runtime_audit["binderlen"],
            runtime_audit["hotspot_0idx"],
            runtime_audit["contigmap_derived_hotspot_0idx"],
            runtime_audit["model_hotspot_tensor_0idx"],
            {
                chain: runtime_audit["chain_idx"].count(chain)
                for chain in sorted(set(runtime_audit["chain_idx"]))
            },
        )

        if sampler.inf_conf.write_trajectory:
            # trajectory pdbs
            traj_prefix = (
                os.path.dirname(out_prefix) + "/traj/" + os.path.basename(out_prefix)
            )
            os.makedirs(os.path.dirname(traj_prefix), exist_ok=True)

            out = f"{traj_prefix}_Xt-1_traj.pdb"
            writepdb_multi(
                out,
                denoised_xyz_stack,
                bfacts,
                final_seq.squeeze(),
                use_hydrogens=False,
                backbone_only=False,
                chain_ids=sampler.chain_idx,
            )

            out = f"{traj_prefix}_pX0_traj.pdb"
            writepdb_multi(
                out,
                px0_xyz_stack,
                bfacts,
                final_seq.squeeze(),
                use_hydrogens=False,
                backbone_only=False,
                chain_ids=sampler.chain_idx,
            )

        if conf.inference.empty_cache_per_design and torch.cuda.is_available():
            torch.cuda.empty_cache()

        log.info(f"Finished design in {(time.time()-start_time)/60:.2f} minutes")


if __name__ == "__main__":
    main()
