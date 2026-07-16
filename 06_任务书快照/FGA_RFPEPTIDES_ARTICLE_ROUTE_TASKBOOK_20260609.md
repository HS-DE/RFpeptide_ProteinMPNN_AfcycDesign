# FGA RFpeptides Article-Style Route Taskbook

Date: 2026-06-09

Status: draft for discussion. This document records a proposed RFpeptides branch
that follows the logic of the RFpeptides macrocycle design paper. It is not yet
an approved execution SOP.

Clean-restart decision:

This RFpeptides article-style branch should be treated as a clean route. Do not
reuse previous ColabDesign candidates, Boltz / AlphaFold-Multimer prediction
results, or previous Patch_A / Patch_B definitions as direct inputs. They may
be inspected as historical context, but the RFpeptides route must independently
redefine target sites, target crops, hotspots, and numbering mappings before
any RFpeptides backbone generation.

Project roots:

```text
Main project:      C:\SH\fga_cyclic_peptide_design
WSL project path:  /mnt/c/SH/fga_cyclic_peptide_design
Upstream clone:    C:\SH\peptide_str\rfd_macro
WSL upstream path: /mnt/c/SH/peptide_str/rfd_macro
```

Boss-facing flowchart:

```text
FGA_RFPEPTIDES_ARTICLE_ROUTE_FLOWCHART_20260616.svg
FGA_RFPEPTIDES_ARTICLE_ROUTE_FLOWCHART_20260616.md
FGA_RFPEPTIDES_SITE2_SELECTION_PITCH_20260616.pptx
```

Primary reference:

- Paper: Accurate de novo design of high-affinity protein-binding macrocycles
  using deep learning.
- Upstream implementation: RFdiffusion / RFpeptides in the RosettaCommons
  RFdiffusion repository.
- Local code currently cloned at `C:\SH\peptide_str\rfd_macro`.

Important local evidence from the cloned repository:

- `rfd_macro/README.md` contains a dedicated section named
  `Macrocyclic peptide design with RFpeptides`.
- `rfd_macro/examples/design_macrocyclic_binder.sh` gives the binder example.
- The binder example uses:
  - `contigmap.contigs=[12-18 A3-117/0]`
  - `inference.cyclic=True`
  - `diffuser.T=50`
  - `inference.cyc_chains='a'`
  - 3-6 style target hotspots.
- The example command is not automatically valid for FGA. Before running any
  FGA job, local example behavior must be checked for contig interpretation,
  generated chain naming, `cyc_chains`, and hotspot residue numbering.

## 1. Branch Decision Boundary

This branch is not a Boltz branch.

For this article-style RFpeptides route, Boltz-2 should not be part of the main
decision logic. Boltz can remain as a separate independent comparison branch,
but it should not define whether an RFpeptides design succeeded.

Reason:

- RFpeptides proposes a target-bound macrocycle backbone.
- ProteinMPNN designs sequences for that backbone.
- The article-style validation question is whether downstream structure
  validation can recover the RFpeptides design pose and target interface.
- If Boltz independently predicts a different binding pose, that may be useful
  information, but it no longer validates the RFpeptides design hypothesis.

## 2. Scientific Goal

Design macrocyclic peptide binders against newly selected native fibrinogen /
FGA exposed target sites, using the RFpeptides macrocycle protocol as the
upstream generation method.

Target-site policy:

- Previous `Patch_A` / `Patch_B` names are legacy labels from the ColabDesign
  branch.
- For this RFpeptides article-style branch, target sites must be selected again
  from the native fibrinogen/FGA structure.
- Until that target-site rediscovery is complete, use provisional labels such
  as `RFpep_Site_1` and `RFpep_Site_2`.
- If the rediscovered sites overlap old Patch_A / Patch_B, record that as a
  convergence observation, not as an assumption.

Patch_C and all old patch labels remain historical ColabDesign branch concepts
unless deliberately reintroduced after RFpeptides-specific target-site review.

The output of this branch is not just a sequence list. The intended output is:

```text
RFpeptides backbone design
+ ProteinMPNN sequence
+ relaxed design model
+ AfCycDesign validation model, or documented RF2/AF2 initial-guess fallback
+ target-site/contact/RMSD/interface scores
+ final ranked candidate table
```

## 3. Core Article-Style Logic

The proposed workflow is:

0. Rediscover RFpeptides-specific FGA target sites from native fibrinogen/FGA
   structure and biological/design rationale.
1. Prepare FGA target crops around the rediscovered RFpeptides target sites.
2. Use RFpeptides / RFdiffusion to generate target-bound macrocyclic peptide
   backbones.
3. Collect generated backbone PDB/TRB files and verify basic geometry.
4. Use ProteinMPNN to design sequences on RFpeptides backbones.
5. Relax and score sequence-designed complexes with Rosetta-style scoring.
6. Validate the sequence-designed complexes with AfCycDesign as the primary
   validation route.
7. If AfCycDesign cannot be installed or cannot run correctly, use an RF2 or
   AF2 initial-guess style validation route as a documented fallback.
8. Compare validation structures back to RFpeptides design poses.
9. Rank only candidates that preserve the selected target site and macrocycle
   geometry.

This differs from the ColabDesign route:

- ColabDesign optimized peptide sequence inside an AlphaFold-derived design
  loop.
- RFpeptides first generates a macrocycle-target structural hypothesis.
- ProteinMPNN and validation models then test whether that hypothesis can be
  supported by an actual sequence.

## 4. Non-Negotiable No-Shortcut Rule

The project-level testing rule remains active:

- Small test = fewer designs/jobs/target sites.
- Bad shortcut = incomplete or scientifically misleading method.

For this branch:

- It is acceptable to run only a few RFpeptides designs.
- It is not acceptable to lower `diffuser.T` below the article/repository
  protocol merely to save time.
- `diffuser.T=15` is only an installation smoke test because RFdiffusion rejects
  smaller values.
- Discussion/pilot/production runs should use `diffuser.T=50` unless a written
  reason is recorded.
- ProteinMPNN, relaxation/scoring, and AfCycDesign-first validation should not
  be skipped if the result is being interpreted as a candidate.
- RF2/AF2 initial-guess style validation is a fallback only after AfCycDesign is
  shown to be unavailable, not an equally preferred first choice.

### Large installation and download rule

To avoid wasting interactive time and context on long transfers or environment
installation, large setup operations are user-run by default.

- Codex should not automatically download or install large model weights,
  containers, datasets, repositories, CUDA/toolchain packages, or other
  time-consuming dependencies.
- Codex should first provide complete WSL/Linux commands, expected paths,
  approximate purpose, verification commands, and any relevant version pins.
- The user runs the large installation or download manually and reports the
  output. Codex then verifies files, versions, imports, GPU access, and protocol
  compatibility, and continues debugging from the concrete result.
- Lightweight read-only checks and small dependency fixes may be performed by
  Codex when useful, but a potentially large or long-running operation requires
  explicit user confirmation first.
- This rule changes who performs environment setup; it does not permit reduced
  scientific parameters or incomplete validation.

## 5. Branch Isolation

Do not overwrite existing branches:

```text
results/raw_designs/colabdesign_outputs/
results/boltz_predictions_*/
results/colabfold_predictions_*/
```

Recommended new result root:

```text
C:\SH\fga_cyclic_peptide_design\results\rfpeptides_article_route_clean_20260612\
```

Recommended subdirectories:

```text
results/rfpeptides_article_route_clean_20260612/
  00_site_discovery/
  00_target_inputs/
  01_rfpeptides_jobs/
  02_rfpeptides_backbones/
  03_backbone_qc/
  04_proteinmpnn_inputs/
  05_proteinmpnn_sequences/
  06_rosetta_relax/
  07_afcycdesign_validation/
  07b_rf2_af2ig_fallback_validation/
  08_validation_parse/
  09_ranked_candidates/
  logs/
  reports/
```

## 6. Planned Project-Side Scripts

These scripts do not need to exist before this taskbook is reviewed. They are
the proposed integration points.

```text
scripts/18_discover_rfpeptides_target_sites.py
scripts/19_prepare_rfpeptides_article_inputs.py
scripts/20_make_rfpeptides_article_jobs.py
scripts/external/run_rfpeptides_article_batch.sh
scripts/21_collect_rfpeptides_backbones.py
scripts/22_prepare_proteinmpnn_jobs.py
scripts/external/run_proteinmpnn_batch.sh
scripts/23_collect_proteinmpnn_sequences.py
scripts/24_prepare_rosetta_relax_jobs.py
scripts/external/run_rosetta_relax_batch.sh
scripts/25_parse_rosetta_scores.py
scripts/26_prepare_afcycdesign_jobs.py
scripts/external/run_afcycdesign_independent_recovery.py
scripts/26b_prepare_rf2_af2ig_fallback_jobs.py
scripts/27_parse_article_route_validation.py
scripts/28_rank_rfpeptides_article_candidates.py
```

The numbering intentionally starts after the existing ColabDesign / Boltz /
ColabFold scripts.

## 7. Stage -1: RFpeptides Target-Site Rediscovery

Goal:

Select RFpeptides-specific target sites from the native fibrinogen/FGA structure
without inheriting the old ColabDesign Patch_A / Patch_B definitions.

Hard rule:

Old patch files are historical references only. They may be opened to understand
what happened previously, but the RFpeptides target-site decision must be made
again and documented from first principles.

Allowed historical references:

```text
data/annotations/FGA_epitope_candidates.csv
data/annotations/FGA_surface_residues.csv
data/annotations/FGA_structure_mapping.csv
results/boltz_predictions_*/
results/colabfold_predictions_*/
ColabDesign candidate reports
```

These are not allowed as direct inputs for target-site selection:

```text
previous Patch_A / Patch_B as fixed design sites
previous ColabDesign candidate positions
previous Boltz or AFM predicted binding poses
```

Primary inputs for rediscovery:

```text
data/structures/raw/3GHG.pdb
data/structures/prepared/fibrinogen_3GHG_clean.pdb
FGA UniProt mapping
native fibrinogen structural context
surface exposure / neighbor-count / SASA-style analysis
manual review of feasible macrocycle binding pockets
```

Target-site selection criteria:

- visible in 3GHG or otherwise structurally justifiable;
- exposed enough for macrocycle access;
- contains a compact set of residues suitable for 3-6 RFpeptides hotspots;
- has enough local target context to preserve structure after crop;
- avoids relying on flexible/unresolved regions unless explicitly marked
  exploratory;
- is not chosen only because old ColabDesign/Boltz/AFM runs produced output
  there.

Practical Stage -1 work order:

1. Rebuild a visible FGA residue inventory from the raw and prepared 3GHG
   structures.
2. Calculate or summarize residue-level exposure and local neighbor context.
3. Cluster accessible FGA surface residues into candidate macrocycle-binding
   target sites without using old Patch_A / Patch_B labels as seeds.
4. Optionally run fpocket on the full native cleaned PDB as supplementary
   pocket/groove evidence. This cross-check is not a hard reject condition and
   does not replace RSA, geometry, occlusion, chemical-anchor, or manual
   review.
5. Manually review candidate sites in the native fibrinogen context for
   macrocycle accessibility and crop feasibility.
6. Pick the first one or two RFpeptides target sites for pilot work and record
   rejected or deferred sites.
7. Assign provisional labels such as `RFpep_Site_1` and `RFpep_Site_2` only
   after the selection rationale is written.

Stage -1 WSL environment after the 2026-06-15 boss-review adjustment:

Use WSL for this step so Python, FreeSASA, and Linux `fpocket` are resolved
inside the same environment. First-time setup only:

```bash
cd /mnt/c/SH/fga_cyclic_peptide_design
source ~/fga_model_envs/miniforge3/etc/profile.d/conda.sh

conda create -n fga_stage1_fpocket -c conda-forge python=3.11 pyyaml freesasa fpocket -y
conda activate fga_stage1_fpocket

which python
which fpocket
```

Routine Stage -1 command:

```bash
cd /mnt/c/SH/fga_cyclic_peptide_design
source ~/fga_model_envs/miniforge3/etc/profile.d/conda.sh
conda activate fga_stage1_fpocket

python scripts/18_discover_rfpeptides_target_sites.py \
  --config config/project.yaml \
  --output-root results/rfpeptides_article_route_clean_20260615_fpocket \
  --max-candidates 100 \
  --propose-sites 10 \
  --hotspots-per-site 4 \
  --rsa-surface-threshold 0.20 \
  --proposal-max-uniprot-overlap-fraction 0.20 \
  --proposal-min-uniprot-center-distance 25 \
  --proposal-min-center-distance 18 \
  --enable-fpocket \
  --fpocket-bin fpocket
```

This command only performs target-site rediscovery and site quality scoring. It
does not prepare RFpeptides target crops and does not run RFpeptides backbone
generation.

If `which fpocket` fails inside the active WSL environment, the script records
`fpocket_status=not_available` and still writes the original Stage -1 outputs.
In that case, FreeSASA/RSA, geometry, occlusion, chemical-anchor, and PyMOL
review remain usable; fpocket can be installed and rerun later as supplementary
evidence.

Expected outputs:

```text
00_site_discovery/FGA_rfpeptides_target_site_candidates.csv
00_site_discovery/FGA_rfpeptides_target_site_selection.md
00_site_discovery/FGA_rfpeptides_target_site_quality_review.md
00_site_discovery/FGA_rfpeptides_rejected_or_deferred_sites.csv
00_site_discovery/RFpep_Site_*_review.pml
00_site_discovery/fpocket_native_context/
```

Suggested target-site candidate columns:

```text
site_id
site_label
original_chain_ids
original_pdb_residue_numbers
uniprot_residue_numbers
center_x
center_y
center_z
exposure_summary
structural_context_summary
macrocycle_accessibility_rationale
initial_hotspot_candidates
selection_status
selection_notes
legacy_patch_overlap
sasa_status
rsa_status
fpocket_status
fpocket_error_summary
nearest_fpocket_pocket_id
nearest_fpocket_distance_to_hotspot
nearest_fpocket_distance_to_site
nearest_fpocket_pocket_score
nearest_fpocket_druggability_score
nearest_fpocket_volume
fpocket_support_status
fpocket_support_reason
site_quality_tier
site_quality_reason
stage0_crop_allowed
```

fpocket support interpretation:

```text
strong_support: pocket distance to hotspot <= 6 A
moderate_support: pocket distance to hotspot <= 10 A, or pocket distance to site <= 6 A
weak_support: pocket distance to site <= 12 A
no_nearby_pocket: no nearby pocket/groove detected by fpocket
not_evaluated: fpocket was not run, unavailable, or failed
```

fpocket status interpretation:

```text
not_run: fpocket was not requested
not_available: fpocket executable was not found; script continues
failed: fpocket ran but failed; script continues and records an error summary
completed: fpocket completed and parsed results were attached to candidates
```

Stage -1 completion gate:

Stage -1 is complete only when the selected RFpeptides target sites have been
named and documented. Only after this point may the project use stable labels
such as:

```text
RFpep_Site_1
RFpep_Site_2
```

If a selected site overlaps old Patch_A or Patch_B, record:

```text
legacy_patch_overlap=Patch_A or Patch_B
```

but do not rename it as Patch_A / Patch_B unless the review explicitly decides
to preserve the old label.

Stage -1 quality gate:

- Only `site_quality_tier=high` or `site_quality_tier=medium` may enter Stage 0
  target crop preparation.
- `site_quality_tier=low` and `site_quality_tier=reject` remain review records
  only and must not be used to make RFpeptides target inputs.
- `legacy_patch_overlap` is historical context only and must not be used as a
  positive or negative quality score.
- If FreeSASA is unavailable, `sasa_status=proxy_only` must be recorded and the
  neighbor-count exposure proxy should be treated as lower confidence.
- fpocket is supplementary pocket/groove evidence only. It may help prioritize
  manual review, but it must not be used as a hard reject, must not change
  `site_quality_tier`, and must not be reported as final validation.

Current Stage -1 selection decision:

2026-06-15 update: after boss review, rerun Stage -1 with optional fpocket
cross-check before finalizing whether the 2026-06-12 Site_2-only pilot remains
the best first Stage 1 target choice.

For the clean 2026-06-12 RFpeptides restart, Stage 0 target inputs were
prepared for the three automatically proposed `high` sites:

```text
RFpep_Site_2 = RFpep_Candidate_004, chain G, hotspots G199,G196,G200,G198
RFpep_Site_3 = RFpep_Candidate_010, chain A, hotspots A30,A32,A29,A28
RFpep_Site_4 = RFpep_Candidate_015, chain J, hotspots J99,J103,J106,J98
```

`RFpep_Site_1`, `RFpep_Site_5`, and `RFpep_Site_6` remain documented
proposed/deferred alternatives for this pilot decision, but they are not the
first Stage 0 target-input set.

Current Stage 1 pilot decision:

```text
Proceed to RFpeptides pilot generation with RFpep_Site_2 only.
RFpep_Site_3 and RFpep_Site_4 have Stage 0 files but are deferred for now.
```

## 8. Stage 0: Target Preparation

Goal:

Create RFpeptides-compatible target PDBs for the newly selected RFpeptides
target sites.

Hard rule:

RFpeptides target sites are target regions, not hotspot lists. A site is the
broader surface area we want the macrocycle to bind. Hotspots are a small,
deliberately chosen subset of residues used to guide RFpeptides toward that
site. Do not pass every site residue as a hotspot, and do not treat a
hotspot-only contact as proof that the whole site objective was satisfied.

Inputs:

```text
data/structures/raw/3GHG.pdb
data/structures/prepared/fibrinogen_3GHG_clean.pdb
00_site_discovery/FGA_rfpeptides_target_site_candidates.csv
00_site_discovery/FGA_rfpeptides_target_site_selection.md
```

Do not use the old patch table as the Stage 0 driver. The old patch table may
be referenced only for historical comparison.

Preparation requirements:

- Preserve the native fibrinogen structural context needed around the selected
  target site.
- Crop target enough to keep runtime manageable.
- Leave approximately 10 A of target context around the intended binding site
  where possible.
- Avoid arbitrary residue renumbering unless the mapping is recorded.
- Export a crop-renumbering mapping table from RFpeptides target residue numbers
  back to original 3GHG residue IDs and UniProt positions.
- Export a hotspot residue selection rationale for each selected target site.
- Explain why each selected hotspot residue is suitable for RFpeptides guidance:
  surface exposure, location inside the intended target site, distance from
  unstable loops or chain breaks when possible, and expected contribution to
  target-site steering.
- Keep unselected target-site residues available for downstream contact
  recovery checks; they are not discarded just because they are not hotspots.
- Confirm chain IDs used by RFpeptides after cropping and renumbering.

Expected outputs:

```text
00_target_inputs/RFpep_Site_2_target.pdb
00_target_inputs/RFpep_Site_2_crop_renumbering_mapping.csv
00_target_inputs/RFpep_Site_2_hotspots.txt
00_target_inputs/RFpep_Site_2_hotspot_selection_rationale.md
00_target_inputs/RFpep_Site_2_site_residues.csv
00_target_inputs/RFpep_Site_3_target.pdb
00_target_inputs/RFpep_Site_3_crop_renumbering_mapping.csv
00_target_inputs/RFpep_Site_3_hotspots.txt
00_target_inputs/RFpep_Site_3_hotspot_selection_rationale.md
00_target_inputs/RFpep_Site_3_site_residues.csv
00_target_inputs/RFpep_Site_4_target.pdb
00_target_inputs/RFpep_Site_4_crop_renumbering_mapping.csv
00_target_inputs/RFpep_Site_4_hotspots.txt
00_target_inputs/RFpep_Site_4_hotspot_selection_rationale.md
00_target_inputs/RFpep_Site_4_site_residues.csv
```

Planned Stage 0 command:

```bash
cd /mnt/c/SH/fga_cyclic_peptide_design
source ~/fga_model_envs/miniforge3/etc/profile.d/conda.sh
conda activate fga_stage1_fpocket

python scripts/19_prepare_rfpeptides_article_inputs.py \
  --config config/project.yaml \
  --input-root results/rfpeptides_article_route_clean_20260615_fpocket \
  --output-root results/rfpeptides_article_route_clean_20260615_fpocket \
  --selected-sites RFpep_Site_2 \
  --crop-context-radius 10 \
  --hotspots-per-site 4
```

Stage 0 command status:

- `scripts/19_prepare_rfpeptides_article_inputs.py` is implemented and this
  command was rerun for `RFpep_Site_2` only on 2026-06-16 after the fpocket
  Stage -1 update.
- Current Stage 0 result: crop residues=86, atoms=712, hotspots
  `A85,A82,A86,A84`.
- The command must read only the clean-restart Stage -1 outputs from
  `00_site_discovery/`.
- It must not read the old Patch_A / Patch_B / Patch_C table as the driver for
  target crop generation.
- It must not launch RFpeptides backbone generation.

Stage 0 decision note:

- We considered returning to Stage -1 for a `RFpep_Site_2` same-chain check,
  meaning a review of chain-G candidates similar to `RFpep_Site_2`.
- The same-chain alternatives were manually reviewed and did not look more
  suitable than the current `RFpep_Site_2` representative.
- No additional same-chain Site_2-family Stage 0 preparation will be run for
  now. Cross-chain family comparisons may be revisited later only if the pilot
  target set needs replacement or expansion.

Minimum mapping table columns:

```text
rfpeptides_chain
rfpeptides_residue_number
original_pdb_id
original_chain_id
original_pdb_residue_number
original_insertion_code
uniprot_accession
uniprot_residue_number
is_target_site_residue
is_selected_hotspot
selection_or_exclusion_note
```

Stage 0 completion gate:

Stage 0 is complete only when all selected RFpeptides target sites have:

- target PDB;
- crop-renumbering mapping table;
- target-site residue table;
- hotspot list;
- hotspot selection rationale;
- reviewed contig target residue range;
- reviewed chain ID that will be used by `cyc_chains` and `ppi.hotspot_res`.

No RFpeptides backbone generation is allowed before this gate is complete.

Open decision:

- Whether RFpeptides should use the full visible FGA/fibrinogen local context
  or a tighter site-centered crop.

## 9. Stage 1: RFpeptides Backbone Generation

Goal:

Generate macrocyclic peptide backbones bound to the selected RFpeptides target
sites.

Stage 0 dependency:

Do not generate RFpeptides backbones until Stage -1 target-site rediscovery and
Stage 0 target inputs, hotspot rationales, and crop-renumbering mappings have
all been completed and reviewed.

Core command shape is adapted from the local RFpeptides binder example:
`C:/SH/peptide_str/rfd_macro/examples/design_macrocyclic_binder.sh`.

Stage 1 job-preparation status:

- `scripts/20_make_rfpeptides_article_jobs.py` is implemented.
- It prepares a command table and run script only; it does not run RFpeptides.
- This preparation step was rerun on 2026-06-16 against
  `results/rfpeptides_article_route_clean_20260615_fpocket` after the Site_2
  Stage 0 refresh.
- The current Stage 1 pilot job table contains only `RFpep_Site_2`.
- Stage 1 generation was run on 2026-06-16 and completed all 10 requested
  RFpep_Site_2 backbone designs: 10 `.pdb` files, 10 `.trb` files, and 20
  trajectory files were written under `02_rfpeptides_backbones/RFpep_Site_2/`.
- No traceback, runtime error, CUDA out-of-memory, or killed-process message
  was found in the Stage 1 log at completion.
- A larger 1000-design RFpep_Site_2 Stage 1 run script was prepared on
  2026-06-16 in a separate output root:
  `results/rfpeptides_article_route_clean_20260615_fpocket_stage1_N1000_no_traj`.
  This keeps the completed 10-design pilot untouched and disables trajectory
  output with `inference.write_trajectory=False` to reduce disk usage.
- The 1000-design RFpep_Site_2 scale-up run completed on 2026-06-17:
  1000 `.pdb` files and 1000 `.trb` files were written, indices 0-999 are
  complete for both file types, no trajectory directory was produced, and no
  traceback, runtime error, CUDA out-of-memory, killed-process, or generic
  error line was found in the Stage 1 scale-up log.
- The generated command uses the real Stage 0 target range and hotspots:
  `contigmap.contigs=[12-18 A1-86/0]` and
  `ppi.hotspot_res=[A85,A82,A86,A84]`.
- The generated macrocycle chain is still expected to be controlled by
  `inference.cyc_chains='a'`, matching the local RFpeptides example; this must
  be confirmed in the first parsed RFpeptides outputs.

Stage 1 job-preparation command:

```bash
cd /mnt/c/SH/fga_cyclic_peptide_design
source ~/fga_model_envs/miniforge3/etc/profile.d/conda.sh
conda activate fga_stage1_fpocket

python scripts/20_make_rfpeptides_article_jobs.py \
  --input-root results/rfpeptides_article_route_clean_20260615_fpocket \
  --output-root results/rfpeptides_article_route_clean_20260615_fpocket \
  --rfpeptides-root /mnt/c/SH/peptide_str/rfd_macro \
  --selected-sites RFpep_Site_2 \
  --num-designs 10 \
  --length-min 12 \
  --length-max 18 \
  --diffuser-t 50
```

Generated run script:

```text
results/rfpeptides_article_route_clean_20260615_fpocket/01_rfpeptides_jobs/run_rfpeptides_stage1_site2.sh
```

Prepared larger run script for overnight scale-up:

```text
results/rfpeptides_article_route_clean_20260615_fpocket_stage1_N1000_no_traj/01_rfpeptides_jobs/run_rfpeptides_stage1_site2_N1000_no_traj.sh
```

Stage 1 N1000 job-preparation command:

```bash
cd /mnt/c/SH/fga_cyclic_peptide_design
source ~/fga_model_envs/miniforge3/etc/profile.d/conda.sh
conda activate fga_stage1_fpocket

python scripts/20_make_rfpeptides_article_jobs.py \
  --input-root results/rfpeptides_article_route_clean_20260615_fpocket \
  --output-root results/rfpeptides_article_route_clean_20260615_fpocket_stage1_N1000_no_traj \
  --rfpeptides-root /mnt/c/SH/peptide_str/rfd_macro \
  --selected-sites RFpep_Site_2 \
  --num-designs 1000 \
  --length-min 12 \
  --length-max 18 \
  --diffuser-t 50 \
  --run-script-name run_rfpeptides_stage1_site2_N1000_no_traj.sh \
  --extra-override inference.write_trajectory=False
```

The generated run script currently contains:

```bash
cd /mnt/c/SH/peptide_str/rfd_macro
source ~/fga_model_envs/miniforge3/etc/profile.d/conda.sh
conda activate SE3nv

./scripts/run_inference.py --config-name base \
  inference.output_prefix=/mnt/c/SH/fga_cyclic_peptide_design/results/rfpeptides_article_route_clean_20260615_fpocket/02_rfpeptides_backbones/RFpep_Site_2/RFpep_Site_2_L12_18 \
  inference.num_designs=10 \
  'contigmap.contigs=[12-18 A1-86/0]' \
  inference.input_pdb=/mnt/c/SH/fga_cyclic_peptide_design/results/rfpeptides_article_route_clean_20260615_fpocket/00_target_inputs/RFpep_Site_2_target.pdb \
  inference.cyclic=True \
  diffuser.T=50 \
  inference.cyc_chains='a' \
  'ppi.hotspot_res=[A85,A82,A86,A84]'
```

Important:

- This is a small design-count pilot, not a reduced-constraint run.
- The command keeps target PDB, target contig, cyclic generation, cyclic chain,
  diffuser timesteps, and hotspot guidance.
- Potentials are off for the first baseline pilot.
- Do not treat any Stage 1 output as a final peptide; Stage 1 generates
  backbone candidates only.

First discussion pilot:

```text
RFpep_Site_2: 10-20 backbones
Length range: 12-18
diffuser.T: 50
potentials: off unless a written reason is recorded
hotspots: 3-6 residues per selected target site
```

Pilot scope decision:

```text
Run RFpep_Site_2 only for the first Stage 1 pilot.
RFpep_Site_3 is deferred because its current Stage 0 crop is very short.
RFpep_Site_4 is deferred as a backup target until Site_2 behavior is reviewed.
```

Why start with potentials off:

The RFdiffusion README notes that potentials can interact oddly with hotspot
PPI settings and recommends cautious case-by-case exploration. For FGA, the
first article-style pilot should establish the baseline without potentials.

Expected outputs:

```text
02_rfpeptides_backbones/RFpep_Site_2/*.pdb
02_rfpeptides_backbones/RFpep_Site_2/*.trb
logs/rfpeptides_stage1_*.log
```

## 10. Stage 2: Backbone QC

Goal:

Reject failed RFpeptides generations before sequence design.

Hard rule:

Stage 2 must not only prove that a macrocyclic peptide touches the Stage 0
target crop. It must prove that the peptide binds near the intended
`RFpep_Site_2` residues and/or the selected hotspot residues. A design that
contacts only a distant part of the A1-A86 target crop is not a target-site
recovery pass.

Minimum QC checks:

- Output PDB exists.
- Output TRB exists.
- Designed chain exists.
- Designed chain length is in the requested range.
- Designed chain is marked/treated as cyclic.
- Macrocycle geometry is evaluated as head-to-tail macrocycle geometry, not as
  Cys-Cys disulfide geometry.
- Peptide is near the target crop, not floating far away.
- Peptide contacts the intended RFpep_Site_2 target-site residues, not only
  arbitrary context residues in the crop.
- Peptide contacts or approaches the selected RFpeptides hotspot residues
  (`A85,A82,A86,A84`) closely enough to support hotspot recovery.
- Peptide-target contacts are checked separately for whole crop, target-site
  residues, and hotspot residues.
- Designs that satisfy only whole-crop contact but fail site/hotspot proximity
  are flagged as `crop_only_contact` and do not pass target-site recovery.
- No severe chain parsing or residue numbering ambiguity.

`macrocycle_geometry_status` definition:

This branch does not use any Cys-Cys / SG-SG logic. Do not look for terminal
Cys, SG atoms, or disulfide distance as a geometry gate. The geometry status
should describe whether the RFpeptides output is consistent with a macrocyclic
peptide generated by the RFpeptides cyclic protocol.

Suggested values:

```text
pass_head_to_tail_macrocycle
warn_chain_or_residue_numbering_unclear
warn_cyclic_metadata_missing_but_geometry_close
fail_open_chain_or_no_cyclic_evidence
fail_parse_error
```

Minimum geometry evidence:

- cyclic chain is the generated peptide chain requested by `cyc_chains`;
- designed chain length matches the requested macrocycle length range;
- N/C terminal backbone atoms are close enough to be consistent with the
  RFpeptides cyclic protocol, or the TRB/cyclic metadata explicitly supports
  cyclic generation;
- no unrelated disulfide-specific checks are used.

Suggested output table:

```text
03_backbone_qc/FGA_rfpeptides_backbones_qc.csv
```

Suggested columns:

```text
design_id
target_site_id
length
target_pdb
rf_pdb
trb
num_target_contacts
num_target_site_contacts
num_hotspot_contacts
peptide_target_min_distance
peptide_site_min_distance
peptide_hotspot_min_distance
target_contact_status
target_site_recovery_status
hotspot_recovery_status
macrocycle_geometry_status
pass_backbone_qc
qc_notes
```

Stage 2 implementation status:

- `scripts/21_collect_rfpeptides_backbones.py` is implemented.
- The script reads the Stage 0 crop-renumbering map, target-site residues, and
  hotspot residues, then scores each Stage 1 RFpeptides backbone against three
  separate contact layers: whole target crop, RFpep_Site_2 target-site
  residues, and selected hotspots.
- The script explicitly flags designs that contact only the crop but miss
  RFpep_Site_2 as `crop_only_contact`.
- The script uses head-to-tail macrocycle geometry, not Cys-Cys / SG-SG
  disulfide geometry.
- The generated PyMOL review script writes Windows-style `C:/...` paths by
  default so it can be opened directly in Windows PyMOL after WSL execution.

Stage 2 command used for the 1000-design RFpep_Site_2 run:

```bash
cd /mnt/c/SH/fga_cyclic_peptide_design
source ~/fga_model_envs/miniforge3/etc/profile.d/conda.sh
conda activate fga_stage1_fpocket

python scripts/21_collect_rfpeptides_backbones.py \
  --stage0-root results/rfpeptides_article_route_clean_20260615_fpocket \
  --stage1-root results/rfpeptides_article_route_clean_20260615_fpocket_stage1_N1000_no_traj \
  --selected-sites RFpep_Site_2
```

Stage 2 result snapshot from the 2026-06-17 run:

```text
total_backbones: 1000
pass_backbone_qc: 2
crop_only_contact: 368
hotspot_contact_pass: 2
hotspot_near_pass: 6
macrocycle_geometry_status: pass_head_to_tail_macrocycle for all 1000
```

Passing backbone IDs:

```text
RFpep_Site_2_0007
RFpep_Site_2_0394
```

Interpretation:

- The RFpeptides cyclic geometry was stable across the 1000-design run.
- Most generated backbones did not recover the intended Site_2 / hotspot
  position; many contacted only the broader crop and are not valid target-site
  recovery passes.
- Downstream work will focus on `RFpep_Site_2_0007` as the primary Stage 2
  backbone. `RFpep_Site_2_0394` is retained only as a weak backup / manual
  review structure.
- The `RFpep_Site_2_0007` peptide chain is currently a 12-residue RFpeptides
  backbone with poly-Gly placeholder sequence (`GGGGGGGGGGGG`). It is not a
  designed peptide sequence yet; sequence design starts in Stage 3.
- The two passing backbones are Stage 2 backbone candidates only. They are not
  final peptide candidates and still require sequence design, structure
  relaxation/scoring, validation, negative screening, and final ranking.

## 11. Stage 3: ProteinMPNN Sequence Design

Goal:

Assign sequences to RFpeptides backbones.

Reason:

RFdiffusion/RFpeptides generates backbone structure. The designed chain may be
poly-Gly or otherwise not a final sequence. ProteinMPNN is needed to create
sequences compatible with each backbone.

Recommended scope for first real pilot:

```text
ProteinMPNN sequences per RFpeptides backbone: 4-8
Fixed target chain sequence: yes
Design peptide chain only: yes
Preserve macrocycle chain length: yes
```

Current Stage 3 implementation status:

- `scripts/22_prepare_proteinmpnn_jobs.py` is implemented.
- ProteinMPNN-FastRelax environment is available as
  `proteinmpnn_binder_design`.
- The `proteinmpnn_binder_design` environment should keep `numpy<2`
  (`numpy==1.26.4` tested), because the current Torch/ProteinMPNN stack is not
  compatible with NumPy 2.x.
- The local `dl_binder_design` checkout is at
  `/mnt/c/SH/peptide_str/dl_binder_design`.
- The first Stage 3 job is scoped to the primary Stage 2 backbone
  `RFpep_Site_2_0007`.
- The input complex keeps peptide chain `B` first and target crop chain `A`
  second, matching `dl_interface_design.py` expectations.
- `scripts/22_prepare_proteinmpnn_jobs.py` supports two Stage 3B modes:
  ProteinMPNN-FastRelax (`relax_cycles > 0`) and ProteinMPNN-only
  (`relax_cycles == 0`).
- In ProteinMPNN-FastRelax mode, `dl_interface_design.py` disallows
  `seqs_per_struct > 1`, so the preparation script duplicates
  `RFpep_Site_2_0007` into 8 independent input tags and runs one
  ProteinMPNN-FastRelax sequence per tag.
- In ProteinMPNN-only mode, the preparation script keeps one input copy of
  `RFpep_Site_2_0007` and requests multiple sequences with
  `seqs_per_struct`.
- The generated 3B run scripts create their checkpoint directories before
  launching ProteinMPNN.

Stage 3A input-preparation command for the original ProteinMPNN-FastRelax
route:

```bash
cd /mnt/c/SH/fga_cyclic_peptide_design
source ~/fga_model_envs/miniforge3/etc/profile.d/conda.sh
conda activate fga_stage1_fpocket

python scripts/22_prepare_proteinmpnn_jobs.py \
  --stage2-root results/rfpeptides_article_route_clean_20260615_fpocket_stage1_N1000_no_traj \
  --selected-backbones RFpep_Site_2_0007 \
  --dl-binder-design-root /mnt/c/SH/peptide_str/dl_binder_design \
  --seqs-per-backbone 8 \
  --relax-cycles 1 \
  --temperature 0.10
```

Stage 3A prepares inputs only. It creates:

```text
04_proteinmpnn_inputs/pdbs/
04_proteinmpnn_inputs/FGA_rfpeptides_stage3B_proteinmpnn_fastrelax_runlist.txt
04_proteinmpnn_inputs/FGA_rfpeptides_stage3B_proteinmpnn_fastrelax_jobs.csv
04_proteinmpnn_inputs/run_stage3B_proteinmpnn_fastrelax_site2_0007.sh
```

Stage 3B ProteinMPNN-FastRelax run command:

```bash
cd /mnt/c/SH/fga_cyclic_peptide_design

bash results/rfpeptides_article_route_clean_20260615_fpocket_stage1_N1000_no_traj/04_proteinmpnn_inputs/run_stage3B_proteinmpnn_fastrelax_site2_0007.sh
```

Note: the older generated script name below is also a Stage 3B run script, not
an input-preparation script:

```text
results/rfpeptides_article_route_clean_20260615_fpocket_stage1_N1000_no_traj/04_proteinmpnn_inputs/run_proteinmpnn_fastrelax_stage3_site2_0007.sh
```

Stage 3B output directory:

```text
05_proteinmpnn_sequences/fastrelax_pdbs/
```

Current Stage 3B result:

```text
RFpep_Site_2_0007 ProteinMPNN-FastRelax completed.
8 input tags were processed.
8 sequence-designed relaxed PDB files were produced.
```

Stage 3C collection/QC command:

```bash
cd /mnt/c/SH/fga_cyclic_peptide_design
source ~/fga_model_envs/miniforge3/etc/profile.d/conda.sh
conda activate fga_stage1_fpocket

python scripts/23_collect_proteinmpnn_sequences.py \
  --stage0-root results/rfpeptides_article_route_clean_20260615_fpocket \
  --stage3-root results/rfpeptides_article_route_clean_20260615_fpocket_stage1_N1000_no_traj \
  --selected-backbones RFpep_Site_2_0007
```

Stage 3C outputs:

```text
05_proteinmpnn_sequences/FGA_rfpeptides_stage3_sequences_qc.csv
05_proteinmpnn_sequences/FGA_rfpeptides_stage3_sequences_qc_pass.csv
05_proteinmpnn_sequences/FGA_rfpeptides_stage3_sequences_qc.md
05_proteinmpnn_sequences/fasta/
05_proteinmpnn_sequences/RFpep_Site_2_stage3C_sequence_qc_review.pml
```

Current Stage 3C result:

```text
Parsed Stage 3B relaxed PDBs: 8
Passed Stage 3C sequence/relax QC: 0
```

Interpretation:

- All 8 ProteinMPNN-FastRelax outputs have valid 12-aa peptide sequences.
- None preserves direct RFpep_Site_2 / hotspot contact after FastRelax.
- None preserves head-to-tail macrocycle geometry under the Stage 3C cutoff.
- Therefore, these 8 Stage 3B outputs should not enter Stage 4.

Stage 3 diagnostic control:

```text
05_proteinmpnn_sequences/diagnostic_fastrelax_control/
```

Control question: does the same Rosetta FastRelax mover damage the original
poly-Gly `RFpep_Site_2_0007` backbone even without ProteinMPNN sequence
redesign?

Result:

```text
before FastRelax:
  hotspot_min_distance: 3.157 A
  num_hotspot_contacts: 3
  terminal_C_N_distance: 1.283 A
  macrocycle_geometry_status: pass_head_to_tail_macrocycle

after same FastRelax on poly-Gly:
  hotspot_min_distance: 1.843 A
  num_hotspot_contacts: 2
  terminal_C_N_distance: 8.655 A
  macrocycle_geometry_status: fail_open_chain_or_no_cyclic_evidence
```

Interpretation:

- Same FastRelax alone preserves RFpep_Site_2 / hotspot contact in this
  poly-Gly control.
- Same FastRelax alone is sufficient to open the head-to-tail macrocycle
  geometry under the current unconstrained setup.
- The detached-and-open Stage 3B outputs therefore point to two problems:
  unconstrained FastRelax opens the ring, while ProteinMPNN sequence changes
  plus FastRelax can also move the peptide away from the intended hotspot.

Open decision:

- How to rerun Stage 3 so sequence generation preserves both hotspot contact
  and head-to-tail macrocycle geometry.

Conservative recommendation:

- Do not continue with unconstrained ProteinMPNN-FastRelax outputs from the
  current 3B run.
- Next Stage 3 rerun should either use ProteinMPNN-only on the RFpeptides
  backbone first, or use a constrained FastRelax route that explicitly preserves
  hotspot contact and head-to-tail macrocycle geometry.
- ProteinMPNN-only outputs are still intermediate sequence proposals, not final
  candidates.

Stage 3 rerun branch: ProteinMPNN-only.

Rationale: first test whether sequence assignment alone can preserve the
RFpeptides backbone pose, RFpep_Site_2 / hotspot proximity, and head-to-tail
macrocycle geometry. This branch intentionally skips FastRelax because the
diagnostic control showed that the current unconstrained FastRelax setup opens
the macrocycle.

Stage 3A ProteinMPNN-only preparation command:

```bash
cd /mnt/c/SH/fga_cyclic_peptide_design
source ~/fga_model_envs/miniforge3/etc/profile.d/conda.sh
conda activate fga_stage1_fpocket

python scripts/22_prepare_proteinmpnn_jobs.py \
  --stage2-root results/rfpeptides_article_route_clean_20260615_fpocket_stage1_N1000_no_traj \
  --output-root results/rfpeptides_article_route_clean_20260615_fpocket_stage1_N1000_no_traj \
  --selected-backbones RFpep_Site_2_0007 \
  --dl-binder-design-root /mnt/c/SH/peptide_str/dl_binder_design \
  --seqs-per-backbone 8 \
  --relax-cycles 0 \
  --temperature 0.10 \
  --omit-aas CX
```

ProteinMPNN-only preparation outputs:

```text
04_proteinmpnn_inputs/pdbs_proteinmpnn_only/
04_proteinmpnn_inputs/FGA_rfpeptides_stage3B_proteinmpnn_only_runlist.txt
04_proteinmpnn_inputs/FGA_rfpeptides_stage3B_proteinmpnn_only_jobs.csv
04_proteinmpnn_inputs/run_stage3B_proteinmpnn_only_site2_0007.sh
05_proteinmpnn_sequences/proteinmpnn_only_pdbs/
```

Stage 3B ProteinMPNN-only run command:

```bash
cd /mnt/c/SH/fga_cyclic_peptide_design

bash results/rfpeptides_article_route_clean_20260615_fpocket_stage1_N1000_no_traj/04_proteinmpnn_inputs/run_stage3B_proteinmpnn_only_site2_0007.sh
```

Stage 3C ProteinMPNN-only collection/QC command:

```bash
cd /mnt/c/SH/fga_cyclic_peptide_design
source ~/fga_model_envs/miniforge3/etc/profile.d/conda.sh
conda activate fga_stage1_fpocket

python scripts/23_collect_proteinmpnn_sequences.py \
  --stage0-root results/rfpeptides_article_route_clean_20260615_fpocket \
  --stage3-root results/rfpeptides_article_route_clean_20260615_fpocket_stage1_N1000_no_traj \
  --stage3-mode proteinmpnn_only \
  --selected-backbones RFpep_Site_2_0007
```

ProteinMPNN-only Stage 3C outputs:

```text
05_proteinmpnn_sequences/FGA_rfpeptides_stage3_proteinmpnn_only_sequences_qc.csv
05_proteinmpnn_sequences/FGA_rfpeptides_stage3_proteinmpnn_only_sequences_qc_pass.csv
05_proteinmpnn_sequences/FGA_rfpeptides_stage3_proteinmpnn_only_sequences_qc.md
05_proteinmpnn_sequences/fasta_proteinmpnn_only/
05_proteinmpnn_sequences/RFpep_Site_2_stage3C_proteinmpnn_only_sequence_qc_review.pml
```

Stage 3D-1 side-chain repack-only cleanup:

Rationale: ProteinMPNN-only preserved RFpep_Site_2 / hotspot contact and
head-to-tail macrocycle geometry, but all direct ProteinMPNN-only outputs
failed Stage 3C because of local severe clashes. Stage 3D-1 keeps the sequence
and backbone fixed, skips FastRelax, and only repacks peptide side chains plus
nearby target side chains.

Stage 3D-1 command:

```bash
cd /mnt/c/SH/fga_cyclic_peptide_design
source ~/fga_model_envs/miniforge3/etc/profile.d/conda.sh
conda activate proteinmpnn_binder_design

python scripts/24_stage3d1_sidechain_repack.py \
  --stage0-root results/rfpeptides_article_route_clean_20260615_fpocket \
  --stage3-root results/rfpeptides_article_route_clean_20260615_fpocket_stage1_N1000_no_traj \
  --selected-backbones RFpep_Site_2_0007
```

Stage 3D-1 outputs:

```text
05_proteinmpnn_sequences/stage3d1_repack_only_pdbs/
05_proteinmpnn_sequences/FGA_rfpeptides_stage3D1_sidechain_repack_qc.csv
05_proteinmpnn_sequences/FGA_rfpeptides_stage3D1_sidechain_repack_qc_pass.csv
05_proteinmpnn_sequences/FGA_rfpeptides_stage3D1_sidechain_repack_qc.md
05_proteinmpnn_sequences/FGA_rfpeptides_stage3D1_sidechain_repack_skipped_inputs.csv
05_proteinmpnn_sequences/RFpep_Site_2_stage3D1_sidechain_repack_review.pml
```

Current Stage 3D-1 result:

```text
unique ProteinMPNN-only sequences repacked: 7
duplicate ProteinMPNN-only inputs skipped: 1
pass_stage3d1_qc: 7
```

Interpretation:

- Repack-only removed the severe clash in all 7 unique ProteinMPNN-only
  sequences.
- RFpep_Site_2 contact, hotspot contact, and head-to-tail macrocycle geometry
  were preserved after repack.
- Stage 3D-1 is still not final validation; the 7 passing structures should
  next enter Stage 4-style scoring / validation.

Stage 3C sequence table columns:

```text
backbone_id
sequence_design_id
peptide_sequence
peptide_site_min_distance
peptide_hotspot_min_distance
macrocycle_terminal_cn_distance
pass_stage3c_qc
qc_failure_reasons
```

## 12. Stage 4: Rosetta Score-Only Interface Scoring

Goal:

Score Stage 3D-1 sequence-designed complexes and compute interface metrics
without ordinary unconstrained FastRelax.

Current route decision:

- Use Stage 3D-1 repack-only PDBs as Stage 4 inputs.
- Do not use ordinary unconstrained FastRelax at this point. The diagnostic
  control showed that the current unconstrained FastRelax mover can open the
  RFpep_Site_2 macrocycle.
- Stage 4 is a score/rank/filter layer, not final validation.

Pilot required metrics:

- `pyrosetta_score_status`
- scored PDB exists and is parseable
- Rosetta total score
- no-repack separation `ddg_proxy_no_repack`
- simple sequence properties and sequence-liability notes
- peptide-target contact count
- selected target-site contact count after scoring
- hotspot contact count after scoring
- peptide-target / site / hotspot minimum distance after scoring
- head-to-tail macrocycle C-N distance
- severe-clash status
- detached/collapsed complex flag
- `stage4_priority_rank`, `stage4_priority_class`, and top validation selection
- parser notes for chain, residue, and macrocycle handling

Article-style target metrics:

- Binding / interface energy, often reported as `ddG`. Current pilot uses a
  no-repack separated-state score proxy, not a mature final ddG.
- Interface contact / contact molecular surface, often reported as `CMS`.
- Aggregation or developability proxy, often reported as `SAP`.

Implementation rule:

- `SAP` and `CMS` are target metrics for the mature article-style route.
- If `SAP` or `CMS` are not implemented during the first pilot, record them as
  `not_available` instead of blocking the pilot.
- Do not promote pilot results to final candidates until the missing
  article-style target metrics are either implemented or explicitly waived in a
  written review.

Expected outputs:

```text
06_rosetta_scoring/scored_pdbs/
06_rosetta_scoring/FGA_rfpeptides_stage4_rosetta_interface_scores.csv
06_rosetta_scoring/FGA_rfpeptides_stage4_rosetta_interface_scores_pass.csv
06_rosetta_scoring/FGA_rfpeptides_stage4_top_validation_candidates.csv
06_rosetta_scoring/FGA_rfpeptides_stage4_rosetta_interface_scores.md
06_rosetta_scoring/RFpep_Site_2_stage4_rosetta_score_only_review.pml
```

Command:

```bash
cd /mnt/c/SH/fga_cyclic_peptide_design
source ~/fga_model_envs/miniforge3/etc/profile.d/conda.sh
conda activate proteinmpnn_binder_design

python scripts/25_stage4_rosetta_interface_scoring.py \
  --stage0-root results/rfpeptides_article_route_clean_20260615_fpocket \
  --stage3-root results/rfpeptides_article_route_clean_20260615_fpocket_stage1_N1000_no_traj \
  --selected-backbones RFpep_Site_2_0007
```

Current Stage 4A-v2 result:

```text
input_rows_scored: 7
pass_stage4_qc: 7
top_validation_candidates: 5
target_site_recovery_status: site_contact_pass for all 7
hotspot_recovery_status: hotspot_contact_pass for all 7
macrocycle_geometry_status: pass_head_to_tail_macrocycle for all 7
clash_status: pass_no_severe_clash for all 7
ddg_proxy_status: unfavorable_interface_proxy for all 7
sequence_liability_notes: none for all 7
```

Current Stage 4A-v2 top validation candidates:

1. `EYLLYGRPSSPL` (`priority_1_validation_ready`, `ddg_proxy_no_repack=0.554`)
2. `DYLLTGRPSSPL` (`priority_1_validation_ready`, `ddg_proxy_no_repack=0.663`)
3. `ELLLTGRPSSPL` (`priority_1_validation_ready`, `ddg_proxy_no_repack=0.854`)
4. `LEELYGLPSSPL` (`priority_2_backup`, `ddg_proxy_no_repack=1.678`)
5. `LYELYGLPSSPL` (`priority_2_backup`, `ddg_proxy_no_repack=1.848`)

Thresholds:

Do not freeze thresholds until the FGA pilot distribution is visible.
Record article-inspired thresholds and local observed distributions side by
side. A candidate should not pass only because one metric is strong.

Interpretation:

- Stage 4 confirms that all 7 Stage 3D-1 structures still satisfy the hard
  geometry/contact checks.
- The Rosetta no-repack separated-state score proxy is positive for all 7, so
  interface-energy evidence is weak or not clearly favorable. The top rows are
  selected for downstream structure prediction validation, not because they are
  confirmed good binders.
- `CMS` and `SAP` remain mature-route target metrics and are recorded as
  `not_available` in the current pilot.

2026-06-29 Stage 1 batch01 L12-24 scale-up update:

- Stage 1 batch:
  `results/rfpeptides_article_route_clean_20260623_stage1_N10000_L12_24_batch01`
- Stage 1 outputs: 10000 `.pdb` and 10000 `.trb`; trajectory output remained
  disabled.
- Stage 2 backbone QC: 10000 parsed, 6 passed. The `hotspot_near_pass`
  backbone `RFpep_Site_2_2330` was not advanced.
- Advanced to ProteinMPNN-only: `RFpep_Site_2_4373`, `RFpep_Site_2_9132`,
  `RFpep_Site_2_2168`, `RFpep_Site_2_3254`, `RFpep_Site_2_2193`.
- Diagnostic ProteinMPNN-FastRelax on these 5 backbones generated 40 structures,
  but Stage 3C passed 0/40; all 40 failed macrocycle geometry after FastRelax.
  This reinforces the decision not to use ordinary unconstrained FastRelax as
  the mainline.
- ProteinMPNN-only generated 40 structures; Stage 3C passed 22/40 directly.
- Stage 3D-1 side-chain repack-only repacked 34 unique/contact-preserving rows
  and passed 34/34.
- Stage 4A-v2 score-only evaluated 34 rows and passed 34/34 hard QC.
- Stage 4A-v2 top validation candidates from batch01:
  1. `GLRNREEVPEELRKIVK` (`ddg_proxy_no_repack=-26.892`)
  2. `GLRNKEEIPPELREIVK` (`ddg_proxy_no_repack=-26.143`)
  3. `GLKNLEEVPPELREIVK` (`ddg_proxy_no_repack=-24.760`)
  4. `GLKNKEEVPEELLKIVK` (`ddg_proxy_no_repack=-24.760`)
  5. `QNILKGEIVLGVSDEEMA` (`ddg_proxy_no_repack=-24.162`)
- These are validation candidates only, not final peptide candidates.

## 13. Stage 5: AfCycDesign-First Validation

Goal:

Test whether the sequence-designed macrocycle-target complex can be recovered
by an independent structure prediction / design-validation model.

Primary validation route:

AfCycDesign is the preferred validation route for this branch because the
branch is explicitly based on the RFpeptides macrocycle design paper. The first
implementation attempt should therefore be AfCycDesign-oriented, not a vague
mix of validation engines.

Fallback rule:

Use RF2 or AF2 initial-guess style validation only if AfCycDesign cannot be
installed, cannot run on the local system, or cannot support the required FGA
macrocycle inputs after a documented attempt. When fallback is used, the report
must label the result as:

```text
validation_route=rf2_or_af2ig_fallback
afcycdesign_status=unavailable_or_failed_with_reason
```

Main validation question:

```text
Does the validation model recover the RFpeptides design pose at the intended
FGA target site?
```

This is different from asking:

```text
Can some model predict some complex for this sequence?
```

Validation inputs:

- Target chain / target crop.
- ProteinMPNN-designed macrocycle sequence.
- Initial design pose if using an initial-guess validation mode.

Validation outputs:

```text
07_structure_validation/
  FGA_rfpeptides_stage5_candidate_manifest.csv
  FGA_rfpeptides_stage5_candidate_manifest.md
  FGA_rfpeptides_stage5_protocol_audit.md
  inputs/
  jobs/
  predictions/  # created only after approved prediction starts
07b_rf2_af2ig_fallback_validation/
08_validation_parse/FGA_rfpeptides_validation_scores.csv
```

Suggested validation columns:

```text
design_id
backbone_id
mpnn_sequence_id
validation_engine
validation_route
afcycdesign_status
validation_pdb
validation_success
interface_pae_or_iPAE
plddt_peptide
plddt_interface
design_to_validation_target_aligned_rmsd
peptide_backbone_rmsd
target_site_contact_recovery
hotspot_contact_count
off_site_contact_count
same_target_site_flag
validation_notes
```

Critical interpretation rule:

- A validation structure that binds elsewhere on FGA is not an RFpeptides
  selected target-site success.
- A validation structure with a strong score but poor agreement to the
  RFpeptides design pose should be treated as a separate hypothesis, not as
  confirmation of the original RFpeptides design.

### Stage 5 preparation status (2026-07-10)

Candidate merge, protocol audit, inputs, and review-gated jobs are prepared.
No Stage 5 structure prediction has been run.

Combined output root:

```text
results/rfpeptides_article_route_clean_20260623_stage5_batch01_batch02/
  07_structure_validation/
```

Prepared candidates:

1. batch01 / `RFpep_Site_2_3254` / `GLRNREEVPEELRKIVK`
2. batch01 / `RFpep_Site_2_2193` / `QNILKGEIVLGVSDEEMA`
3. batch01 / `RFpep_Site_2_9132` / `DLEKGTKWSPDASLE`
4. batch02 / `RFpep_Site_2_1518` / `KGDAENPEVDALLG`
5. batch02 / `RFpep_Site_2_5529` / `LLLGLGTGTLEE`

This first round covers five distinct backbone families. No backbone receives
a second sequence in the first round; same-backbone sequence-sensitivity jobs
can be added later if recovery differs among backbone families.

`RFpep_Site_2_5529` is a valid Stage 4 hard-QC backbone family: six of its
sequence structures passed Stage 4. Its selected row has four target-site and
three hotspot contacts, but a weak `ddg_proxy_no_repack=-2.918`. It is included
for backbone/contact diversity and as a validation control, not as an
energy-ranked lead.

The manifest reports raw `ddg_proxy_no_repack` together with:

```text
ddg_proxy_per_peptide_residue
ddg_proxy_per_target_contact
ddg_proxy_per_site_contact
```

These normalizations help compare different peptide lengths/contact counts,
but none is an experimental binding energy.

Legacy Stage 5A prediction protocol used for the completed 125-model run:

```text
engine: ColabDesign gamma AfCycDesign prediction
gamma_commit: 5ab4efaba2321a6c3c314b82d2fff8e0241f5c2d
validation_test_type: independent_recovery
input: target sequence + cyclic peptide sequence
cyclic_chain: peptide only
cyclic_encoding: relative-position cyclic offset
template_mode: none
use_initial_guess: false
target_msa_mode: single_sequence
peptide_msa_mode: single_sequence
msa_rows_input: 1
use_mlm: true
mlm_replace_fraction: 0.15
seeds_per_candidate: 5
models_per_seed: 5
recycles: 6
seed_jobs: 25
model_evaluations: 125
```

AfCycDesign supports head-to-tail cyclic geometry through the cyclic relative
position encoding. This is not an explicit terminal C-N chemical bond record,
so terminal C-N distance and cyclic topology must still be checked in every
prediction.

The existing local `ColabDesign-cyclic-binder` environment supports cyclic
binder design but lacks the gamma `contrib.predict` and `contrib.cyclic`
prediction modules. The generated setup script prepares a separate pinned
gamma source checkout and leaves the existing environment unchanged.

Current environment checkpoint (2026-07-14):

```text
pinned_gamma_source: prepared
gamma_commit: 5ab4efaba2321a6c3c314b82d2fff8e0241f5c2d
contrib.predict: present
contrib.cyclic: present
stage5_python_overlay: ready_IPython_installed
stage5_protocol_preflight: passed
current_preflight_blocker: none
stage5_predictions_started: true
smoke_job_completed: S5_04_batch02_RFpep_Site_2_1518_seed00
exploratory_prediction_jobs_completed: 25/25
prediction_models_completed: 125/125
stage5_validation_claim_status: no_candidate_passed_independent_recovery
```

The pinned gamma modules, cyclic offset helper, Stage 5 Python overlay, and AF
parameter directory have passed the generated protocol preflight. The next
action is one full-parameter smoke prediction only:

```text
candidate: batch02 / RFpep_Site_2_1518
seed: 00
models: all 5
recycles: 6
template_mode: none
use_initial_guess: false
validation_test_type: independent_recovery
```

Do not start the remaining 24 prediction jobs until this smoke output has been
parsed and checked for chain identity, confidence outputs, target-site recovery,
and terminal C-N cyclic-topology evidence.

Preparation command:

```bash
cd /mnt/c/SH/fga_cyclic_peptide_design
source ~/fga_model_envs/miniforge3/etc/profile.d/conda.sh
conda activate fga_stage1_fpocket

python scripts/26_prepare_afcycdesign_jobs.py
```

Generated scripts are review-gated. The master runner exits unless
`RUN_STAGE5_PREDICTIONS=YES` is explicitly set after protocol review.

Run the single approved smoke job only:

```bash
cd /mnt/c/SH/fga_cyclic_peptide_design

RUN_STAGE5_PREDICTIONS=YES bash \
  results/rfpeptides_article_route_clean_20260623_stage5_batch01_batch02/07_structure_validation/jobs/run_S5_04_batch02_RFpep_Site_2_1518_seed00.sh
```

This one-command environment assignment does not unlock later jobs. Do not run
`run_stage5_afcycdesign_independent_recovery_all.sh` at this checkpoint.

Smoke result checkpoint (2026-07-13):

```text
technical_execution: pass
models_generated: 5/5
models_with_head_to_tail_CN_geometry_pass: 5/5
models_with_Site_2_contact: 0/5
models_with_hotspot_contact: 0/5
models_with_geometry_recovery: 0/5
models_with_strict_recovery: 0/5
target_aligned_peptide_backbone_RMSD_A: 26.896-43.649
peptide_pLDDT_0_100: 46.44-56.63
interface_PAE_A: 26.219-27.665
iPTM: 0.03701-0.05684
remaining_24_jobs_release_status: hold
```

Interpretation: the gamma runner, all five AF parameter sets, output writing,
and relative-position cyclic offset are technically compatible. However, this
seed did not independently recover the designed target site or peptide pose,
and its interface confidence is weak. This is a failed independent-recovery
result, not evidence that the peptide is a validated binder.

Parse completed or partial Stage 5 outputs with:

```bash
cd /mnt/c/SH/fga_cyclic_peptide_design

~/fga_model_envs/colabdesign-py310/.pixi/envs/default/bin/python \
  scripts/27_collect_afcycdesign_validation.py \
  --stage5-root results/rfpeptides_article_route_clean_20260623_stage5_batch01_batch02/07_structure_validation \
  --stage0-root results/rfpeptides_article_route_clean_20260615_fpocket
```

The parser reports model-level target-aligned peptide backbone RMSD, Site_2 and
hotspot recovery, off-site contacts, same-site status, terminal C-N geometry,
severe clashes, peptide pLDDT, interface PAE, and seed-level recovery counts.

Full five-candidate exploratory result checkpoint (2026-07-14):

```text
candidates: 5
seeds: 25/25
models: 125/125
target_sequence_match: 125/125
peptide_sequence_match: 125/125
head_to_tail_CN_geometry_pass: 125/125
models_with_any_Site_2_contact: 31/125
models_with_hotspot_contact: 0/125
same_target_site_flag: 0/125
models_with_severe_interchain_clash: 16/125
geometry_recovery_pass: 0/125
confidence_pass: 0/125
strict_independent_recovery_pass: 0/125
target_backbone_RMSD_A: 7.026-36.949
target_aligned_peptide_backbone_RMSD_A: 8.599-59.289
minimum_hotspot_distance_A: 10.834
peptide_pLDDT_0_100_overall_range: 38.0-69.4
interface_PAE_A_overall_range: 23.0-28.2
maximum_iPTM: 0.127
```

The 31 broad Site_2 contacts involved A9, A13, or A79. None involved the
selected hotspots A82, A84, A85, or A86. The result therefore does not satisfy
the Stage 5 requirement that recovery occur at RFpep_Site_2 near its hotspots.

Interpretation: all jobs and protocol metadata are complete and internally
consistent, so this is not a missing-output or sequence-mix-up problem. The
sequence-only, no-template protocol preserved terminal cyclic geometry but did
not reliably recover the isolated target crop conformation or the designed
peptide pose. These results are protocol-exploration evidence only. They do not
prove that the five peptides cannot bind experimentally, and none may advance
to final ranking as an independently validated candidate.

### Stage 5A-v2 identity and metric corrections (2026-07-14)

The completed 125-model output above is retained as:

```text
Stage 5A legacy single-sequence + MLM 0.15 independent-recovery run
```

It remains valid as protocol-exploration evidence, but its output directory is
not reused for Stage 5A-v2. The old candidate/job IDs did not bind the peptide
sequence or complete protocol identity strongly enough to the cache key.

Updated scripts:

```text
scripts/26_prepare_afcycdesign_jobs.py
scripts/external/run_afcycdesign_independent_recovery.py
scripts/27_collect_afcycdesign_validation.py
```

Stage 5A-v2 now applies the following safeguards:

- Candidate IDs contain the first 8 characters of peptide-sequence SHA1.
- Job IDs additionally contain an 8-character protocol hash.
- A completed output is reused only if metadata matches target sequence,
  peptide sequence, both sequence hashes, seed, requested recycle count,
  forward-pass count, fixed ColabDesign commit, cyclic chain index, template
  mode, initial-guess mode, target/peptide MSA modes, and MLM settings.
- The runner verifies that the imported ColabDesign module is under the pinned
  source directory and that its commit marker equals the job specification.
- Stage 4 CSV peptide sequence must equal design-PDB chain B; Stage 0 target
  sequence must equal every design-PDB chain A and the job target sequence.
- The target must be 86 aa, all target/peptide sequences must use legal amino
  acids, and recorded peptide lengths must agree with sequence and PDB.

New Stage 5A-v2 output root:

```text
results/rfpeptides_article_route_clean_20260623_stage5A_v2_batch01_batch02/
  07_structure_validation/
```

The 5 candidates and 25 seed jobs were regenerated there, but no Stage 5A-v2
complex prediction has been started.

Metric corrections:

```text
plddt_mean_fraction
plddt_peptide_mean_fraction
plddt_mean_100
plddt_peptide_mean_100
interface_pae_global_mean_A
interface_pae_site2_mean_A
interface_pae_site2_median_A
interface_pae_hotspot_mean_A
```

The pinned gamma implementation stores `model.aux["plddt"]` on a 0-1 scale;
reports use the explicit 0-100 fields. PAE remains in angstrom. The old field
named `interface_pae_normalized` was actually the bidirectional mean over all
target-peptide residue pairs and was not normalized; Stage 5A-v2 no longer uses
that name.

Recycle accounting is explicit:

```text
requested_recycles: 6
forward_passes: 7
```

### Stage 5 target-only controls (2026-07-14)

Before any Stage 5B target-structure-conditioned test, isolate whether the
86-aa target crop itself can be recovered from sequence. Four target-only
groups are prepared:

```text
A: single sequence + MLM 0.15
B: single sequence + no MLM
C: real unpaired target homolog MSA + MLM 0.15
D: real unpaired target homolog MSA + no MLM
```

All groups use three seeds and all five model parameter sets, for 12 planned
seed jobs and 60 planned model evaluations. Every group remains:

```text
template_mode: none
use_initial_guess: false
peptide_included: false
validation_test_type: sequence_based_target_recovery_control
```

Prepared scripts:

```text
scripts/28_prepare_stage5_target_controls.py
scripts/external/run_afcycdesign_target_recovery_control.py
scripts/29_collect_stage5_target_controls.py
```

Control output root:

```text
results/rfpeptides_article_route_clean_20260623_stage5_target_controls_v1/
  07_structure_validation_target_controls/
```

Current preparation result:

```text
groups_planned: 4
seed_jobs_planned: 12
A/B_seed_jobs_completed: 6/6
A/B_model_evaluations_completed: 30/30
C/D_seed_jobs_blocked_missing_real_homolog_MSA: 6
target_control_models_passing_recovery: 0/30
Stage5B_started: false
```

No real FGA homolog A3M was found locally. C/D must not be approximated with a
fake paired MSA. The user should obtain a real unpaired homolog A3M whose query
matches `data/input/FGA_full_length_1_866.fasta`; the preparation script then
verifies the query and crops UniProt columns 134-219 to the current 86-aa target.
The peptide remains single-sequence in Stage 5A and is absent from these
target-only controls.

Target-control parsing compares:

```text
target_global_CA_RMSD_A
Site_2_local_CA_RMSD_A
hotspot_local_CA_RMSD_A
target_mean_pLDDT_100
Site_2_mean_pLDDT_100
hotspot_mean_pLDDT_100
```

Site_2 and hotspot RMSD are measured after one global target CA alignment.
This control must be reviewed before deciding whether to revisit sequence-only
Stage 5A, use a larger target context, or define a Stage 5B pose-stability test.

Available target-only result checkpoint:

```text
A_single_sequence_MLM_0.15:
  seeds_completed: 3/3
  models_completed: 15/15
  recovery_pass: 0/15
  median_target_global_CA_RMSD_A: 36.525
  median_Site_2_local_CA_RMSD_A: 65.135
  median_hotspot_local_CA_RMSD_A: 81.704
  median_target_pLDDT_100: 87.866
  median_Site_2_pLDDT_100: 67.684

B_single_sequence_no_MLM:
  seeds_completed: 3/3
  models_completed: 15/15
  recovery_pass: 0/15
  median_target_global_CA_RMSD_A: 10.984
  median_Site_2_local_CA_RMSD_A: 18.817
  median_hotspot_local_CA_RMSD_A: 26.015
  median_target_pLDDT_100: 78.138
  median_Site_2_pLDDT_100: 64.885
```

Turning MLM off improves the target-shape metrics substantially in this small
control, but does not recover the reference target or Site_2. High global
pLDDT in group A is therefore not evidence of correct target conformation.

Under deterministic `dropout=false` inference, group B has no remaining random
MSA-mask source, and its three nominal seeds produce identical model outputs.
The job metadata still records all seeds, but they must not be described as
three independent structural samples. Group A varies slightly because MLM
masking supplies a stochastic component. C/D remain required before deciding
whether real homolog information can recover the 86-aa crop.

### Stage 5B target-structure-conditioned recovery (prepared 2026-07-14)

The user elected to prepare the complete Stage 5B matrix without a separate
prediction pilot because the expected runtime is manageable. A strict static
preflight remains mandatory and is executed again by the master runner before
the first model prediction.

Stage 5B is defined as:

```text
validation_test_type: target_structure_conditioned_recovery
interpretation: target-conditioned peptide pose-recovery test
independent_recovery: false
experimental_binding_validation: false
final_candidate_decision: false
```

The model receives the Stage 0 `RFpep_Site_2` target-only structure as a
template for chain A. It receives the candidate peptide sequence and peptide
cyclic relative-position offset, but no peptide template, peptide coordinates,
initial guess, or complete Stage 4 complex. The Stage 4 PDB is staged separately
and is read only after prediction for target alignment and pose/contact recovery
analysis.

Prepared scripts:

```text
scripts/30_prepare_stage5b_target_conditioned_jobs.py
scripts/external/run_afcycdesign_target_conditioned_recovery.py
scripts/31_collect_stage5b_validation.py
```

Independent output root:

```text
results/rfpeptides_article_route_clean_20260623_stage5B_batch01_batch02/
  07_structure_validation_target_conditioned/
```

Prepared candidate set:

1. batch01 / `RFpep_Site_2_3254` / `GLRNREEVPEELRKIVK`
2. batch01 / `RFpep_Site_2_2193` / `QNILKGEIVLGVSDEEMA`
3. batch01 / `RFpep_Site_2_9132` / `DLEKGTKWSPDASLE`
4. batch02 / `RFpep_Site_2_1518` / `KGDAENPEVDALLG`
5. batch02 / `RFpep_Site_2_5529` / `LLLGLGTGTLEE`

Run matrix and fixed protocol identity:

```text
candidates: 5
seeds_per_candidate: 5
model_parameter_sets_per_seed: 5
seed_jobs: 25
model_predictions: 125
target_template_coverage_required: 86/86
peptide_template_coverage_required: 0
template_mode: target_only
use_initial_guess: false
target_msa_mode: single_sequence
peptide_msa_mode: single_sequence
use_mlm: false
use_dropout: true
requested_recycles: 6
forward_passes: 7
cyclic_chain_index: 1
cyclic_offset_applied_to: peptide_only
gamma_commit: 5ab4efaba2321a6c3c314b82d2fff8e0241f5c2d
protocol_hash: 0e506709
```

Candidate IDs include peptide-sequence SHA1. Job IDs additionally include the
target-template SHA1, protocol hash, and seed. A completed cache entry is reused
only when its metadata matches target and peptide sequences/hashes, seed,
recycles, fixed commit, cyclic chain, template mode, initial-guess mode, target
template hash and coverage, MSA modes, MLM/dropout settings, and all five
expected PDB/NPZ model outputs. A mismatch forces recomputation.

Static preflight result:

```text
candidate_specs_checked: 5/5
target_template_chain: A only
target_template_coverage: 86/86 for all candidates
peptide_template_coverage: 0 for all candidates
use_initial_guess: false
complete_Stage4_complex_loaded_as_model_input: false
pinned_gamma_commit_check: pass
shell_syntax_check: pass
model_forward_pass_started_by_preflight: false
full_predictions_started: false
```

The target-only homolog-MSA controls C/D are still blocked because no real FGA
homolog A3M is available locally. Proceeding with Stage 5B does not resolve that
missing sequence-only control and must not be reported as if it did.

Regenerate the manifest and jobs only when the candidate set or protocol changes:

```bash
cd /mnt/c/SH/fga_cyclic_peptide_design
source ~/fga_model_envs/miniforge3/etc/profile.d/conda.sh
conda activate fga_stage1_fpocket

python scripts/30_prepare_stage5b_target_conditioned_jobs.py
```

Run the non-predictive preflight:

```bash
cd /mnt/c/SH/fga_cyclic_peptide_design

bash \
  results/rfpeptides_article_route_clean_20260623_stage5B_batch01_batch02/07_structure_validation_target_conditioned/jobs/check_stage5B_target_conditioned_protocol.sh
```

After reviewing the manifest and preflight, the user may start the complete
25-job/125-model run manually:

```bash
cd /mnt/c/SH/fga_cyclic_peptide_design

RUN_STAGE5B_PREDICTIONS=YES bash \
  results/rfpeptides_article_route_clean_20260623_stage5B_batch01_batch02/07_structure_validation_target_conditioned/jobs/run_stage5B_target_conditioned_recovery_all.sh
```

Collect complete or partial outputs with the AfCycDesign Python environment:

```bash
cd /mnt/c/SH/fga_cyclic_peptide_design

~/fga_model_envs/colabdesign-py310/.pixi/envs/default/bin/python \
  scripts/31_collect_stage5b_validation.py \
  --stage5b-root results/rfpeptides_article_route_clean_20260623_stage5B_batch01_batch02/07_structure_validation_target_conditioned \
  --stage0-root results/rfpeptides_article_route_clean_20260615_fpocket
```

The collector writes model-level metrics, candidate summaries, seed-diversity
hashes, and a Markdown report. Candidate classes are limited to:

```text
stage5B_strong_support
stage5B_partial_support
stage5B_not_recovered
stage5B_protocol_failure
```

These classes are internal computational evidence tiers. Even
`stage5B_strong_support` does not create a final peptide candidate or establish
experimental binding.

### Stage 5B execution and upstream hotspot-mapping audit (2026-07-15)

The complete Stage 5B matrix was executed and parsed:

```text
seed_jobs_completed: 25/25
model_predictions_parsed: 125/125
protocol_identity_valid: 125/125
target_template_input_verified: 125/125
peptide_template_coverage_zero: 125/125
macrocycle_terminal_CN_geometry_pass: 125/125
effective_unique_predictions: 125/125
seed_diversity_status: seed_diverse for all 5 candidates
confidence_supportive: 0/125
severe_interchain_clash: 26/125
target_template_recovery_pass: 0/125
reference_designs_valid_for_corrected_Site_2: 0/5
stage5B_strong_or_partial_support: 0/5
```

The gamma model log did not populate a `multi` ranking field. Because pTM and
ipTM were complete, `ranking_confidence` was recovered without rerunning models
as `0.8 * ipTM + 0.2 * pTM` and explicitly tagged
`derived_0.8_iPTM_plus_0.2_pTM` for all 125 rows (range 0.101-0.213).

The target/template RMSD distributions were:

```text
target_global_CA_RMSD_A: 3.533 / 5.812 / 36.504
Site_2_local_CA_RMSD_A: 5.729 / 9.976 / 65.796
hotspot_local_CA_RMSD_A: 8.136 / 14.755 / 83.615
```

The first collector pass exposed a post-processing frame bug: Stage 4 reference
complexes and the Stage 0 target template are in different global coordinate
frames. Predicted peptides had been aligned to the Stage 0 target but compared
against untransformed Stage 4 peptide coordinates, creating artificial peptide
RMSD values of approximately 52-96 A. `scripts/31_collect_stage5b_validation.py`
now aligns both predicted and reference targets to the Stage 0 frame before
comparing peptide backbones. No AfCycDesign prediction needed to be rerun.

The corrected peptide-pose analysis then revealed a more important upstream
problem. Stage 0 hotspot values `82,84,85,86` are crop sequence positions. In an
RFpeptides output with a 17-residue binder, the target chain is numbered
`18-103`; the intended hotspot residues therefore appear as output PDB residues
`99,101,102,103`.

Two old assumptions were incorrect:

1. RFpeptides `get_idx0_hotspots()` used receptor-local indices without adding
   the binder length. Input hotspots were therefore applied to complex-global
   positions `82,84,85,86`, not `99,101,102,103` for the 17-residue example.
2. Stage 2 and downstream QC matched Stage 0 crop positions directly against
   raw RFpeptides/Rosetta PDB residue numbers. This reproduced the same shifted
   hotspot and generated false-positive Site_2/hotspot pass labels.

#### 故障起点与传播链（必须保留的路线审计记录）

这条路线**最早从 Stage 1 RFpeptides hotspot conditioning 开始失效**，不是
从 ProteinMPNN、Rosetta 或 Stage 5 AfCycDesign 才开始出错。各阶段边界如下：

1. **Stage -1 site rediscovery：不是本次错误起点。**
   `RFpep_Site_2` 的表面暴露、RSA、几何、遮挡、chemical anchor 和 fpocket
   旁证仍是对 native target context 的位点筛选；本次没有发现这些指标使用了
   RFpeptides 输出复合物的偏移编号。
2. **Stage 0 target preparation：位点定义本身有效。**
   crop 内 hotspot `82,84,85,86` 表示 86-aa target crop 的序列位置。Stage 0
   target-only PDB、序列和 crop mapping 可以继续作为修复后路线的输入。
3. **Stage 1 RFpeptides generation：第一个实质错误。**
   RFpeptides 将长度为 `L` 的 peptide 放在 target 前面组装复合物，但旧
   `get_idx0_hotspots()` 把 target-local hotspot index 直接当作 complex-global
   index，没有加 `L`。对 crop position `p`，正确的 complex 0-based index
   应为 `L + (p - 1)`，对应输出 PDB residue number `L + p`。例如 `L=17`
   时，正确热点 PDB 编号是 `99,101,102,103`，旧流程却把条件施加在
   `82,84,85,86`。因此模型从生成开始就被引导到错误的 target 区域。
4. **Stage 2 backbone QC：第二个错误，并把 Stage 1 问题掩盖成假阳性。**
   旧 QC 又把 crop position `82,84,85,86` 直接与输出 PDB residue number
   比较，没有按 target-chain residue order 映射。生成和 QC 恰好共享同一种
   偏移错误，所以若 peptide 接近错误区域，仍会被标记为 Site_2/hotspot pass。
5. **Stage 3：没有制造最初的位点偏移，但继承了错误 backbone。**
   ProteinMPNN-only 只为既有 backbone 生成序列；3D-1 只做 side-chain repack，
   不移动 backbone，因此都不可能把 off-target backbone 自动移回正确 Site_2。
   旧 Stage 3C/3D-1 contact QC 同样继承了错误 residue mapping。
6. **Stage 4：能量计算本身可以运行，但排序对象是错误位点上的复合物。**
   `ddg_proxy_no_repack` 等数值只能描述这些旧 pose 的 score-only 界面，不能
   作为正确 Site_2 binder 的排序证据。
7. **Stage 5A/5B：协议没有发现链顺序、cyclic offset、peptide template 或
   initial guess 误用，但输入候选前提已经失效。** Stage 5B 的 125 个模型可
   保留为 protocol engineering record，不能解读为对 5 个正确 Site_2 候选的
   independent/conditioned recovery 结论。

Stage 5B collector 还单独发现并修复了一个**后处理坐标系错误**：预测 peptide
已经对齐到 Stage 0 target frame，但旧代码曾将其与未变换的 Stage 4 peptide
坐标比较，导致约 `52-96 A` 的人工 RMSD。这个错误会污染 RMSD 报告，但它
不是最初的 off-target 根因；修复坐标系后，Stage 1 hotspot 偏移仍然成立。

受影响范围必须按下面方式隔离：

- 所有使用修复前 `rfd_macro` 和该 Site_2 hotspot conditioning 生成的旧
  RFpeptides backbones 都视为 **pre-hotspot-fix / suspect off-target**。
- 旧 batch01/batch02 的 Stage 2-5 pass、ProteinMPNN 序列和 Stage 4 排名不得
  作为正确 Site_2 候选证据；其他修复前批次在进入下游前也必须重新审计。
- Stage -1 位点发现和 Stage 0 target crop 不因本次问题作废，可由修复后的
  Stage 1 重新开始。
- 旧文件保留用于追溯，不覆盖、不删除，也不与修复后生产结果混合汇总。

Independent checks of the five Stage 1 backbones selected for Stage 5 showed:

```text
RFpep_Site_2_3254: corrected hotspot contacts=0, hotspot min distance=23.308 A
RFpep_Site_2_2193: corrected hotspot contacts=0, hotspot min distance=19.954 A
RFpep_Site_2_9132: corrected hotspot contacts=0, hotspot min distance=16.985 A
RFpep_Site_2_1518: corrected hotspot contacts=0, hotspot min distance=24.560 A
RFpep_Site_2_5529: corrected hotspot contacts=0, hotspot min distance=21.360 A
```

Their staged Stage 4 references likewise have zero corrected Site_2/hotspot
contacts. The Stage 5B collector now checks this reference premise explicitly;
all five rows are classified as `stage5B_protocol_failure`. The 125 predictions
remain useful as an engineering record, but they are not a valid negative test
of five correctly targeted Site_2 candidates.

Corrected files:

```text
rfd_macro/rfdiffusion/inference/utils.py
rfd_macro/rfdiffusion/inference/model_runners.py
scripts/20_make_rfpeptides_article_jobs.py
scripts/21_collect_rfpeptides_backbones.py
scripts/23_collect_proteinmpnn_sequences.py
scripts/31_collect_stage5b_validation.py
```

The RFpeptides runner now aborts unless a static preflight proves that hotspot
indices include the binder-length offset. Stage 2 maps crop positions through
target-chain residue order and records the actual PDB residue numbers per row.
Stage 3C/3D-1/4 inherit the corrected mapping, including object-specific PyMOL
selections.

The old batch01/batch02 Stage 2-5 results must be retained as legacy/off-target
protocol-exploration outputs and must not be used for final ranking. They may be
re-screened for accidental correct-site contacts, but they were not generated
with correct intended-hotspot conditioning.

A new 20-design, no-trajectory smoke job is prepared but not run:

```text
results/rfpeptides_article_route_clean_20260715_hotspotfix_smoke/
  01_rfpeptides_jobs/run_rfpeptides_stage1_site2_hotspotfix_smoke_N20.sh
```

Run it manually:

```bash
cd /mnt/c/SH/fga_cyclic_peptide_design

bash \
  results/rfpeptides_article_route_clean_20260715_hotspotfix_smoke/01_rfpeptides_jobs/run_rfpeptides_stage1_site2_hotspotfix_smoke_N20.sh
```

Then run corrected Stage 2 QC into a versioned subdirectory, preserving old
results:

```bash
cd /mnt/c/SH/fga_cyclic_peptide_design
source ~/fga_model_envs/miniforge3/etc/profile.d/conda.sh
conda activate fga_stage1_fpocket

python scripts/21_collect_rfpeptides_backbones.py \
  --stage0-root results/rfpeptides_article_route_clean_20260615_fpocket \
  --stage1-root results/rfpeptides_article_route_clean_20260715_hotspotfix_smoke \
  --output-subdir 03_backbone_qc_hotspot_mapping_v2 \
  --selected-sites RFpep_Site_2
```

Do not restart ProteinMPNN or Stage 5 until the smoke TRB confirms that input
crop hotspots map to the correct complex-global target positions and corrected
Stage 2 finds genuine Site_2/hotspot contacts.

### Pre-hotspot-fix output quarantine (2026-07-16)

All identified Stage 1-5 outputs that used or inherited the incorrect Site_2
hotspot mapping were moved into a dedicated quarantine directory:

```text
C:/SH/fga_cyclic_peptide_design/results/
  _archived_invalid_site2_hotspot_mapping_20260716/
    README_错误路线归档.md
    DO_NOT_USE_FOR_SITE2_DESIGN.txt
    archive_manifest.csv
    archive_manifest_pre_move.csv
    affected_results/
```

The move was performed within the same filesystem. No affected output was
deleted and no second multi-GiB copy was created. Every entry was checked
against its pre-move file count and byte count:

```text
archived_entries: 13
archived_files: 99,147
archived_bytes: 3,812,075,662
archived_size_GiB: 3.55
inventory_status: moved_verified for 13/13 entries
```

The quarantine contains:

```text
20260612/01_rfpeptides_jobs
20260615_fpocket/01_rfpeptides_jobs
20260615_fpocket/02_rfpeptides_backbones
20260615_fpocket/logs
20260615_fpocket_stage1_N1000_no_traj
20260623_stage1_N10000_L12_24_batch01 through batch05
20260623_stage5_batch01_batch02
20260623_stage5A_v2_batch01_batch02
20260623_stage5B_batch01_batch02
```

The following remain outside quarantine because they do not depend on the old
peptide hotspot pose premise or are part of the corrected route:

```text
20260612/00_site_discovery and 00_target_inputs
20260615_fpocket/00_site_discovery and 00_target_inputs
20260623_stage5_target_controls_v1
20260715_hotspotfix_smoke
```

The corrected N20 smoke produced 20 backbones. Its first versioned Stage 2 QC
table contained six pass rows, but three had zero direct hotspot contacts and
passed only through the old `hotspot_near_pass` allowance. Stage 2 now rejects
near-only site/hotspot states by default; `--allow-near-pass` exists only as an
explicit compatibility option. The strict rerun was written to:

```text
results/rfpeptides_article_route_clean_20260715_hotspotfix_smoke/
  03_backbone_qc_hotspot_mapping_v3_strict/
```

Strict result:

```text
parsed_backbones: 20
strict_pass: 3
strict_pass_ids: RFpep_Site_2_0002, RFpep_Site_2_0017, RFpep_Site_2_0018
direct_hotspot_contacts: 2, 1, 1 respectively
severe_clash: none for all 3
macrocycle_geometry: pass for all 3
```

These three remain smoke candidates pending manual structure review; they are
not production candidates and do not yet justify a large rerun.

The active restart boundary remains corrected Stage 1 followed by strict Stage
2. No quarantined Stage 3, Stage 4, Stage 5A, or Stage 5B output may be restored
to candidate ranking.

## 14. Stage 6: Final Ranking

Goal:

Rank only designs that pass all required stages.

Required pass gates:

```text
RFpeptides generation pass
Backbone QC pass
ProteinMPNN sequence exists
Rosetta relax / scoring pass
AfCycDesign validation pass, or documented RF2/AF2 initial-guess fallback pass
Target-site contact recovery pass
macrocycle_geometry_status is acceptable
No severe off-site binding
```

Suggested output:

```text
09_ranked_candidates/FGA_rfpeptides_article_route_ranked.csv
reports/FGA_RFpeptides_Article_Route_Report_YYYYMMDD.md
```

Suggested final table columns:

```text
rank
design_id
target_site_id
length
peptide_sequence
rfpeptides_pdb
relaxed_pdb
validation_pdb
ddG
CMS
SAP
interface_pae_or_iPAE
design_to_validation_rmsd
target_site_contact_recovery
macrocycle_geometry_status
same_target_site_flag
final_status
final_notes
```

## 15. Proposed Pilot Before Full Campaign

This is a discussion proposal, not a command to run yet.

Pilot A:

```text
Target-site scope: RFpep_Site_2 only for the first Stage 1 pilot
Length scope: 12-18 mixed range
RFpeptides backbones: 10-20 for RFpep_Site_2
diffuser.T: 50
hotspots: 3-6 per selected target site
ProteinMPNN sequences: 4-8 per backbone
Rosetta relax: yes if environment is ready
Validation: AfCycDesign first; RF2/AF2 initial-guess style only as documented fallback
```

Success criteria for deciding whether to scale:

- RFpeptides outputs parse reliably.
- A meaningful fraction of generated macrocycles contact the intended target
  site.
- ProteinMPNN sequences are diverse and valid.
- Relaxed complexes do not collapse or detach.
- Validation recovers at least some RFpeptides design poses.
- RFpep_Site_2 behavior is interpretable enough to justify scaling or adding
  deferred targets.

## 16. Full Campaign Concept

Only after pilot review:

```text
RFpep_Site_2: hundreds to thousands of RFpeptides backbones
RFpep_Site_3 / RFpep_Site_4: deferred; revisit only after Site_2 pilot review
Length: 12-18, optionally split into fixed length bins later
ProteinMPNN: multiple sequences per backbone
Rosetta pilot-required metrics: all sequence-designed candidates
Rosetta article-style target metrics: implement and apply as the route matures
Validation: AfCycDesign-first on the chosen subset or all candidates, depending on compute
```

The full campaign count should be decided from pilot pass rates and available
GPU/CPU time.

## 17. Open Decisions Before Any New Scripts

These must be resolved before implementation:

1. Exact RFpeptides-specific target-site definitions from clean site discovery.
2. Exact target crop definitions for each selected target site.
3. Exact hotspot residue lists after target crop renumbering.
4. Whether to use one mixed `12-18` contig or separate fixed length bins.
5. Whether Rosetta is installed and usable for ProteinMPNN-FastRelax.
6. Whether AfCycDesign can be installed and used as the primary validation
   route, or whether a documented RF2/AF2 initial-guess fallback is necessary.
7. How to display head-to-tail macrocycle geometry in reports without using any
   Cys-Cys / SG-SG terminology from the old ColabDesign branch.
8. Whether to add RFpep_Site_3 / RFpep_Site_4 or other exploratory extra sites
   after the RFpep_Site_2 pilot, based on Stage -1 and Stage 1 evidence.

## 18. What This Branch Must Not Do

Do not:

- Mix outputs into ColabDesign raw design directories.
- Use old Patch_A / Patch_B / Patch_C definitions as direct RFpeptides inputs.
- Use Boltz as the main validation gate for this branch.
- Treat RFpeptides backbone-only outputs as final peptide candidates.
- Treat ProteinMPNN sequences as final without structural validation.
- Use `diffuser.T=15` except for installation smoke tests.
- Remove Rosetta scoring or AfCycDesign-first validation to make the route
  faster.
- Use Cys-Cys / SG-SG geometry checks for this branch.
- Report any peptide as final before all gates pass.

## 19. Current Recovery Step (updated 2026-07-16)

Current state:

```text
Stage -1 site discovery: retained
Stage 0 target crop and mapping: retained
pre-hotspot-fix Stage 1-5 outputs: quarantined
corrected Stage 1 N20 smoke: completed
strict Stage 2 direct-contact QC: 3/20 pass
ProteinMPNN / Stage 4 / Stage 5 restart: blocked pending manual review
```

The three strict smoke pass IDs are:

```text
RFpep_Site_2_0002
RFpep_Site_2_0017
RFpep_Site_2_0018
```

Review script:

```text
results/rfpeptides_article_route_clean_20260715_hotspotfix_smoke/
  03_backbone_qc_hotspot_mapping_v3_strict/
  RFpep_Site_2_stage2_top_pass_review.pml
```

Before any larger generation run:

1. Confirm the runner's static hotspot-offset preflight passes.
2. Inspect at least one TRB mapping from each peptide length represented.
3. Confirm mapped hotspot PDB residue numbers equal `peptide_length + crop_position`.
4. Open all three strict pass structures in PyMOL and verify genuine Site_2 and
   hotspot contacts rather than contact with another crop surface.
5. Record manual accept/reject decisions before selecting a new small scale.

Do not restart ProteinMPNN, Rosetta scoring, or Stage 5 until this manual review
is complete. The next scale must remain small enough to audit before any
production-size RFpeptides run.
