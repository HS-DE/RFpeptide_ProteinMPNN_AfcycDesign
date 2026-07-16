# FGA RFpeptides Stage 2 Backbone QC

Status: RFpeptides backbone outputs parsed and checked before sequence design.

Important rule: Stage 2 does not treat whole-crop contact as sufficient.
A design must make direct contacts within `contact_cutoff_A` to the intended
`RFpep_Site_2` target-site residues and selected hotspots. Near-only states are
reported but fail by default. Designs that contact only a distant part of the
target crop are flagged as `crop_only_contact`.

Output directory:

```text
/mnt/c/SH/fga_cyclic_peptide_design/results/rfpeptides_article_route_clean_20260716_runtimefix_smoke_N1_v3/03_backbone_qc
```

Parameters:

```text
stage0_root: results/rfpeptides_article_route_clean_20260615_fpocket
stage1_root: results/rfpeptides_article_route_clean_20260716_runtimefix_smoke_N1_v3
output_subdir: 03_backbone_qc
selected_sites: RFpep_Site_2
contact_cutoff_A: 5.0
site_near_distance_A: 6.0
hotspot_near_distance_A: 8.0
severe_clash_distance_A: 1.2
min_target_contacts: 3
min_site_contacts: 1
min_hotspot_contacts: 1
allow_near_pass: False
allow_missing_runtime_audit: False
```

Target-site residue numbers:

```text
6,9,10,13,79,80,81,82,83,84,85,86
```

These are Stage 0 crop sequence positions, not raw residue numbers in an
RFpeptides output PDB. Each design maps them by target-chain residue order;
the mapped output-PDB residue numbers are recorded per row.

Production QC requires a full Stage 1 -> Stage 2 provenance closure. The final
runtime-audit JSON must exactly equal the audit embedded in the TRB; locked
Stage 0 file hashes must match; the TRB ContigMap, helper indices, actual model
hotspot tensor, writer `chain_idx`/`idx_pdb`, and output PDB identities must all
agree for every hotspot. Any mismatch is a hard failure and cannot enter Stage
3. Structural contact, macrocycle, and clash QC remain separate requirements.

Hotspot residue numbers:

```text
82,84,85,86
```

## Counts

```text
total_backbones: 1
pass_backbone_qc: 0
crop_only_contact: 0
```

## Target-Site Recovery Status

- site_missed: 1

## Runtime Audit Status

- pass: 1

## JSON/TRB Audit Identity

- pass_exact_match: 1

## Hotspot Provenance Closure

- pass_stage0_trb_tensor_pdb_closed: 1

## Model Hotspot Tensor

- pass_helper_contigmap_tensor_equal: 1

## Target Sequence Identity Status

- pass_exact_stage0_target_sequence: 1

## Hotspot Recovery Status

- hotspot_missed: 1

## Macrocycle Geometry Status

- pass_head_to_tail_macrocycle: 1

## Top Passing Backbones

No backbones passed Stage 2 QC.
