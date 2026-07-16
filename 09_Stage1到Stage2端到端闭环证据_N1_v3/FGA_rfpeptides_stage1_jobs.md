# FGA RFpeptides Stage 1 Jobs

Status: RFpeptides command table and run script prepared. No backbone
generation was run by this preparation script.

Run script:

```text
/mnt/c/SH/fga_cyclic_peptide_design/results/rfpeptides_article_route_clean_20260716_runtimefix_smoke_N1_v3/01_rfpeptides_jobs/run_rfpeptides_stage1_site2_runtimefix_smoke_N1_v3.sh
```

| rfpeptides_job_id | site_label | num_designs | contigmap_contigs | hotspot_res | output_prefix | status |
| --- | --- | --- | --- | --- | --- | --- |
| RFpep_Site_2_L17_17_N1 | RFpep_Site_2 | 1 | [17-17 A1-86/0] | ['A82','A84','A85','A86'] | /mnt/c/SH/fga_cyclic_peptide_design/results/rfpeptides_article_route_clean_20260716_runtimefix_smoke_N1_v3/02_rfpeptides_backbones/RFpep_Site_2/RFpep_Site_2_L17_17 | pending_manual_execution |

Important:

- The runner pins one RFpeptides runtime root, Git commit, and four source-file
  hashes. It aborts if preflight and inference would import different code.
- This is a small pilot in design count, not a reduced-constraint run.
- The command keeps the RFpeptides binder requirements: target PDB, target
  contig, cyclic generation, cyclic chain, diffuser timesteps, and hotspot
  residues.
- RFpep_Site_3 and RFpep_Site_4 remain deferred and are not included in this
  Stage 1 pilot job table.
