#!/bin/bash
set -eo pipefail

mkdir -p '/mnt/c/SH/fga_cyclic_peptide_design/results/rfpeptides_article_route_clean_20260615_fpocket_stage1_N1000_no_traj/04_proteinmpnn_inputs/work'
mkdir -p '/mnt/c/SH/fga_cyclic_peptide_design/results/rfpeptides_article_route_clean_20260615_fpocket_stage1_N1000_no_traj/05_proteinmpnn_sequences/fastrelax_pdbs'
mkdir -p '/mnt/c/SH/fga_cyclic_peptide_design/results/rfpeptides_article_route_clean_20260615_fpocket_stage1_N1000_no_traj/logs'
mkdir -p '/mnt/c/SH/fga_cyclic_peptide_design/results/rfpeptides_article_route_clean_20260615_fpocket_stage1_N1000_no_traj/04_proteinmpnn_inputs/checkpoints'
cd '/mnt/c/SH/fga_cyclic_peptide_design/results/rfpeptides_article_route_clean_20260615_fpocket_stage1_N1000_no_traj/04_proteinmpnn_inputs/work'
rm -f '/mnt/c/SH/fga_cyclic_peptide_design/results/rfpeptides_article_route_clean_20260615_fpocket_stage1_N1000_no_traj/04_proteinmpnn_inputs/checkpoints/proteinmpnn_fastrelax_stage3.checkpoint'
source ~/fga_model_envs/miniforge3/etc/profile.d/conda.sh
conda activate 'proteinmpnn_binder_design'
set -u

echo '[ProteinMPNN-FastRelax] starting Stage 3'
{
python '/mnt/c/SH/peptide_str/dl_binder_design/mpnn_fr/dl_interface_design.py' \
-pdbdir '/mnt/c/SH/fga_cyclic_peptide_design/results/rfpeptides_article_route_clean_20260615_fpocket_stage1_N1000_no_traj/04_proteinmpnn_inputs/pdbs' \
-outpdbdir '/mnt/c/SH/fga_cyclic_peptide_design/results/rfpeptides_article_route_clean_20260615_fpocket_stage1_N1000_no_traj/05_proteinmpnn_sequences/fastrelax_pdbs' \
-runlist '/mnt/c/SH/fga_cyclic_peptide_design/results/rfpeptides_article_route_clean_20260615_fpocket_stage1_N1000_no_traj/04_proteinmpnn_inputs/FGA_rfpeptides_stage3_runlist.txt' \
-checkpoint_name '/mnt/c/SH/fga_cyclic_peptide_design/results/rfpeptides_article_route_clean_20260615_fpocket_stage1_N1000_no_traj/04_proteinmpnn_inputs/checkpoints/proteinmpnn_fastrelax_stage3.checkpoint' \
-relax_cycles 1 \
-seqs_per_struct 1 \
-temperature 0.1 \
-omit_AAs 'CX'
} 2>&1 | tee '/mnt/c/SH/fga_cyclic_peptide_design/results/rfpeptides_article_route_clean_20260615_fpocket_stage1_N1000_no_traj/logs/proteinmpnn_fastrelax_stage3_RFpep_Site_2_0007.log'
echo '[ProteinMPNN-FastRelax] finished Stage 3'
