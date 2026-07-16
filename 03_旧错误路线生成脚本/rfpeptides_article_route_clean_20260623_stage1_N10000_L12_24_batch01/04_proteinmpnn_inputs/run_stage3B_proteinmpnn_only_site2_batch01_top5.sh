#!/bin/bash
set -eo pipefail

mkdir -p '/mnt/c/SH/fga_cyclic_peptide_design/results/rfpeptides_article_route_clean_20260623_stage1_N10000_L12_24_batch01/04_proteinmpnn_inputs/work'
mkdir -p '/mnt/c/SH/fga_cyclic_peptide_design/results/rfpeptides_article_route_clean_20260623_stage1_N10000_L12_24_batch01/05_proteinmpnn_sequences/proteinmpnn_only_pdbs'
mkdir -p '/mnt/c/SH/fga_cyclic_peptide_design/results/rfpeptides_article_route_clean_20260623_stage1_N10000_L12_24_batch01/logs'
mkdir -p '/mnt/c/SH/fga_cyclic_peptide_design/results/rfpeptides_article_route_clean_20260623_stage1_N10000_L12_24_batch01/04_proteinmpnn_inputs/checkpoints'
cd '/mnt/c/SH/fga_cyclic_peptide_design/results/rfpeptides_article_route_clean_20260623_stage1_N10000_L12_24_batch01/04_proteinmpnn_inputs/work'
rm -f '/mnt/c/SH/fga_cyclic_peptide_design/results/rfpeptides_article_route_clean_20260623_stage1_N10000_L12_24_batch01/04_proteinmpnn_inputs/checkpoints/proteinmpnn_only_stage3B.checkpoint'
source ~/fga_model_envs/miniforge3/etc/profile.d/conda.sh
conda activate 'proteinmpnn_binder_design'
set -u

echo '[proteinmpnn_only] starting Stage 3'
{
python '/mnt/c/SH/peptide_str/dl_binder_design/mpnn_fr/dl_interface_design.py' \
-pdbdir '/mnt/c/SH/fga_cyclic_peptide_design/results/rfpeptides_article_route_clean_20260623_stage1_N10000_L12_24_batch01/04_proteinmpnn_inputs/pdbs_proteinmpnn_only' \
-outpdbdir '/mnt/c/SH/fga_cyclic_peptide_design/results/rfpeptides_article_route_clean_20260623_stage1_N10000_L12_24_batch01/05_proteinmpnn_sequences/proteinmpnn_only_pdbs' \
-runlist '/mnt/c/SH/fga_cyclic_peptide_design/results/rfpeptides_article_route_clean_20260623_stage1_N10000_L12_24_batch01/04_proteinmpnn_inputs/FGA_rfpeptides_stage3B_proteinmpnn_only_runlist.txt' \
-checkpoint_name '/mnt/c/SH/fga_cyclic_peptide_design/results/rfpeptides_article_route_clean_20260623_stage1_N10000_L12_24_batch01/04_proteinmpnn_inputs/checkpoints/proteinmpnn_only_stage3B.checkpoint' \
-relax_cycles 0 \
-seqs_per_struct 8 \
-temperature 0.1 \
-omit_AAs 'CX'
} 2>&1 | tee '/mnt/c/SH/fga_cyclic_peptide_design/results/rfpeptides_article_route_clean_20260623_stage1_N10000_L12_24_batch01/logs/proteinmpnn_only_stage3B_RFpep_Site_2_0007.log'
echo '[proteinmpnn_only] finished Stage 3'
