#!/bin/bash
set -eo pipefail

cd '/mnt/c/SH/peptide_str/rfd_macro'
source ~/fga_model_envs/miniforge3/etc/profile.d/conda.sh
conda activate 'SE3nv'
set -u

mkdir -p '/mnt/c/SH/fga_cyclic_peptide_design/results/rfpeptides_article_route_clean_20260623_stage1_N10000_L12_24_batch01/02_rfpeptides_backbones/RFpep_Site_2'
mkdir -p '/mnt/c/SH/fga_cyclic_peptide_design/results/rfpeptides_article_route_clean_20260623_stage1_N10000_L12_24_batch01/logs'

echo '[RFpeptides] starting RFpep_Site_2_L12_24_N10000'
{
./scripts/run_inference.py --config-name base \
inference.output_prefix=/mnt/c/SH/fga_cyclic_peptide_design/results/rfpeptides_article_route_clean_20260623_stage1_N10000_L12_24_batch01/02_rfpeptides_backbones/RFpep_Site_2/RFpep_Site_2_L12_24 \
inference.num_designs=10000 \
'contigmap.contigs=[12-24 A1-86/0]' \
inference.input_pdb=/mnt/c/SH/fga_cyclic_peptide_design/results/rfpeptides_article_route_clean_20260615_fpocket/00_target_inputs/RFpep_Site_2_target.pdb \
inference.cyclic=True \
diffuser.T=50 \
inference.cyc_chains='a' \
'ppi.hotspot_res=[A85,A82,A86,A84]' \
inference.write_trajectory=False
} 2>&1 | tee '/mnt/c/SH/fga_cyclic_peptide_design/results/rfpeptides_article_route_clean_20260623_stage1_N10000_L12_24_batch01/logs/rfpeptides_stage1_RFpep_Site_2_L12_24_N10000.log'
echo '[RFpeptides] finished RFpep_Site_2_L12_24_N10000'
