#!/bin/sh

# # check Roulette MR range
# sbatch check_roulette_mr_range.sh

# gnomAD MPAC SNP analyses
sbatch gnomAD_selection_snp_per_chr.sh
Rscript gnomAD_selection_snp_combine_chr.R

# gnomAD VEP SNP analyses
sbatch gnomAD_selection_vep_per_chr.sh
Rscript gnomAD_selection_vep_combine_chr.R

# gnomAD MPAD indels anlayses
sbatch gnomAD_selection_indels_per_chr.sh
Rscript gnomAD_selection_indels_combine_chr.R
