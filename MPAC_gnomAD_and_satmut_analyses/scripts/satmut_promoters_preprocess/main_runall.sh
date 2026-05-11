#!/bin/bash
#SBATCH --job-name=satmut_preprocess
#SBATCH --mem=64G
#SBATCH --cpus-per-task=1
#SBATCH --time=1-00:00:00
#SBATCH --output=satmut_preprocess_%j.out
#SBATCH --error=satmut_preprocess_%j.err

# conda activate mpac

mkdir -p ../../data/satmut_promoters_preprocess/
mkdir -p ../../data/satmut_promoters_preprocess/satmut_promoters_regions
mkdir -p ../../data/satmut_promoters_preprocess/satmut_promoters_pred
mkdir -p ../../data/satmut_promoters_preprocess/satmut_promoters_pred_phyloP
mkdir -p ../../data/satmut_promoters_preprocess/satmut_promoters_final
mkdir -p ../../data/satmut_promoters_preprocess/satmut_promoters_final_base
mkdir -p ../../data/satmut_promoters_preprocess/satmut_promoters_final_prom
mkdir -p ../../data/satmut_promoters_preprocess/satmut_promoters_final_dist
mkdir -p ../../data/satmut_promoters_preprocess/satmut_promoters_gene_metadata

# Preprocess regions and annotations
Rscript preprocess_promoter_regions.R

# Annotate variants by annotations
Rscript preprocess_promoter_pred.R
Rscript preprocess_promoter_phyloP.R
Rscript preprocess_promoter_final.R

# Annotate summaries at base, promoter, and distance level
Rscript preprocess_promoter_base.R
Rscript preprocess_promoter_prom.R
Rscript preprocess_promoter_dist.R

# Annotate gene-level coding constraint and gene expression
Rscript preprocess_gene_metadata.R
