#!/bin/bash
#SBATCH --job-name=indel_summaries
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --time=4:00:00
#SBATCH --output=indel_summaries_chr%a_%j.out
#SBATCH --error=indel_summaries_chr%a_%j.err
#SBATCH --array=1-22

module load miniconda
conda activate mpac

Rscript gnomAD_selection_indels_per_chr.R ${SLURM_ARRAY_TASK_ID}
