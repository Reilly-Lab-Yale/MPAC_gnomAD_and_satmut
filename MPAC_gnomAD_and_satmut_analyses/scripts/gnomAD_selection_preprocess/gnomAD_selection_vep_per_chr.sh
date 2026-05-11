#!/bin/bash
#SBATCH --job-name=vep_summaries
#SBATCH --mem=64G
#SBATCH --cpus-per-task=1
#SBATCH --time=1-00:00:00
#SBATCH --output=vep_summaries_chr%a_%j.out
#SBATCH --error=vep_summaries_chr%a_%j.err
#SBATCH --array=1-22

module load miniconda
conda activate mpac

Rscript gnomAD_selection_vep_per_chr.R ${SLURM_ARRAY_TASK_ID}
