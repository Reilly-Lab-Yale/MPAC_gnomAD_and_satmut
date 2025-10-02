#!/bin/bash
#SBATCH --job-name=merge_vcfs
#SBATCH --output=merge_vcfs_%A_%a.out
#SBATCH --error=merge_vcfs_%A_%a.err
#SBATCH --array=1-22%22
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --time=08:00:00
#SBATCH -p ycga

module load miniconda
conda activate mcn_varef

# Define the directories
DIR1="/home/mcn26/varef/data/Malinois/gnomAD_variants/"
DIR2="/gpfs/gibbs/pi/reilly/VariantEffects/data/gnomAD/gnomAD_genomes_v3.1.2/"
OUT_DIR="/gpfs/gibbs/pi/reilly/VariantEffects/scripts/noon_data/0.merge/out/" # make sure this directory exists

# Define the chromosome based on the Slurm array task ID
CHR="chr${SLURM_ARRAY_TASK_ID}"

bcftools merge  \
   "${DIR1}gnomad.genomes.v3.1.2.sites.${CHR}.converted.vcf.gz"\
   "${DIR2}gnomad.genomes.v3.1.2.sites.${CHR}.subinfo.vcf.gz"\
  -Oz -o "${OUT_DIR}combined.${CHR}.vcf.gz"


# Index the combined VCF file
bcftools index -f -t --threads ${SLURM_CPUS_PER_TASK} "${OUT_DIR}combined.${CHR}.vcf.gz"
