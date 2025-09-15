#!/bin/bash
#SBATCH -p ycga
#SBATCH --job-name=dump_indels
#SBATCH --array=1-22
#SBATCH --mem-per-cpu=1G
#SBATCH --cpus-per-task=4

DATA_ROOT="/gpfs/gibbs/pi/reilly/VariantEffects/data/gnomAD/gnomAD_genomes_v3.1.2"

DATA_OUTPUT_ROOT="/home/mcn26/varef/scripts/noon_data/split_off_indel"

CHR="chr${SLURM_ARRAY_TASK_ID}"

module load VCFtools

vcftools --gzvcf ${DATA_ROOT}/gnomad.genomes.v3.1.2.sites.${CHR}.subinfo.vcf.gz --keep-only-indels --recode --recode-INFO-all  --stdout | bgzip -c > ${DATA_OUTPUT_ROOT}/gnomad_${CHR}_indels_only.vcf.gz



