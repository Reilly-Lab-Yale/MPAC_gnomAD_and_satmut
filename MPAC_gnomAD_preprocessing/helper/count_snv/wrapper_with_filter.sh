#!/bin/bash
#SBATCH -p ycga
#SBATCH -t 6:00:00
#SBATCH -c 1
#SBATCH --mem-per-cpu=4G
#SBATCH -J count_snv_filter
#SBATCH -o /vast/palmer/pi/reilly/VariantEffects/scripts/noon_scripts/MPAC_gnomAD_preprocessing/helper/count_snv/filter_job_%j.out

module load miniconda
conda activate speedracer

root_dir="/gpfs/gibbs/pi/reilly/VariantEffects/data/gnomAD/gnomAD_genomes_v3.1.2"
out_dir="/vast/palmer/pi/reilly/VariantEffects/scripts/noon_scripts/MPAC_gnomAD_preprocessing/helper/count_snv"

zcat ${root_dir}/gnomad.genomes.v3.1.2.sites.chr?.subinfo.vcf.gz \
     ${root_dir}/gnomad.genomes.v3.1.2.sites.chr??.subinfo.vcf.gz \
    | grep -v "^#" \
    | pypy3 ${out_dir}/fastprocess_with_filter.py \
    | sort | uniq -c > ${out_dir}/output_with_filter.txt

echo "DONE"
