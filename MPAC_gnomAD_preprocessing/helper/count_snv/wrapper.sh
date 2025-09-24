#!/bin/bash
#SBATCH -p ycga
#SBATCH -t 6:00:00
#SBATCH -c 1
#SBATCH --mem-per-cpu=4G

module load miniconda
conda activate speedracer

root_dir="/gpfs/gibbs/pi/reilly/VariantEffects/data/gnomAD/gnomAD_genomes_v3.1.2"

zcat ${root_dir}/gnomad.genomes.v3.1.2.sites.chr?.subinfo.vcf.gz \
     ${root_dir}/gnomad.genomes.v3.1.2.sites.chr??.subinfo.vcf.gz | grep -v "^#" | pypy3 fastprocess.py | sort | uniq -c > output.txt
