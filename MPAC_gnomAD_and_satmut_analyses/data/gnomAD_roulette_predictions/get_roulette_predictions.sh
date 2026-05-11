#!/bin/sh
#SBATCH --time=1-00:00:00

for i in $(seq 1 22); do
  wget http://genetics.bwh.harvard.edu/downloads/Vova/Roulette/${i}_rate_v5.2_TFBS_correction_all.vcf.bgz
  wget http://genetics.bwh.harvard.edu/downloads/Vova/Roulette/${i}_rate_v5.2_TFBS_correction_all.vcf.bgz.csi
done
