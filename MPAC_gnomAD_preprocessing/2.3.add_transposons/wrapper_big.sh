#!/bin/bash
#SBATCH --array=1-2,4
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=64G
#SBATCH -t 14-00:00:00
#SBATCH -p ycga_long
#SBATCH --mail-user mackenzie.noon@yale.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -J big_addrep

module load miniconda
conda activate mcn_varef

jupyter nbconvert --to script add_rep.ipynb


export which_chr="chr${SLURM_ARRAY_TASK_ID}"



# Execute the converted Python script
spark-submit --executor-memory 100g --driver-memory 20g add_rep.py
