#!/bin/bash
#SBATCH --array=4-22
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=10G
#SBATCH -t 2-00:00:00
#SBATCH -p ycga
#SBATCH --mail-user mackenzie.noon@yale.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -J addrep

module load miniconda
conda activate mcn_varef

jupyter nbconvert --to script add_rep.ipynb


export which_chr="chr${SLURM_ARRAY_TASK_ID}"



# Execute the converted Python script
spark-submit --executor-memory 15g --driver-memory 4g add_rep.py
