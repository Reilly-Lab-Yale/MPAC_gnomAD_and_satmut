#!/bin/bash
#SBATCH --array=0-9
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=16G
#SBATCH -t 14-00:00:00
#SBATCH -p ycga_long
#SBATCH --mail-user mackenzie.noon@yale.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -J 1addrep

module load miniconda
conda activate mcn_varef

jupyter nbconvert --to script add_rep_broken.ipynb


export which_chr="chr2"
export part=${SLURM_ARRAY_TASK_ID}


# Execute the converted Python script
spark-submit --executor-memory 54g --driver-memory 9g add_rep_broken.py
