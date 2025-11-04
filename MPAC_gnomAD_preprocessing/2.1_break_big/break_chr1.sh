#!/bin/bash
#SBATCH -p ycga
#SBATCH -t 8:00:00
#SBATCH -J break_chr1
#SBATCH -c 1
#SBATCH --mem-per-cpu=12G
INFILE="/home/mcn26/varef/scripts/noon_data/indel/2.0.annotate/annotated_output_chr1.csv.gz/part-00001-15599e78-07fb-4200-bc47-e0375a2f4961-c000.csv.gz"
#last step consolidated, so just one input file

module load miniconda
conda activate speedracer
echo "computing total"
TOTAL=$(zcat $INFILE | wc -l)
echo "total: ${TOTAL}"
echo "split!"
N=10

OUTDIR="/home/mcn26/varef/scripts/noon_data/indel/2.1_break/chr1"

# Run the splitter
zcat "$INFILE" | pypy3 break.py "$N" "$TOTAL" part_ "$OUTDIR"
echo "done"