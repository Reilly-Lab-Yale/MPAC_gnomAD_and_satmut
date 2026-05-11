#!/bin/bash

# get unique regions
gunzip -cd collapsed_motifs_overlapping_consensus_footprints_hg38_GRanges.txt.gz | awk 'NR > 1 {print $1 "\t" $2-1 "\t" $3}' | sort -k1,1 -k2,2n | bedtools merge -i - > tf.footprints.regions.bed

# # get complement regions
# bedtools complement -i tf.footprints.regions.bed -g hg38.chrom.sizes.sorted > tf.footprints.complement.bed
