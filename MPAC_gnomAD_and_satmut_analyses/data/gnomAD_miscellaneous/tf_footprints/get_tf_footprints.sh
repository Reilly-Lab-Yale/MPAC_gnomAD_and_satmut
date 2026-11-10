#!/bin/bash

# get from zenodo
wget -nc https://zenodo.org/records/3905306/files/collapsed_motifs_overlapping_consensus_footprints_hg38.bed.gz

# get unique regions
gunzip -cd collapsed_motifs_overlapping_consensus_footprints_hg38.bed.gz | grep '^chr' | cut -f1-3 | LC_ALL=C sort -k1,1 -k2,2n | bedtools merge -i - | gzip > tf.footprints.regions.bed.gz

# # get complement regions
# bedtools complement -i <(gunzip -cd tf.footprints.regions.bed.gz) -g hg38.chrom.sizes.sorted > tf.footprints.complement.bed
