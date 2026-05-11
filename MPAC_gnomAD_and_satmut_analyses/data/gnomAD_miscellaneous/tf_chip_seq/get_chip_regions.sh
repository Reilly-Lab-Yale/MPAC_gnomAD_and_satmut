#!/bin/bash

# get unique regions
gunzip -c chip.for.steve.txt.gz | awk 'NR > 1 {print $7 "\t" $8-1 "\t" $9}' | sort -k1,1 -k2,2n | bedtools merge -i - > chip.for.steve.regions.bed

# lift over to hg38
# used https://genome.ucsc.edu/cgi-bin/hgLiftOver service to get hglft_genome_19baa_20b520.bed
sort -k1,1 -k2,2n hglft_genome_19baa_20b520.bed | bedtools merge -i - > chip.for.steve.regions.hg38.bed

# # get complement regions
# sort -k1,1 -k2,2n hg38.chrom.sizes > hg38.chrom.sizes.sorted
# bedtools complement -i chip.for.steve.regions.hg38.bed -g hg38.chrom.sizes.sorted > chip.for.steve.complement.hg38.bed
