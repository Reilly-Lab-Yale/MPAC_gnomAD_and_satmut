#!/bin/sh

bedtools complement -i GRCh38-cCREs.V4.bed.gz -g <(sort -k1,1 -k2,2n ../reference_genomes/Homo_sapiens-hg38-GRCh38/UCSC_hg38/hg38.chrom.sizes) | bgzip > GRCh38-DNase-low.V4.bed.gz 
