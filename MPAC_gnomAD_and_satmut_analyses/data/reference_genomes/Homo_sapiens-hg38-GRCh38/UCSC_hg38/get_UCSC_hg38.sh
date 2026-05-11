#!/bin/sh
wget -nc https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
wget -nc https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes
wget -nc https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chromAlias.txt
wget -nc https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/md5sum.txt

# bgzip and index
gzip -d hg38.fa.gz
bgzip hg38.fa
samtools faidx hg38.fa.gz
