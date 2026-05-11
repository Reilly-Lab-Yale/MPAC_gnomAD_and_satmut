#!/bin/sh
wget -nc http://ftp.ensembl.org/pub/release-107/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
wget -nc http://ftp.ensembl.org/pub/release-107/fasta/homo_sapiens/dna/CHECKSUMS

# bgzip and index
gunzip -d Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
bgzip Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa
samtools faidx Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
