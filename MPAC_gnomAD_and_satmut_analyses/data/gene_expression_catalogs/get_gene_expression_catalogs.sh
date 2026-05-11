#!/bin/bash

# Download gene expression catalogs from the Human Protein Atlas (version 25.0)
# https://www.proteinatlas.org/about/download

# RNA expression in cell lines (1206 cell lines)
wget "https://www.proteinatlas.org/download/tsv/rna_celline.tsv.zip"
unzip rna_celline.tsv.zip
gzip rna_celline.tsv
rm rna_celline.tsv.zip

# RNA expression per single cell type (154 cell types)
wget "https://www.proteinatlas.org/download/tsv/rna_single_cell_type.tsv.zip"
unzip rna_single_cell_type.tsv.zip
gzip rna_single_cell_type.tsv
rm rna_single_cell_type.tsv.zip

# RNA expression consensus tissue (51 tissues)
wget "https://www.proteinatlas.org/download/tsv/rna_tissue_consensus.tsv.zip"
unzip rna_tissue_consensus.tsv.zip
gzip rna_tissue_consensus.tsv
rm rna_tissue_consensus.tsv.zip
