# Identifying non-coding variant effects at scale via machine learning models of cis-regulatory reporter assays

## Introduction

Non-coding variants are a major contributor to human disease, yet interpreting their functional impact at scale remains a challenge. MPAC is an ensemble of machine-learning models trained on Massively Parallel Reporter Assay (MPRA) data that predicts the cis-regulatory impact of non-coding variants. This repository contains code to analyze and visualize two MPAC applications: (1) population-level purifying selection analyses using gnomAD variants, quantifying the relationship between predicted regulatory function and evolutionary constraint across ENCODE cis-regulatory elements; and (2) in-silico saturation mutagenesis across human promoters, identifying positions under widespread selection against variants predicted to disrupt promoter activity.

Code used to analyze and visualize MPAC gnomAD and saturation mutagenesis predictions. MPAC predictions available at https://doi.org/10.5281/zenodo.15178434. See complementary repository https://github.com/john-c-butts/MPAC/.

**Citation**: John C. Butts, Stephen Rong, Sager J. Gosai, Rodrigo I. Castro, Mackenzie Noon, Kehinde Adeniran, Rohit Ghosh, Pardis C. Sabeti, Ryan Tewhey, Steven K. Reilly. Identifying non-coding variant effects at scale via machine learning models of cis-regulatory reporter assays. *bioRxiv* 2025. https://doi.org/10.1101/2025.04.16.648420

Last updated April 25th 2026 by Stephen Rong (current: srong AT ic DOT ac DOT uk, previous: stephen DOT rong AT yale DOT edu). Contact corresponding authors or us with questions.

## Content descriptions

### Scripts
- **scripts/gnomAD_selection_preprocess**: Preprocessing scripts for MPAC gnomAD predictions. Per-chromosome scripts for SNPs, indels, and VEP consequence analyses (VEP annotations are extracted from gnomAD VCFs) are run as Slurm array jobs, merging gnomAD variant calls with MPAC prediction files and additional annotations (Roulette mutation rates, phyloP conservation scores, TF ChIP-seq peaks, TF footprints, ENCODE cCREs), then combined into genome-wide summary tables.

- **scripts/gnomAD_selection_analysis**: Analyses and visualizations for MPAC gnomAD predictions (RMarkdown). Includes purifying selection analyses for SNPs, indels, and VEP consequences, with comparisons across ENCODE cCRE classes using allelic skew, activity, emVar specificity, phyloP constraint, mutation rates, and singleton vs. common allele frequency tests.

- **scripts/satmut_promoters_preprocess**: Preprocessing scripts for MPAC saturation mutagenesis predictions. Per-chromosome scripts annotate variants with promoter regions, phyloP conservation scores, exon/splice masks, and emVar classifications, then summarize at base, promoter, and TSS-distance levels with gene-level constraint and expression metadata.

- **scripts/satmut_promoters_analysis**: Analyses and visualizations for MPAC saturation mutagenesis predictions (RMarkdown). Includes meta-promoter distance profiles, genome-wide heatmaps, cross-cell-line correlation analyses, stratified promoter analyses by gene constraint, individual promoter plots, and gene set enrichment analyses.

### Results
- **results/satmut_promoters_analysis**: Corresponding results folder for scripts/satmut_promoters_analysis.

- **results/gnomAD_purifying_selection**: Corresponding results folder for scripts/gnomAD_selection_analysis.

### Data
- **data/gencode_filtered_regions**: Code and BED files for masking exonic and splice regions.

- **data/gene_constraint_metrics**: GeneBayes s_het, gnomAD LOEUF/MOEUF, AlphaMissense gene constraint, and ABC.

- **data/gene_expression_catalogs**: Human Protein Atlas gene expression annotations for K562, HepG2, and SK-N-SH.

- **data/gene_regulatory_elements**: ENCODE cCRE annotations and cell-line-specific BigWig signal tracks.

- **data/gnomAD_genomes_v3**: gnomAD v3.1.2 VCFs and subsetted annotations.

- **data/gnomAD_indels_predictions**: MPAC gnomAD indel predictions.

- **data/gnomAD_indels_summaries**: Per-chromosome and combined summary tables for gnomAD indel analyses.

- **data/gnomAD_miscellaneous**: TF ChIP-seq peak and TF footprint annotation files.

- **data/gnomAD_roulette_predictions**: Roulette mutation rate predictions.

- **data/gnomAD_snp_predictions**: MPAC gnomAD SNP predictions.

- **data/gnomAD_snp_summaries**: Per-chromosome and combined summary tables for gnomAD SNP analyses.

- **data/gnomAD_vep_summaries**: Per-chromosome and combined summary tables for gnomAD VEP consequence analyses.

- **data/reference_genomes**: Reference genome sequences (Ensembl GRCh38.p13 and UCSC hg38).

- **data/satmut_promoters_meme_pwms**: PWM meme files for saturation mutagenesis motif analysis.

- **data/satmut_promoters_predictions**: Copies of MPAC saturation mutagenesis variant predictions.

- **data/satmut_promoters_preprocess**: Preprocessed intermediate files for saturation mutagenesis analyses.

- **data/zoonomia_phylop**: Zoonomia phyloP base-level annotations.

### Misc
- **misc**: Supplementary Tables.
  - `Supplementary_Data_X1_Promoter-level_function_x_constraint_correlations_and_gene_metadata_revised.xlsx`: Promoter-level function x constraint correlations and gene metadata.
  - `Supplementary_Data_X2_Enrichr_results_for_promoter-level_function_x_constraint_correlation_classes_revised.xlsx`: Enrichr results for promoter-level function x constraint correlation classes.

## Software dependencies
Analyses performed on Yale University HPC as Slurm scripts or in RStudio. RMarkdown analysis scripts (`.Rmd`) in `scripts/satmut_promoters_analysis` and `scripts/gnomAD_selection_analysis` were rendered with RStudio.

### Environment setup
All dependencies are managed via conda, with two additional R packages installed after environment creation. To set up:
```bash
conda env create -f mpac.yaml
conda activate mpac
R -e 'install.packages(c("enrichR", "remotes")); remotes::install_github("snystrom/memes"); BiocManager::install("universalmotif")'
```

Note: The `memes` R package additionally requires the [MEME Suite](https://meme-suite.org/) command-line tools to be installed and available on `$PATH`.

### Software versions
- htslib=1.21
- bcftools=1.21
- bedtools=v2.31.1
- R=4.4.1 (2024-06-14)
- MEME Suite=5.5.5 (required by memes R package)
