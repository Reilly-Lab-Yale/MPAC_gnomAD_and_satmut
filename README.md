# Identifying non-coding variant effects at scale via machine learning models of cis-regulatory reporter assays

## Introduction

Non-coding variants are a major contributor to human disease, yet interpreting their functional impact at scale remains a challenge. MPAC is an ensemble of machine-learning models trained on Massively Parallel Reporter Assay (MPRA) data that predicts the cis-regulatory impact of non-coding variants. This repository contains code to analyze and visualize two MPAC applications: (1) population-level purifying selection analyses using gnomAD variants, quantifying the relationship between predicted regulatory function and evolutionary constraint across ENCODE cis-regulatory elements; and (2) in-silico saturation mutagenesis across human promoters, identifying positions under widespread selection against variants predicted to disrupt promoter activity.

Preprint available at https://www.biorxiv.org/content/10.1101/2025.04.16.648420v1. MPAC predictions available at https://doi.org/10.5281/zenodo.15178434. See https://github.com/Reilly-Lab-Yale/coda_mpac for MPAC itself. See https://github.com/john-c-butts/MPAC/ for code which produces figures (maintained by John Butts). See https://github.com/Reilly-Lab-Yale/predictions_plusplus for supplemental predictions.

**Citation**: John C. Butts, Stephen Rong, Sager J. Gosai, Rodrigo I. Castro, Mackenzie Noon, Kehinde Adeniran, Rohit Ghosh, Pardis C. Sabeti, Ryan Tewhey, Steven K. Reilly. Identifying non-coding variant effects at scale via machine learning models of cis-regulatory reporter assays. *bioRxiv* 2025. https://doi.org/10.1101/2025.04.16.648420

Everything necessary is in the main branch *except* for indel prediction filtering, which is performed with an alternate version of the MPAC_gnomAD_preprocessing pipeline in the `corrected_indel` branch.

See subdirectory READMEs for detailed content descriptions and software dependencies / setup instructions.

**Last updated**: April 25th 2026 by Stephen Rong (current: srong AT ic DOT ac DOT uk, previous: stephen DOT rong AT yale DOT edu) and Mackenzie Noon (mackenzie DOT noon AT yale DOT edu). Contact corresponding authors or us with questions.

## MPAC_gnomAD_preprocessing
Contains code used to preprocess gnomAD predictions into intermediate files and summary tables. See README there for details. The "Helper" directory also contains analyses for many small tasks referenced elsewhere. Maintained by Mackenzie Noon.

## MPAC_gnomAD_and_satmut_analyses
Contains code used to analyze and visualize data related to the gnomAD and promoter saturation mutagenesis work. See README therein for details. Maintained by Stephen Rong.
