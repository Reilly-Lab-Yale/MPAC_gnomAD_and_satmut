#!/bin/sh

# conda activate mpac

# gnomAD MPAC snp analyses
Rscript -e "rmarkdown::render('gnomAD_selection_snp.Rmd')"

# gnomAD MPAC vep analyses
Rscript -e "rmarkdown::render('gnomAD_selection_vep.Rmd')"

# gnomAD MPAC indels analyses
Rscript -e "rmarkdown::render('gnomAD_selection_indels.Rmd')"
