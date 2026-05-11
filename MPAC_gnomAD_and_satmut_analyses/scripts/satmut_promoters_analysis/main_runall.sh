#!/bin/sh

# conda activate mpac

# Motif analyses
Rscript -e "rmarkdown::render('satmut_promoters_meme.Rmd')"

# Meta-promoter analyses
Rscript -e "rmarkdown::render('satmut_promoters_overall_dist.Rmd')"

# Big heatmap analyses
Rscript -e "rmarkdown::render('satmut_promoters_overall_bigheat.Rmd')"

# Overall correlation analyses
Rscript -e "rmarkdown::render('satmut_promoters_overall_corr.Rmd')"

# Stratified promoter analyses
Rscript -e "rmarkdown::render('satmut_promoters_strat_promoters.Rmd')"

# Overview of individual promoters analyses (must run before individual_promoter)
Rscript -e "rmarkdown::render('satmut_promoters_overall_promoter.Rmd')"

# Individual promoter analyses (depends on overall_promoter output)
Rscript -e "rmarkdown::render('satmut_promoters_individual_promoter.Rmd')"

# Gene set enrichment analysis using enrichR (depends on overall_promoter output)
Rscript -e "rmarkdown::render('satmut_promoter_enrichr_gene_sets.Rmd')"
