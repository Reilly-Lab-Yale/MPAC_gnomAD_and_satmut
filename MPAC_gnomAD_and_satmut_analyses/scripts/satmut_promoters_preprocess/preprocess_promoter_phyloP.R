#!/bin/R

# Annotate phyloP scores from Zoonomia mammal and Kuderna primate and mammal

# load libraries
library(data.table)
library(tidyverse)

# load bioconductor
library(rtracklayer)

# loop chromosomes
for (chr in paste0("chr", c(22:1))) {
	print(chr)

	# load predications
	satmut_promoters_pred <- as_tibble(fread(paste0("../../data/satmut_promoters_preprocess/satmut_promoters_pred/satmut_promoters_pred_", chr, ".txt.gz"))) %>% 
		dplyr::select(chrom, pos, ref, alt, id, strand, tss_dist)  # , ref_orig, alt_orig)

	# load phyloP from hg38
	satmut_promoters_phyloP_mam241 <- import("../../data/zoonomia_phylop/data_download/241-mammalian-2020v2.phylop-Homo_sapiens.bigWig", 
		"bigWig", selection=BigWigSelection(import(paste0("../../data/satmut_promoters_preprocess/satmut_promoters_regions/satmut_promoters_promoter_regions_", chr, ".bed"))))

	# join preds and phyloP
	satmut_promoters_phyloP_mam241 <- satmut_promoters_phyloP_mam241 %>% 
		as_tibble() %>% 
		mutate(
			chrom = seqnames, 
			pos = start, 
			phyloP_mam241 = score, 
			phyloP_mam241_pos = pmax(0, score), 
			phyloP_mam241_neg = pmin(0, score), 
			phyloP_mam241_cons = (score >= 2.27)
		) %>% 
		dplyr::select(chrom, pos, phyloP_mam241, phyloP_mam241_pos, phyloP_mam241_neg, phyloP_mam241_cons)

	satmut_promoters_pred_phyloP <- satmut_promoters_pred %>% 
		left_join(satmut_promoters_phyloP_mam241)

	# save preds and phyloP
	write_tsv(satmut_promoters_pred_phyloP, gzfile(paste0("../../data/satmut_promoters_preprocess/satmut_promoters_pred_phyloP/satmut_promoters_pred_phyloP_", chr, ".txt.gz")))
}
