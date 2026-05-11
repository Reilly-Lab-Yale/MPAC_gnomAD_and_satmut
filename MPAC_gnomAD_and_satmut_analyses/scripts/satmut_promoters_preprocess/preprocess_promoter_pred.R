#!/bin/R

# Process predictions and add in emVar annotations and promoter/masked regions

# load libraries
library(data.table)
library(tidyverse)
library(stringi)

# load bioconductor
library(plyranges)
library(rtracklayer)

# load promoter/exon data
satmut_promoters_promoter_regions_1kb <- readRDS("../../data/satmut_promoters_preprocess/satmut_promoters_regions/satmut_promoters_promoter_regions_1kb.rds")
satmut_promoters_promoter_regions_750bp <- readRDS("../../data/satmut_promoters_preprocess/satmut_promoters_regions/satmut_promoters_promoter_regions_750bp.rds")
satmut_promoters_promoter_regions_500bp <- readRDS("../../data/satmut_promoters_preprocess/satmut_promoters_regions/satmut_promoters_promoter_regions_500bp.rds")
satmut_promoters_promoter_regions_250bp <- readRDS("../../data/satmut_promoters_preprocess/satmut_promoters_regions/satmut_promoters_promoter_regions_250bp.rds")
satmut_promoters_promoter_regions_tss <- readRDS("../../data/satmut_promoters_preprocess/satmut_promoters_regions/satmut_promoters_promoter_regions_tss.rds")
satmut_promoters_all_exon_regions <- rtracklayer::import("../../data/gencode_filtered_regions/gencode.v44.basic.annotation.exons.autosomes.v2.bed", extraCols = c(id = "character"))
satmut_promoters_all_exon_splice_regions <- rtracklayer::import("../../data/gencode_filtered_regions/gencode.v44.basic.annotation.exons.splice.autosomes.v2.bed", extraCols = c(id = "character"))

# loop chromosomes
for (chr in paste0("chr", c(22:1))) {
	print(chr)

	# load predictions
	satmut_promoters_pred <- as_tibble(fread(paste0("../../data/satmut_promoters_predictions/gencode.v44.canonical.protein.coding.1kb.promoters.sat.mut.updated.pos.", chr, ".vcf.gz")))
	mpac_cols <- c("K562_ref", "HepG2_ref", "SKNSH_ref", "K562_alt", "HepG2_alt", "SKNSH_alt", "K562_skew", "HepG2_skew", "SKNSH_skew")
	satmut_promoters_pred <- satmut_promoters_pred %>% 
		mutate(INFO = stri_split_fixed(INFO, ";")) %>% unnest(INFO) %>% mutate(INFO = as.numeric(gsub(".*=", "", INFO))) %>% 
		(function(x) {x$TEMP <- rep(mpac_cols, nrow(x)/length(mpac_cols)); return(x)})(.) %>% pivot_wider(names_from="TEMP", values_from=INFO)

	# clean up formating
	satmut_promoters_pred <- satmut_promoters_pred %>% 
		mutate(id = gsub("\\.\\..*", "", id))
	# satmut_promoters_pred <- satmut_promoters_pred %>% 
	# 	distinct()  # gets rid of repeats daue to original bed file used for prediction containing multiple promoters of some genes

	# get strand  # and flip ref and alt
	satmut_promoters_promoter_regions_1kb_strand <- satmut_promoters_promoter_regions_1kb %>% 
		as_tibble() %>% 
		dplyr::select(id, strand)

	satmut_promoters_pred <- satmut_promoters_pred %>% 
		left_join(satmut_promoters_promoter_regions_1kb_strand)

	# get distance to tss for meta promoters
	satmut_promoters_promoter_regions_tss_pos <- satmut_promoters_promoter_regions_tss %>% 
		as_tibble() %>% 
		mutate(tss = ifelse(strand == "+", start, end)) %>% 
		dplyr::select(id, tss)

	satmut_promoters_pred <- satmut_promoters_pred %>% 
		left_join(satmut_promoters_promoter_regions_tss_pos) %>% 
		mutate(tss_dist = ifelse(strand == "+", pos - tss, tss - pos)) %>% 
		dplyr::select(-tss)

	# get final columns
	satmut_promoters_pred <- satmut_promoters_pred %>% 
		dplyr::select(chrom, pos, ref, alt, id, strand, tss_dist, everything())

	# get emVar summaries
	satmut_promoters_pred <- satmut_promoters_pred %>% 
		# activity - temp
		mutate(
			K562_activity = pmax(K562_ref, K562_alt),
			HepG2_activity = pmax(HepG2_ref, HepG2_alt),
			SKNSH_activity = pmax(SKNSH_ref, SKNSH_alt)
		) %>% 
		# abs skew
		mutate(
			K562_abs_skew = abs(K562_skew),
			HepG2_abs_skew = abs(HepG2_skew),
			SKNSH_abs_skew = abs(SKNSH_skew)
		) %>% 
		# mean activity
		mutate(
			avgKHS_ref = (K562_ref+HepG2_ref+SKNSH_ref)/3, 
			avgKHS_alt = (K562_alt+HepG2_alt+SKNSH_alt)/3, 
			avgKHS_skew = (K562_skew+HepG2_skew+SKNSH_skew)/3, 
			avgKHS_activity = (K562_activity+HepG2_activity+SKNSH_activity)/3, 
			avgKHS_abs_skew = (K562_abs_skew+HepG2_abs_skew+SKNSH_abs_skew)/3, 
		) %>% 
		# emVars
		mutate(
			K562_emVar = (K562_abs_skew > 0.5),
			HepG2_emVar = (HepG2_abs_skew > 0.5),
			SKNSH_emVar = (SKNSH_abs_skew > 0.5),
			avgKHS_emVar = (avgKHS_abs_skew > 0.5),
		) %>% 
		mutate(
			anyKHS_emVar = (K562_emVar | HepG2_emVar | SKNSH_emVar)
		) %>% 
		# emVars pos
		mutate(
			K562_emVar_pos = (K562_skew > 0.5),
			HepG2_emVar_pos = (HepG2_skew > 0.5),
			SKNSH_emVar_pos = (SKNSH_skew > 0.5),
			avgKHS_emVar_pos = (avgKHS_skew > 0.5),
		) %>% 
		mutate(
			anyKHS_emVar_pos = (K562_emVar_pos | HepG2_emVar_pos | SKNSH_emVar_pos),
		) %>% 
		# emVars neg
		mutate(
			K562_emVar_neg = (K562_skew < -0.5),
			HepG2_emVar_neg = (HepG2_skew < -0.5),
			SKNSH_emVar_neg = (SKNSH_skew < -0.5),
			avgKHS_emVar_neg = (avgKHS_skew < -0.5),
		) %>% 
		mutate(
			anyKHS_emVar_neg = (K562_emVar_neg | HepG2_emVar_neg | SKNSH_emVar_neg)
		) %>% 
		# pleiotropy
		mutate(emVar_category = 
			case_when(
				K562_emVar & SKNSH_emVar & HepG2_emVar ~ "K562+SKNSH+HepG2 emVar",
				K562_emVar & SKNSH_emVar ~ "K562+SKNSH emVar",
				K562_emVar & HepG2_emVar ~ "K562+HepG2 emVar",
				SKNSH_emVar & HepG2_emVar ~ "SKNSH+HepG2 emVar",
				K562_emVar ~ "K562 emVar",
				HepG2_emVar ~ "HepG2 emVar",
				SKNSH_emVar ~ "SKNSH emVar",
				.default = "none"
			)
		) %>% 
		mutate(
			emVar_count = (K562_emVar + HepG2_emVar + SKNSH_emVar)
		)

	# get promoter overlap
	overlap_promoter_1kb <- as_tibble(satmut_promoters_promoter_regions_1kb) %>% 
		mutate(chrom = seqnames) %>%  # , id = name) %>% 
		dplyr::select(chrom, start, end, id)
	overlap_promoter_regions_750bp <- as_tibble(satmut_promoters_promoter_regions_750bp) %>% 
		mutate(chrom = seqnames) %>%  # , id = name) %>% 
		dplyr::select(chrom, start, end, id)
	overlap_promoter_regions_500bp <- as_tibble(satmut_promoters_promoter_regions_500bp) %>% 
		mutate(chrom = seqnames) %>%  # , id = name) %>% 
		dplyr::select(chrom, start, end, id)
	overlap_promoter_regions_250bp <- as_tibble(satmut_promoters_promoter_regions_250bp) %>% 
		mutate(chrom = seqnames) %>%  # , id = name) %>% 
		dplyr::select(chrom, start, end, id)

	satmut_promoters_pred <- satmut_promoters_pred %>% 
		left_join(overlap_promoter_1kb) %>% 
		mutate(in_promoter_1kb = ((pos >= start) & (pos <= end))) %>% 
		dplyr::select(-start, -end)

	satmut_promoters_pred <- satmut_promoters_pred %>% 
		left_join(overlap_promoter_regions_750bp) %>% 
		mutate(in_promoter_750bp = ((pos >= start) & (pos <= end))) %>% 
		dplyr::select(-start, -end)

	satmut_promoters_pred <- satmut_promoters_pred %>% 
		left_join(overlap_promoter_regions_500bp) %>% 
		mutate(in_promoter_500bp = ((pos >= start) & (pos <= end))) %>% 
		dplyr::select(-start, -end)

	satmut_promoters_pred <- satmut_promoters_pred %>% 
		left_join(overlap_promoter_regions_250bp) %>% 
		mutate(in_promoter_250bp = ((pos >= start) & (pos <= end))) %>% 
		dplyr::select(-start, -end)

	# remove first base of exon
	satmut_promoters_pred <- satmut_promoters_pred %>% 
		filter(in_promoter_1kb)

	# get promoter quartile
	satmut_promoters_pred <- satmut_promoters_pred %>% 
		mutate(in_promoter_quartile = 
			ifelse(in_promoter_250bp, "in_promoter_250bp", 
				ifelse(in_promoter_500bp, "in_promoter_500bp", 
					ifelse(in_promoter_750bp, "in_promoter_750bp", 
						ifelse(in_promoter_1kb, "in_promoter_1kb", NA)))))

	# get exon overlap
	satmut_promoters_pred_gr <- satmut_promoters_pred %>% 
		mutate(start = pos, end = pos) %>% 
		GRanges()

	overlap_all_exon_regions <- queryHits(findOverlaps(satmut_promoters_pred_gr, satmut_promoters_all_exon_regions))
	satmut_promoters_pred <- satmut_promoters_pred %>% 
		mutate(in_all_exon = (row_number() %in% overlap_all_exon_regions))

	overlap_all_exon_splice_regions <- queryHits(findOverlaps(satmut_promoters_pred_gr, satmut_promoters_all_exon_splice_regions))
	satmut_promoters_pred <- satmut_promoters_pred %>% 
		mutate(in_all_exon_splice = (row_number() %in% overlap_all_exon_splice_regions))

	# rearrange
	satmut_promoters_pred <- satmut_promoters_pred %>% 
		dplyr::select(
			chrom, pos, ref, alt, id, strand, tss_dist, starts_with("in_"), 
			starts_with("K562_"), starts_with("HepG2_"), starts_with("SKNSH_"), 
			starts_with("avgKHS_"), everything()
		)

	# save predictions
	write_tsv(satmut_promoters_pred, gzfile(paste0("../../data/satmut_promoters_preprocess/satmut_promoters_pred/satmut_promoters_pred_", chr, ".txt.gz")))
}
