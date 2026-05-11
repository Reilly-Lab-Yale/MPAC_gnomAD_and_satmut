#!/bin/R

# Get promoter level predictions

library(tidyverse)
library(data.table)

set.seed(123456789)

# promoter-level data
for (chr in paste0("chr", c(22:1))) {
	print(chr)

	# load base-level data
	satmut_promoters_final_base <- as_tibble(fread(paste0("../../data/satmut_promoters_preprocess/satmut_promoters_final_base/satmut_promoters_final_base_", chr, ".txt.gz")))

	# summ promoter-level data
	prom_level_shared <- function(x) {
		x %>% 
		group_by(chrom, id, strand) %>% 
		summarise(
			# what is the number of variants of each promoter
			n_in_prom = dplyr::n(),

			# what is the average ref pred of each promoter in each cell line
			K562_ref_prom = mean(K562_ref_base, na.rm=T), 
			HepG2_ref_prom = mean(HepG2_ref_base, na.rm=T), 
			SKNSH_ref_prom = mean(SKNSH_ref_base, na.rm=T), 
			avgKHS_ref_prom = mean(avgKHS_ref_base, na.rm=T), 

			# what is the average alt pred of each promoter in each cell line
			K562_alt_prom = mean(K562_alt_base, na.rm=T), 
			HepG2_alt_prom = mean(HepG2_alt_base, na.rm=T), 
			SKNSH_alt_prom = mean(SKNSH_alt_base, na.rm=T), 
			avgKHS_alt_prom = mean(avgKHS_alt_base, na.rm=T), 

			# what is the average skew skew of each promoter in each cell line
			K562_skew_prom = mean(K562_skew_base, na.rm=T), 
			HepG2_skew_prom = mean(HepG2_skew_base, na.rm=T), 
			SKNSH_skew_prom = mean(SKNSH_skew_base, na.rm=T), 
			avgKHS_skew_prom = mean(avgKHS_skew_base, na.rm=T), 

			# what is the average activity of each promoter in each cell line
			K562_activity_prom = mean(K562_activity_base, na.rm=T), 
			HepG2_activity_prom = mean(HepG2_activity_base, na.rm=T), 
			SKNSH_activity_prom = mean(SKNSH_activity_base, na.rm=T), 
			avgKHS_activity_prom = mean(avgKHS_activity_base, na.rm=T), 

			# what is the average abs skew of each promoter in each cell line
			K562_abs_skew_prom = mean(K562_abs_skew_base, na.rm=T), 
			HepG2_abs_skew_prom = mean(HepG2_abs_skew_base, na.rm=T), 
			SKNSH_abs_skew_prom = mean(SKNSH_abs_skew_base, na.rm=T), 
			avgKHS_abs_skew_prom = mean(avgKHS_abs_skew_base, na.rm=T), 

			# what is the number of emVar of each promoter in each cell line
			K562_emVar_prom = sum(K562_emVar_base, na.rm=T), 
			HepG2_emVar_prom = sum(HepG2_emVar_base, na.rm=T), 
			SKNSH_emVar_prom = sum(SKNSH_emVar_base, na.rm=T), 
			avgKHS_emVar_prom = sum(avgKHS_emVar_base, na.rm=T), 
			anyKHS_emVar_prom = sum(anyKHS_emVar_base, na.rm=T), 

			# what is the number of pos emVar of each promoter in each cell line
			K562_emVar_pos_prom = sum(K562_emVar_pos_base, na.rm=T), 
			HepG2_emVar_pos_prom = sum(HepG2_emVar_pos_base, na.rm=T), 
			SKNSH_emVar_pos_prom = sum(SKNSH_emVar_pos_base, na.rm=T), 
			avgKHS_emVar_pos_prom = sum(avgKHS_emVar_pos_base, na.rm=T), 
			anyKHS_emVar_pos_prom = sum(anyKHS_emVar_pos_base, na.rm=T), 

			# what is the number of pos emVar of each promoter in each cell line
			K562_emVar_neg_prom = sum(K562_emVar_neg_base, na.rm=T), 
			HepG2_emVar_neg_prom = sum(HepG2_emVar_neg_base, na.rm=T), 
			SKNSH_emVar_neg_prom = sum(SKNSH_emVar_neg_base, na.rm=T), 
			avgKHS_emVar_neg_prom = sum(avgKHS_emVar_neg_base, na.rm=T), 
			anyKHS_emVar_neg_prom = sum(anyKHS_emVar_neg_base, na.rm=T), 

			# what is the average constraint of each promoter
			# 	only use mam241 for downstream analyses
			phyloP_mam241_prom = mean(phyloP_mam241_base, na.rm=T), 
			phyloP_mam241_pos_prom = mean(phyloP_mam241_pos_base, na.rm=T), 
			phyloP_mam241_neg_prom = mean(phyloP_mam241_neg_base, na.rm=T), 
			phyloP_mam241_cons_prom = mean(phyloP_mam241_cons_base, na.rm=T),

			# what is the correlation between phyloP pos and abs skew
			# 	flag is to avoid na's and regions where phyloP is exactly zero
			flag_correlation_mam241_vs_K562_n_abs = sum(!((is.na(phyloP_mam241_base) | (phyloP_mam241_base == 0)) | is.na(K562_abs_skew_base))), 
			correlation_mam241_vs_K562_rho_est_abs = ifelse(flag_correlation_mam241_vs_K562_n_abs >= 50, cor.test(phyloP_mam241_pos_base, K562_abs_skew_base, na.action = na.omit, method = "spearman")$estimate, NA), 
			correlation_mam241_vs_K562_rho_pval_abs = ifelse(flag_correlation_mam241_vs_K562_n_abs >= 50, cor.test(phyloP_mam241_pos_base, K562_abs_skew_base, na.action = na.omit, method = "spearman")$p.value, NA), 

			flag_correlation_mam241_vs_HepG2_n_abs = sum(!((is.na(phyloP_mam241_base) | (phyloP_mam241_base == 0)) | is.na(HepG2_abs_skew_base))), 
			correlation_mam241_vs_HepG2_rho_est_abs = ifelse(flag_correlation_mam241_vs_HepG2_n_abs >= 50, cor.test(phyloP_mam241_pos_base, HepG2_abs_skew_base, na.action = na.omit, method = "spearman")$estimate, NA), 
			correlation_mam241_vs_HepG2_rho_pval_abs = ifelse(flag_correlation_mam241_vs_HepG2_n_abs >= 50, cor.test(phyloP_mam241_pos_base, HepG2_abs_skew_base, na.action = na.omit, method = "spearman")$p.value, NA), 

			flag_correlation_mam241_vs_SKNSH_n_abs = sum(!((is.na(phyloP_mam241_base) | (phyloP_mam241_base == 0)) | is.na(SKNSH_abs_skew_base))), 
			correlation_mam241_vs_SKNSH_rho_est_abs = ifelse(flag_correlation_mam241_vs_SKNSH_n_abs >= 50, cor.test(phyloP_mam241_pos_base, SKNSH_abs_skew_base, na.action = na.omit, method = "spearman")$estimate, NA), 
			correlation_mam241_vs_SKNSH_rho_pval_abs = ifelse(flag_correlation_mam241_vs_SKNSH_n_abs >= 50, cor.test(phyloP_mam241_pos_base, SKNSH_abs_skew_base, na.action = na.omit, method = "spearman")$p.value, NA), 

			flag_correlation_mam241_vs_avgKHS_n_abs = sum(!((is.na(phyloP_mam241_base) | (phyloP_mam241_base == 0)) | is.na(avgKHS_abs_skew_base))), 
			correlation_mam241_vs_avgKHS_rho_est_abs = ifelse(flag_correlation_mam241_vs_avgKHS_n_abs >= 50, cor.test(phyloP_mam241_pos_base, avgKHS_abs_skew_base, na.action = na.omit, method = "spearman")$estimate, NA), 
			correlation_mam241_vs_avgKHS_rho_pval_abs = ifelse(flag_correlation_mam241_vs_avgKHS_n_abs >= 50, cor.test(phyloP_mam241_pos_base, avgKHS_abs_skew_base, na.action = na.omit, method = "spearman")$p.value, NA), 

			# what is the correlation between phyloP pos and pos skew only
			flag_correlation_mam241_vs_K562_n_pos = sum(!((is.na(phyloP_mam241_base) | (phyloP_mam241_base == 0)) | is.na(K562_pos_skew_base))), 
			correlation_mam241_vs_K562_rho_est_pos = ifelse(flag_correlation_mam241_vs_K562_n_pos >= 50, cor.test(phyloP_mam241_pos_base, K562_pos_skew_base, na.action = na.omit, method = "spearman")$estimate, NA), 
			correlation_mam241_vs_K562_rho_pval_pos = ifelse(flag_correlation_mam241_vs_K562_n_pos >= 50, cor.test(phyloP_mam241_pos_base, K562_pos_skew_base, na.action = na.omit, method = "spearman")$p.value, NA), 

			flag_correlation_mam241_vs_HepG2_n_pos = sum(!((is.na(phyloP_mam241_base) | (phyloP_mam241_base == 0)) | is.na(HepG2_pos_skew_base))), 
			correlation_mam241_vs_HepG2_rho_est_pos = ifelse(flag_correlation_mam241_vs_HepG2_n_pos >= 50, cor.test(phyloP_mam241_pos_base, HepG2_pos_skew_base, na.action = na.omit, method = "spearman")$estimate, NA), 
			correlation_mam241_vs_HepG2_rho_pval_pos = ifelse(flag_correlation_mam241_vs_HepG2_n_pos >= 50, cor.test(phyloP_mam241_pos_base, HepG2_pos_skew_base, na.action = na.omit, method = "spearman")$p.value, NA), 

			flag_correlation_mam241_vs_SKNSH_n_pos = sum(!((is.na(phyloP_mam241_base) | (phyloP_mam241_base == 0)) | is.na(SKNSH_pos_skew_base))), 
			correlation_mam241_vs_SKNSH_rho_est_pos = ifelse(flag_correlation_mam241_vs_SKNSH_n_pos >= 50, cor.test(phyloP_mam241_pos_base, SKNSH_pos_skew_base, na.action = na.omit, method = "spearman")$estimate, NA), 
			correlation_mam241_vs_SKNSH_rho_pval_pos = ifelse(flag_correlation_mam241_vs_SKNSH_n_pos >= 50, cor.test(phyloP_mam241_pos_base, SKNSH_pos_skew_base, na.action = na.omit, method = "spearman")$p.value, NA), 

			flag_correlation_mam241_vs_avgKHS_n_pos = sum(!((is.na(phyloP_mam241_base) | (phyloP_mam241_base == 0)) | is.na(avgKHS_pos_skew_base))), 
			correlation_mam241_vs_avgKHS_rho_est_pos = ifelse(flag_correlation_mam241_vs_avgKHS_n_pos >= 50, cor.test(phyloP_mam241_pos_base, avgKHS_pos_skew_base, na.action = na.omit, method = "spearman")$estimate, NA), 
			correlation_mam241_vs_avgKHS_rho_pval_pos = ifelse(flag_correlation_mam241_vs_avgKHS_n_pos >= 50, cor.test(phyloP_mam241_pos_base, avgKHS_pos_skew_base, na.action = na.omit, method = "spearman")$p.value, NA), 

			# what is the correlation between phyloP pos and neg skew only
			flag_correlation_mam241_vs_K562_n_neg = sum(!((is.na(phyloP_mam241_base) | (phyloP_mam241_base == 0)) | is.na(K562_neg_skew_base))), 
			correlation_mam241_vs_K562_rho_est_neg = ifelse(flag_correlation_mam241_vs_K562_n_neg >= 50, cor.test(phyloP_mam241_pos_base, K562_neg_skew_base, na.action = na.omit, method = "spearman")$estimate, NA), 
			correlation_mam241_vs_K562_rho_pval_neg = ifelse(flag_correlation_mam241_vs_K562_n_neg >= 50, cor.test(phyloP_mam241_pos_base, K562_neg_skew_base, na.action = na.omit, method = "spearman")$p.value, NA), 

			flag_correlation_mam241_vs_HepG2_n_neg = sum(!((is.na(phyloP_mam241_base) | (phyloP_mam241_base == 0)) | is.na(HepG2_neg_skew_base))), 
			correlation_mam241_vs_HepG2_rho_est_neg = ifelse(flag_correlation_mam241_vs_HepG2_n_neg >= 50, cor.test(phyloP_mam241_pos_base, HepG2_neg_skew_base, na.action = na.omit, method = "spearman")$estimate, NA), 
			correlation_mam241_vs_HepG2_rho_pval_neg = ifelse(flag_correlation_mam241_vs_HepG2_n_neg >= 50, cor.test(phyloP_mam241_pos_base, HepG2_neg_skew_base, na.action = na.omit, method = "spearman")$p.value, NA), 

			flag_correlation_mam241_vs_SKNSH_n_neg = sum(!((is.na(phyloP_mam241_base) | (phyloP_mam241_base == 0)) | is.na(SKNSH_neg_skew_base))), 
			correlation_mam241_vs_SKNSH_rho_est_neg = ifelse(flag_correlation_mam241_vs_SKNSH_n_neg >= 50, cor.test(phyloP_mam241_pos_base, SKNSH_neg_skew_base, na.action = na.omit, method = "spearman")$estimate, NA), 
			correlation_mam241_vs_SKNSH_rho_pval_neg = ifelse(flag_correlation_mam241_vs_SKNSH_n_neg >= 50, cor.test(phyloP_mam241_pos_base, SKNSH_neg_skew_base, na.action = na.omit, method = "spearman")$p.value, NA), 

			flag_correlation_mam241_vs_avgKHS_n_neg = sum(!((is.na(phyloP_mam241_base) | (phyloP_mam241_base == 0)) | is.na(avgKHS_neg_skew_base))), 
			correlation_mam241_vs_avgKHS_rho_est_neg = ifelse(flag_correlation_mam241_vs_avgKHS_n_neg >= 50, cor.test(phyloP_mam241_pos_base, avgKHS_neg_skew_base, na.action = na.omit, method = "spearman")$estimate, NA), 
			correlation_mam241_vs_avgKHS_rho_pval_neg = ifelse(flag_correlation_mam241_vs_avgKHS_n_neg >= 50, cor.test(phyloP_mam241_pos_base, avgKHS_neg_skew_base, na.action = na.omit, method = "spearman")$p.value, NA), 
		) %>% 
		ungroup() %>% 
		mutate(
			# what is the tissue-specifity of each promoter in each cell line
			# 	note, highly correlated between mingap of mean and mean of mingap
			K562_activity_mingap_prom_alt = K562_activity_prom - pmax(HepG2_activity_prom, SKNSH_activity_prom), 
			HepG2_activity_mingap_prom_alt = HepG2_activity_prom - pmax(K562_activity_prom, SKNSH_activity_prom), 
			SKNSH_activity_mingap_prom_alt = SKNSH_activity_prom - pmax(K562_activity_prom, HepG2_activity_prom), 

			K562_abs_skew_mingap_prom_alt = K562_abs_skew_prom - pmax(HepG2_abs_skew_prom, SKNSH_abs_skew_prom), 
			HepG2_abs_skew_mingap_prom_alt = HepG2_abs_skew_prom - pmax(K562_abs_skew_prom, SKNSH_abs_skew_prom), 
			SKNSH_abs_skew_mingap_prom_alt = SKNSH_abs_skew_prom - pmax(K562_abs_skew_prom, HepG2_abs_skew_prom), 
		)
	}

	shuffle_shared <- function(x) {
		x %>% 
			group_by(chrom, id, strand, 
				# in_promoter_250bp  # keep for convenience
			) %>% 
			mutate(
				phyloP_mam241_pos_base = sample(phyloP_mam241_pos_base)
			) %>% 
			ungroup()
	}

	# summ promoter-level data, masked version
	satmut_promoters_final_prom_250bp_level <- satmut_promoters_final_base %>% 
		filter(in_promoter_250bp, !in_all_exon_splice) %>% 
		prom_level_shared() %>% 
		mutate(analysis_version = "promoter_250bp")

	# summ promoter-level data, masked version
	satmut_promoters_final_prom_500bp_level <- satmut_promoters_final_base %>% 
		filter(in_promoter_500bp, !in_all_exon_splice) %>% 
		prom_level_shared() %>% 
		mutate(analysis_version = "promoter_500bp")

	# summ promoter-level data, masked version
	satmut_promoters_final_prom_750bp_level <- satmut_promoters_final_base %>% 
		filter(in_promoter_750bp, !in_all_exon_splice) %>% 
		prom_level_shared() %>% 
		mutate(analysis_version = "promoter_750bp")

	# summ promoter-level data, masked version
	satmut_promoters_final_prom_1kb_level <- satmut_promoters_final_base %>% 
		filter(in_promoter_1kb, !in_all_exon_splice) %>% 
		prom_level_shared() %>% 
		mutate(analysis_version = "promoter_1kb")

	# summ promoter-level data, masked shuffled version
	satmut_promoters_final_prom_250bp_level_shuffled <- satmut_promoters_final_base %>% 
		shuffle_shared() %>% 
		filter(in_promoter_250bp, !in_all_exon_splice) %>% 
		prom_level_shared() %>% 
		mutate(analysis_version = "shuffled_promoter_250bp")

	# summ promoter-level data, masked shuffled version
	satmut_promoters_final_prom_500bp_level_shuffled <- satmut_promoters_final_base %>% 
		shuffle_shared() %>% 
		filter(in_promoter_500bp, !in_all_exon_splice) %>% 
		prom_level_shared() %>% 
		mutate(analysis_version = "shuffled_promoter_500bp")

	# summ promoter-level data, masked shuffled version
	satmut_promoters_final_prom_750bp_level_shuffled <- satmut_promoters_final_base %>% 
		shuffle_shared() %>% 
		filter(in_promoter_750bp, !in_all_exon_splice) %>% 
		prom_level_shared() %>% 
		mutate(analysis_version = "shuffled_promoter_750bp")

	# summ promoter-level data, masked shuffled version
	satmut_promoters_final_prom_1kb_level_shuffled <- satmut_promoters_final_base %>% 
		shuffle_shared() %>% 
		filter(in_promoter_1kb, !in_all_exon_splice) %>% 
		prom_level_shared() %>% 
		mutate(analysis_version = "shuffled_promoter_1kb")

	# bind them all
	satmut_promoters_final_prom <- 
		satmut_promoters_final_prom_250bp_level %>% 
		bind_rows(satmut_promoters_final_prom_500bp_level) %>% 
		bind_rows(satmut_promoters_final_prom_750bp_level) %>% 
		bind_rows(satmut_promoters_final_prom_1kb_level) %>% 
		bind_rows(satmut_promoters_final_prom_250bp_level_shuffled) %>% 
		bind_rows(satmut_promoters_final_prom_500bp_level_shuffled) %>% 
		bind_rows(satmut_promoters_final_prom_750bp_level_shuffled) %>% 
		bind_rows(satmut_promoters_final_prom_1kb_level_shuffled)

	# save all promoter-level
	write_tsv(satmut_promoters_final_prom, gzfile(paste0("../../data/satmut_promoters_preprocess/satmut_promoters_final_prom/satmut_promoters_final_prom_", chr, ".txt.gz")))
}

# bind promoter get all
satmut_promoters_final_prom_all <- list()
for (chr in paste0("chr", c(22:1))) {
	print(chr)
	satmut_promoters_final_prom_all[[chr]] <- as_tibble(fread(paste0("../../data/satmut_promoters_preprocess/satmut_promoters_final_prom/satmut_promoters_final_prom_", chr, ".txt.gz")))
}

satmut_promoters_final_prom_all <- satmut_promoters_final_prom_all %>% bind_rows()
write_tsv(satmut_promoters_final_prom_all, gzfile(paste0("../../data/satmut_promoters_preprocess/satmut_promoters_final_prom/satmut_promoters_final_prom_all.txt.gz")))

# filter for supplement
satmut_promoters_final_prom_all <- as_tibble(fread("../../data/satmut_promoters_preprocess/satmut_promoters_final_prom/satmut_promoters_final_prom_all.txt.gz"))

satmut_promoters_final_prom_supp <- satmut_promoters_final_prom_all %>% 
	filter(analysis_version %in% c("promoter_250bp", "shuffled_promoter_250bp")) %>% 
	mutate(chrom = as.integer(gsub("chr", "", chrom))) %>% 
	arrange(chrom) %>% 
	mutate(chrom = paste0("chr", chrom)) %>% 
	dplyr::select(id, n_in_prom, avgKHS_ref_prom, avgKHS_alt_prom, avgKHS_skew_prom, avgKHS_activity_prom, avgKHS_abs_skew_prom, avgKHS_emVar_prom, avgKHS_emVar_pos_prom, avgKHS_emVar_neg_prom, phyloP_mam241_prom, phyloP_mam241_pos_prom, phyloP_mam241_neg_prom, phyloP_mam241_cons_prom, flag_correlation_mam241_vs_avgKHS_n_abs, correlation_mam241_vs_avgKHS_rho_est_abs, correlation_mam241_vs_avgKHS_rho_pval_abs, flag_correlation_mam241_vs_avgKHS_n_pos, correlation_mam241_vs_avgKHS_rho_est_pos, correlation_mam241_vs_avgKHS_rho_pval_pos, flag_correlation_mam241_vs_avgKHS_n_neg, correlation_mam241_vs_avgKHS_rho_est_neg, correlation_mam241_vs_avgKHS_rho_pval_neg, analysis_version)

write_tsv(satmut_promoters_final_prom_supp, gzfile(paste0("../../data/satmut_promoters_preprocess/satmut_promoters_final_prom/satmut_promoters_final_prom_supp.txt.gz")))
