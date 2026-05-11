#!/bin/R

# Get base level and promoter level predictions

library(tidyverse)
library(data.table)

# base-level summaries
for (chr in paste0("chr", c(22:1))) {
	print(chr)
	
	# load variant data
	satmut_promoters_final <- as_tibble(fread(paste0("../../data/satmut_promoters_preprocess/satmut_promoters_final/satmut_promoters_final_", chr, ".txt.gz")))

	# summ base-level data
	satmut_promoters_final_base <- satmut_promoters_final %>% 
		group_by(
			chrom, pos, id, strand, tss_dist, 
			in_promoter_1kb, in_promoter_750bp, in_promoter_500bp, in_promoter_250bp, 
			in_promoter_quartile, in_all_exon, in_all_exon_splice  # , in_TFfoot
		) %>% 
		summarise(
			# basic summaries
			K562_ref_base = mean(K562_ref, na.rm=T), 
			HepG2_ref_base = mean(HepG2_ref, na.rm=T), 
			SKNSH_ref_base = mean(SKNSH_ref, na.rm=T), 
			avgKHS_ref_base = mean(avgKHS_ref, na.rm=T), 

			K562_alt_base = mean(K562_alt, na.rm=T), 
			HepG2_alt_base = mean(HepG2_alt, na.rm=T), 
			SKNSH_alt_base = mean(SKNSH_alt, na.rm=T), 
			avgKHS_alt_base = mean(avgKHS_alt, na.rm=T), 
		
			K562_skew_base = mean(K562_skew, na.rm=T), 
			HepG2_skew_base = mean(HepG2_skew, na.rm=T), 
			SKNSH_skew_base = mean(SKNSH_skew, na.rm=T), 
			avgKHS_skew_base = mean(avgKHS_skew, na.rm=T), 

			K562_activity_base = mean(K562_activity, na.rm=T), 
			HepG2_activity_base = mean(HepG2_activity, na.rm=T), 
			SKNSH_activity_base = mean(SKNSH_activity, na.rm=T), 
			avgKHS_activity_base = mean(avgKHS_activity, na.rm=T), 
		
			K562_abs_skew_base = mean(K562_abs_skew, na.rm=T), 
			HepG2_abs_skew_base = mean(HepG2_abs_skew, na.rm=T), 
			SKNSH_abs_skew_base = mean(SKNSH_abs_skew, na.rm=T), 
			avgKHS_abs_skew_base = mean(avgKHS_abs_skew, na.rm=T), 

			# how many emVars are at this base
			K562_emVar_base = sum(K562_emVar, na.rm=T), 
			HepG2_emVar_base = sum(HepG2_emVar, na.rm=T), 
			SKNSH_emVar_base = sum(SKNSH_emVar, na.rm=T), 
			avgKHS_emVar_base = sum(avgKHS_emVar, na.rm=T), 
			anyKHS_emVar_base = sum(anyKHS_emVar, na.rm=T), 

			# how many pos emVars are at this base
			K562_emVar_pos_base = sum(K562_emVar_pos, na.rm=T), 
			HepG2_emVar_pos_base = sum(HepG2_emVar_pos, na.rm=T), 
			SKNSH_emVar_pos_base = sum(SKNSH_emVar_pos, na.rm=T), 
			avgKHS_emVar_pos_base = sum(avgKHS_emVar_pos, na.rm=T), 
			anyKHS_emVar_pos_base = sum(anyKHS_emVar_pos, na.rm=T), 

			# how many neg emVars are at this base
			K562_emVar_neg_base = sum(K562_emVar_neg, na.rm=T), 
			HepG2_emVar_neg_base = sum(HepG2_emVar_neg, na.rm=T), 
			SKNSH_emVar_neg_base = sum(SKNSH_emVar_neg, na.rm=T), 
			avgKHS_emVar_neg_base = sum(avgKHS_emVar_neg, na.rm=T), 
			anyKHS_emVar_neg_base = sum(anyKHS_emVar_neg, na.rm=T), 

			# how much conservation
			phyloP_mam241_base = mean(phyloP_mam241, na.rm=T), 
			phyloP_mam241_pos_base = mean(phyloP_mam241_pos, na.rm=T), 
			phyloP_mam241_neg_base = mean(phyloP_mam241_neg, na.rm=T), 
			phyloP_mam241_cons_base = mean(phyloP_mam241_cons, na.rm=T)
		) %>% 
		# additional features used for signed correlations
		mutate(
			sign_K562_skew_base = ifelse(K562_skew_base > 0, "pos", "neg"),
			sign_HepG2_skew_base = ifelse(HepG2_skew_base > 0, "pos", "neg"),
			sign_SKNSH_skew_base = ifelse(SKNSH_skew_base > 0, "pos", "neg"),
			sign_avgKHS_skew_base = ifelse(avgKHS_skew_base > 0, "pos", "neg")
		) %>% 
		mutate(
			K562_pos_skew_base = ifelse(sign_K562_skew_base == "pos", K562_skew_base, NA), 
			HepG2_pos_skew_base = ifelse(sign_HepG2_skew_base == "pos", HepG2_skew_base, NA), 
			SKNSH_pos_skew_base = ifelse(sign_SKNSH_skew_base == "pos", SKNSH_skew_base, NA), 
			avgKHS_pos_skew_base = ifelse(sign_avgKHS_skew_base == "pos", avgKHS_skew_base, NA), 

			K562_neg_skew_base = ifelse(sign_K562_skew_base == "neg", K562_skew_base, NA), 
			HepG2_neg_skew_base = ifelse(sign_HepG2_skew_base == "neg", HepG2_skew_base, NA), 
			SKNSH_neg_skew_base = ifelse(sign_SKNSH_skew_base == "neg", SKNSH_skew_base, NA), 
			avgKHS_neg_skew_base = ifelse(sign_avgKHS_skew_base == "neg", avgKHS_skew_base, NA), 
		) %>% 
		ungroup() %>% 
		dplyr::select(-sign_K562_skew_base, -sign_HepG2_skew_base, -sign_SKNSH_skew_base, -sign_avgKHS_skew_base)

	# save all base-level
	write_tsv(satmut_promoters_final_base, gzfile(paste0("../../data/satmut_promoters_preprocess/satmut_promoters_final_base/satmut_promoters_final_base_", chr, ".txt.gz")))
}
