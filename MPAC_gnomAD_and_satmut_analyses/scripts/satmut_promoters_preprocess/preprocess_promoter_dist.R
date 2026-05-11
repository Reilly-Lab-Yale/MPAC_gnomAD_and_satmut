#!/bin/R

# Get distance level predictions

library(tidyverse)
library(data.table)

# dist-level data
satmut_promoters_final_dist_list <- list()
for (chr in paste0("chr", c(22:1))) {
	print(chr)
	satmut_promoters_final_temp <- as_tibble(fread(paste0("../../data/satmut_promoters_preprocess/satmut_promoters_final/satmut_promoters_final_", chr, ".txt.gz")))
	satmut_promoters_final_dist_temp <- satmut_promoters_final_temp %>% 
		filter(!in_all_exon_splice) %>% 
		group_by(
			tss_dist, 
			in_promoter_1kb, 
			in_promoter_750bp, 
			in_promoter_500bp, 
			in_promoter_250bp,
			) %>% 
			summarise(
				n_in_tss_dist = n(),

				K562_ref = sum(K562_ref, na.rm=T), 
				HepG2_ref = sum(HepG2_ref, na.rm=T), 
				SKNSH_ref = sum(SKNSH_ref, na.rm=T), 
				avgKHS_ref = sum(avgKHS_ref, na.rm=T), 

				K562_alt = sum(K562_alt, na.rm=T), 
				HepG2_alt = sum(HepG2_alt, na.rm=T), 
				SKNSH_alt = sum(SKNSH_alt, na.rm=T), 
				avgKHS_alt = sum(avgKHS_alt, na.rm=T), 
			
				K562_skew = sum(K562_skew, na.rm=T), 
				HepG2_skew = sum(HepG2_skew, na.rm=T), 
				SKNSH_skew = sum(SKNSH_skew, na.rm=T), 
				avgKHS_skew = sum(avgKHS_skew, na.rm=T), 

				K562_activity = sum(K562_activity, na.rm=T), 
				HepG2_activity = sum(HepG2_activity, na.rm=T), 
				SKNSH_activity = sum(SKNSH_activity, na.rm=T), 
				avgKHS_activity = sum(avgKHS_activity, na.rm=T), 
			
				K562_abs_skew = sum(K562_abs_skew, na.rm=T), 
				HepG2_abs_skew = sum(HepG2_abs_skew, na.rm=T), 
				SKNSH_abs_skew = sum(SKNSH_abs_skew, na.rm=T), 
				avgKHS_abs_skew = sum(avgKHS_abs_skew, na.rm=T), 

				phyloP_mam241 = sum(phyloP_mam241, na.rm=T)
			) %>% 
			ungroup()

	satmut_promoters_final_dist_list[[chr]] <- satmut_promoters_final_dist_temp
}

satmut_promoters_final_dist_all <- satmut_promoters_final_dist_list %>% bind_rows() %>% 
	group_by(
		tss_dist, 
		in_promoter_1kb, 
		in_promoter_750bp, 
		in_promoter_500bp, 
		in_promoter_250bp
		) %>% 
		summarise(
			n_in_tss_dist = sum(n_in_tss_dist, na.rm=T),

			K562_ref_dist = sum(K562_ref, na.rm=T), 
			HepG2_ref_dist = sum(HepG2_ref, na.rm=T), 
			SKNSH_ref_dist = sum(SKNSH_ref, na.rm=T), 
			avgKHS_ref_dist = sum(avgKHS_ref, na.rm=T), 

			K562_alt_dist = sum(K562_alt, na.rm=T), 
			HepG2_alt_dist = sum(HepG2_alt, na.rm=T), 
			SKNSH_alt_dist = sum(SKNSH_alt, na.rm=T), 
			avgKHS_alt_dist = sum(avgKHS_alt, na.rm=T), 
		
			K562_skew_dist = sum(K562_skew, na.rm=T), 
			HepG2_skew_dist = sum(HepG2_skew, na.rm=T), 
			SKNSH_skew_dist = sum(SKNSH_skew, na.rm=T), 
			avgKHS_skew_dist = sum(avgKHS_skew, na.rm=T), 

			K562_activity_dist = sum(K562_activity, na.rm=T), 
			HepG2_activity_dist = sum(HepG2_activity, na.rm=T), 
			SKNSH_activity_dist = sum(SKNSH_activity, na.rm=T), 
			avgKHS_activity_dist = sum(avgKHS_activity, na.rm=T), 
		
			K562_abs_skew_dist = sum(K562_abs_skew, na.rm=T), 
			HepG2_abs_skew_dist = sum(HepG2_abs_skew, na.rm=T), 
			SKNSH_abs_skew_dist = sum(SKNSH_abs_skew, na.rm=T), 
			avgKHS_abs_skew_dist = sum(avgKHS_abs_skew, na.rm=T), 

			phyloP_mam241_dist = sum(phyloP_mam241, na.rm=T)
		) %>% 
		mutate(
			K562_ref_dist = K562_ref_dist/n_in_tss_dist,
			HepG2_ref_dist = HepG2_ref_dist/n_in_tss_dist,
			SKNSH_ref_dist = SKNSH_ref_dist/n_in_tss_dist,
			avgKHS_ref_dist = avgKHS_ref_dist/n_in_tss_dist,

			K562_alt_dist = K562_alt_dist/n_in_tss_dist,
			HepG2_alt_dist = HepG2_alt_dist/n_in_tss_dist,
			SKNSH_alt_dist = SKNSH_alt_dist/n_in_tss_dist,
			avgKHS_alt_dist = avgKHS_alt_dist/n_in_tss_dist,
		
			K562_skew_dist = K562_skew_dist/n_in_tss_dist,
			HepG2_skew_dist = HepG2_skew_dist/n_in_tss_dist,
			SKNSH_skew_dist = SKNSH_skew_dist/n_in_tss_dist,
			avgKHS_skew_dist = avgKHS_skew_dist/n_in_tss_dist,

			K562_activity_dist = K562_activity_dist/n_in_tss_dist,
			HepG2_activity_dist = HepG2_activity_dist/n_in_tss_dist,
			SKNSH_activity_dist = SKNSH_activity_dist/n_in_tss_dist,
			avgKHS_activity_dist = avgKHS_activity_dist/n_in_tss_dist,
		
			K562_abs_skew_dist = K562_abs_skew_dist/n_in_tss_dist,
			HepG2_abs_skew_dist = HepG2_abs_skew_dist/n_in_tss_dist,
			SKNSH_abs_skew_dist = SKNSH_abs_skew_dist/n_in_tss_dist,
			avgKHS_abs_skew_dist = avgKHS_abs_skew_dist/n_in_tss_dist,

			phyloP_mam241_dist = phyloP_mam241_dist/n_in_tss_dist
		) %>% 
		ungroup()	

write_tsv(satmut_promoters_final_dist_all, gzfile(paste0("../../data/satmut_promoters_preprocess/satmut_promoters_final_dist/satmut_promoters_final_dist_all.txt.gz")))
