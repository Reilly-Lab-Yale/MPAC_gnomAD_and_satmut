#!/bin/R

# Created helper BED files for VCF region subsetting

# load libraries
library(data.table)
library(tidyverse)
library(stringi)

# load bioconductor
library(plyranges)
library(rtracklayer)

# load promoter/exon data
# satmut_promoters_promoter_regions_1kb <- rtracklayer::import("../../data/gencode_filtered_regions/gencode.v44.protein.coding.canonical.autosomes.0.based.bed", extraCols = c(id = "character"))
satmut_promoters_promoter_regions_1kb <- rtracklayer::import("../../data/gencode_filtered_regions/gencode.v44.protein.coding.1kb.promoters.autosomes.v2.bed", extraCols = c(id = "character"))

# save regions
for (chr in paste0("chr", c(22:1))) {
	# save with "chr"
	satmut_promoters_promoter_regions_1kb_temp <- filter(satmut_promoters_promoter_regions_1kb, seqnames==chr)
	rtracklayer::export(satmut_promoters_promoter_regions_1kb_temp, paste0("../../data/satmut_promoters_preprocess/satmut_promoters_regions/satmut_promoters_promoter_regions_", chr, ".bed"), "bed")

	# save without "chr"
	seqlevels(satmut_promoters_promoter_regions_1kb_temp) <- gsub("chr", "", seqlevels(satmut_promoters_promoter_regions_1kb_temp))
	rtracklayer::export(satmut_promoters_promoter_regions_1kb_temp, paste0("../../data/satmut_promoters_preprocess/satmut_promoters_regions/satmut_promoters_promoter_regions_", gsub("chr", "", chr), ".bed"), "bed")
}

# helper functions
extend_granges <- function(x, upstream=0, downstream=0) {
	# based on https://support.bioconductor.org/p/78652/
	if (any(strand(x) == "*")) {
		warning("'*' ranges were treated as '+'")
	}
	on_plus <- strand(x) == "+" | strand(x) == "*"
	new_start <- start(x) - ifelse(on_plus, upstream, downstream)
	new_end <- end(x) + ifelse(on_plus, downstream, upstream)
	ranges(x) <- IRanges(new_start, new_end)
	return(trim(x))
}

# pull subregions
satmut_promoters_promoter_regions_1kb <- satmut_promoters_promoter_regions_1kb %>% extend_granges(0, 0)
satmut_promoters_promoter_regions_750bp <- satmut_promoters_promoter_regions_1kb %>% extend_granges(-250, 0)
satmut_promoters_promoter_regions_500bp <- satmut_promoters_promoter_regions_1kb %>% extend_granges(-500, 0)
satmut_promoters_promoter_regions_250bp <- satmut_promoters_promoter_regions_1kb %>% extend_granges(-750, 0)
satmut_promoters_promoter_regions_tss <- satmut_promoters_promoter_regions_1kb %>% extend_granges(-1000, 0)

saveRDS(satmut_promoters_promoter_regions_1kb, "../../data/satmut_promoters_preprocess/satmut_promoters_regions/satmut_promoters_promoter_regions_1kb.rds")
saveRDS(satmut_promoters_promoter_regions_750bp, "../../data/satmut_promoters_preprocess/satmut_promoters_regions/satmut_promoters_promoter_regions_750bp.rds")
saveRDS(satmut_promoters_promoter_regions_500bp, "../../data/satmut_promoters_preprocess/satmut_promoters_regions/satmut_promoters_promoter_regions_500bp.rds")
saveRDS(satmut_promoters_promoter_regions_250bp, "../../data/satmut_promoters_preprocess/satmut_promoters_regions/satmut_promoters_promoter_regions_250bp.rds")
saveRDS(satmut_promoters_promoter_regions_tss, "../../data/satmut_promoters_preprocess/satmut_promoters_regions/satmut_promoters_promoter_regions_tss.rds")
