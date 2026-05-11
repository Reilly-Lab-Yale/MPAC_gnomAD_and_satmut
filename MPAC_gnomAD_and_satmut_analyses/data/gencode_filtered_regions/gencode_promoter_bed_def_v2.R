#!/usr/bin/env Rscript

# load libraries
library(tidyverse)
library(data.table)
library(rtracklayer)
library(plyranges)

# helper functions
extend_granges <- function(x, upstream = 0, downstream = 0) {
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

write_bed7 <- function(gr, file) {
	df <- data.frame(
		chrom = as.character(seqnames(gr)),
		start = start(gr) - 1L,  # convert to 0-based
		end = end(gr),
		name = gr$gene_name,
		score = 0L,
		strand = as.character(strand(gr)),
		label = paste(gr$gene_id, gr$transcript_id, gr$gene_name, sep = "_")
	)
	fwrite(df, file, sep = "\t", col.names = FALSE)
}

# load gencode gff3
gff_full <- rtracklayer::import("gencode.v44.basic.annotation.gff3.gz")

# filter canonical protein-coding transcripts on autosomes per spec:
#   - type == "transcript"
#   - transcript_type == "protein_coding"
#   - chr1-22
#   - non-missing hgnc_id
#   - tag matches "Ensembl_canonical" but not "readthrough"
canonical_transcripts <- gff_full %>%
	filter(type == "transcript") %>%
	filter(transcript_type == "protein_coding") %>%
	filter(seqnames %in% paste0("chr", 1:22)) %>%
	as_tibble() %>%
	filter(!is.na(hgnc_id)) %>%
	filter(grepl("Ensembl_canonical", tag)) %>%
	filter(!grepl("readthrough", tag)) %>%
	GRanges()

# sanity check: no duplicate gene_names
dupes <- canonical_transcripts %>%
	as_tibble() %>%
	count(gene_name) %>%
	filter(n > 1)
if (nrow(dupes) > 0) {
	warning("Duplicate gene_names found: ", paste(dupes$gene_name, collapse = ", "))
	print(canonical_transcripts %>% filter(gene_name %in% dupes$gene_name))
	stop("Resolve duplicate gene_names before proceeding")
}

canonical_tx_ids <- canonical_transcripts$transcript_id

# exons from canonical transcripts only
canonical_exons <- gff_full %>%
	filter(type == "exon") %>%
	filter(transcript_id %in% canonical_tx_ids)

# all exons on autosomes
all_exons <- gff_full %>%
	filter(type == "exon") %>%
	filter(seqnames %in% paste0("chr", 1:22))

# save genes, canonical exons, all exons
write_bed7(canonical_transcripts, "gencode.v44.protein.coding.genes.autosomes.v2.bed")
write_bed7(canonical_exons, "gencode.v44.protein.coding.exons.autosomes.v2.bed")
write_bed7(all_exons, "gencode.v44.basic.annotation.exons.autosomes.v2.bed")

# get TSS and promoters
canonical_tss <- canonical_transcripts %>%
	anchor_5p() %>%
	mutate(width = 1)

canonical_promoters_1kb <- canonical_transcripts %>%
	flank_upstream(width = 1000)
canonical_promoters_750bp <- canonical_transcripts %>%
	flank_upstream(width = 750)
canonical_promoters_500bp <- canonical_transcripts %>%
	flank_upstream(width = 500)
canonical_promoters_250bp <- canonical_transcripts %>%
	flank_upstream(width = 250)

# save TSS and promoters
write_bed7(canonical_tss, "gencode.v44.protein.coding.tss.autosomes.v2.bed")
write_bed7(canonical_promoters_1kb, "gencode.v44.protein.coding.1kb.promoters.autosomes.v2.bed")
write_bed7(canonical_promoters_750bp, "gencode.v44.protein.coding.750bp.promoters.autosomes.v2.bed")
write_bed7(canonical_promoters_500bp, "gencode.v44.protein.coding.500bp.promoters.autosomes.v2.bed")
write_bed7(canonical_promoters_250bp, "gencode.v44.protein.coding.250bp.promoters.autosomes.v2.bed")

# identify first, internal, last exons for canonical protein-coding exons
canonical_exons_splice <- canonical_exons %>%
	as_tibble() %>%
	mutate(exon_number = as.integer(exon_number)) %>%
	group_by(transcript_name) %>%
	mutate(first = ((exon_number == 1) & (max(exon_number) > 1))) %>%
	mutate(internal = !(exon_number %in% c(1, max(exon_number)))) %>%
	mutate(last = ((exon_number == max(exon_number)) & (max(exon_number) > 1))) %>%
	ungroup()

# extend exons to splice regions (based on maxentscan)
acceptor <- 20
donor <- 6
canonical_exons_splice <- canonical_exons_splice %>%
	mutate(start = start - ifelse(strand == "+", 0, donor) * first, end = end + ifelse(strand == "+", donor, 0) * first) %>%
	mutate(start = start - ifelse(strand == "+", acceptor, donor) * internal, end = end + ifelse(strand == "+", donor, acceptor) * internal) %>%
	mutate(start = start - ifelse(strand == "+", acceptor, 0) * last, end = end + ifelse(strand == "+", 0, acceptor) * last) %>%
	GRanges()

# save canonical splice exons
write_bed7(canonical_exons_splice, "gencode.v44.protein.coding.exons.splice.autosomes.v2.bed")

# identify first, internal, last exons for all exons
all_exons_splice <- all_exons %>%
	as_tibble() %>%
	mutate(exon_number = as.integer(exon_number)) %>%
	group_by(transcript_name) %>%
	mutate(first = ((exon_number == 1) & (max(exon_number) > 1))) %>%
	mutate(internal = !(exon_number %in% c(1, max(exon_number)))) %>%
	mutate(last = ((exon_number == max(exon_number)) & (max(exon_number) > 1))) %>%
	ungroup()

# extend exons to splice regions (based on maxentscan)
all_exons_splice <- all_exons_splice %>%
	mutate(start = start - ifelse(strand == "+", 0, donor) * first, end = end + ifelse(strand == "+", donor, 0) * first) %>%
	mutate(start = start - ifelse(strand == "+", acceptor, donor) * internal, end = end + ifelse(strand == "+", donor, acceptor) * internal) %>%
	mutate(start = start - ifelse(strand == "+", acceptor, 0) * last, end = end + ifelse(strand == "+", 0, acceptor) * last) %>%
	GRanges()

# save all splice exons
write_bed7(all_exons_splice, "gencode.v44.basic.annotation.exons.splice.autosomes.v2.bed")
