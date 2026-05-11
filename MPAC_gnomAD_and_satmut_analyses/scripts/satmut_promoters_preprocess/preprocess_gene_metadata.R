#!/bin/R

# Get gene level constraint metrics and other annotations

# Load libraries
library(tidyverse)
library(data.table)
library(rtracklayer)

# Load gene names
# gencode_gene_names <- as_tibble(rtracklayer::import("../../data/gencode_filtered_regions/gencode.v44.protein.coding.canonical.autosomes.0.based.bed", extraCols = c(id = "character"))) %>% 
gencode_gene_names <- as_tibble(rtracklayer::import("../../data/gencode_filtered_regions/gencode.v44.protein.coding.1kb.promoters.autosomes.v2.bed", extraCols = c(id = "character"))) %>% 
	separate(id, sep="_", into=c("ensembl_gene", "ensembl_tx", "gene_name"), remove=F) %>% 
	mutate(ensembl_gene = gsub("\\..*", "", ensembl_gene)) %>% 
	mutate(ensembl_tx = gsub("\\..*", "", ensembl_tx)) %>% 
	dplyr::select(-name, -score)

# Load LOEUF scores
loeuf <- as_tibble(fread("../../data/gene_constraint_metrics/loeuf_moeuf/gnomad.v4.0.constraint_metrics.tsv.gz"))
loeuf <- distinct(filter(loeuf, mane_select)[c("transcript", "lof.oe_ci.upper")])
names(loeuf) <- c("ensembl_tx", "loeuf")

# Load MOEUF scores
moeuf <- as_tibble(fread("../../data/gene_constraint_metrics/loeuf_moeuf/gnomad.v4.0.constraint_metrics.tsv.gz"))
moeuf <- distinct(filter(moeuf, mane_select)[c("transcript", "mis.oe_ci.upper")])
names(moeuf) <- c("ensembl_tx", "moeuf")

# Load s_het scores
s_het <- as_tibble(fread("../../data/gene_constraint_metrics/s_het/media-1.tsv.gz"))
s_het <- distinct(s_het[c("ensg", "post_mean")])
names(s_het) <- c("ensembl_gene", "s_het")

# Load AlphaMissense
alphmis <- as_tibble(fread("../../data/gene_constraint_metrics/alpha_missense/Supplementary_Table_S4.txt.gz"))
alphmis <- distinct(alphmis[c("gene", "mean_am_pathogenicity")])
names(alphmis) <- c("gene_name", "alphmis")

# Gene Expression
gene_expression <- as_tibble(fread("../../data/gene_expression_catalogs/rna_celline_filtered.tsv.gz"))
gene_expression <- gene_expression %>% 
	mutate(`Cell line` = c("K-562" = "K562", "Hep-G2" = "HepG2", "SK-N-SH" = "SKNSH")[`Cell line`])

gene_expression <- gene_expression %>% 
	pivot_wider(names_from = `Cell line`, values_from = c(TPM, pTPM, nTPM))

gene_expression <- gene_expression %>% 
	dplyr::select(Gene, TPM_K562, TPM_HepG2, TPM_SKNSH) %>% 
	dplyr::rename(ensembl_gene = Gene) %>% 
	rowwise() %>% 
	mutate(
		logTPM_K562 = log2(TPM_K562+1), 
		logTPM_HepG2 = log2(TPM_HepG2+1), 
		logTPM_SKNSH = log2(TPM_SKNSH+1), 
		logTPM_max = pmax(logTPM_K562, logTPM_HepG2, logTPM_SKNSH)
	) %>% 
	ungroup()

# Tau score function
compute_tau <- function(x) {
	# Remove NA and require at least 2 tissues
	x <- x[!is.na(x)]
	if (length(x) < 2) return(NA_real_)
	x_max <- max(x)
	if (x_max == 0) return(NA_real_)  # not expressed anywhere
	x_hat <- x / x_max
	sum(1 - x_hat) / (length(x) - 1)
}

# Tissue-level tau (bulk RNA-seq, nTPM)
rna_tissue <- as_tibble(fread("../../data/gene_expression_catalogs/rna_tissue_consensus.tsv.gz"))
tau_tissue <- rna_tissue %>%
	group_by(Gene) %>%
	summarise(tau_tissue = compute_tau(nTPM), .groups = "drop") %>%
	dplyr::rename(ensembl_gene = Gene)

# Cell-type-level tau (single-cell RNA-seq, nCPM)
rna_celltype <- as_tibble(fread("../../data/gene_expression_catalogs/rna_single_cell_type.tsv.gz"))
tau_celltype <- rna_celltype %>%
	group_by(Gene) %>%
	summarise(tau_celltype = compute_tau(nCPM), .groups = "drop") %>%
	dplyr::rename(ensembl_gene = Gene)

# Median expression across tissues (bulk RNA-seq)
expr_tissue <- rna_tissue %>%
	group_by(Gene) %>%
	summarise(median_nTPM_tissue = median(nTPM, na.rm = TRUE), .groups = "drop") %>%
	mutate(log_median_nTPM_tissue = log2(median_nTPM_tissue + 1)) %>%
	dplyr::rename(ensembl_gene = Gene)

# Median expression across cell types (single-cell RNA-seq)
expr_celltype <- rna_celltype %>%
	group_by(Gene) %>%
	summarise(median_nCPM_celltype = median(nCPM, na.rm = TRUE), .groups = "drop") %>%
	mutate(log_median_nCPM_celltype = log2(median_nCPM_celltype + 1)) %>%
	dplyr::rename(ensembl_gene = Gene)

# Load EDS scores
eds <- as_tibble(fread("../../data/gene_constraint_metrics/eds/enhancer_domain_score.txt.gz"))
names(eds) <- c("ensembl_gene", "eds")
eds <- distinct(eds)

# Load ABC predictions (only needed columns to save memory)
abc <- as_tibble(fread(
	"../../data/gene_constraint_metrics/abc/AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz",
	select = c("chr", "start", "end", "TargetGene", "ABC.Score", "CellType", "class")
))

# Create unique enhancer ID from coordinates
abc <- abc %>% mutate(enhancer_id = paste0(chr, ":", start, "-", end))

# Per-gene ABC metrics (conditioned on active cell types to avoid breadth confounding)
abc_per_gene <- abc %>%
	group_by(TargetGene, CellType) %>%
	summarise(
		n_distal_enhancers = n_distinct(enhancer_id[class != "promoter"]),
		n_total_enhancers = n_distinct(enhancer_id),
		frac_distal_count = n_distal_enhancers / n_total_enhancers,
		distal_abc_score = sum(ABC.Score[class != "promoter"], na.rm = TRUE),
		total_abc_score = sum(ABC.Score, na.rm = TRUE),
		frac_distal_score = distal_abc_score / total_abc_score,
		.groups = "drop"
	) %>%
	group_by(TargetGene) %>%
	summarise(
		mean_n_enhancers_when_active = mean(n_distal_enhancers),
		mean_abc_score_when_active = mean(distal_abc_score),
		mean_frac_distal_count = mean(frac_distal_count, na.rm = TRUE),
		mean_frac_distal_score = mean(frac_distal_score, na.rm = TRUE),
		.groups = "drop"
	) %>%
	dplyr::rename(gene_name = TargetGene)

# BigWig signal over promoter regions
bigwig_dir <- "../../data/gene_regulatory_elements"

bigwig_files <- list(
	DNase_K562    = file.path(bigwig_dir, "ENCFF414OGC.bigWig"),
	H3K4me3_K562  = file.path(bigwig_dir, "ENCFF806YEZ.bigWig"),
	H3K27ac_K562  = file.path(bigwig_dir, "ENCFF849TDM.bigWig"),
	DNase_HepG2   = file.path(bigwig_dir, "ENCFF546MZK.bigWig"),
	H3K4me3_HepG2 = file.path(bigwig_dir, "ENCFF732PJK.bigWig"),
	H3K27ac_HepG2 = file.path(bigwig_dir, "ENCFF795ONN.bigWig"),
	DNase_SKNSH   = file.path(bigwig_dir, "ENCFF280RMA.bigWig"),
	H3K4me3_SKNSH = file.path(bigwig_dir, "ENCFF651WOM.bigWig"),
	H3K27ac_SKNSH = file.path(bigwig_dir, "ENCFF262UEH.bigWig")
)

# Build 1kb and 250bp promoter GRanges
gr_1000 <- GRanges(
	seqnames = gencode_gene_names$seqnames,
	ranges = IRanges(start = gencode_gene_names$start, end = gencode_gene_names$end),
	strand = gencode_gene_names$strand
)

gr_250 <- gr_1000
# + strand: TSS = end, 250bp = [end-249, end]
# - strand: TSS = start, 250bp = [start, start+249]
plus_idx <- as.character(strand(gr_1000)) == "+"
minus_idx <- as.character(strand(gr_1000)) == "-"
start(gr_250)[plus_idx] <- end(gr_1000)[plus_idx] - 249
end(gr_250)[minus_idx] <- start(gr_1000)[minus_idx] + 249

# Compute mean signal for each BigWig over each region
bigwig_signal <- tibble(row_idx = seq_len(nrow(gencode_gene_names)))

for (bw_name in names(bigwig_files)) {
	message("Processing ", bw_name, "...")
	bw <- BigWigFile(bigwig_files[[bw_name]])

	sig_1000 <- summary(bw, gr_1000, type = "mean")
	sig_250  <- summary(bw, gr_250, type = "mean")

	bigwig_signal[[paste0(bw_name, "_1kb")]] <- unlist(lapply(sig_1000, function(x) x$score))
	bigwig_signal[[paste0(bw_name, "_250bp")]]  <- unlist(lapply(sig_250, function(x) x$score))
}

bigwig_signal <- bigwig_signal %>% dplyr::select(-row_idx)

# Collate and save
satmut_promoters_gene_metadata <- gencode_gene_names %>% 
	left_join(loeuf) %>% 
	left_join(moeuf) %>% 
	left_join(s_het) %>% 
	left_join(alphmis) %>% 
	left_join(gene_expression) %>% 
	left_join(tau_tissue) %>%
	left_join(tau_celltype) %>%
	left_join(expr_tissue) %>%
	left_join(expr_celltype) %>%
	left_join(eds) %>%
	left_join(abc_per_gene) %>%
	bind_cols(bigwig_signal)

colSums(is.na(satmut_promoters_gene_metadata))

write_tsv(satmut_promoters_gene_metadata, gzfile("../../data/satmut_promoters_preprocess/satmut_promoters_gene_metadata/satmut_promoters_gene_metadata.txt.gz"))

# Filter for supplement
satmut_promoters_gene_metadata_supp <- satmut_promoters_gene_metadata %>% 
	dplyr::select(
		seqnames, start, end, width, strand, 
		id, ensembl_gene, ensembl_tx, gene_name, 
		loeuf, moeuf, s_het, alphmis, 
		TPM_K562, TPM_HepG2, TPM_SKNSH, 
		logTPM_K562, logTPM_HepG2, logTPM_SKNSH, 
		tau_tissue, tau_celltype,
		median_nTPM_tissue, log_median_nTPM_tissue, 
		median_nCPM_celltype, log_median_nCPM_celltype, 
		eds, mean_n_enhancers_when_active, mean_abc_score_when_active, 
		mean_frac_distal_count, mean_frac_distal_score, 
		ends_with("_250bp")
	)

write_tsv(satmut_promoters_gene_metadata_supp, gzfile("../../data/satmut_promoters_preprocess/satmut_promoters_gene_metadata/satmut_promoters_gene_metadata_supp.txt.gz"))
