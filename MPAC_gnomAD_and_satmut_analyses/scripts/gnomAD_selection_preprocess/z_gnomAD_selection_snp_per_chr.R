#!/usr/bin/env Rscript

library(data.table)
library(rtracklayer)
library(GenomicRanges)

# Parse chromosome number from command line
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript gnomAD_selection_snps_per_chr.R <chr_num>")
chr_num <- as.integer(args[1])
chr_str <- paste0("chr", chr_num)
cat(paste0("Processing ", chr_str, "\n"))

# Paths
pred_path     <- "../../data/gnomAD_snp_predictions/"
gnomad_path   <- "../../data/gnomAD_genomes_v3/"
roulette_path <- "../../data/gnomAD_roulette_predictions/"
ccre_path     <- "../../data/gene_regulatory_elements/"
mask_path     <- "../../data/gencode_filtered_regions/gencode.v44.basic.annotation.exons.splice.autosomes.v2.bed"
phylop_path   <- "../../data/zoonomia_phylop/data_download/241-mammalian-2020v2.phylop-Homo_sapiens.bigWig"
tf_chip_path  <- "../../data/gnomAD_miscellaneous/tf_chip_seq/chip.for.steve.regions.hg38.bed.gz"
tf_foot_path  <- "../../data/gnomAD_miscellaneous/tf_footprints/tf.footprints.regions.bed.gz"
output_path   <- "../../data/gnomAD_snp_summaries/"

dir.create(output_path, showWarnings = FALSE, recursive = TRUE)

# Funnel tracking
funnel <- data.table(step = integer(), description = character(), n = integer())
add_funnel <- function(step, desc, n) {
	funnel <<- rbind(funnel, data.table(step = step, description = desc, n = as.integer(n)))
	cat(paste0("  [Step ", step, "] ", desc, ": ", n, "\n"))
}

# Mutation classification
classify_mutation <- function(ref, alt, pn) {
	is_ti <- (ref == "C" & alt == "T") | (ref == "T" & alt == "C") |
	         (ref == "A" & alt == "G") | (ref == "G" & alt == "A")
	is_cpg <- (ref == "C" & alt == "T" & substr(pn, 4, 4) == "G") |
	          (ref == "G" & alt == "A" & substr(pn, 2, 2) == "C")
	has_pn <- !is.na(pn)
	fcase(
		!has_pn,           "unknown",
		is_cpg,            "CpG",
		is_ti,             "non-CpG-ti",
		default =          "non-CpG-tv"
	)
}

# Skew bins
skew_breaks <- c(-Inf, -1.5, -1.0, -0.5, -0.2, -0.05, 0.05, 0.2, 0.5, 1.0, 1.5, Inf)
skew_labels <- c("(-Inf,-1.5)", "[-1.5,-1.0)", "[-1.0,-0.5)", "[-0.5,-0.2)",
                 "[-0.2,-0.05)", "[-0.05,0.05)", "[0.05,0.2)", "[0.2,0.5)",
                 "[0.5,1.0)", "[1.0,1.5)", "[1.5,Inf)")

# Activity bins
activity_breaks <- c(-Inf, 1, 2, 3, 4, 5, Inf)
activity_labels <- c("[-Inf,1)", "[1,2)", "[2,3)", "[3,4)", "[4,5)", "[5,Inf)")

# Load cCRE annotations
cat("Loading cCRE annotations\n")
ccre <- fread(
	cmd = paste0("gunzip -cd ", ccre_path, "GRCh38-cCREs.V4.bed.gz"),
	header = FALSE,
	col.names = c("chr", "start", "end", "id1", "id2", "class")
)
ccre[, start := start + 1L]
ccre[, ccre_class := fifelse(
	class %in% c("PLS", "pELS", "dELS"),
	class,
	"Other cCREs"
)]
ccre[, priority := fcase(
	ccre_class == "PLS",   1L,
	ccre_class == "pELS",  2L,
	ccre_class == "dELS",  3L,
	default = 4L
)]
ccre <- ccre[chr == chr_str, .(chr, start, end, ccre_class, priority)]
setkey(ccre, chr, start, end)
cat(paste0("  cCRE intervals: ", nrow(ccre), "\n"))

# Load exon/splice mask
cat("Loading exon/splice mask\n")
mask <- fread(
	mask_path,
	header = FALSE,
	select = 1:3,
	col.names = c("chr", "start", "end")
)
mask[, start := start + 1L]
mask <- mask[chr == chr_str]
mask[, mask_flag := TRUE]
setkey(mask, chr, start, end)
cat(paste0("  Mask intervals: ", nrow(mask), "\n"))

# Load TF ChIP-seq peaks
cat("Loading TF ChIP-seq peaks\n")
tf_chip <- fread(
	cmd = paste0("gunzip -cd ", tf_chip_path),
	header = FALSE,
	select = 1:3,
	col.names = c("chr", "start", "end")
)
tf_chip[, start := start + 1L]
tf_chip <- tf_chip[chr == chr_str]
tf_chip[, tf_chip_flag := TRUE]
setkey(tf_chip, chr, start, end)
cat(paste0("  TF ChIP-seq intervals: ", nrow(tf_chip), "\n"))

# Load TF footprints
cat("Loading TF footprints\n")
tf_foot <- fread(
	cmd = paste0("gunzip -cd ", tf_foot_path),
	header = FALSE,
	select = 1:3,
	col.names = c("chr", "start", "end")
)
tf_foot[, start := start + 1L]
tf_foot <- tf_foot[chr == chr_str]
tf_foot[, tf_foot_flag := TRUE]
setkey(tf_foot, chr, start, end)
cat(paste0("  TF footprint intervals: ", nrow(tf_foot), "\n"))

# Load gnomAD subinfo
cat("Loading gnomAD subinfo\n")
gnomad <- fread(
	cmd = paste0(
		"gunzip -cd ", gnomad_path,
		"gnomad.genomes.v3.1.2.sites.", chr_str, ".subinfo.vcf.gz",
		" | grep -v '^#'",
		" | awk -F'\\t' 'BEGIN{OFS=\"\\t\"}{",
		"n=split($8,a,\";\");ac=\".\";an=\".\";af=\".\";",
		"for(i=1;i<=n;i++){split(a[i],kv,\"=\");",
		"if(kv[1]==\"AC\")ac=kv[2];",
		"if(kv[1]==\"AN\")an=kv[2];",
		"if(kv[1]==\"AF\")af=kv[2]}",
		"print $1,$2,$4,$5,$7,ac,an,af}'"
	),
	header = FALSE,
	col.names = c("chr", "pos", "ref", "alt", "FILTER", "AC", "AN", "AF")
)
gnomad[, `:=`(AC = as.integer(AC), AN = as.integer(AN), AF = as.numeric(AF))]

# Filter to SNVs only
gnomad <- gnomad[nchar(ref) == 1L & nchar(alt) == 1L]
add_funnel(1L, "gnomAD SNVs loaded", nrow(gnomad))

# QC filter
gnomad <- gnomad[FILTER == "PASS" & AC > 0L & AN >= 76156L]
add_funnel(2L, "gnomAD after QC filter", nrow(gnomad))

# Apply exon/splice mask
gnomad[, `:=`(start = pos, end = pos)]
setkey(gnomad, chr, start, end)
ov <- foverlaps(gnomad, mask, type = "any", nomatch = NA)
gnomad <- ov[is.na(mask_flag)]
gnomad[, `:=`(start = NULL, end = NULL, mask_flag = NULL, i.start = NULL, i.end = NULL)]
add_funnel(3L, "gnomAD after exon/splice mask", nrow(gnomad))
rm(ov)

# Annotate gnomAD with cCREs for pre-merge counts
cat("Computing gnomAD cCRE totals (pre-merge)\n")
gnomad[, `:=`(start = pos, end = pos)]
setkey(gnomad, chr, start, end)
gnomad_ov <- foverlaps(gnomad, ccre, type = "any", nomatch = NA)
gnomad_ov[is.na(ccre_class), `:=`(ccre_class = "non-cCRE", priority = 5L)]
setorder(gnomad_ov, i.start, ref, alt, priority)
gnomad_dedup <- unique(gnomad_ov, by = c("i.start", "ref", "alt"))
summary_gnomad <- gnomad_dedup[, .(n = .N), by = .(ccre_class)]
summary_gnomad <- rbind(summary_gnomad, data.table(ccre_class = "TOTAL", n = nrow(gnomad_dedup)))
rm(gnomad_ov, gnomad_dedup)
gnomad[, `:=`(start = NULL, end = NULL)]

# Load MPAC predictions
cat("Loading MPAC predictions\n")
mpac <- fread(
	cmd = paste0(
		"gunzip -cd ", pred_path,
		"gnomad.genomes.v3.1.2.sites.", chr_str, ".vcf.gz",
		" | tail -n +2",
		" | awk -F'\\t' 'BEGIN{OFS=\"\\t\"}{",
		"n=split($6,a,\";\");",
		"kr=\".\";hr=\".\";sr=\".\";ks=\".\";hs=\".\";ss=\".\";",
		"for(i=1;i<=n;i++){split(a[i],kv,\"=\");",
		"if(kv[1]==\"K562__ref\")kr=kv[2];",
		"if(kv[1]==\"HepG2__ref\")hr=kv[2];",
		"if(kv[1]==\"SKNSH__ref\")sr=kv[2];",
		"if(kv[1]==\"K562__skew\")ks=kv[2];",
		"if(kv[1]==\"HepG2__skew\")hs=kv[2];",
		"if(kv[1]==\"SKNSH__skew\")ss=kv[2]}",
		"print $1,$2,$4,$5,kr,hr,sr,ks,hs,ss}'"
	),
	header = FALSE,
	col.names = c("chr", "pos", "ref", "alt",
	              "K562_ref", "HepG2_ref", "SKNSH_ref",
	              "K562_skew", "HepG2_skew", "SKNSH_skew")
)
cat(paste0("  MPAC rows loaded: ", nrow(mpac), "\n"))

# Validate parsing
num_cols <- c("K562_ref", "HepG2_ref", "SKNSH_ref", "K562_skew", "HepG2_skew", "SKNSH_skew")
bad_mask <- Reduce(`|`, lapply(num_cols, function(col) is.na(mpac[[col]]) | !is.finite(mpac[[col]])))
if (sum(bad_mask) > 0) {
	cat(paste0("  WARNING: ", sum(bad_mask), " rows with parsing failures:\n"))
	print(head(mpac[bad_mask], 10))
	bad_file <- paste0(output_path, "snp_parsing_errors_", chr_str, ".tsv")
	fwrite(mpac[bad_mask], bad_file, sep = "\t")
	cat(paste0("  Bad rows written to: ", bad_file, "\n"))
	mpac <- mpac[!bad_mask]
}
add_funnel(4L, "MPAC predictions (good parsing)", nrow(mpac))

# Compute activity (max of ref and alt = max of ref and ref+skew) and mean skew
mpac[, `:=`(
	K562_activity  = pmax(K562_ref, K562_ref + K562_skew),
	HepG2_activity = pmax(HepG2_ref, HepG2_ref + HepG2_skew),
	SKNSH_activity = pmax(SKNSH_ref, SKNSH_ref + SKNSH_skew),
	mean_skew      = (K562_skew + HepG2_skew + SKNSH_skew) / 3
)]
mpac[, mean_activity := (K562_activity + HepG2_activity + SKNSH_activity) / 3]

# Drop raw ref columns (no longer needed)
mpac[, `:=`(K562_ref = NULL, HepG2_ref = NULL, SKNSH_ref = NULL)]

# Merge gnomAD + MPAC
cat("Merging gnomAD + MPAC\n")
merged <- merge(gnomad, mpac, by = c("chr", "pos", "ref", "alt"), all = FALSE)
add_funnel(5L, "Merged gnomAD + MPAC", nrow(merged))
rm(gnomad, mpac)
gc(verbose = FALSE)

# Compute MAF and MAC
merged[, `:=`(
  MAF = pmin(AF, 1 - AF),
  MAC = pmin(AC, AN - AC)
)]

# Classify on MAF/MAC
merged[, af_class := fcase(
  MAC == 1L,      "SINGLETON",
  MAF < 0.0001,   "ULTRARARE",
  MAF < 0.001,    "RARE",
  MAF < 0.01,     "LOW_FREQ",
  default = "COMMON"
)]

# Load Roulette
cat("Loading Roulette\n")
roulette <- fread(
	cmd = paste0(
		"gunzip -cd ", roulette_path, chr_num,
		"_rate_v5.2_TFBS_correction_all.vcf.bgz",
		" | grep -v '^#'",
		" | awk -F'\\t' 'BEGIN{OFS=\"\\t\"}{",
		"n=split($8,a,\";\");pn=\".\";mr=\".\";",
		"for(i=1;i<=n;i++){split(a[i],kv,\"=\");",
		"if(kv[1]==\"PN\")pn=kv[2];",
		"if(kv[1]==\"MR\")mr=kv[2]}",
		"print $1,$2,$4,$5,pn,mr}'"
	),
	header = FALSE,
	col.names = c("chr", "pos", "ref", "alt", "PN", "MR")
)
roulette[, MR := as.numeric(MR)]
roulette[, chr := paste0("chr", chr)]
add_funnel(6L, "Roulette loaded", nrow(roulette))

# Merge with Roulette
cat("Merging with Roulette\n")
merged <- merge(merged, roulette, by = c("chr", "pos", "ref", "alt"), all.x = TRUE)
add_funnel(7L, "Roulette coverage", sum(!is.na(merged$MR)))
rm(roulette)
gc(verbose = FALSE)

# Query phyloP scores from bigWig
cat("Querying phyloP scores\n")
unique_pos <- sort(unique(merged$pos))
query_gr <- GRanges(
	seqnames = chr_str,
	ranges = IRanges(start = unique_pos, width = 1)
)
bw_gr <- import.bw(phylop_path, which = query_gr)

hits <- findOverlaps(query_gr, bw_gr)
phylop_dt <- data.table(
	pos = unique_pos[queryHits(hits)],
	phyloP_score = bw_gr$score[subjectHits(hits)]
)
phylop_dt <- unique(phylop_dt, by = "pos")

merged <- merge(merged, phylop_dt, by = "pos", all.x = TRUE)
merged[, is_conserved := !is.na(phyloP_score) & phyloP_score > 2.27]

n_with_phylop <- sum(!is.na(merged$phyloP_score))
cat(paste0("  Variants with phyloP scores: ", n_with_phylop, " / ", nrow(merged), "\n"))
rm(unique_pos, query_gr, bw_gr, hits, phylop_dt)
gc(verbose = FALSE)

# Annotate TF ChIP-seq peak overlap
cat("Annotating TF ChIP-seq peak overlap\n")
merged[, `:=`(start = pos, end = pos)]
setkey(merged, chr, start, end)
ov_chip <- foverlaps(merged, tf_chip, type = "any", nomatch = NA)
merged[, is_tf_chip_peak := !is.na(ov_chip$tf_chip_flag)]
merged[, `:=`(start = NULL, end = NULL)]
rm(ov_chip)

# Annotate TF footprint overlap
cat("Annotating TF footprint overlap\n")
merged[, `:=`(start = pos, end = pos)]
setkey(merged, chr, start, end)
ov_foot <- foverlaps(merged, tf_foot, type = "any", nomatch = NA)
merged[, is_tf_footprint := !is.na(ov_foot$tf_foot_flag)]
merged[, `:=`(start = NULL, end = NULL)]
rm(ov_foot, tf_chip, tf_foot)

cat(paste0("  Variants in TF ChIP-seq peaks: ", sum(merged$is_tf_chip_peak), " / ", nrow(merged), "\n"))
cat(paste0("  Variants in TF footprints: ", sum(merged$is_tf_footprint), " / ", nrow(merged), "\n"))

# Classify mutation type
merged[, mut_class := classify_mutation(ref, alt, PN)]

# Annotate cCRE class via foverlaps
merged[, `:=`(start = pos, end = pos)]
setkey(merged, chr, start, end)

ov <- foverlaps(merged, ccre, type = "any", nomatch = NA)
ov[is.na(ccre_class), `:=`(ccre_class = "non-cCRE", priority = 5L)]

setorder(ov, i.start, ref, alt, priority)
annotated <- unique(ov, by = c("i.start", "ref", "alt"))
add_funnel(8L, "Annotated with cCREs", nrow(annotated))
rm(merged, ov)
gc(verbose = FALSE)

# Bin skew and activity values
annotated[, `:=`(
	mean_skew_bin     = cut(mean_skew,     breaks = skew_breaks,     labels = skew_labels,     right = FALSE),
	K562_skew_bin     = cut(K562_skew,     breaks = skew_breaks,     labels = skew_labels,     right = FALSE),
	HepG2_skew_bin    = cut(HepG2_skew,    breaks = skew_breaks,     labels = skew_labels,     right = FALSE),
	SKNSH_skew_bin    = cut(SKNSH_skew,    breaks = skew_breaks,     labels = skew_labels,     right = FALSE),
	mean_activity_bin  = cut(mean_activity,  breaks = activity_breaks, labels = activity_labels, right = FALSE),
	K562_activity_bin  = cut(K562_activity,  breaks = activity_breaks, labels = activity_labels, right = FALSE),
	HepG2_activity_bin = cut(HepG2_activity, breaks = activity_breaks, labels = activity_labels, right = FALSE),
	SKNSH_activity_bin = cut(SKNSH_activity, breaks = activity_breaks, labels = activity_labels, right = FALSE)
)]

skew_types        <- c("mean_skew", "K562_skew", "HepG2_skew", "SKNSH_skew")
skew_bin_cols     <- c("mean_skew_bin", "K562_skew_bin", "HepG2_skew_bin", "SKNSH_skew_bin")
activity_types    <- c("mean_activity", "K562_activity", "HepG2_activity", "SKNSH_activity")
activity_bin_cols <- c("mean_activity_bin", "K562_activity_bin", "HepG2_activity_bin", "SKNSH_activity_bin")

# emVar classification
annotated[, `:=`(
	emvar_K562  = abs(K562_skew) > 0.5,
	emvar_HepG2 = abs(HepG2_skew) > 0.5,
	emvar_SKNSH = abs(SKNSH_skew) > 0.5
)]
annotated[, emvar_class := fcase(
	emvar_K562 & emvar_HepG2 & emvar_SKNSH, "K562+HepG2+SKNSH",
	emvar_K562 & emvar_HepG2,                "K562+HepG2",
	emvar_K562 & emvar_SKNSH,                "K562+SKNSH",
	emvar_HepG2 & emvar_SKNSH,               "HepG2+SKNSH",
	emvar_K562,                               "K562",
	emvar_HepG2,                              "HepG2",
	emvar_SKNSH,                              "SKNSH",
	default =                                 "none"
)]

# Active flag: active (log2FC > 1) in at least one cell type
annotated[, is_active := (K562_activity > 1) | (HepG2_activity > 1) | (SKNSH_activity > 1)]

# Summary 1a: skew x cCRE x mut_class
cat("Computing summary 1a (skew x cCRE x mut_class)\n")
summary1a <- rbindlist(lapply(seq_along(skew_types), function(j) {
	annotated[, .(
		n              = .N,
		n_MR           = sum(!is.na(MR)),
		sum_MR         = sum(MR, na.rm = TRUE),
		sum_MR_sq      = sum(MR^2, na.rm = TRUE),
		n_phyloP       = sum(!is.na(phyloP_score)),
		sum_phyloP     = sum(phyloP_score, na.rm = TRUE),
		sum_phyloP_sq  = sum(phyloP_score^2, na.rm = TRUE),
		n_conserved    = sum(is_conserved),
		n_tf_chip_peak = sum(is_tf_chip_peak),
		n_tf_footprint = sum(is_tf_footprint)
	), by = .(skew_bin = get(skew_bin_cols[j]), ccre_class, mut_class)][
		, skew_type := skew_types[j]]
}))

# Summary 1b: skew x cCRE
cat("Computing summary 1b (skew x cCRE)\n")
summary1b <- rbindlist(lapply(seq_along(skew_types), function(j) {
	annotated[, .(
		n              = .N,
		n_MR           = sum(!is.na(MR)),
		sum_MR         = sum(MR, na.rm = TRUE),
		sum_MR_sq      = sum(MR^2, na.rm = TRUE),
		n_phyloP       = sum(!is.na(phyloP_score)),
		sum_phyloP     = sum(phyloP_score, na.rm = TRUE),
		sum_phyloP_sq  = sum(phyloP_score^2, na.rm = TRUE),
		n_conserved    = sum(is_conserved),
		n_tf_chip_peak = sum(is_tf_chip_peak),
		n_tf_footprint = sum(is_tf_footprint)
	), by = .(skew_bin = get(skew_bin_cols[j]), ccre_class)][
		, skew_type := skew_types[j]]
}))

# Summary 2: activity x cCRE
cat("Computing summary 2 (activity x cCRE)\n")
summary2 <- rbindlist(lapply(seq_along(activity_types), function(j) {
	annotated[, .(
		n              = .N,
		n_MR           = sum(!is.na(MR)),
		sum_MR         = sum(MR, na.rm = TRUE),
		sum_MR_sq      = sum(MR^2, na.rm = TRUE),
		n_phyloP       = sum(!is.na(phyloP_score)),
		sum_phyloP     = sum(phyloP_score, na.rm = TRUE),
		sum_phyloP_sq  = sum(phyloP_score^2, na.rm = TRUE),
		n_conserved    = sum(is_conserved),
		n_tf_chip_peak = sum(is_tf_chip_peak),
		n_tf_footprint = sum(is_tf_footprint)
	), by = .(activity_bin = get(activity_bin_cols[j]), ccre_class)][
		, activity_type := activity_types[j]]
}))

# Summary 3: emVar x cCRE
cat("Computing summary 3 (emVar specificity x cCRE)\n")
summary3 <- annotated[, .(
	n              = .N,
	n_MR           = sum(!is.na(MR)),
	sum_MR         = sum(MR, na.rm = TRUE),
	sum_MR_sq      = sum(MR^2, na.rm = TRUE),
	n_phyloP       = sum(!is.na(phyloP_score)),
	sum_phyloP     = sum(phyloP_score, na.rm = TRUE),
	sum_phyloP_sq  = sum(phyloP_score^2, na.rm = TRUE),
	n_conserved    = sum(is_conserved),
	n_tf_chip_peak = sum(is_tf_chip_peak),
	n_tf_footprint = sum(is_tf_footprint)
), by = .(emvar_class, ccre_class)]

# Summary 4a: skew x cCRE x AF x mut_class
cat("Computing summary 4a (skew x cCRE x AF x mut_class)\n")
summary4a <- rbindlist(lapply(seq_along(skew_types), function(j) {
	annotated[, .(
		n              = .N,
		n_MR           = sum(!is.na(MR)),
		sum_MR         = sum(MR, na.rm = TRUE),
		sum_MR_sq      = sum(MR^2, na.rm = TRUE),
		n_phyloP       = sum(!is.na(phyloP_score)),
		sum_phyloP     = sum(phyloP_score, na.rm = TRUE),
		sum_phyloP_sq  = sum(phyloP_score^2, na.rm = TRUE),
		n_conserved    = sum(is_conserved),
		n_tf_chip_peak = sum(is_tf_chip_peak),
		n_tf_footprint = sum(is_tf_footprint)
	), by = .(skew_bin = get(skew_bin_cols[j]), ccre_class, af_class, mut_class)][
		, skew_type := skew_types[j]]
}))

# Summary 4b: skew x cCRE x AF
cat("Computing summary 4b (skew x cCRE x AF)\n")
summary4b <- rbindlist(lapply(seq_along(skew_types), function(j) {
	annotated[, .(
		n              = .N,
		n_MR           = sum(!is.na(MR)),
		sum_MR         = sum(MR, na.rm = TRUE),
		sum_MR_sq      = sum(MR^2, na.rm = TRUE),
		n_phyloP       = sum(!is.na(phyloP_score)),
		sum_phyloP     = sum(phyloP_score, na.rm = TRUE),
		sum_phyloP_sq  = sum(phyloP_score^2, na.rm = TRUE),
		n_conserved    = sum(is_conserved),
		n_tf_chip_peak = sum(is_tf_chip_peak),
		n_tf_footprint = sum(is_tf_footprint)
	), by = .(skew_bin = get(skew_bin_cols[j]), ccre_class, af_class)][
		, skew_type := skew_types[j]]
}))

# Summary 5: totals by cCRE (with activity and skew sufficient statistics)
cat("Computing summary 5 (totals by cCRE)\n")
summary5 <- annotated[, .(
	n                    = .N,
	n_MR                 = sum(!is.na(MR)),
	sum_MR               = sum(MR, na.rm = TRUE),
	sum_MR_sq            = sum(MR^2, na.rm = TRUE),
	n_phyloP             = sum(!is.na(phyloP_score)),
	sum_phyloP           = sum(phyloP_score, na.rm = TRUE),
	sum_phyloP_sq        = sum(phyloP_score^2, na.rm = TRUE),
	n_conserved          = sum(is_conserved),
	n_tf_chip_peak       = sum(is_tf_chip_peak),
	n_tf_footprint       = sum(is_tf_footprint),
	n_active             = sum(is_active),
	sum_mean_activity    = sum(mean_activity),
	sum_mean_activity_sq = sum(mean_activity^2),
	sum_abs_mean_skew    = sum(abs(mean_skew)),
	sum_abs_mean_skew_sq = sum(abs(mean_skew)^2)
), by = .(ccre_class)]
summary5 <- rbind(
	summary5,
	annotated[, .(
		ccre_class           = "TOTAL",
		n                    = .N,
		n_MR                 = sum(!is.na(MR)),
		sum_MR               = sum(MR, na.rm = TRUE),
		sum_MR_sq            = sum(MR^2, na.rm = TRUE),
		n_phyloP             = sum(!is.na(phyloP_score)),
		sum_phyloP           = sum(phyloP_score, na.rm = TRUE),
		sum_phyloP_sq        = sum(phyloP_score^2, na.rm = TRUE),
		n_conserved          = sum(is_conserved),
		n_tf_chip_peak       = sum(is_tf_chip_peak),
		n_tf_footprint       = sum(is_tf_footprint),
		n_active             = sum(is_active),
		sum_mean_activity    = sum(mean_activity),
		sum_mean_activity_sq = sum(mean_activity^2),
		sum_abs_mean_skew    = sum(abs(mean_skew)),
		sum_abs_mean_skew_sq = sum(abs(mean_skew)^2)
	)]
)

# Write per-chromosome outputs
fwrite(summary1a, paste0(output_path, "snp_skew_by_ccre_mutclass_", chr_str, ".tsv"), sep = "\t")
fwrite(summary1b, paste0(output_path, "snp_skew_by_ccre_", chr_str, ".tsv"), sep = "\t")
fwrite(summary2, paste0(output_path, "snp_activity_by_ccre_", chr_str, ".tsv"), sep = "\t")
fwrite(summary3, paste0(output_path, "snp_emvar_by_ccre_", chr_str, ".tsv"), sep = "\t")
fwrite(summary4a, paste0(output_path, "snp_skew_by_ccre_af_mutclass_", chr_str, ".tsv"), sep = "\t")
fwrite(summary4b, paste0(output_path, "snp_skew_by_ccre_af_", chr_str, ".tsv"), sep = "\t")
fwrite(summary5, paste0(output_path, "snp_totals_by_ccre_", chr_str, ".tsv"), sep = "\t")
fwrite(summary_gnomad, paste0(output_path, "snp_gnomad_totals_by_ccre_", chr_str, ".tsv"), sep = "\t")
fwrite(funnel, paste0(output_path, "funnel_snp_", chr_str, ".tsv"), sep = "\t")

cat(paste0("Done: ", chr_str, "\n"))
cat(paste0("  Summary 1a rows: ", nrow(summary1a), "\n"))
cat(paste0("  Summary 1b rows: ", nrow(summary1b), "\n"))
cat(paste0("  Summary 2 rows:  ", nrow(summary2), "\n"))
cat(paste0("  Summary 3 rows:  ", nrow(summary3), "\n"))
cat(paste0("  Summary 4a rows: ", nrow(summary4a), "\n"))
cat(paste0("  Summary 4b rows: ", nrow(summary4b), "\n"))
cat(paste0("  Summary 5 rows:  ", nrow(summary5), "\n"))
