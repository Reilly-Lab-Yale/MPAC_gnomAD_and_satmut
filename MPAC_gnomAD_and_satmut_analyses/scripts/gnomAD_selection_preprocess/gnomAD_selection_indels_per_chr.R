#!/usr/bin/env Rscript

library(data.table)

# Parse chromosome number from command line
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript gnomAD_selection_indels_per_chr.R <chr_num>")
chr_num <- as.integer(args[1])
chr_str <- paste0("chr", chr_num)
cat(paste0("Processing ", chr_str, "\n"))

# Paths
pred_path   <- "../../data/gnomAD_indels_predictions/"
gnomad_path <- "../../data/gnomAD_genomes_v3/"
ccre_path   <- "../../data/gene_regulatory_elements/"
mask_path   <- "../../data/gencode_filtered_regions/gencode.v44.basic.annotation.exons.splice.autosomes.v2.bed"
output_path <- "../../data/gnomAD_indels_summaries/"

dir.create(output_path, showWarnings = FALSE, recursive = TRUE)

# Funnel tracking
funnel <- data.table(step = integer(), description = character(), n = integer())
add_funnel <- function(step, desc, n) {
	funnel <<- rbind(funnel, data.table(step = step, description = desc, n = as.integer(n)))
	cat(paste0("  [Step ", step, "] ", desc, ": ", n, "\n"))
}

# Skew bins
skew_breaks <- c(-Inf, -1.5, -1.0, -0.5, -0.2, -0.05, 0.05, 0.2, 0.5, 1.0, 1.5, Inf)
skew_labels <- c("(-Inf,-1.5)", "[-1.5,-1.0)", "[-1.0,-0.5)", "[-0.5,-0.2)",
                 "[-0.2,-0.05)", "[-0.05,0.05)", "[0.05,0.2)", "[0.2,0.5)",
                 "[0.5,1.0)", "[1.0,1.5)", "[1.5,Inf)")

# Activity bins
activity_breaks <- c(-Inf, 1, 2, 3, 4, 5, Inf)
activity_labels <- c("[-Inf,1)", "[1,2)", "[2,3)", "[3,4)", "[4,5)", "[5,Inf)")

# Helper: compute indel interval (full REF span for deletions, point for insertions)
indel_end <- function(pos, ref, alt) {
	fifelse(nchar(ref) > nchar(alt), pos + nchar(ref) - 1L, pos)
}

# Classify indel type and length
classify_indel <- function(ref, alt) {
	ref_len <- nchar(ref)
	alt_len <- nchar(alt)
	indel_len <- abs(ref_len - alt_len)
	indel_type <- fifelse(alt_len > ref_len, "ins", "del")
	len_label <- fifelse(indel_len >= 10L, "10plus", as.character(indel_len))
	paste0(indel_type, "_", len_label)
}

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

# Load gnomAD indels
cat("Loading gnomAD indels\n")
gnomad <- fread(
	cmd = paste0(
		"gunzip -cd ", gnomad_path,
		"gnomad.genomes.v3.1.2.sites.", chr_str, ".subinfo.vcf.gz",
		" | grep -v '^#'",
		" | awk -F'\\t' 'BEGIN{OFS=\"\\t\"}{",
		"if(length($4)==length($5))next;",
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
add_funnel(1L, "gnomAD indels loaded", nrow(gnomad))

# QC filter
gnomad <- gnomad[FILTER == "PASS" & AC > 0L & AN >= 76156L]
add_funnel(2L, "gnomAD after QC filter", nrow(gnomad))

# Apply exon/splice mask (using full deletion span)
gnomad[, `:=`(start = pos, end = indel_end(pos, ref, alt))]
setkey(gnomad, chr, start, end)
ov <- foverlaps(gnomad, mask, type = "any", nomatch = NA)
gnomad <- ov[is.na(mask_flag)]
gnomad[, `:=`(start = NULL, end = NULL, mask_flag = NULL, i.start = NULL, i.end = NULL)]
add_funnel(3L, "gnomAD after exon/splice mask", nrow(gnomad))
rm(ov)

# Annotate gnomAD with cCREs for pre-merge counts (using full deletion span)
cat("Computing gnomAD cCRE totals (pre-merge)\n")
gnomad[, `:=`(start = pos, end = indel_end(pos, ref, alt))]
setkey(gnomad, chr, start, end)
gnomad_ov <- foverlaps(gnomad, ccre, type = "any", nomatch = NA)
gnomad_ov[is.na(ccre_class), `:=`(ccre_class = "non-cCRE", priority = 5L)]
setorder(gnomad_ov, i.start, ref, alt, priority)
gnomad_dedup <- unique(gnomad_ov, by = c("i.start", "ref", "alt"))
summary_gnomad <- gnomad_dedup[, .(n = .N), by = .(ccre_class)]
summary_gnomad <- rbind(summary_gnomad, data.table(ccre_class = "TOTAL", n = nrow(gnomad_dedup)))
rm(gnomad_ov, gnomad_dedup)
gnomad[, `:=`(start = NULL, end = NULL)]

# Load MPAC indel predictions
cat("Loading MPAC indel predictions\n")
pred <- fread(
	cmd = paste0(
		"gunzip -cd ", pred_path, chr_str, ".sorted.vcf.gz",
		" | grep -v '^#'",
		" | awk -F'\\t' 'BEGIN{OFS=\"\\t\"}{",
		"n=split($8,a,\";\");",
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
cat(paste0("  Prediction rows loaded: ", nrow(pred), "\n"))

# Validate parsing
num_cols <- c("K562_ref", "HepG2_ref", "SKNSH_ref", "K562_skew", "HepG2_skew", "SKNSH_skew")
bad_mask <- Reduce(`|`, lapply(num_cols, function(col) is.na(pred[[col]]) | !is.finite(pred[[col]])))
if (sum(bad_mask) > 0) {
	cat(paste0("  WARNING: ", sum(bad_mask), " rows with parsing failures:\n"))
	print(head(pred[bad_mask], 10))
	bad_file <- paste0(output_path, "indel_parsing_errors_", chr_str, ".tsv")
	fwrite(pred[bad_mask], bad_file, sep = "\t")
	cat(paste0("  Bad rows written to: ", bad_file, "\n"))
	pred <- pred[!bad_mask]
}
add_funnel(4L, "Predictions (good parsing)", nrow(pred))

# Compute activity (max of ref and alt = max of ref and ref+skew) and mean skew
pred[, `:=`(
	K562_activity  = pmax(K562_ref, K562_ref + K562_skew),
	HepG2_activity = pmax(HepG2_ref, HepG2_ref + HepG2_skew),
	SKNSH_activity = pmax(SKNSH_ref, SKNSH_ref + SKNSH_skew),
	mean_skew      = (K562_skew + HepG2_skew + SKNSH_skew) / 3
)]
pred[, mean_activity := (K562_activity + HepG2_activity + SKNSH_activity) / 3]

# Drop raw ref columns (no longer needed)
pred[, `:=`(K562_ref = NULL, HepG2_ref = NULL, SKNSH_ref = NULL)]

# Merge gnomAD + predictions
cat("Merging gnomAD + predictions\n")
merged <- merge(gnomad, pred, by = c("chr", "pos", "ref", "alt"), all = FALSE)
add_funnel(5L, "Merged gnomAD + predictions", nrow(merged))
rm(gnomad, pred)
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

# Classify indel type and length
merged[, indel_class := classify_indel(ref, alt)]

# Annotate cCRE class via foverlaps (using full deletion span)
merged[, `:=`(start = pos, end = indel_end(pos, ref, alt))]
setkey(merged, chr, start, end)

ov <- foverlaps(merged, ccre, type = "any", nomatch = NA)
ov[is.na(ccre_class), `:=`(ccre_class = "non-cCRE", priority = 5L)]

setorder(ov, i.start, ref, alt, priority)
annotated <- unique(ov, by = c("i.start", "ref", "alt"))
add_funnel(6L, "Annotated with cCREs", nrow(annotated))
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

# Summary 1a: skew x cCRE x indel_class
cat("Computing summary 1a (skew x cCRE x indel_class)\n")
summary1a <- rbindlist(lapply(seq_along(skew_types), function(j) {
	annotated[, .(n = .N), by = .(skew_bin = get(skew_bin_cols[j]), ccre_class, indel_class)][
		, skew_type := skew_types[j]]
}))

# Summary 1b: skew x cCRE
cat("Computing summary 1b (skew x cCRE)\n")
summary1b <- rbindlist(lapply(seq_along(skew_types), function(j) {
	annotated[, .(n = .N), by = .(skew_bin = get(skew_bin_cols[j]), ccre_class)][
		, skew_type := skew_types[j]]
}))

# Summary 2: activity x cCRE
cat("Computing summary 2 (activity x cCRE)\n")
summary2 <- rbindlist(lapply(seq_along(activity_types), function(j) {
	annotated[, .(n = .N), by = .(activity_bin = get(activity_bin_cols[j]), ccre_class)][
		, activity_type := activity_types[j]]
}))

# Summary 3: emVar x cCRE
cat("Computing summary 3 (emVar specificity x cCRE)\n")
summary3 <- annotated[, .(n = .N), by = .(emvar_class, ccre_class)]

# Summary 4a: skew x cCRE x AF x indel_class
cat("Computing summary 4a (skew x cCRE x AF x indel_class)\n")
summary4a <- rbindlist(lapply(seq_along(skew_types), function(j) {
	annotated[, .(n = .N), by = .(skew_bin = get(skew_bin_cols[j]), ccre_class, af_class, indel_class)][
		, skew_type := skew_types[j]]
}))

# Summary 4b: skew x cCRE x AF
cat("Computing summary 4b (skew x cCRE x AF)\n")
summary4b <- rbindlist(lapply(seq_along(skew_types), function(j) {
	annotated[, .(n = .N), by = .(skew_bin = get(skew_bin_cols[j]), ccre_class, af_class)][
		, skew_type := skew_types[j]]
}))

# Summary 5: totals by cCRE (with activity and skew sufficient statistics)
cat("Computing summary 5 (totals by cCRE)\n")
summary5 <- annotated[, .(
	n                    = .N,
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
		n_active             = sum(is_active),
		sum_mean_activity    = sum(mean_activity),
		sum_mean_activity_sq = sum(mean_activity^2),
		sum_abs_mean_skew    = sum(abs(mean_skew)),
		sum_abs_mean_skew_sq = sum(abs(mean_skew)^2)
	)]
)

# Write per-chromosome outputs
fwrite(summary1a, paste0(output_path, "indel_skew_by_ccre_indelclass_", chr_str, ".tsv"), sep = "\t")
fwrite(summary1b, paste0(output_path, "indel_skew_by_ccre_", chr_str, ".tsv"), sep = "\t")
fwrite(summary2, paste0(output_path, "indel_activity_by_ccre_", chr_str, ".tsv"), sep = "\t")
fwrite(summary3, paste0(output_path, "indel_emvar_by_ccre_", chr_str, ".tsv"), sep = "\t")
fwrite(summary4a, paste0(output_path, "indel_skew_by_ccre_af_indelclass_", chr_str, ".tsv"), sep = "\t")
fwrite(summary4b, paste0(output_path, "indel_skew_by_ccre_af_", chr_str, ".tsv"), sep = "\t")
fwrite(summary5, paste0(output_path, "indel_totals_by_ccre_", chr_str, ".tsv"), sep = "\t")
fwrite(summary_gnomad, paste0(output_path, "indel_gnomad_totals_by_ccre_", chr_str, ".tsv"), sep = "\t")
fwrite(funnel, paste0(output_path, "funnel_indel_", chr_str, ".tsv"), sep = "\t")

cat(paste0("Done: ", chr_str, "\n"))
cat(paste0("  Summary 1a rows: ", nrow(summary1a), "\n"))
cat(paste0("  Summary 1b rows: ", nrow(summary1b), "\n"))
cat(paste0("  Summary 2 rows:  ", nrow(summary2), "\n"))
cat(paste0("  Summary 3 rows:  ", nrow(summary3), "\n"))
cat(paste0("  Summary 4a rows: ", nrow(summary4a), "\n"))
cat(paste0("  Summary 4b rows: ", nrow(summary4b), "\n"))
cat(paste0("  Summary 5 rows:  ", nrow(summary5), "\n"))
