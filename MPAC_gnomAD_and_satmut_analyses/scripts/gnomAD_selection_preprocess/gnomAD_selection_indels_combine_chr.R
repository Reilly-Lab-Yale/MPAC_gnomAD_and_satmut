#!/usr/bin/env Rscript

library(tidyverse)

path <- "../../data/gnomAD_indels_summaries/"

# Check all chromosomes present for each summary type
prefixes <- c(
	"indel_skew_by_ccre_indelclass",
	"indel_skew_by_ccre",
	"indel_activity_by_ccre",
	"indel_emvar_by_ccre",
	"indel_skew_by_ccre_af_indelclass",
	"indel_skew_by_ccre_af",
	"indel_totals_by_ccre",
	"indel_gnomad_totals_by_ccre",
	"funnel_indel"
)

for (prefix in prefixes) {
	files <- paste0(path, prefix, "_chr", 1:22, ".tsv")
	missing <- !file.exists(files)
	if (any(missing)) {
		cat(paste0("Missing files for ", prefix, ":\n"))
		cat(paste0("  ", files[missing], "\n"))
		stop("Not all chromosomes complete.")
	}
}

read_chr_files <- function(prefix) {
	map_dfr(1:22, ~ read_tsv(
		paste0(path, prefix, "_chr", .x, ".tsv"),
		show_col_types = FALSE
	))
}

# Summary 1a: skew x cCRE x indel_class
cat("Combining summary 1a (skew x cCRE x indel_class)\n")
s1a <- read_chr_files("indel_skew_by_ccre_indelclass") |>
	group_by(skew_type, skew_bin, ccre_class, indel_class) |>
	summarise(n = sum(n), .groups = "drop")
write_tsv(s1a, paste0(path, "gnomAD_indel_skew_by_ccre_indelclass.tsv"))
cat(paste0("  Rows: ", nrow(s1a), "\n"))

# Summary 1b: skew x cCRE
cat("Combining summary 1b (skew x cCRE)\n")
s1b <- read_chr_files("indel_skew_by_ccre") |>
	group_by(skew_type, skew_bin, ccre_class) |>
	summarise(n = sum(n), .groups = "drop")
write_tsv(s1b, paste0(path, "gnomAD_indel_skew_by_ccre.tsv"))
cat(paste0("  Rows: ", nrow(s1b), "\n"))

# Summary 2: activity x cCRE
cat("Combining summary 2 (activity x cCRE)\n")
s2 <- read_chr_files("indel_activity_by_ccre") |>
	group_by(activity_type, activity_bin, ccre_class) |>
	summarise(n = sum(n), .groups = "drop")
write_tsv(s2, paste0(path, "gnomAD_indel_activity_by_ccre.tsv"))
cat(paste0("  Rows: ", nrow(s2), "\n"))

# Summary 3: emVar specificity x cCRE
cat("Combining summary 3 (emVar specificity x cCRE)\n")
s3 <- read_chr_files("indel_emvar_by_ccre") |>
	group_by(emvar_class, ccre_class) |>
	summarise(n = sum(n), .groups = "drop")
write_tsv(s3, paste0(path, "gnomAD_indel_emvar_by_ccre.tsv"))
cat(paste0("  Rows: ", nrow(s3), "\n"))

# Summary 4a: skew x cCRE x AF x indel_class
cat("Combining summary 4a (skew x cCRE x AF x indel_class)\n")
s4a <- read_chr_files("indel_skew_by_ccre_af_indelclass") |>
	group_by(skew_type, skew_bin, ccre_class, af_class, indel_class) |>
	summarise(n = sum(n), .groups = "drop")
write_tsv(s4a, paste0(path, "gnomAD_indel_skew_by_ccre_af_indelclass.tsv"))
cat(paste0("  Rows: ", nrow(s4a), "\n"))

# Summary 4b: skew x cCRE x AF
cat("Combining summary 4b (skew x cCRE x AF)\n")
s4b <- read_chr_files("indel_skew_by_ccre_af") |>
	group_by(skew_type, skew_bin, ccre_class, af_class) |>
	summarise(n = sum(n), .groups = "drop")
write_tsv(s4b, paste0(path, "gnomAD_indel_skew_by_ccre_af.tsv"))
cat(paste0("  Rows: ", nrow(s4b), "\n"))

# Summary 5: totals by cCRE
cat("Combining summary 5 (totals by cCRE)\n")
s5 <- read_chr_files("indel_totals_by_ccre") |>
	group_by(ccre_class) |>
	summarise(
		n                    = sum(n),
		n_active             = sum(n_active),
		sum_mean_activity    = sum(sum_mean_activity),
		sum_mean_activity_sq = sum(sum_mean_activity_sq),
		sum_abs_mean_skew    = sum(sum_abs_mean_skew),
		sum_abs_mean_skew_sq = sum(sum_abs_mean_skew_sq),
		.groups = "drop"
	)
write_tsv(s5, paste0(path, "gnomAD_indel_totals_by_ccre.tsv"))
cat(paste0("  Rows: ", nrow(s5), "\n"))

# gnomAD totals by cCRE (pre-merge)
cat("Combining gnomAD totals by cCRE\n")
s6 <- read_chr_files("indel_gnomad_totals_by_ccre") |>
	group_by(ccre_class) |>
	summarise(n = sum(n), .groups = "drop")
write_tsv(s6, paste0(path, "gnomAD_indel_gnomad_totals_by_ccre.tsv"))
cat(paste0("  Rows: ", nrow(s6), "\n"))

# Merge loss
cat("\nMerge loss (gnomAD indels without predictions):\n")
comparison <- left_join(s6, s5 |> select(ccre_class, n), by = "ccre_class", suffix = c("_gnomad", "_merged")) |>
	mutate(lost = n_gnomad - n_merged)
print(as.data.frame(comparison))

# Funnel
cat("\nCombining funnel\n")
funnel <- read_chr_files("funnel_indel") |>
	group_by(step, description) |>
	summarise(n = sum(n), .groups = "drop") |>
	arrange(step)

cat("\nGenome-wide funnel:\n")
print(as.data.frame(funnel))
write_tsv(funnel, paste0(path, "gnomAD_indel_funnel.tsv"))

# Parsing errors
cat("\nChecking parsing errors\n")
error_files <- paste0(path, "indel_parsing_errors_chr", 1:22, ".tsv")
error_exists <- file.exists(error_files)

if (any(error_exists)) {
	cat("Parsing errors found for chromosomes:\n")
	errors_all <- map_dfr(which(error_exists), function(i) {
		df <- read_tsv(error_files[i], show_col_types = FALSE)
		cat(paste0("  chr", i, ": ", nrow(df), " bad rows\n"))
		df
	})
	write_tsv(errors_all, paste0(path, "indel_parsing_errors_all.tsv"))
	cat(paste0("Total bad rows: ", nrow(errors_all), "\n"))
	cat(paste0("Written to: ", path, "indel_parsing_errors_all.tsv\n"))
} else {
	cat("No parsing errors found.\n")
}

cat("Done\n")
