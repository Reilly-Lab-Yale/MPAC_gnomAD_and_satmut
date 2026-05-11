#!/usr/bin/env Rscript

library(tidyverse)

path <- "../../data/gnomAD_vep_summaries/"

# Check all chromosomes present for each summary type
prefixes <- c(
	"vep_consequence_by_mutclass",
	"vep_consequence",
	"vep_consequence_by_af_mutclass",
	"vep_consequence_by_af",
	"vep_totals",
	"funnel_vep"
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

# Helper: sum all sufficient stat columns
sum_suf_stats <- function(df) {
	df |> summarise(
		n              = sum(n),
		n_MR           = sum(n_MR), 
		sum_MR         = sum(sum_MR),
		sum_MR_sq      = sum(sum_MR_sq),
		n_phyloP       = sum(n_phyloP),
		sum_phyloP     = sum(sum_phyloP),
		sum_phyloP_sq  = sum(sum_phyloP_sq),
		n_conserved    = sum(n_conserved),
		.groups = "drop"
	)
}

# Summary 1a: vep_consequence x mut_class
cat("Combining summary 1a (vep_consequence x mut_class)\n")
s1a <- read_chr_files("vep_consequence_by_mutclass") |>
	group_by(vep_consequence, mut_class) |>
	sum_suf_stats()
write_tsv(s1a, paste0(path, "gnomAD_vep_consequence_by_mutclass.tsv"))
cat(paste0("  Rows: ", nrow(s1a), "\n"))

# Summary 1b: vep_consequence
cat("Combining summary 1b (vep_consequence)\n")
s1b <- read_chr_files("vep_consequence") |>
	group_by(vep_consequence) |>
	sum_suf_stats()
write_tsv(s1b, paste0(path, "gnomAD_vep_consequence.tsv"))
cat(paste0("  Rows: ", nrow(s1b), "\n"))

cat("\nVEP consequence totals:\n")
print(as.data.frame(s1b |> select(vep_consequence, n) |> arrange(desc(n))))

# Summary 2a: vep_consequence x AF x mut_class
cat("\nCombining summary 2a (vep_consequence x AF x mut_class)\n")
s2a <- read_chr_files("vep_consequence_by_af_mutclass") |>
	group_by(vep_consequence, af_class, mut_class) |>
	sum_suf_stats()
write_tsv(s2a, paste0(path, "gnomAD_vep_consequence_by_af_mutclass.tsv"))
cat(paste0("  Rows: ", nrow(s2a), "\n"))

# Summary 2b: vep_consequence x AF
cat("Combining summary 2b (vep_consequence x AF)\n")
s2b <- read_chr_files("vep_consequence_by_af") |>
	group_by(vep_consequence, af_class) |>
	sum_suf_stats()
write_tsv(s2b, paste0(path, "gnomAD_vep_consequence_by_af.tsv"))
cat(paste0("  Rows: ", nrow(s2b), "\n"))

# Summary 3: totals
cat("Combining summary 3 (totals)\n")
s3 <- read_chr_files("vep_totals") |>
	sum_suf_stats()
write_tsv(s3, paste0(path, "gnomAD_vep_totals.tsv"))
cat(paste0("  Rows: ", nrow(s3), "\n"))
print(as.data.frame(s3))

# Funnel
cat("\nCombining funnel\n")
funnel <- read_chr_files("funnel_vep") |>
	group_by(step, description) |>
	summarise(n = sum(n), .groups = "drop") |>
	arrange(step)

cat("\nGenome-wide funnel:\n")
print(as.data.frame(funnel))
write_tsv(funnel, paste0(path, "gnomAD_vep_funnel.tsv"))

cat("Done\n")
