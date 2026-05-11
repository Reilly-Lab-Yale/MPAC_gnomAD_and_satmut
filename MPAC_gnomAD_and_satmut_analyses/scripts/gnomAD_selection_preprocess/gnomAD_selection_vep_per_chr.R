#!/usr/bin/env Rscript

library(data.table)
library(rtracklayer)
library(GenomicRanges)

# Parse chromosome number from command line
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript gnomAD_selection_vep_per_chr.R <chr_num>")
chr_num <- as.integer(args[1])
chr_str <- paste0("chr", chr_num)
cat(paste0("Processing ", chr_str, "\n"))

# Paths
gnomad_path   <- "../../data/gnomAD_genomes_v3/"
roulette_path <- "../../data/gnomAD_roulette_predictions/"
phylop_path   <- "../../data/zoonomia_phylop/data_download/241-mammalian-2020v2.phylop-Homo_sapiens.bigWig"
output_path   <- "../../data/gnomAD_vep_summaries/"

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

# VEP consequence severity ranking (1 = most severe)
vep_severity <- setNames(1L:22L, c(
	"splice_acceptor_variant", "splice_donor_variant", "stop_gained",
	"stop_lost", "start_lost", "missense_variant", "splice_region_variant",
	"incomplete_terminal_codon_variant", "synonymous_variant",
	"stop_retained_variant", "coding_sequence_variant", "mature_miRNA_variant",
	"5_prime_UTR_variant", "3_prime_UTR_variant",
	"non_coding_transcript_exon_variant", "intron_variant",
	"non_coding_transcript_variant", "upstream_gene_variant",
	"downstream_gene_variant", "TF_binding_site_variant",
	"regulatory_region_variant", "intergenic_variant"
))

vep_display <- c(
	"splice_acceptor_variant" = "splice acceptor variant",
	"splice_donor_variant" = "splice donor variant",
	"stop_gained" = "stop gained",
	"stop_lost" = "stop lost",
	"start_lost" = "start lost",
	"missense_variant" = "missense variant",
	"splice_region_variant" = "splice region variant",
	"incomplete_terminal_codon_variant" = "incomplete terminal codon variant",
	"synonymous_variant" = "synonymous variant",
	"stop_retained_variant" = "stop retained variant",
	"coding_sequence_variant" = "coding sequence variant",
	"mature_miRNA_variant" = "mature miRNA variant",
	"5_prime_UTR_variant" = "5 prime UTR variant",
	"3_prime_UTR_variant" = "3 prime UTR variant",
	"non_coding_transcript_exon_variant" = "non coding transcript exon variant",
	"intron_variant" = "intron variant",
	"non_coding_transcript_variant" = "non coding transcript variant",
	"upstream_gene_variant" = "upstream gene variant",
	"downstream_gene_variant" = "downstream gene variant",
	"TF_binding_site_variant" = "TF binding site variant",
	"regulatory_region_variant" = "regulatory region variant",
	"intergenic_variant" = "intergenic variant"
)

# Determine worst VEP consequence per variant
worst_consequence <- function(cons_strings) {
	result <- character(length(cons_strings))

	# Handle empty/missing
	empty <- is.na(cons_strings) | cons_strings == "." | cons_strings == ""
	result[empty] <- "other consequence"

	# Fast path: single term (no comma)
	valid <- !empty
	single <- valid & !grepl(",", cons_strings, fixed = TRUE)
	if (any(single)) {
		display <- vep_display[cons_strings[single]]
		result[single] <- fifelse(is.na(display), "other consequence", as.character(display))
	}

	# Slow path: multiple terms — find worst
	multi <- valid & !single
	if (any(multi)) {
		result[multi] <- vapply(strsplit(cons_strings[multi], ",", fixed = TRUE), function(terms) {
			ranks <- vep_severity[terms]
			ranks[is.na(ranks)] <- 999L
			best_term <- terms[which.min(ranks)]
			display <- vep_display[best_term]
			if (is.na(display)) "other consequence" else as.character(display)
		}, character(1))
	}

	result
}

# Load gnomAD subinfo with VEP consequences
# awk filters to SNVs and extracts unique consequence terms per variant
cat("Loading gnomAD subinfo with VEP\n")
gnomad <- fread(
	cmd = paste0(
		"gunzip -cd ", gnomad_path,
		"gnomad.genomes.v3.1.2.sites.", chr_str, ".subinfo.vcf.gz",
		" | grep -v '^#'",
		" | awk -F'\\t' 'BEGIN{OFS=\"\\t\"}{",
		"if(length($4)!=1||length($5)!=1)next;",
		"n=split($8,a,\";\");",
		"ac=\".\";an=\".\";af=\".\";vep=\".\";",
		"for(i=1;i<=n;i++){",
		"nkv=split(a[i],kv,\"=\");",
		"if(kv[1]==\"AC\")ac=kv[2];",
		"if(kv[1]==\"AN\")an=kv[2];",
		"if(kv[1]==\"AF\")af=kv[2];",
		"if(kv[1]==\"vep\"){vep=kv[2];for(k=3;k<=nkv;k++)vep=vep\"=\"kv[k]}}",
		"nt=split(vep,transcripts,\",\");",
		"delete seen;cons=\"\";",
		"for(j=1;j<=nt;j++){",
		"split(transcripts[j],fields,\"|\");",
		"ns=split(fields[2],csq,\"&\");",
		"for(k=1;k<=ns;k++){",
		"if(csq[k]!=\"\"&&!(csq[k] in seen)){",
		"seen[csq[k]]=1;",
		"if(cons==\"\")cons=csq[k];else cons=cons\",\"csq[k]}}}",
		"print $1,$2,$4,$5,$7,ac,an,af,cons}'"
	),
	header = FALSE,
	col.names = c("chr", "pos", "ref", "alt", "FILTER", "AC", "AN", "AF", "vep_consequences")
)
gnomad[, `:=`(AC = as.integer(AC), AN = as.integer(AN), AF = as.numeric(AF))]
add_funnel(1L, "gnomAD SNVs loaded", nrow(gnomad))

# QC filter
gnomad <- gnomad[FILTER == "PASS" & AC > 0L & AN >= 76156L]
add_funnel(2L, "gnomAD after QC filter", nrow(gnomad))

# Parse worst VEP consequence
cat("Parsing VEP worst consequence\n")
gnomad[, vep_consequence := worst_consequence(vep_consequences)]
gnomad[, vep_consequences := NULL]

n_other <- sum(gnomad$vep_consequence == "other consequence")
cat(paste0("  Variants with 'other consequence': ", n_other, " / ", nrow(gnomad), "\n"))
cat("  VEP consequence distribution:\n")
print(gnomad[, .N, by = vep_consequence][order(-N)])

# Compute MAF and MAC
gnomad[, `:=`(
  MAF = pmin(AF, 1 - AF),
  MAC = pmin(AC, AN - AC)
)]

# Classify on MAF/MAC
gnomad[, af_class := fcase(
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
add_funnel(3L, "Roulette loaded", nrow(roulette))

# Merge with Roulette
cat("Merging with Roulette\n")
merged <- merge(gnomad, roulette, by = c("chr", "pos", "ref", "alt"), all.x = TRUE)
add_funnel(4L, "Roulette coverage", sum(!is.na(merged$MR)))
rm(gnomad, roulette)
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

# Classify mutation type
merged[, mut_class := classify_mutation(ref, alt, PN)]

# Summary 1a: vep_consequence x mut_class
cat("Computing summary 1a (vep_consequence x mut_class)\n")
summary1a <- merged[, .(
	n              = .N,
	n_MR           = sum(!is.na(MR)),
	sum_MR         = sum(MR, na.rm = TRUE),
	sum_MR_sq      = sum(MR^2, na.rm = TRUE),
	n_phyloP       = sum(!is.na(phyloP_score)),
	sum_phyloP     = sum(phyloP_score, na.rm = TRUE),
	sum_phyloP_sq  = sum(phyloP_score^2, na.rm = TRUE),
	n_conserved    = sum(is_conserved)
), by = .(vep_consequence, mut_class)]

# Summary 1b: vep_consequence
cat("Computing summary 1b (vep_consequence)\n")
summary1b <- merged[, .(
	n              = .N,
	n_MR           = sum(!is.na(MR)),
	sum_MR         = sum(MR, na.rm = TRUE),
	sum_MR_sq      = sum(MR^2, na.rm = TRUE),
	n_phyloP       = sum(!is.na(phyloP_score)),
	sum_phyloP     = sum(phyloP_score, na.rm = TRUE),
	sum_phyloP_sq  = sum(phyloP_score^2, na.rm = TRUE),
	n_conserved    = sum(is_conserved)
), by = .(vep_consequence)]

# Summary 2a: vep_consequence x AF x mut_class
cat("Computing summary 2a (vep_consequence x AF x mut_class)\n")
summary2a <- merged[, .(
	n              = .N,
	n_MR           = sum(!is.na(MR)),
	sum_MR         = sum(MR, na.rm = TRUE),
	sum_MR_sq      = sum(MR^2, na.rm = TRUE),
	n_phyloP       = sum(!is.na(phyloP_score)),
	sum_phyloP     = sum(phyloP_score, na.rm = TRUE),
	sum_phyloP_sq  = sum(phyloP_score^2, na.rm = TRUE),
	n_conserved    = sum(is_conserved)
), by = .(vep_consequence, af_class, mut_class)]

# Summary 2b: vep_consequence x AF
cat("Computing summary 2b (vep_consequence x AF)\n")
summary2b <- merged[, .(
	n              = .N,
	n_MR           = sum(!is.na(MR)),
	sum_MR         = sum(MR, na.rm = TRUE),
	sum_MR_sq      = sum(MR^2, na.rm = TRUE),
	n_phyloP       = sum(!is.na(phyloP_score)),
	sum_phyloP     = sum(phyloP_score, na.rm = TRUE),
	sum_phyloP_sq  = sum(phyloP_score^2, na.rm = TRUE),
	n_conserved    = sum(is_conserved)
), by = .(vep_consequence, af_class)]

# Summary 3: totals
cat("Computing summary 3 (totals)\n")
summary3 <- merged[, .(
	n              = .N,
	n_MR           = sum(!is.na(MR)),
	sum_MR         = sum(MR, na.rm = TRUE),
	sum_MR_sq      = sum(MR^2, na.rm = TRUE),
	n_phyloP       = sum(!is.na(phyloP_score)),
	sum_phyloP     = sum(phyloP_score, na.rm = TRUE),
	sum_phyloP_sq  = sum(phyloP_score^2, na.rm = TRUE),
	n_conserved    = sum(is_conserved)
)]

# Write per-chromosome outputs
fwrite(summary1a, paste0(output_path, "vep_consequence_by_mutclass_", chr_str, ".tsv"), sep = "\t")
fwrite(summary1b, paste0(output_path, "vep_consequence_", chr_str, ".tsv"), sep = "\t")
fwrite(summary2a, paste0(output_path, "vep_consequence_by_af_mutclass_", chr_str, ".tsv"), sep = "\t")
fwrite(summary2b, paste0(output_path, "vep_consequence_by_af_", chr_str, ".tsv"), sep = "\t")
fwrite(summary3, paste0(output_path, "vep_totals_", chr_str, ".tsv"), sep = "\t")
fwrite(funnel, paste0(output_path, "funnel_vep_", chr_str, ".tsv"), sep = "\t")

cat(paste0("Done: ", chr_str, "\n"))
cat(paste0("  Summary 1a rows: ", nrow(summary1a), "\n"))
cat(paste0("  Summary 1b rows: ", nrow(summary1b), "\n"))
cat(paste0("  Summary 2a rows: ", nrow(summary2a), "\n"))
cat(paste0("  Summary 2b rows: ", nrow(summary2b), "\n"))
cat(paste0("  Summary 3 rows:  ", nrow(summary3), "\n"))
