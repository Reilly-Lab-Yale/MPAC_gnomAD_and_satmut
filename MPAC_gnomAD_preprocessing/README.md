# Identifying non-coding variant effects at scale via machine learning models of cis-regulatory reporter assays

Code used to preprocess gnomAD predictions into intermediate files and summary tables. Preprint available at X. Publication available at Y.

Contact corresponding authors or me (mackenzie.noon@yale.edu) with questions about these scripts. 

# Purpose

In sum, these scripts take genomic data from a variety of sources, use them to annotate a set of genetic variants, and produce a "count table" which records the number of times every posssible combination of categorical values occurs. Quantative annotations are thresholded to produce categorical variables or averaged.

The generation of the "count table" facilitates rapid exploratory analysis and graphing in downstream analysis. You can find the code for the downstream analysis in the other half of this repository (`MPAC_gnomAD_and_satmut_analyses`).

# Preliminary notes


Manual checks
- This analysis uses many genomic data formats, which inconsistently use 0 and 1 based coordinates. In all nontrivial integrations of multiple formats, manual verification is performed.
- Manual, genome-browser verification of several steps of filtering procedure is included in a few "manual_verification" directories, to avoid any off-by-one errors or similar.

Misc.
- Sometimes presence of variants within a bin is stored in a one-hot fashion : multiple true/false columns, one for each bin. This is not generally desirable, so I convert away with "find_true_column" other times, this is desirable as it is not actually one hot, but multi-hot. Such as presence in genomic region (i.e. some variants might be both in a promoter like sequence & some other category)
- Variants present on non-autosomes are not carried through the analysis.
- This code could be optimized a lot with the use of intervals trees or a similar structure for range-joins & lookups.
- pyspark is used, despite the fact that the computing cluster these analyses were performed on does not support Yarn. As a compromise, Slurm arrays are used throughout, often creating multiple spark instacnces (usually one for each chromosome). This is not terribly efficient : if you are running the code yourself I would ideally recommend using Yarn.
- Multiple slurm "wrapper" scripts are included, each creating separate arrays of separate sets of chromosomes. For example, the large chromosomes (chromosomes 1-3) usually required additional resources, and so needed their own Slurm job. Adjust to your cluster's resources. 

# Steps

## 0. merge 

Uses `bedtools merge`[^d] to combine gnomad variants annotated with Malinois predictions with variants annotated with a subset of gnomad annotations as both sets of annotations are useful for subsequent steps. 

## 1.0.format_zoonomia_phylop

Pre-processes phylo P scores[^b] into a format usable by `bcftools annotate`[^c].

## 2.0.annotate

Removes variants which do **not** meet one or more of the following criteria:
- Variant has non-null AF, AC, AN values
    - (See section "Gnomad" below, or the gnomad website for an explanation of these abbriviations.)
- Variant AN>=76156
    - (See aforementioned "Gnomad" section)
- Has a non-null CADD score
    - (See section "CADD" below).
- Passes all of gnomad's own filters. 

Annotates all variants with
- Which (if any) genomic regions they fall into ("genomic regions" here are retrieved from Encode, see Encode below) in a multi-hot format.
- PhyloP scores[^b] formatted in the previous step.

Adds several columns for donwstream use, computed from existing columns:
- Mean malinois reference activity and mean malinois skew
- MAF : Minor Allele Frequency
    - (See aforementioned "Gnomad" section)
- Rarity codes (simply called "category") and based on AC and MAF.


Note on rarity code `MAF_OR_AC_IS_ZERO` : all variants bearing this code are filtered later steps. 

## 2.2.add_roulette

Adds roulette mutation rate[^r]. 

## 2.3.add_transposons

Annotates all variants with a boolean value that indicates whether they are in a repetitive element or not. 

## 2.5.filter

Removes variants present in exons and splice-sites. 

## 3.0.pleio_and_filter

Computes emVar status and "pleiotropy" : That is, how many cell types (of the three with MPRA results predicted by malinois) is the variant predicted to be an emVar in. Takes values 0-3.

Removes variants called as "MAF_OR_AC_IS_ZERO"

## 3.5.add_TF_footprints

Annotates variants with TF footprints[^tf]

## 3.6.remove_non_SNP

Makes sure all variants are clear AGTC->AGTC SNVs.

## helper

A few additional scripts not directly involved in the generation of the main count table. See the readme files in their respective directories for more information.


# How to reproduce 

- All dependencies are either standard GNU command-line utilities (e.g., zcat, sed, awk, grep...) commonly available on Linux systems, or managed via the Conda environment (including packages installed via pip within the environment).
- The exact environment is recorded in environment.yml (generated via conda env export > environment.yml).
- For convenient installation, packages explicitly installed by the user are listed in environment_minimal.yml and requirements.txt. These are generated with:
    - `conda env export --from-history > environment_minimal.yml`
    - `$CONDA_PREFIX/bin/pip freeze > requirements.txt`
- You can see `conda activate mcn_varef` in the batch files used to submit the relevant Slurm jobs.


- These scripts uses file generated by other analysis... HERE
    - Uses gnomad VCF annotated with malinois activity predictions ...
    - gnomad from ...
    - malinous activity predictions produced in ...
    - The filtering step used masking produced in ...
  



# Helpful excerpts

(The excerpts below help explain some abbriviations used by the various imported data). 


## Gnomad

(from tooltips on gnomad website)

AC : Allele Count : "Alternate allele count in high-quality genotypes"
AN : Allele Number : "Total high quality genotypes"
AF : Allele Frequency : "Alternate Allele Frequency in high quality genotypes"

AF=AC/AN

I define : MAF (minor allele frequency) = min(AF,1-AF)

gnomad filters 1) filter = pass, 2) maf != 0, 3) AN >= 76,156
- 76,156 is [the number of high-quality genomes in 3.1](https://gnomad.broadinstitute.org/news/2020-10-gnomad-v3-1/)
- Therefore, if a site were called in all individuals it would be 76,156\*2 . Gnomad sometimes issues the following warning

> Warning this variant is covered in fewer than 50% of individuals in gnomAD v3.1.2 genomes. This may indicate a low-quality site.

Example : [rs113653250](https://gnomad.broadinstitute.org/variant/1-434284-T-G?dataset=gnomad_r3)


## CADD

>To simplify interpretation in some contexts, we also defined phred-like scores (scaled C  scores) on the basis of the rank of the C score of each variant relative  to all 8.6 billion possible SNVs, ranging from 1 to 99 (Supplementary  Note). For example, substitutions with the highest 10% (10−1) of all  scores—that is, those least likely to be observed human alleles under  our model—were assigned values of 10 or greater (‘≥C10’), whereas  variants  in  the  highest  1%  (10^−2),  0.1%  (10^−3),  etc.  were  assigned  scores ‘≥C20’, ‘≥C30’, etc

[^a]


## ensembl VEP

Consequence codes taken from [here](https://useast.ensembl.org/info/genome/variation/prediction/predicted_data.html) on 2024-06-10.


## Encode
https://screen.encodeproject.org/

>1. cCREs with promoter-like signatures (cCRE-PLS) fall within 200 bp (center to center) of an annotated GENCODE TSS and have high DNase and H3K4me3 signals (evaluated as DNase and H3K4me3 max-Z scores, defined as the maximal DNase or H3K4me3 Z scores across all biosamples with data; see Methods).
>2. cCREs with enhancer-like signatures (cCRE-ELS) have high DNase and H3K27ac max-Z scores and must additionally have a low H3K4me3 max-Z score if they are within 200 bp of an annotated TSS. The subset of cCREs-ELS within 2 kb of a TSS is denoted proximal (cCRE-pELS), while the remaining subset is denoted distal (cCRE-dELS).
>3. DNase-H3K4me3 cCREs have high H3K4me3 max-Z scores but low H3K27ac max-Z scores and do not fall within 200 bp of a TSS.
>4. CTCF-only cCREs have high DNase and CTCF max-Z scores and low H3K4me3 and H3K27ac max-Z scores. 




[^b]: Evolutionary constraint and innovation across hundreds of placental mammals. Science 380, eabn3943 (2023).
[^c]: Twelve years of SAMtools and BCFtools Petr Danecek, James K Bonfield, Jennifer Liddle, John Marshall, Valeriu Ohan, Martin O Pollard, Andrew Whitwham, Thomas Keane, Shane A McCarthy, Robert M Davies, Heng Li GigaScience, Volume 10, Issue 2, February 2021, giab008, https://doi.org/10.1093/gigascience/giab008
[^a]: A general framework for estimating the relative pathogenicity of human genetic variants, Nature Genetics, 2014
[^d]: BEDTools: a flexible suite of utilities for comparing genomic features https://doi.org/10.1093/bioinformatics/btq033
[^r]: A mutation rate model at the basepair resolution identifies the mutagenic effect of polymerase III transcription. Nat Genet 55, 2235–2242 (2023).
[^tf]: Global reference mapping of human transcription factor footprints. Nature 583, 729–736 (2020).

