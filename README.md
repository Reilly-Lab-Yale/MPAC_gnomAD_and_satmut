# Identifying non-coding variant effects at scale via machine learning models of cis-regulatory reporter assays

Preprint available at https://www.biorxiv.org/content/10.1101/2025.04.16.648420v1. MPAC predictions available at https://doi.org/10.5281/zenodo.15178434. See https://github.com/Reilly-Lab-Yale/coda_mpac for MPAC itself. See https://github.com/john-c-butts/MPAC/ for code which produces figures (Maintained by John Butts). See https://github.com/Reilly-Lab-Yale/predictions_plusplus for supplemental predictions.

Everything necessary is in the main branch *except* for indel prediction filtering, which is performed with an alternate version of the  MPAC_gnomAD_preprocessing pipeline in the `corrected_indel` branch.

### MPAC_gnomAD_preprocessing
Contains contains code used to preprocess gnomAD predictions into intermediate files and summary tables. See README there for details. "Helper" directory also contains analyses for many small tasks reference elsehwere. Maintained by Mackenzie Noon.

### MPAC_gnomAD_and_satmut_analyses
Contains code used to analyze and visualize data related to the gnomAD and promoter saturation mutagenesis work. See README therein for details. Maintained by Stephen Rong.
