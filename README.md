# C-A-Landers-Thesis

## Description

This repository contains only commands and scripts modified from or not provided by the original pipeline, which can be found [here](https://github.com/CGSbioinfo/RNASeq_pipeline) 

Its purpose is to provide commands and scripts not available elsewhere and to highlight modifications made by the author \(C A Landers\). All additional, explanatory code comments have been added either during construction of code or to aid in its modification - large sections of some scripts have thus been left uncommented. 

## Comments
Scripts and commands related to the same step of analysis are provided in the same document

Each set of scripts/commands is provided with the following information
1. Designer/s : Individual/s who designed the original scripts/commands
2. Description : Basic description of the purpose of the scripts/commands
3. Indications of the purpose of individual scripts/commands where more than one is provided in a single document

All other information is provided within code comments, including which sections/lines of code were modified by the author \(C A Landers\)

All code is within the Python environment unless otherwise indicated

## SessionInfo

R
```R
R version 3.4.1 (2017-06-30)
Platform: x86_64-redhat-linux-gnu (64-bit)
Running under: CentOS release 6.8 (Final)

Matrix products: default
BLAS: /usr/lib64/R/lib/libRblas.so
LAPACK: /usr/lib64/R/lib/libRlapack.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C
 [9] LC_ADDRESS=C               LC_TELEPHONE=C
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C

attached base packages:
[1] parallel  stats4    stats     graphics  grDevices utils     datasets
[8] methods   base

other attached packages:
 [1] calibrate_1.7.2            MASS_7.3-47
 [3] RColorBrewer_1.1-2         rtracklayer_1.36.4
 [5] gplots_3.0.1               DESeq2_1.16.1
 [7] SummarizedExperiment_1.6.5 DelayedArray_0.2.7
 [9] matrixStats_0.52.2         Biobase_2.36.2
[11] GenomicRanges_1.28.5       GenomeInfoDb_1.12.2
[13] IRanges_2.10.3             S4Vectors_0.14.5
[15] BiocGenerics_0.22.0        ggplot2_2.2.1

loaded via a namespace (and not attached):
 [1] bit64_0.9-7              splines_3.4.1            gtools_3.5.0
 [4] Formula_1.2-2            latticeExtra_0.6-28      blob_1.1.0
 [7] GenomeInfoDbData_0.99.0  Rsamtools_1.28.0         RSQLite_2.0
[10] backports_1.1.1          lattice_0.20-35          digest_0.6.12
[13] XVector_0.16.0           checkmate_1.8.4          colorspace_1.3-2
[16] htmltools_0.3.6          Matrix_1.2-11            plyr_1.8.4
[19] XML_3.98-1.9             genefilter_1.58.1        zlibbioc_1.22.0
[22] xtable_1.8-2             scales_0.5.0             gdata_2.18.0
[25] BiocParallel_1.10.1      htmlTable_1.9            tibble_1.3.4
[28] annotate_1.54.0          nnet_7.3-12              lazyeval_0.2.0
[31] survival_2.41-3          magrittr_1.5             memoise_1.1.0
[34] foreign_0.8-69           tools_3.4.1              data.table_1.10.4
[37] stringr_1.2.0            munsell_0.4.3            locfit_1.5-9.1
[40] cluster_2.0.6            AnnotationDbi_1.38.2     Biostrings_2.44.2
[43] compiler_3.4.1           caTools_1.17.1           rlang_0.1.2
[46] grid_3.4.1               RCurl_1.95-4.8           htmlwidgets_0.9
[49] bitops_1.0-6             base64enc_0.1-3          labeling_0.3
[52] gtable_0.2.0             DBI_0.7                  GenomicAlignments_1.12.2
[55] gridExtra_2.3            knitr_1.17               bit_1.1-12
[58] Hmisc_4.0-3              KernSmooth_2.23-15       stringi_1.1.5
[61] Rcpp_0.12.13             geneplotter_1.54.0       rpart_4.1-11
[64] acepack_1.4.1
>

```


