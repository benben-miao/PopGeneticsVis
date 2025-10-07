#!/bin/bash

cat > install_packages.R << 'EOF'
cran_pkgs <- c("argparse", "CMplot", "dplyr", "tidyr", "stringr", "data.table", "ggplot2", "ggpubr")
new_cran <- cran_pkgs[!(cran_pkgs %in% installed.packages()[,"Package"])]
if (length(new_cran) > 0) install.packages(new_cran, repos = "https://cran.rstudio.com/")

if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
bio_pkgs <- c("methylKit", "GenomicFeatures", "rtracklayer", "GenomicRanges")
new_bio <- bio_pkgs[!(bio_pkgs %in% installed.packages()[,"Package"])]
if (length(new_bio) > 0) BiocManager::install(new_bio)
EOF

R --slave < install_packages.R
rm install_packages.R
