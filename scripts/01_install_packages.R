# scripts/01_install_packages.R

# ================================
# Install required packages
# ================================

options(timeout = 600)

cran_pkgs <- c(
  "tidyverse",
  "data.table",
  "janitor",
  "pheatmap",
  "ggrepel",
  "pROC"
)

bioc_pkgs <- c(
  "GEOquery",
  "limma",
  "edgeR",
  "DESeq2",
  "AnnotationDbi",
  "org.Hs.eg.db",
  "clusterProfiler",
  "enrichplot"
)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

cat("Package installation/check completed.\n")