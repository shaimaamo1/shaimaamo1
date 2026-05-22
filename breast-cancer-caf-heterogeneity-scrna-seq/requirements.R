# ============================================================
# Package requirements for the CAF scRNA-seq project
# Run this file once before running the main analysis script.
# ============================================================

cran_packages <- c(
  "Seurat",
  "Matrix",
  "ggplot2",
  "patchwork",
  "dplyr",
  "harmony",
  "pheatmap",
  "plyr"
)

bioc_packages <- c(
  "clusterProfiler",
  "org.Hs.eg.db",
  "enrichplot"
)

optional_packages <- c(
  "CellChat",
  "NMF",
  "ggalluvial"
)

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

for (pkg in cran_packages) {
  install_if_missing(pkg)
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

message("Core package installation finished.")
message("CellChat is optional and may require separate installation depending on your R version.")
