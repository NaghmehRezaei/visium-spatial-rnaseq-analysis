# ============================================================
# Setup Environment for Visium Spatial RNA-seq Analysis
# ============================================================

# This script documents the R package environment used for
# sequencing-based spatial transcriptomics analysis with Seurat.
# It is intended to be run once per environment.

# ------------------------------------------------------------
# CRAN packages
# ------------------------------------------------------------
install.packages(c(
  "Seurat",
  "tidyverse",
  "patchwork",
  "ggplot2",
  "cowplot"
))

# ------------------------------------------------------------
# Bioconductor packages
# ------------------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c(
  "SpatialExperiment",
  "SingleCellExperiment",
  "scater",
  "scran"
))

# ------------------------------------------------------------
# Version check (recommended)
# ------------------------------------------------------------
sessionInfo()
