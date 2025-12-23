# ============================================================
# Quality Control for Visium Spatial RNA-seq Data
# ============================================================

# This script performs spot-level quality control on a
# Seurat spatial object, including inspection of library
# size and detected feature distributions.

# ------------------------------------------------------------
# Load libraries
# ------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(patchwork)

# ------------------------------------------------------------
# Load raw spatial object
# ------------------------------------------------------------
spatial_obj <- readRDS(
  "01_data_import/visium_seurat_object_raw.rds"
)

# ------------------------------------------------------------
# Calculate QC metrics
# ------------------------------------------------------------
# Percentage of mitochondrial transcripts
# (mitochondrial genes typically start with "MT-")

spatial_obj[["percent.mt"]] <- PercentageFeatureSet(
  spatial_obj,
  pattern = "^MT-"
)

# ------------------------------------------------------------
# Visualize QC metrics
# ------------------------------------------------------------
# Examine distributions to guide filtering thresholds

qc_plots <- VlnPlot(
  spatial_obj,
  features = c("nCount_Spatial", "nFeature_Spatial", "percent.mt"),
  ncol = 3,
  pt.size = 0
)

qc_plots

# ------------------------------------------------------------
# Spot filtering
# ------------------------------------------------------------
# Filtering thresholds should be guided by distributions
# rather than fixed values; conservative defaults are shown.

spatial_obj_filtered <- subset(
  spatial_obj,
  subset =
    nCount_Spatial > 500 &
    nFeature_Spatial > 200 &
    percent.mt < 20
)

# ------------------------------------------------------------
# Summary after filtering
# ------------------------------------------------------------
cat("Spots before QC:", ncol(spatial_obj), "\n")
cat("Spots after QC:", ncol(spatial_obj_filtered), "\n")

# ------------------------------------------------------------
# Save filtered object
# ------------------------------------------------------------
saveRDS(
  spatial_obj_filtered,
  file = "02_quality_control/visium_seurat_object_qc.rds"
)
