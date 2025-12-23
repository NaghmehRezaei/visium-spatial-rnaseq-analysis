# ============================================================
# Normalization and Feature Selection for Visium Spatial RNA-seq
# ============================================================

# This script performs normalization and identification of
# variable features on a quality-controlled Visium Seurat object.

# ------------------------------------------------------------
# Load libraries
# ------------------------------------------------------------
library(Seurat)
library(tidyverse)

# ------------------------------------------------------------
# Load QC-filtered spatial object
# ------------------------------------------------------------
spatial_obj <- readRDS(
  "02_quality_control/visium_seurat_object_qc.rds"
)

# ------------------------------------------------------------
# Normalize expression data
# ------------------------------------------------------------
# Log-normalization is commonly used for Visium spatial data
# to stabilize variance across spots.

spatial_obj <- NormalizeData(
  spatial_obj,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

# ------------------------------------------------------------
# Identify highly variable features
# ------------------------------------------------------------
# Variable genes are used for dimensionality reduction
# and clustering.

spatial_obj <- FindVariableFeatures(
  spatial_obj,
  selection.method = "vst",
  nfeatures = 2000
)

# ------------------------------------------------------------
# Inspect variable features
# ------------------------------------------------------------
top_variable_genes <- head(
  VariableFeatures(spatial_obj),
  20
)

print(top_variable_genes)

# ------------------------------------------------------------
# Save normalized object
# ------------------------------------------------------------
saveRDS(
  spatial_obj,
  file = "03_normalization/visium_seurat_object_norm.rds"
)
