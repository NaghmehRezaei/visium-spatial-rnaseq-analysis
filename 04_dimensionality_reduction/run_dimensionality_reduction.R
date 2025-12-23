# ============================================================
# Dimensionality Reduction for Visium Spatial RNA-seq
# ============================================================

# This script performs scaling, PCA, and UMAP on a normalized
# Visium Seurat object to enable clustering and visualization.

# ------------------------------------------------------------
# Load libraries
# ------------------------------------------------------------
library(Seurat)
library(tidyverse)

# ------------------------------------------------------------
# Load normalized spatial object
# ------------------------------------------------------------
spatial_obj <- readRDS(
  "03_normalization/visium_seurat_object_norm.rds"
)

# ------------------------------------------------------------
# Scale data
# ------------------------------------------------------------
# Scale only variable features to reduce memory usage
# and focus on informative genes.

spatial_obj <- ScaleData(
  spatial_obj,
  features = VariableFeatures(spatial_obj),
  verbose = FALSE
)

# ------------------------------------------------------------
# Principal Component Analysis (PCA)
# ------------------------------------------------------------
spatial_obj <- RunPCA(
  spatial_obj,
  features = VariableFeatures(spatial_obj),
  verbose = FALSE
)

# Inspect PCA results
print(spatial_obj[["pca"]])

# ------------------------------------------------------------
# Determine dimensionality
# ------------------------------------------------------------
# Visual inspection of variance explained is commonly used.

ElbowPlot(spatial_obj, ndims = 50)

# ------------------------------------------------------------
# Run UMAP
# ------------------------------------------------------------
# Use a conservative number of PCs based on PCA inspection.

spatial_obj <- RunUMAP(
  spatial_obj,
  dims = 1:30,
  verbose = FALSE
)

# ------------------------------------------------------------
# Save object
# ------------------------------------------------------------
saveRDS(
  spatial_obj,
  file = "04_dimensionality_reduction/visium_seurat_object_dr.rds"
)
