# ============================================================
# Clustering for Visium Spatial RNA-seq
# ============================================================

# This script performs graph-based clustering on a Visium
# spatial Seurat object using PCA-derived embeddings.

# ------------------------------------------------------------
# Load libraries
# ------------------------------------------------------------
library(Seurat)
library(tidyverse)

# ------------------------------------------------------------
# Load dimensionality-reduced object
# ------------------------------------------------------------
spatial_obj <- readRDS(
  "04_dimensionality_reduction/visium_seurat_object_dr.rds"
)

# ------------------------------------------------------------
# Construct nearest-neighbor graph
# ------------------------------------------------------------
spatial_obj <- FindNeighbors(
  spatial_obj,
  dims = 1:30,
  verbose = FALSE
)

# ------------------------------------------------------------
# Perform clustering
# ------------------------------------------------------------
# Resolution controls cluster granularity; values between
# 0.2 and 0.8 are commonly explored for Visium data.

spatial_obj <- FindClusters(
  spatial_obj,
  resolution = 0.4,
  verbose = FALSE
)

# ------------------------------------------------------------
# Inspect clustering results
# ------------------------------------------------------------
table(Idents(spatial_obj))

# ------------------------------------------------------------
# Save clustered object
# ------------------------------------------------------------
saveRDS(
  spatial_obj,
  file = "05_clustering/visium_seurat_object_clustered.rds"
)
