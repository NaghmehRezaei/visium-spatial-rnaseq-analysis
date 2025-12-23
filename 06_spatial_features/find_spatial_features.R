# ============================================================
# Identification of Spatially Variable Features (Visium)
# ============================================================

# This script identifies genes with spatially structured
# expression patterns across tissue using Seurat's
# spatial feature detection methods.

# ------------------------------------------------------------
# Load libraries
# ------------------------------------------------------------
library(Seurat)
library(tidyverse)

# ------------------------------------------------------------
# Load clustered spatial object
# ------------------------------------------------------------
spatial_obj <- readRDS(
  "05_clustering/visium_seurat_object_clustered.rds"
)

# ------------------------------------------------------------
# Identify spatially variable features
# ------------------------------------------------------------
# This method accounts for spatial autocorrelation
# across neighboring spots.

spatial_obj <- FindSpatiallyVariableFeatures(
  spatial_obj,
  assay = "Spatial",
  selection.method = "markvariogram"
)

# ------------------------------------------------------------
# Inspect top spatially variable genes
# ------------------------------------------------------------
top_spatial_genes <- head(
  SpatiallyVariableFeatures(spatial_obj),
  20
)

print(top_spatial_genes)

# ------------------------------------------------------------
# Save updated object
# ------------------------------------------------------------
saveRDS(
  spatial_obj,
  file = "06_spatial_features/visium_seurat_object_spatial_features.rds"
)
