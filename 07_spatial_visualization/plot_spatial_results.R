# ============================================================
# Spatial Visualization for Visium Spatial RNA-seq
# ============================================================

# This script visualizes clustering results and spatially
# variable genes in the context of tissue images.

# ------------------------------------------------------------
# Load libraries
# ------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(patchwork)

# ------------------------------------------------------------
# Load spatial feature–annotated object
# ------------------------------------------------------------
spatial_obj <- readRDS(
  "06_spatial_features/visium_seurat_object_spatial_features.rds"
)

# ------------------------------------------------------------
# Spatial visualization of clusters
# ------------------------------------------------------------
# Display clusters overlaid on histology image

p_clusters <- SpatialDimPlot(
  spatial_obj,
  label = TRUE,
  label.size = 3
)

p_clusters

# ------------------------------------------------------------
# Spatial visualization of selected genes
# ------------------------------------------------------------
# Visualize spatial expression of representative genes

genes_to_plot <- head(
  SpatiallyVariableFeatures(spatial_obj),
  4
)

p_genes <- SpatialFeaturePlot(
  spatial_obj,
  features = genes_to_plot,
  ncol = 2
)

p_genes

# ------------------------------------------------------------
# Save plots (optional)
# ------------------------------------------------------------
ggsave(
  filename = "07_spatial_visualization/spatial_clusters.png",
  plot = p_clusters,
  width = 6,
  height = 6
)

ggsave(
  filename = "07_spatial_visualization/spatial_genes.png",
  plot = p_genes,
  width = 8,
  height = 6
)
