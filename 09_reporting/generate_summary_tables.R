# ============================================================
# Reporting & Summary Tables for Visium Spatial RNA-seq
# ============================================================

# This script generates summary tables suitable for
# downstream visualization, interpretation, and reporting.

# ------------------------------------------------------------
# Load libraries
# ------------------------------------------------------------
library(Seurat)
library(tidyverse)

# ------------------------------------------------------------
# Load final spatial object
# ------------------------------------------------------------
spatial_obj <- readRDS(
  "08_integration_scRNA/visium_seurat_object_with_reference.rds"
)

# ------------------------------------------------------------
# Cluster composition summary
# ------------------------------------------------------------
cluster_summary <- spatial_obj@meta.data %>%
  count(seurat_clusters) %>%
  rename(
    cluster = seurat_clusters,
    n_spots = n
  )

write_csv(
  cluster_summary,
  "09_reporting/cluster_spot_counts.csv"
)

# ------------------------------------------------------------
# Marker gene identification (spatial context)
# ------------------------------------------------------------
markers <- FindAllMarkers(
  spatial_obj,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

write_csv(
  markers,
  "09_reporting/spatial_cluster_markers.csv"
)

# ------------------------------------------------------------
# Export Seurat metadata for record keeping
# ------------------------------------------------------------
metadata_out <- spatial_obj@meta.data %>%
  rownames_to_column("spot_id")

write_csv(
  metadata_out,
  "09_reporting/spatial_metadata.csv"
)

# ------------------------------------------------------------
# Notes
# ------------------------------------------------------------
# Outputs generated here are intended for:
# - figure annotation
# - supplementary tables
# - downstream biological interpretation
