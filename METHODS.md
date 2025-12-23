# Methods Appendix — Visium Spatial RNA-seq Analysis

## Overview

This document describes the analytical workflow used to process and analyze
10x Genomics Visium spatial transcriptomics data. The methods emphasize
reproducibility, analytical transparency, and biological interpretability,
and are written in a manuscript-style format rather than as a software manual.

---

## 1. Data Import

Raw Visium outputs, including filtered feature-barcode matrices and
associated histology images, were imported into R using Seurat’s spatial
data infrastructure. Spatial coordinates and imaging metadata were retained
for downstream tissue-aware analyses.

Scripts: `01_data_import/import_visium_data.R`

---

## 2. Quality Control

Quality control was performed at the **spot level**. Low-quality spots were
identified based on total transcript counts, detected gene numbers, and
mitochondrial read proportions. Thresholds were selected conservatively to
minimize removal of biologically relevant tissue regions.

Scripts: `02_quality_control/spatial_qc.R`

---

## 3. Normalization and Feature Selection

Data were normalized using log-normalization. Highly variable genes were
identified to capture dominant sources of biological variation while reducing
technical noise. Normalization was applied uniformly across all spatial spots.

Scripts: `03_normalization/normalize_spatial_data.R`

---

## 4. Dimensionality Reduction

Principal component analysis (PCA) was performed using variable genes, followed
by non-linear dimensionality reduction (UMAP) for visualization. The number of
principal components retained was selected based on variance explained and
biological interpretability.

Scripts: `04_dimensionality_reduction/run_umap.R`

---

## 5. Clustering of Spatial Spots

Unsupervised clustering was performed in reduced dimensional space to identify
transcriptionally distinct spatial domains. Clustering resolution was chosen
to balance granularity with interpretability.

Scripts: `05_clustering_annotation/cluster_spatial_domains.R`

---

## 6. Identification of Spatially Variable Genes

Spatially variable genes were identified using spatial autocorrelation-based
methods implemented in Seurat. These genes highlight transcriptional programs
with non-random spatial organization across the tissue section.

Scripts: `06_spatial_features/find_spatial_features.R`

---

## 7. Spatial Visualization

Clusters and gene expression patterns were projected back onto tissue images
to enable spatial interpretation. Visualization focused on identifying
anatomical structure, regional specialization, and spatial gradients.

Scripts: `07_spatial_visualization/plot_spatial_results.R`

---

## 8. Integration with scRNA-seq Reference (Optional)

To support biological interpretation, a processed scRNA-seq reference dataset
was optionally used to map transcriptional signatures onto spatial spots using
anchor-based transfer. Transferred labels were interpreted as **probabilistic
enrichment**, not single-cell assignments.

Scripts: `08_integration_scRNA/map_scRNA_reference.R`

---

## 9. Reporting and Output Generation

Final outputs included cluster composition summaries, spatial marker gene
tables, and exported metadata suitable for figure generation and supplementary
materials. All results were saved in spreadsheet-compatible formats.

Scripts: `09_reporting/generate_summary_tables.R`

---

## Notes on Interpretation and Reproducibility

Spatial transcriptomics data represent aggregated transcriptional signals from
multiple cells within each spot. All analyses were designed to reflect this
resolution. Scripts may require parameter tuning for different tissues,
platform versions, or experimental designs.
