# 🧬 Visium Spatial RNA-seq Analysis

**End-to-end spatial transcriptomics analysis using 10x Genomics Visium data
and Seurat (R).**

This repository documents a **complete, reproducible workflow** for processing,
analyzing, and interpreting Visium spatial RNA-sequencing datasets, from raw
outputs to biologically interpretable spatial patterns.

The emphasis is on **best practices used in research laboratories**, with
transparent analytical rationale rather than black-box automation.

---

## 📌 Scope of This Repository

This workflow covers:

- Quality control of spatial spots
- Normalization and feature selection
- Dimensionality reduction and clustering
- Identification of spatially variable genes
- Visualization of clusters and gene expression in tissue context
- Optional integration with scRNA-seq references
- Generation of summary tables for figures and manuscripts

📄 **Detailed analytical methods are provided in [`METHODS.md`](METHODS.md).**

---

## 🧪 Technologies & Tools

- **Seurat (R)** — spatial data processing and visualization  
- **10x Genomics Visium** — spatial transcriptomics platform  
- **tidyverse** — data manipulation and reporting  
- **patchwork / ggplot2** — visualization  

---

## 📂 Repository Structure

```text
visium-spatial-rnaseq-analysis/
├── 00_setup/                      # Environment & package setup
├── 01_data_import/                # Import Visium data
├── 02_quality_control/             # Spot-level QC
├── 03_normalization/               # Normalization & feature selection
├── 04_dimensionality_reduction/    # PCA / UMAP
├── 05_clustering_annotation/       # Clustering of spatial spots
├── 06_spatial_features/            # Spatially variable gene detection
├── 07_spatial_visualization/       # Tissue-level plots
├── 08_integration_scRNA/           # scRNA-seq reference mapping (optional)
├── 09_reporting/                   # Summary tables & exports
├── METHODS.md
├── README.md
└── LICENSE

---

🧠 Interpretation Notes

Spatial transcriptomics data represent spot-level mixtures, not single cells.

Cluster and reference-mapping results should be interpreted as enriched
transcriptional programs, not definitive cell identities.

Parameters may require tuning for different tissues or experimental designs.

---

🔁 Reproducibility

This repository reflects real analysis pipelines used in research settings.
Raw sequencing data and sensitive metadata are intentionally excluded.

---

👩‍🔬 Author

Naghmeh Rezaei
Computational Biology · Spatial & Single-Cell Genomics
GitHub: https://github.com/NaghmehRezaei
