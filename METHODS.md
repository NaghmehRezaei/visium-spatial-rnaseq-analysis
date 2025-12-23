\# Methods Appendix: Spatial Transcriptomics Analysis (Seurat)



\## Overview



This document describes analytical workflows used for the processing

and interpretation of sequencing-based spatial transcriptomics data,

such as 10x Genomics Visium. The methods emphasize analytical rationale,

reproducibility, and biological interpretation.



---



\## 1. Data Import and Object Initialization



Spatial transcriptomics data were loaded into Seurat objects using

platform-compatible input formats. Spot-level gene expression matrices,

spatial coordinates, and associated histology images were incorporated

to enable joint analysis of transcriptional profiles and tissue context.



Initial inspection of spot-level metadata was performed to verify data

integrity prior to downstream analysis.



---



\## 2. Quality Control and Spot Filtering



Quality control metrics were evaluated at the spatial spot level,

including library size and detected feature counts. Low-quality spots

were filtered using conservative thresholds to minimize technical noise

while retaining biologically relevant spatial information.



Quality control decisions were guided by metric distributions rather

than fixed cutoffs.



---



\## 3. Normalization and Variance Stabilization



Expression values were normalized to account for differences in

sequencing depth across spots. Log-normalization or variance-stabilizing

approaches were applied to enable meaningful comparison of expression

patterns across the tissue section.



Highly variable features were identified for downstream dimensionality

reduction and clustering.



---



\## 4. Dimensionality Reduction and Clustering



Principal component analysis was performed on variable genes to reduce

dimensionality. Lower-dimensional embeddings were used for clustering

of spatial spots, with clustering resolution selected based on cluster

stability and interpretability.



Dimensionality reduction results were evaluated alongside spatial

coordinates to assess correspondence between transcriptional and

anatomical structure.



---



\## 5. Identification of Spatially Variable Features



Spatially variable genes were identified using methods implemented in

Seurat that account for spatial autocorrelation across neighboring

spots. Genes exhibiting regionally restricted expression patterns were

prioritized for biological interpretation.



Spatial feature analysis enabled the identification of tissue domains

with distinct transcriptional programs.



---



\## 6. Spatial Visualization



Gene expression patterns were visualized in the context of tissue

images using spatial feature plots. Both individual marker genes and

gene sets were examined to assess spatial organization of biological

processes.



Visualization parameters were adjusted to balance clarity and spatial

resolution.



---



\## 7. Integration with Single-Cell RNA-seq (Optional)



When appropriate, spatial data were compared with single-cell RNA-seq

references to aid cell-type interpretation. Integration was performed

using anchor-based approaches to map single-cell transcriptional

signatures onto spatial locations.



This step was used to support biological interpretation rather than to

assign definitive cell identities.



---



\## 8. Reporting and Reproducibility



Final outputs included spatial expression plots, clustering summaries,

and annotated result tables. Scripts are organized to enable stepwise

execution and reproducibility across datasets and computing

environments.



