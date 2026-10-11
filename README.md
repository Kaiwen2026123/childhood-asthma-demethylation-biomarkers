# childhood-asthma-demethylation-biomarkers
R code for integrative bioinformatics and single-cell RNA-seq analysis of demethylation-associated transcriptomic biomarkers in childhood asthma

## Project Overview
This repository contains reproducible R scripts for integrative analysis combining single-cell RNA-seq and DNA methylation data from peripheral blood mononuclear cells (PBMC) of children with asthma and control subjects.
The workflow covers raw data preprocessing, quality control, cell annotation, differential expression, differentially methylated region (DMR) analysis, functional enrichment, visualization, machine-learning biomarker screening and molecular docking.

> ⚠️ Note: Raw sequencing files, large Seurat `.rds` objects, methylation matrices and intermediate large files are **not uploaded to this repository** due to file size limits. Users need to prepare raw data locally before running these scripts.

## Script Workflow
| Script | Description |
|---|---|
| `01_data_preprocess.R` | scRNA-seq raw matrix import, quality control, doublet filtering, Harmony batch correction, cell type annotation, monocyte subclustering and SingleR reference-based annotation |
| `02_DE_analysis.R` | Cell-type-specific differential gene expression analysis (Asthma vs Control) |
| `03_dmr_analysis.R` | Differentially methylated region (DMR) identification from methylation array data |
| `04_go_kegg_enrich.R` | GO / KEGG functional enrichment analysis of DEGs and DMR-associated genes |
| `05_visualization.R` | Figure generation: UMAP, dot plot, heatmap, volcano plot, FeaturePlot with publication-style theme |
| `06_machine_learning.R` | Construction of machine-learning model to screen biomarker gene signature |
| `07_molecular_docking.R` | Preparation and analysis for molecular docking |

## R Environment & Required Packages
R ≥ 4.3
```r
# Core packages
library(Seurat)
library(tidyverse)
library(patchwork)
library(harmony)
library(scDblFinder)
library(SingleR)
library(celldex)
library(limma)
library(minfi)
library(clusterProfiler)
library(ggsci)
library(glmnet)
