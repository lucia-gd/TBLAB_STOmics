# TBLAB_STOmics

## Introduction
TBLAB_STOMics is a bioinformatics R script for the analysis, visualization and annotation of spatial transcriptomics data using Seurat. It allows its visualization, clustering and annotation of clusters to cell types using a mannually curated dataset.

## Pipeline summary
The script performs the following major steps:
1. Data loading: load 10X Genomics spatial stranscriptomics data.
2. Preprocessing: normalize, scale, and reduce dimensionality of data
3. Clustering: identify clusters.
4. Marker discovery: detect differentially expressed genes per cluster.
5. Annotation: assign cell types using reference datasets with decreasing specificity (D1 - D4).
6. Visualization: generate annotated UMAP and spatial plots.
7. Output: save results including plots, top markers, and cell type assignments.


## Get started
### Prerequisites
Have R installed, along with the following packages:
- library(Seurat)
- library(SeuratData)
- library(ggplot2)
- library(patchwork)
- library(dplyr)
- library(presto)
- library(devtools)


### Required Data
- 10x spatial transcriptomics data
- Annotation reference CSVs: D1.csv, D2.csv, D3.csv, D4.csv; each containing at least:
  - symbol (gene symbol)
  - cell_name (cell type)
  - Additional metadata (optional but recommended for extended annotation):
    - tissue_type, cancer_type and tissue_class


 ## Usage
 Rscript seurat_integrated_all.R <sample> <data_dir> <n.markers> <max.dimensions> <clustering.resolution>

 ## Example:
 Rscript seurat_integrated_all.R A694Tumor /path/to/outs 10 10 0.5

 ## Script arguments:


 ## Citation
 


    
