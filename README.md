# TBLAB_STOmics

## Introduction
TBLAB_STOMics is a bioinformatics R script for the analysis, visualization and annotation of spatial transcriptomics data using Seurat. It allows its visualization, clustering and annotation of clusters to cell types using a mannually curated dataset.

## Pipeline summary
The script performs the following major steps:
1. Data loading: load 10X Genomics spatial stranscriptomics data.
2. Preprocessing: normalize, scale, and reduce dimensionality of data
3. Clustering: identify clusters.
4. Marker discovery: detect differentially expressed genes per cluster.
5. Annotation: assign cell types using reference datasets with increasing specificity (D1 - D4).
6. Visualization: generate annotated UMAP and spatial plots.
7. Output: save results including plots, top markers, and cell type assignments.


## Prerequisites
Have R installed, along with the following packages:
- library(Seurat)
- library(SeuratData)
- library(ggplot2)
- library(patchwork)
- library(dplyr)
- library(presto)
- library(devtools)


## Usage
