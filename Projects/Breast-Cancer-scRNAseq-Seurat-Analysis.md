# Single-Cell RNA-seq Analysis of Breast Cancer Tumor Microenvironment

## Project Overview

This project focuses on analyzing single-cell RNA sequencing data from human breast cancer samples to explore cellular heterogeneity within the tumor microenvironment. The analysis was performed using the **Seurat R package** and aimed to identify major cell populations, compare tumor subtypes, and investigate differentially expressed genes between breast cancer groups such as **TNBC** and **ER+** samples.

The workflow included quality control, normalization, feature selection, dimensionality reduction, clustering, cell-type annotation, and differential gene expression analysis.

## Objectives

- Analyze breast cancer single-cell RNA-seq data using Seurat.
- Perform quality control and preprocessing of raw gene expression data.
- Identify highly variable genes and reduce dimensionality using PCA and UMAP.
- Cluster cells based on transcriptomic similarity.
- Annotate major cell populations using known marker genes.
- Compare gene expression patterns between breast cancer subtypes.
- Identify differentially expressed genes that may contribute to tumor progression and subtype-specific biology.

## Dataset

The project used publicly available breast cancer single-cell RNA-seq data from GEO. The dataset contains gene expression profiles from individual cells, allowing detailed investigation of tumor heterogeneity and the tumor microenvironment.

## Workflow

### 1. Data Loading

- Imported 10X-compatible gene expression matrices.
- Created a Seurat object for downstream analysis.
- Added sample metadata to the Seurat object.

### 2. Quality Control

- Filtered low-quality cells based on the number of detected genes.
- Calculated mitochondrial gene percentage.
- Removed poor-quality cells and potential outliers.

### 3. Normalization and Feature Selection

- Normalized gene expression data.
- Identified highly variable genes for downstream analysis.

### 4. Scaling and Dimensionality Reduction

- Scaled gene expression values.
- Performed PCA to summarize major sources of variation.
- Used UMAP for 2D visualization of cellular populations.

### 5. Clustering

- Constructed a nearest-neighbor graph.
- Clustered cells using Seurat clustering methods.
- Tested clustering resolution to obtain biologically meaningful groups.

### 6. Cell-Type Annotation

- Annotated clusters using known marker genes.
- Identified major cell populations within the tumor microenvironment.

### 7. Differential Expression Analysis

- Compared gene expression between selected groups.
- Identified marker genes and subtype-associated differentially expressed genes.
- Visualized results using volcano plots, feature plots, violin plots, and dot plots.

## Tools and Technologies

- **R**
- **Seurat**
- **Single-cell RNA-seq analysis**
- **PCA**
- **UMAP**
- **Clustering**
- **Differential gene expression analysis**
- **Data visualization**
- **GEO datasets**

## Key Skills Demonstrated

- Bioinformatics data analysis
- Single-cell transcriptomics
- Tumor microenvironment analysis
- Breast cancer subtype comparison
- R programming
- Data preprocessing and quality control
- Statistical analysis of gene expression data
- Biological interpretation of computational results
- Scientific visualization and reporting

## Project Outcome

This project provided insight into the cellular composition of breast cancer samples and highlighted transcriptional differences between tumor subtypes. The analysis demonstrated how single-cell RNA-seq can be used to study tumor heterogeneity, identify cell populations, and explore genes associated with cancer progression and subtype-specific behavior.

## Repository Structure

```text
├── data/                  # Input data and metadata
├── scripts/               # R scripts for preprocessing and analysis
├── results/               # Output files, marker genes, and figures
├── figures/               # UMAPs, volcano plots, violin plots, and other visualizations
└── README.md              # Project documentation
```

## Author

**Shaimaa Mohamed El Haddad**  
Biomedical Science Student | Computational Biology & Genomics Track

---

## GitHub Repository Description

Single-cell RNA-seq analysis of breast cancer tumor microenvironment using Seurat in R, including quality control, normalization, PCA, UMAP, clustering, cell-type annotation, and differential gene expression analysis between breast cancer subtypes.
