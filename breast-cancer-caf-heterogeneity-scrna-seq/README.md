# Decoding CAF Heterogeneity in Breast Cancer Metastasis Using scRNA-seq

[![R](https://img.shields.io/badge/R-4.x-blue)](https://www.r-project.org/)
[![Seurat](https://img.shields.io/badge/Seurat-scRNA--seq-6A5ACD)](https://satijalab.org/seurat/)
[![Project](https://img.shields.io/badge/Project-Bioinformatics-green)]()
[![Status](https://img.shields.io/badge/Status-Completed-brightgreen)]()

## Project Summary

This project investigates **Cancer-Associated Fibroblast (CAF) heterogeneity** in breast cancer using **single-cell RNA sequencing (scRNA-seq)** data. The analysis focuses on identifying major tumor microenvironment cell populations, extracting fibroblast/CAF populations, re-clustering CAFs, and studying transcriptional differences between **TNBC** and **ER-positive** breast cancer groups.

CAF populations are important components of the tumor microenvironment because they can contribute to extracellular matrix remodeling, inflammatory signaling, immune regulation, tumor progression, and metastatic behavior.

## Dataset

| Item | Description |
|---|---|
| Dataset | GSE176078 |
| Disease | Breast cancer |
| Data type | Single-cell RNA sequencing |
| Main tool | Seurat in R |
| Main comparison | TNBC vs ER-positive |

## Research Aim

To decode CAF heterogeneity in breast cancer and identify gene-expression patterns that may distinguish aggressive tumor microenvironment states, especially between TNBC and ER-positive groups.

## Analysis Workflow

```text
Raw scRNA-seq matrix
        ↓
Create Seurat object
        ↓
Quality control and filtering
        ↓
Normalization and highly variable gene selection
        ↓
Scaling and PCA
        ↓
Harmony batch correction
        ↓
UMAP visualization and clustering
        ↓
Cluster marker identification
        ↓
Cell type annotation
        ↓
Fibroblast / CAF extraction
        ↓
CAF re-clustering and subtype marker analysis
        ↓
TNBC vs ER-positive differential expression
        ↓
GO / KEGG pathway enrichment
        ↓
Optional CellChat communication analysis
```

## Repository Structure

```text
breast-cancer-caf-heterogeneity-scrna-seq/
├── README.md
├── requirements.R
├── .gitignore
├── scripts/
│   └── scRNAseq_GSE176078_full_analysis.R
└── docs/
    └── project_notes.md
```

## Main Script

The full analysis workflow is available here:

```text
scripts/scRNAseq_GSE176078_full_analysis.R
```

The script performs:

- Data loading using `Read10X()`
- Seurat object creation
- Metadata integration
- Quality control
- Normalization
- Highly variable gene selection
- PCA and UMAP
- Harmony batch correction
- Clustering
- Marker gene detection
- Cell type annotation
- CAF extraction and re-clustering
- Differential expression analysis
- Volcano plot generation
- Heatmap generation
- GO and KEGG enrichment
- Optional CellChat analysis

## How to Run

### 1. Install packages

Open RStudio and run:

```r
source("requirements.R")
```

### 2. Prepare input data

Place the GSE176078 count matrix and metadata in one local folder. The script supports common 10X-style folders such as:

```text
tenx_for_Read10X_final/
filtered_feature_bc_matrix/
filtered_gene_bc_matrices/
```

### 3. Edit the data path

Open the main script and change only this line:

```r
data_dir <- "D:/Bioinformatics_Project/GSE176078"
```

### 4. Run the analysis

Run the full script in RStudio:

```r
source("scripts/scRNAseq_GSE176078_full_analysis.R")
```

All outputs will be saved automatically in:

```text
analysis_outputs/
```

## Expected Outputs

| Output type | Examples |
|---|---|
| QC figures | Violin plots, scatter plots |
| Dimensionality reduction | PCA, elbow plot, UMAP |
| Clustering | Clustered UMAP, marker genes |
| Annotation | Annotated UMAP, marker DotPlot |
| CAF analysis | CAF UMAP, CAF marker tables |
| Differential expression | TNBC vs ER-positive DE table |
| Visualization | Volcano plot, DE heatmap |
| Pathway analysis | GO and KEGG enrichment tables and dotplots |
| Communication analysis | Optional CellChat network plots |

## Important Marker Genes

| Cell population | Marker genes |
|---|---|
| Epithelial / tumor cells | EPCAM, KRT8, KRT18, KRT19, MUC1 |
| Fibroblasts / CAFs | COL1A1, COL1A2, DCN, LUM, FAP, ACTA2 |
| Immune cells | PTPRC, CD3D, CD3E, MS4A1, LYZ, NKG7 |
| Endothelial cells | PECAM1, VWF, KDR, CLDN5 |
| Pericytes | RGS5, PDGFRB, MCAM, CSPG4 |
| Proliferating cells | MKI67, TOP2A, STMN1 |

## Differential Expression

The main differential expression comparison is:

```text
TNBC vs ER-positive
```

This version uses molecular subtype comparison rather than primary vs metastatic comparison. If metadata labels differ, update these variables at the top of the script:

```r
de_group_col <- "subtype"
de_group_1 <- "TNBC"
de_group_2 <- "ER+"
```

## Team Contributions

### Shimaa Mohamed

- Data loading and preprocessing
- Quality control
- Normalization and highly variable gene selection
- PCA, UMAP, clustering
- Differential expression analysis
- Volcano plot and DE heatmap generation

### Esraa Mohamed

- Cluster annotation
- Fibroblast identification using marker genes
- CAF extraction and re-clustering
- CAF subtype marker interpretation
- Annotated UMAP and marker DotPlot interpretation

### Menatalla Essam

- Cell-cell communication analysis using CellChat
- Ligand-receptor interaction interpretation
- Pathway enrichment interpretation
- Biological interpretation of DE and communication outputs

## Tools and Packages

- R
- Seurat
- Harmony
- ggplot2
- dplyr
- patchwork
- pheatmap
- clusterProfiler
- org.Hs.eg.db
- enrichplot
- CellChat optional

## Notes

Raw sequencing data, large expression matrices, generated `.rds` files, and full output folders are intentionally excluded from GitHub. The repository contains the analysis code, documentation, and reproducibility instructions.

## Project Status

Completed as a bioinformatics course project focused on breast cancer tumor microenvironment analysis using single-cell RNA sequencing.
