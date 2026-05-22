# Decoding Cancer-Associated Fibroblast (CAF) Heterogeneity in Breast Cancer Metastasis Using Single-Cell RNA Sequencing

This project analyzes single-cell RNA sequencing data from breast cancer samples to investigate tumor microenvironment heterogeneity, with a focus on Cancer-Associated Fibroblasts (CAFs) and their potential role in breast cancer progression and metastasis.

## Project Overview

Cancer-Associated Fibroblasts are important stromal cells within the tumor microenvironment. They can contribute to extracellular matrix remodeling, inflammatory signaling, immune modulation, and tumor progression. This project uses a single-cell RNA-seq workflow to process breast cancer data, identify major cell populations, extract fibroblast/CAF populations, re-cluster CAFs, and perform differential expression analysis.

## Dataset

- Dataset: GSE176078
- Disease context: Breast cancer
- Data type: Single-cell RNA sequencing
- Analysis platform: R / Seurat

## Main Workflow

1. Load 10X-format count matrix using `Read10X()`
2. Create a Seurat object
3. Add cell-level metadata
4. Perform quality control and filtering
5. Normalize gene expression data
6. Identify highly variable genes
7. Run PCA for dimensionality reduction
8. Apply Harmony for batch correction
9. Run UMAP and clustering
10. Identify cluster marker genes
11. Annotate major cell types using known marker genes
12. Extract fibroblast/CAF cells
13. Re-cluster CAF populations
14. Identify CAF subtype marker genes
15. Run differential expression analysis comparing TNBC vs ER-positive groups
16. Perform GO and KEGG pathway enrichment
17. Generate CellChat-ready outputs and optional cell-cell communication analysis

## Team Contributions

### Shimaa Mohamed
- Dataset loading and preprocessing
- Quality control
- Normalization and highly variable gene selection
- PCA, UMAP, clustering
- Differential expression analysis
- Output figures: QC plots, PCA plot, UMAP plots, clustered UMAP, volcano plot, and DE heatmap

### Esraa Mohamed
- Cluster annotation
- Fibroblast identification using marker genes
- CAF extraction and re-clustering
- CAF subtype identification
- Biological interpretation of DE results
- Output figures: annotated UMAP, marker DotPlot, fibroblast UMAP, CAF subtype UMAP, and marker heatmap

### Menatalla Essam
- Cell-cell communication analysis using CellChat
- Ligand-receptor interaction analysis
- Pathway enrichment interpretation
- Linking DE genes to biological pathways
- Output figures: network plots, bubble plots, and pathway enrichment plots

## Repository Structure

```text
breast-cancer-caf-heterogeneity-scrna-seq/
├── README.md
├── scripts/
│   └── scRNAseq_GSE176078_full_analysis.R
└── docs/
    └── project_notes.md
```

## How to Run

1. Download the GSE176078 expression matrix and metadata.
2. Place the data in one local folder.
3. Open the R script in RStudio.
4. Change only the `data_dir` path at the top of the script:

```r
data_dir <- "D:/Bioinformatics_Project/GSE176078"
```

5. Run the full script.
6. Output files will be saved automatically in:

```text
analysis_outputs/
```

## Main Outputs

- QC violin plots
- QC scatter plots
- Highly variable gene plot
- PCA plot
- Elbow plot
- UMAP by sample
- Clustered UMAP
- Annotated UMAP
- Marker DotPlot
- Fibroblast/CAF UMAP
- CAF marker tables
- TNBC vs ER-positive DE table
- Volcano plot
- DE heatmap
- GO enrichment results
- KEGG enrichment results
- Optional CellChat communication outputs

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

The main differential expression comparison in this version is TNBC vs ER-positive breast cancer groups. If the metadata uses different column names or labels, update the `de_group_col`, `de_group_1`, and `de_group_2` settings at the top of the script.
