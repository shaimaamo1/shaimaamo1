# Decoding Cancer-Associated Fibroblast Heterogeneity in the Breast Cancer Microenvironment Using Single-Cell RNA Sequencing

## Project Overview

This project investigates **Cancer-Associated Fibroblast (CAF) heterogeneity** in the breast cancer tumor microenvironment using single-cell RNA sequencing (scRNA-seq). Breast tumors contain multiple interacting cell populations, not only malignant epithelial cells, and scRNA-seq enables transcriptomic analysis at single-cell resolution.

The project focused on identifying major cell populations, isolating CAFs, re-clustering CAF populations, annotating CAF subtypes, and exploring CAF functional diversity through differential expression analysis and cell–cell communication analysis.

The analysis was performed using **R**, **Seurat**, **Harmony**, and **CellChat**.

## Research Aim

The aim of this project was to investigate breast cancer CAF heterogeneity using scRNA-seq data and analyze CAF function through:

- Major cell-type annotation
- CAF isolation and re-clustering
- CAF subtype annotation
- Differential gene expression analysis
- Cell–cell communication analysis
- Biological interpretation of CAF-related signaling pathways

## Dataset

The project used the publicly available human breast cancer scRNA-seq dataset **GSE176078**.

Dataset characteristics:

- Human breast cancer single-cell RNA-seq data
- 26 primary breast tumor samples
- Breast cancer subtypes including **ER+**, **HER2+**, and **TNBC**
- Input files included sparse count matrix, barcodes, genes, and metadata
- Initial loaded dataset contained **100,064 cells**

## Analysis Workflow

```text
Raw Data
↓
Quality Control
↓
Filtering
↓
Normalization
↓
Highly Variable Gene Selection
↓
Scaling
↓
PCA
↓
Harmony Batch Correction
↓
UMAP Visualization
↓
Clustering
↓
Cluster Marker Identification
↓
Major Cell Type Annotation
↓
CAF Subsetting
↓
CAF Re-clustering
↓
CAF Subtype Annotation
↓
Differential Expression Analysis
↓
CellChat Communication Analysis
↓
Biological Interpretation
```

## Methods

### 1. Data Loading

The scRNA-seq count matrix was imported into R and used to create a Seurat object. Sample metadata were added to support downstream subtype comparison and biological interpretation.

### 2. Quality Control and Filtering

Quality control was performed to remove low-quality cells and possible technical noise.

QC metrics included:

- `nFeature_RNA`
- `nCount_RNA`
- `percent.mito`

Filtering criteria:

- `nFeature_RNA > 300`
- `nFeature_RNA < 6000`
- `percent.mito < 10%`

After filtering, **81,820 cells** were retained from the original **100,064 cells**.

### 3. Normalization, Highly Variable Genes, and PCA

Gene expression data were normalized using log normalization. Highly variable genes were selected to capture genes with high cell-to-cell expression variability. The normalized and scaled data were then used for principal component analysis to summarize major sources of transcriptomic variation.

### 4. Batch Correction, UMAP, and Clustering

Harmony batch correction was included in the workflow to reduce sample-level technical variation. UMAP was used to visualize cells in two-dimensional space, and graph-based clustering identified **12 transcriptionally distinct clusters** in the breast cancer scRNA-seq dataset.

### 5. Differential Expression Analysis: TNBC vs ER+

Differential expression analysis was performed to compare gene expression patterns between **TNBC** and **ER+** cells.

Subtype-associated differentially expressed genes included:

- Higher in TNBC: `CYBA`, `PFN1`, `CORO1A`, `RAC2`
- Higher in ER+: `AGR3`, `MALAT1`

Results were visualized using volcano plots and heatmaps to highlight significant subtype-associated expression patterns.

### 6. Major Cell-Type Annotation

Cluster-specific marker genes were identified using `FindAllMarkers()`. Canonical marker genes were then used to annotate major cell populations in the tumor microenvironment.

Annotated cell populations included:

- T cells
- Epithelial cells
- Myeloid cells
- Endothelial cells
- CAFs
- Pericytes
- Plasma cells
- B cells
- Proliferating cells
- Basal epithelial cells
- pDCs
- Unknown cells

Marker expression was validated using DotPlot and FeaturePlot visualizations.

### 7. CAF Subsetting and Re-clustering

CAF populations were isolated from the annotated Seurat object and reanalyzed independently to investigate fibroblast heterogeneity in greater detail.

CAF-specific workflow:

```text
Annotated Seurat Object
↓
CAF Subsetting
↓
Re-normalization
↓
PCA
↓
UMAP
↓
CAF Re-clustering
↓
CAF Subtype Annotation
```

CAF re-clustering revealed **5 transcriptionally distinct CAF populations**.

### 8. CAF Subtype Annotation

CAF subtypes were annotated using known CAF marker genes.

| CAF Subtype | Marker Genes |
|---|---|
| myCAF | ACTA2, TAGLN |
| iCAF | CXCL12, IL6 |
| matrixCAF | COL1A1, FN1 |
| apCAF | HLA-DRA, CD74 |
| proliferativeCAF | MKI67, TOP2A |

DotPlot and FeaturePlot visualizations were used to validate CAF subtype marker expression across re-clustered CAF populations.

### 9. Cell–Cell Communication Analysis Using CellChat

Cell–cell communication analysis was performed using **CellChat** to investigate signaling interactions between CAF subtypes.

The analysis included:

- Creating a CellChat object from the CAF dataset
- Loading the human CellChat ligand–receptor database
- Detecting overexpressed communication genes
- Identifying overexpressed ligand–receptor interactions
- Computing communication probabilities
- Aggregating interaction networks
- Interpreting pathway-level communication patterns

## Key Results

- The project identified major cell populations within the breast cancer tumor microenvironment.
- The initial dataset contained **100,064 cells**, and **81,820 cells** were retained after quality control.
- UMAP clustering identified **12 transcriptionally distinct clusters**.
- Differential expression analysis highlighted subtype-associated genes between **TNBC** and **ER+** cells.
- CAFs were successfully isolated and re-clustered into **5 CAF subtypes**.
- CAF subtype annotation identified **myCAF**, **iCAF**, **matrixCAF**, **apCAF**, and **proliferativeCAF** populations.
- CellChat analysis revealed active communication between CAF subtypes.
- **COLLAGEN**, **MK**, and **FN1** signaling pathways showed strong information flow.
- **myCAF** and **apCAF** displayed stronger outgoing signaling activity compared with other CAF subtypes.
- These findings suggest that CAF subtypes may contribute differently to extracellular matrix remodeling and tumor microenvironment regulation.

## Key Figures and Visual Outputs

The project generated and interpreted several biological visualizations, including:

- QC violin plots before and after filtering
- Highly variable gene plot
- UMAP clustering plot
- TNBC vs ER+ volcano plot
- Top DE gene heatmap
- Marker-gene DotPlots
- Marker-gene FeaturePlots
- Annotated major cell-type UMAP
- CAF re-clustering UMAP
- CAF subtype UMAP
- CAF marker DotPlot and FeaturePlot
- CellChat signaling and communication network outputs

## Tools and Technologies

- R
- Seurat
- Harmony
- CellChat
- GEO dataset analysis
- Single-cell RNA sequencing analysis
- PCA and UMAP
- Graph-based clustering
- Differential gene expression analysis
- Marker gene-based annotation
- Ligand–receptor communication analysis
- Scientific visualization

## Skills Demonstrated

- Single-cell RNA-seq data analysis
- Breast cancer tumor microenvironment analysis
- Cancer-associated fibroblast subtype identification
- R programming for bioinformatics
- Quality control and preprocessing of scRNA-seq data
- Dimensionality reduction and clustering
- Marker gene-based cell-type annotation
- Differential gene expression analysis
- Cell–cell communication analysis using CellChat
- Biological interpretation of transcriptomic and signaling results
- Scientific visualization and presentation

## Project Limitation

One limitation of this study is that the dataset included **primary breast cancer samples only**. Therefore, metastatic CAF comparison could not be performed.

## Conclusion

This project demonstrates a complete single-cell transcriptomics workflow for investigating CAF heterogeneity in breast cancer. By combining Seurat-based clustering, CAF subtype annotation, differential expression analysis, and CellChat communication analysis, the study provides insight into the functional diversity of CAF populations and their potential roles in regulating the breast cancer tumor microenvironment.

## Repository Structure

```text
├── data/                  # Input data and metadata
├── scripts/               # R scripts for preprocessing and analysis
├── results/               # Output files, marker genes, and CellChat results
├── figures/               # UMAPs, volcano plots, heatmaps, DotPlots, FeaturePlots
├── presentation/          # Final project presentation
└── README.md              # Project documentation
```

## Author

**Shaimaa Mohamed El Haddad**  
Biomedical Science Student | Computational Biology & Genomics Track

---

## GitHub Repository Description

Single-cell RNA-seq analysis of breast cancer CAF heterogeneity using Seurat, Harmony, and CellChat to annotate CAF subtypes, perform differential expression analysis, and explore cell–cell communication within the tumor microenvironment.
