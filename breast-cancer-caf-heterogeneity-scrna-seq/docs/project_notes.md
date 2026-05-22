# Project Notes

## Project Title

Decoding Cancer-Associated Fibroblast (CAF) Heterogeneity in Breast Cancer Metastasis Using Single-Cell RNA Sequencing

## Biological Question

How do Cancer-Associated Fibroblast populations vary across breast cancer samples, and which molecular programs may be associated with aggressive tumor phenotypes?

## Main Analysis Focus

The workflow focuses on:

- Preprocessing breast cancer scRNA-seq data
- Identifying major tumor microenvironment cell populations
- Extracting fibroblast/CAF cells
- Re-clustering CAFs to explore CAF heterogeneity
- Comparing TNBC and ER-positive groups using differential expression analysis
- Linking significant genes to biological pathways
- Preparing outputs for cell-cell communication interpretation

## Differential Expression Comparison

The selected comparison in this version is:

```text
TNBC vs ER-positive
```

This comparison is used instead of primary vs metastatic comparison.

## Important Marker Genes

### Fibroblast / CAF markers

- COL1A1
- COL1A2
- DCN
- LUM
- FAP
- ACTA2

### Epithelial / tumor markers

- EPCAM
- KRT8
- KRT18
- KRT19
- MUC1

### Immune markers

- PTPRC
- CD3D
- CD3E
- MS4A1
- LYZ
- NKG7

### Endothelial markers

- PECAM1
- VWF
- KDR
- CLDN5

## Expected Figures

- QC violin plot before filtering
- QC violin plot after filtering
- PCA plot
- Elbow plot
- UMAP by sample
- Clustered UMAP
- Annotated UMAP
- Marker DotPlot
- Fibroblast/CAF UMAP
- CAF subtype marker DotPlot
- Volcano plot for TNBC vs ER-positive comparison
- DE heatmap
- GO/KEGG enrichment dotplots
- CellChat network plots if CellChat is installed

## Notes for Reproducibility

The script saves all outputs automatically inside an `analysis_outputs` folder. Large raw data files and generated RDS files should usually not be uploaded to GitHub unless required, because they can be too large.
