# ============================================================
# Breast Cancer scRNA-seq Project - GSE176078
# Full workflow: preprocessing, clustering, annotation, CAF analysis,
# differential expression, pathway enrichment, and CellChat-ready outputs.
# ============================================================

# -----------------------------
# 0) User settings
# -----------------------------
data_dir <- "D:/Bioinformatics_Project/GSE176078"   # CHANGE THIS ONLY
out_dir  <- file.path(data_dir, "analysis_outputs")

# Differential expression comparison
# We changed the DE comparison to molecular subtype comparison.
# Main comparison used here: TNBC vs ER-positive.
de_group_col <- "subtype"        # common alternatives: "Subtype", "molecular_subtype", "ER_status"
de_group_1   <- "TNBC"
de_group_2   <- "ER+"

# If your metadata writes ER-positive differently, add allowed names here
tnbc_names <- c("TNBC", "Triple-negative", "Triple negative", "Triple_Negative")
erpos_names <- c("ER+", "ER-positive", "ER positive", "ER_pos", "ERpos", "Luminal", "Luminal A", "Luminal B")

set.seed(1234)

# -----------------------------
# 1) Packages
# -----------------------------
cran_packages <- c("Seurat", "Matrix", "ggplot2", "patchwork", "dplyr", "harmony", "pheatmap")
bioc_packages <- c("clusterProfiler", "org.Hs.eg.db", "enrichplot")
optional_packages <- c("CellChat", "NMF", "ggalluvial")

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

for (pkg in cran_packages) install_if_missing(pkg)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(harmony)
  library(pheatmap)
})

has_clusterprofiler <- requireNamespace("clusterProfiler", quietly = TRUE) &&
  requireNamespace("org.Hs.eg.db", quietly = TRUE)

has_cellchat <- requireNamespace("CellChat", quietly = TRUE)

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
fig_dir <- file.path(out_dir, "figures")
marker_dir <- file.path(out_dir, "markers")
rds_dir <- file.path(out_dir, "rds")
pathway_dir <- file.path(out_dir, "pathway_enrichment")
cellchat_dir <- file.path(out_dir, "cellchat_outputs")

dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(marker_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(rds_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(pathway_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(cellchat_dir, showWarnings = FALSE, recursive = TRUE)

save_plot <- function(plot_obj, filename, width = 8, height = 6, dpi = 300) {
  ggsave(file.path(fig_dir, filename), plot = plot_obj, width = width, height = height, dpi = dpi)
}

find_first_existing <- function(paths) {
  paths[file.exists(paths)][1]
}

# -----------------------------
# 2) Load 10X dataset
# -----------------------------
tenx_dir <- find_first_existing(c(
  file.path(data_dir, "tenx_for_Read10X_final"),
  file.path(data_dir, "filtered_feature_bc_matrix"),
  file.path(data_dir, "filtered_gene_bc_matrices"),
  data_dir
))

if (is.na(tenx_dir)) stop("No valid 10X folder found. Check data_dir.")

counts_raw <- Read10X(data.dir = tenx_dir)

counts <- if (is.list(counts_raw)) {
  if ("Gene Expression" %in% names(counts_raw)) {
    counts_raw[["Gene Expression"]]
  } else {
    counts_raw[[1]]
  }
} else {
  counts_raw
}

seurat_obj <- CreateSeuratObject(
  counts = counts,
  project = "GSE176078_BRCA",
  min.cells = 3,
  min.features = 200
)

# -----------------------------
# 3) Add metadata
# -----------------------------
metadata_file <- find_first_existing(c(
  file.path(data_dir, "metadata.csv"),
  file.path(data_dir, "meta_data.csv"),
  file.path(data_dir, "cell_metadata.csv")
))

if (!is.na(metadata_file)) {
  metadata <- read.csv(metadata_file, row.names = 1, check.names = FALSE)
  shared_cells <- intersect(colnames(seurat_obj), rownames(metadata))
  
  if (length(shared_cells) > 0) {
    seurat_obj <- subset(seurat_obj, cells = shared_cells)
    metadata <- metadata[colnames(seurat_obj), , drop = FALSE]
    seurat_obj <- AddMetaData(seurat_obj, metadata)
  } else {
    warning("Metadata was found, but cell names did not match the Seurat object.")
  }
} else {
  warning("No metadata.csv file found. The workflow will continue, but group comparisons may need metadata.")
}

# Standardize some common metadata column names
meta_cols <- colnames(seurat_obj@meta.data)

if (!"sample" %in% meta_cols) {
  possible_sample_cols <- c("Sample", "sample_id", "Patient", "patient", "orig.ident")
  found_sample <- intersect(possible_sample_cols, meta_cols)[1]
  seurat_obj$sample <- seurat_obj@meta.data[[found_sample]]
}

if (!de_group_col %in% colnames(seurat_obj@meta.data)) {
  possible_group_cols <- c("subtype", "Subtype", "molecular_subtype", "Molecular_Subtype",
                           "ER_status", "ER.Status", "condition", "group", "diagnosis")
  found_group <- intersect(possible_group_cols, colnames(seurat_obj@meta.data))[1]
  if (!is.na(found_group)) {
    de_group_col <- found_group
  }
}

# Clean TNBC / ER+ labels if the subtype column exists
if (de_group_col %in% colnames(seurat_obj@meta.data)) {
  group_values <- as.character(seurat_obj@meta.data[[de_group_col]])
  group_values[group_values %in% tnbc_names] <- "TNBC"
  group_values[group_values %in% erpos_names] <- "ER+"
  seurat_obj@meta.data[[de_group_col]] <- group_values
}

# -----------------------------
# 4) Quality control
# -----------------------------
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
seurat_obj[["percent.ribo"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RP[SL]")

qc_before <- VlnPlot(
  seurat_obj,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3,
  pt.size = 0.05
)
save_plot(qc_before, "01_QC_before_filtering.png", width = 12, height = 5)

qc_scatter1 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
qc_scatter2 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
save_plot(qc_scatter1 + qc_scatter2, "02_QC_scatterplots.png", width = 12, height = 5)

seurat_obj_filtered <- subset(
  seurat_obj,
  subset = nFeature_RNA >= 300 &
    nFeature_RNA <= 6000 &
    percent.mt <= 20
)

qc_after <- VlnPlot(
  seurat_obj_filtered,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3,
  pt.size = 0.05
)
save_plot(qc_after, "03_QC_after_filtering.png", width = 12, height = 5)

# -----------------------------
# 5) Normalization and HVGs
# -----------------------------
seurat_obj_filtered <- NormalizeData(
  seurat_obj_filtered,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

seurat_obj_filtered <- FindVariableFeatures(
  seurat_obj_filtered,
  selection.method = "vst",
  nfeatures = 3000
)

top10_hvg <- head(VariableFeatures(seurat_obj_filtered), 10)

hvg_plot <- VariableFeaturePlot(seurat_obj_filtered)
hvg_plot <- LabelPoints(plot = hvg_plot, points = top10_hvg, repel = TRUE)
save_plot(hvg_plot, "04_highly_variable_genes.png", width = 8, height = 6)

# -----------------------------
# 6) Scaling and PCA
# -----------------------------
seurat_obj_filtered <- ScaleData(
  seurat_obj_filtered,
  features = VariableFeatures(seurat_obj_filtered)
)

seurat_obj_filtered <- RunPCA(
  seurat_obj_filtered,
  features = VariableFeatures(seurat_obj_filtered),
  npcs = 50
)

pca_plot <- DimPlot(seurat_obj_filtered, reduction = "pca", group.by = "sample")
save_plot(pca_plot, "05_PCA_by_sample.png", width = 8, height = 6)

elbow_plot <- ElbowPlot(seurat_obj_filtered, ndims = 50)
save_plot(elbow_plot, "06_elbow_plot.png", width = 8, height = 6)

# -----------------------------
# 7) Batch correction with Harmony
# -----------------------------
batch_col <- if ("sample" %in% colnames(seurat_obj_filtered@meta.data)) "sample" else "orig.ident"

seurat_obj_filtered <- RunHarmony(
  object = seurat_obj_filtered,
  group.by.vars = batch_col,
  reduction = "pca",
  dims.use = 1:30,
  assay.use = "RNA"
)

# -----------------------------
# 8) UMAP and clustering
# -----------------------------
seurat_obj_filtered <- RunUMAP(
  seurat_obj_filtered,
  reduction = "harmony",
  dims = 1:30,
  min.dist = 0.3,
  spread = 1
)

seurat_obj_filtered <- FindNeighbors(
  seurat_obj_filtered,
  reduction = "harmony",
  dims = 1:30
)

seurat_obj_filtered <- FindClusters(
  seurat_obj_filtered,
  resolution = 0.5
)

umap_all <- DimPlot(seurat_obj_filtered, reduction = "umap", group.by = "sample")
save_plot(umap_all, "07_UMAP_all_cells_by_sample.png", width = 8, height = 6)

umap_clusters <- DimPlot(seurat_obj_filtered, reduction = "umap", label = TRUE)
save_plot(umap_clusters, "08_clustered_UMAP_before_annotation.png", width = 8, height = 6)

if (de_group_col %in% colnames(seurat_obj_filtered@meta.data)) {
  umap_group <- DimPlot(seurat_obj_filtered, reduction = "umap", group.by = de_group_col)
  save_plot(umap_group, "09_UMAP_by_DE_group.png", width = 8, height = 6)
}

saveRDS(seurat_obj_filtered, file.path(rds_dir, "seurat_after_preprocessing_clustering.rds"))

# -----------------------------
# 9) Marker genes before annotation
# -----------------------------
cluster_markers <- FindAllMarkers(
  seurat_obj_filtered,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

write.csv(cluster_markers, file.path(marker_dir, "01_all_cluster_markers_before_annotation.csv"), row.names = FALSE)

fc_col <- if ("avg_log2FC" %in% colnames(cluster_markers)) "avg_log2FC" else "avg_logFC"

top10_markers <- cluster_markers %>%
  group_by(cluster) %>%
  arrange(desc(.data[[fc_col]]), .by_group = TRUE) %>%
  slice_head(n = 10) %>%
  ungroup()

write.csv(top10_markers, file.path(marker_dir, "02_top10_markers_per_cluster.csv"), row.names = FALSE)

top_marker_genes <- unique(top10_markers$gene)
top_marker_genes <- top_marker_genes[top_marker_genes %in% rownames(seurat_obj_filtered)]

if (length(top_marker_genes) > 2) {
  heatmap_clusters <- DoHeatmap(seurat_obj_filtered, features = top_marker_genes, size = 3) + NoLegend()
  save_plot(heatmap_clusters, "10_marker_heatmap_before_annotation.png", width = 12, height = 8)
}

# -----------------------------
# 10) Cluster annotation using marker genes
# -----------------------------
marker_list <- list(
  Epithelial_Tumor = c("EPCAM", "KRT8", "KRT18", "KRT19", "MUC1"),
  Fibroblast_CAF  = c("COL1A1", "COL1A2", "DCN", "LUM", "FAP", "ACTA2"),
  Immune          = c("PTPRC", "CD3D", "CD3E", "MS4A1", "LYZ", "NKG7"),
  Endothelial     = c("PECAM1", "VWF", "KDR", "CLDN5"),
  Pericyte        = c("RGS5", "PDGFRB", "MCAM", "CSPG4"),
  Proliferating   = c("MKI67", "TOP2A", "STMN1")
)

for (nm in names(marker_list)) {
  marker_list[[nm]] <- marker_list[[nm]][marker_list[[nm]] %in% rownames(seurat_obj_filtered)]
}

valid_markers <- marker_list[lengths(marker_list) > 0]

seurat_obj_filtered <- AddModuleScore(
  seurat_obj_filtered,
  features = valid_markers,
  name = names(valid_markers)
)

score_cols <- grep(paste0("^", names(valid_markers), collapse = "|"), colnames(seurat_obj_filtered@meta.data), value = TRUE)

cluster_scores <- seurat_obj_filtered@meta.data %>%
  mutate(cluster = as.character(seurat_clusters)) %>%
  group_by(cluster) %>%
  summarise(across(all_of(score_cols), mean), .groups = "drop")

annotation_map <- data.frame(
  cluster = cluster_scores$cluster,
  predicted_cell_type = apply(cluster_scores[, score_cols, drop = FALSE], 1, function(x) {
    best <- score_cols[which.max(x)]
    gsub("[0-9]+$", "", best)
  })
)

write.csv(annotation_map, file.path(marker_dir, "03_cluster_annotation_map.csv"), row.names = FALSE)

seurat_obj_filtered$cell_type <- plyr::mapvalues(
  x = as.character(seurat_obj_filtered$seurat_clusters),
  from = annotation_map$cluster,
  to = annotation_map$predicted_cell_type,
  warn_missing = FALSE
)

annotated_umap <- DimPlot(seurat_obj_filtered, reduction = "umap", group.by = "cell_type", label = TRUE, repel = TRUE)
save_plot(annotated_umap, "11_annotated_UMAP_cell_types.png", width = 9, height = 7)

feature_markers <- unique(unlist(marker_list))
feature_markers <- feature_markers[feature_markers %in% rownames(seurat_obj_filtered)]

dot_plot <- DotPlot(seurat_obj_filtered, features = feature_markers, group.by = "cell_type") + RotatedAxis()
save_plot(dot_plot, "12_cell_type_marker_DotPlot.png", width = 12, height = 6)

saveRDS(seurat_obj_filtered, file.path(rds_dir, "seurat_annotated_all_cells.rds"))

# -----------------------------
# 11) Extract fibroblasts / CAFs
# -----------------------------
fibroblast_cells <- colnames(seurat_obj_filtered)[seurat_obj_filtered$cell_type == "Fibroblast_CAF"]

if (length(fibroblast_cells) < 50) {
  warning("Few fibroblast cells were identified by automatic annotation. Check marker expression manually.")
}

caf_obj <- subset(seurat_obj_filtered, cells = fibroblast_cells)

caf_obj <- NormalizeData(caf_obj)
caf_obj <- FindVariableFeatures(caf_obj, selection.method = "vst", nfeatures = 2000)
caf_obj <- ScaleData(caf_obj, features = VariableFeatures(caf_obj))
caf_obj <- RunPCA(caf_obj, features = VariableFeatures(caf_obj), npcs = 30)
caf_obj <- RunUMAP(caf_obj, dims = 1:20)
caf_obj <- FindNeighbors(caf_obj, dims = 1:20)
caf_obj <- FindClusters(caf_obj, resolution = 0.4)

fib_umap <- DimPlot(caf_obj, reduction = "umap", label = TRUE)
save_plot(fib_umap, "13_fibroblast_CAF_UMAP_clusters.png", width = 8, height = 6)

caf_markers <- list(
  myCAF = c("ACTA2", "TAGLN", "MYL9", "TPM2"),
  iCAF  = c("IL6", "CXCL12", "CXCL14", "LIF"),
  apCAF = c("HLA-DRA", "HLA-DRB1", "CD74"),
  ECM_CAF = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM")
)

caf_markers <- lapply(caf_markers, function(x) x[x %in% rownames(caf_obj)])

caf_features <- unique(unlist(caf_markers))
if (length(caf_features) > 1) {
  caf_dot <- DotPlot(caf_obj, features = caf_features) + RotatedAxis()
  save_plot(caf_dot, "14_CAF_subtype_marker_DotPlot.png", width = 10, height = 5)
}

caf_cluster_markers <- FindAllMarkers(
  caf_obj,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

write.csv(caf_cluster_markers, file.path(marker_dir, "04_CAF_cluster_markers.csv"), row.names = FALSE)

saveRDS(caf_obj, file.path(rds_dir, "caf_reclustered_object.rds"))

# -----------------------------
# 12) Differential expression: TNBC vs ER+
# -----------------------------
if (de_group_col %in% colnames(seurat_obj_filtered@meta.data)) {
  
  available_groups <- unique(as.character(seurat_obj_filtered@meta.data[[de_group_col]]))
  
  if (all(c(de_group_1, de_group_2) %in% available_groups)) {
    
    Idents(seurat_obj_filtered) <- de_group_col
    
    de_results <- FindMarkers(
      seurat_obj_filtered,
      ident.1 = de_group_1,
      ident.2 = de_group_2,
      min.pct = 0.10,
      logfc.threshold = 0.25,
      test.use = "wilcox"
    )
    
    de_results$gene <- rownames(de_results)
    write.csv(de_results, file.path(marker_dir, "05_DE_TNBC_vs_ERpositive_all_cells.csv"), row.names = FALSE)
    
    de_results$minus_log10_padj <- -log10(de_results$p_val_adj + 1e-300)
    de_results$direction <- "Not significant"
    de_results$direction[de_results$p_val_adj < 0.05 & de_results[[fc_col]] > 0.25] <- paste0("Higher in ", de_group_1)
    de_results$direction[de_results$p_val_adj < 0.05 & de_results[[fc_col]] < -0.25] <- paste0("Higher in ", de_group_2)
    
    volcano <- ggplot(de_results, aes(x = .data[[fc_col]], y = minus_log10_padj, color = direction)) +
      geom_point(alpha = 0.75, size = 1.2) +
      geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed") +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
      theme_classic() +
      labs(
        title = paste0("Differential expression: ", de_group_1, " vs ", de_group_2),
        x = "Average log2 fold-change",
        y = "-log10 adjusted p-value"
      )
    
    save_plot(volcano, "15_volcano_TNBC_vs_ERpositive.png", width = 8, height = 6)
    
    sig_genes <- de_results %>%
      filter(p_val_adj < 0.05, abs(.data[[fc_col]]) > 0.25) %>%
      arrange(p_val_adj)
    
    top_de_genes <- head(sig_genes$gene, 30)
    top_de_genes <- top_de_genes[top_de_genes %in% rownames(seurat_obj_filtered)]
    
    if (length(top_de_genes) > 2) {
      de_heatmap <- DoHeatmap(seurat_obj_filtered, features = top_de_genes, group.by = de_group_col, size = 3) + NoLegend()
      save_plot(de_heatmap, "16_DE_heatmap_TNBC_vs_ERpositive.png", width = 10, height = 7)
    }
    
  } else {
    warning("TNBC and ER+ labels were not both found in the selected metadata column. Check de_group_col and label names.")
  }
  
} else {
  warning("No valid DE metadata column found. Differential expression step was skipped.")
}

# Optional DE inside CAFs only
if (de_group_col %in% colnames(caf_obj@meta.data)) {
  caf_groups <- unique(as.character(caf_obj@meta.data[[de_group_col]]))
  
  if (all(c(de_group_1, de_group_2) %in% caf_groups)) {
    Idents(caf_obj) <- de_group_col
    
    caf_de <- FindMarkers(
      caf_obj,
      ident.1 = de_group_1,
      ident.2 = de_group_2,
      min.pct = 0.10,
      logfc.threshold = 0.25,
      test.use = "wilcox"
    )
    
    caf_de$gene <- rownames(caf_de)
    write.csv(caf_de, file.path(marker_dir, "06_DE_TNBC_vs_ERpositive_CAFs_only.csv"), row.names = FALSE)
  }
}

# -----------------------------
# 13) Pathway enrichment using significant DE genes
# -----------------------------
if (exists("sig_genes") && nrow(sig_genes) > 5 && has_clusterprofiler) {
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  
  genes_for_enrichment <- sig_genes$gene[1:min(300, nrow(sig_genes))]
  
  gene_map <- bitr(
    genes_for_enrichment,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )
  
  if (nrow(gene_map) > 5) {
    ego <- enrichGO(
      gene = unique(gene_map$ENTREZID),
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.20,
      readable = TRUE
    )
    
    write.csv(as.data.frame(ego), file.path(pathway_dir, "01_GO_BP_enrichment_DE_genes.csv"), row.names = FALSE)
    
    if (nrow(as.data.frame(ego)) > 0) {
      go_plot <- dotplot(ego, showCategory = 15) + ggtitle("GO Biological Process enrichment")
      ggsave(file.path(pathway_dir, "01_GO_BP_dotplot.png"), go_plot, width = 9, height = 6, dpi = 300)
    }
    
    ekegg <- enrichKEGG(
      gene = unique(gene_map$ENTREZID),
      organism = "hsa",
      pvalueCutoff = 0.05
    )
    
    write.csv(as.data.frame(ekegg), file.path(pathway_dir, "02_KEGG_enrichment_DE_genes.csv"), row.names = FALSE)
    
    if (nrow(as.data.frame(ekegg)) > 0) {
      kegg_plot <- dotplot(ekegg, showCategory = 15) + ggtitle("KEGG pathway enrichment")
      ggsave(file.path(pathway_dir, "02_KEGG_dotplot.png"), kegg_plot, width = 9, height = 6, dpi = 300)
    }
  }
} else {
  message("Pathway enrichment skipped: no significant DE genes or required Bioconductor packages missing.")
}

# -----------------------------
# 14) CellChat analysis
# -----------------------------
# CellChat installation can be system-dependent.
# If it is installed, this block will run. If not, the project still produces all Seurat outputs.

if (has_cellchat) {
  library(CellChat)
  
  data.input <- GetAssayData(seurat_obj_filtered, assay = "RNA", slot = "data")
  meta <- seurat_obj_filtered@meta.data
  
  if (!"cell_type" %in% colnames(meta)) {
    stop("CellChat needs cell_type annotation.")
  }
  
  cellchat <- createCellChat(object = data.input, meta = meta, group.by = "cell_type")
  CellChatDB <- CellChatDB.human
  cellchat@DB <- CellChatDB
  
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat)
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  
  saveRDS(cellchat, file.path(cellchat_dir, "cellchat_all_cell_types.rds"))
  
  png(file.path(cellchat_dir, "01_cellchat_interaction_number_network.png"), width = 1600, height = 1400, res = 200)
  groupSize <- as.numeric(table(cellchat@idents))
  netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = TRUE,
                   label.edge = FALSE, title.name = "Number of interactions")
  dev.off()
  
  png(file.path(cellchat_dir, "02_cellchat_interaction_strength_network.png"), width = 1600, height = 1400, res = 200)
  netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = TRUE,
                   label.edge = FALSE, title.name = "Interaction strength")
  dev.off()
  
  comm_table <- subsetCommunication(cellchat)
  write.csv(comm_table, file.path(cellchat_dir, "03_ligand_receptor_communication_table.csv"), row.names = FALSE)
  
} else {
  message("CellChat package is not installed. Skipping CellChat block.")
}

# -----------------------------
# 15) Save final objects and session info
# -----------------------------
saveRDS(seurat_obj_filtered, file.path(rds_dir, "final_annotated_seurat_object.rds"))

if (exists("caf_obj")) {
  saveRDS(caf_obj, file.path(rds_dir, "final_CAF_object.rds"))
}

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))

message("Analysis finished. Check the analysis_outputs folder.")
