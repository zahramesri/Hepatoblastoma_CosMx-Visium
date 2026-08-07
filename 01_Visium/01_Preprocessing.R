###############################################################################
# Spatial Transcriptomics - Hepatoblastoma (Sindhi Cohort)
# Script 01: Preprocessing&UMAP — 11 Samples (HB_1364 Excluded)
#
# Description:
#   Merges SCT-normalized Visium objects from 11 samples (HB_1364 excluded
#   due to poor spatial mapping quality), performs dimensionality reduction,
#   clustering, and cluster marker identification.
#
# Input:  00_SCTNormalized_perSample.RData
#           - liver       : named list of per-sample SCT Seurat objects
#           - colMeta     : sample-level metadata data frame
#           - st.features : variable features selected jointly across samples
# Output: Downstream Seurat objects and cluster marker tables
#
# References:
#   https://satijalab.org/seurat/articles/merge_vignette.html
#   https://github.com/satijalab/seurat/issues/5761
#   https://github.com/satijalab/seurat/issues/2814
###############################################################################

# ── Libraries ─────────────────────────────────────────────────────────────────
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyverse)
library(Matrix)
library(future)
library(xlsx)
library(scales)
library(enrichR)

# ── Paths (edit to match your environment) ────────────────────────────────────
data_dir   <- "path/to/data"
output_dir <- "path/to/output"

# ── Load data ─────────────────────────────────────────────────────────────────
load(file.path(data_dir, "00_SCTNormalized_PerPatient.RData"))

options(future.globals.maxSize = 8000 * 1024^2)

# ── Merge samples ──────────────────────────────────────────────────────────────
sample_ids <- colMeta$SampleID   # HB_ IDs, 11 samples (HB_1364 absent)

mergedLiverlist <- merge(
  liver[[1]],
  y            = liver[2:11],
  add.cell.ids = sample_ids,
  merge.data   = TRUE
)

# ── QC ────────────────────────────────────────────────────────────────────────
mergedLiverlist[["percent.mt"]] <- PercentageFeatureSet(mergedLiverlist, pattern = "^MT-")

VlnPlot(mergedLiverlist,
        features = c("nFeature_Spatial", "nCount_Spatial"),
        pt.size  = 0, ncol = 2) &
  theme(axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank())

SpatialFeaturePlot(mergedLiverlist, features = "nCount_Spatial", ncol = 3) +
  theme(legend.position = "right")

plot_list <- lapply(seq_along(liver), function(i) {
  FeatureScatter(liver[[i]], feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial") + NoLegend()
})

pdf(file.path(output_dir, "QC/EachSample_corrValue.pdf"))
for (p in plot_list) print(p)
dev.off()

mergedMultiData <- subset(mergedLiverlist,
                          subset = nCount_Spatial > 500 & nCount_Spatial < 50000 &
                            nFeature_Spatial > 500)

SpatialFeaturePlot(mergedMultiData, features = "nCount_Spatial", ncol = 3) +
  theme(legend.position = "right")

# ── Metadata ───────────────────────────────────────────────────────────────────
colMeta$TransplantRelapse <- paste(colMeta$`Transplant /Resection`, colMeta$`Tumor relapse`)

VariableFeatures(mergedMultiData) <- st.features

seurat_metadata             <- mergedMultiData@meta.data
seurat_metadata$SampleID    <- seurat_metadata$orig.ident
metadata_combined           <- merge(seurat_metadata, colMeta, by = "SampleID")
rownames(metadata_combined) <- rownames(mergedMultiData@meta.data)
mergedMultiData@meta.data   <- metadata_combined

save(file = file.path(output_dir, "00_MergedMultiData_sctNormalized_excludingHB_1364.RData"),
     mergedMultiData)
#Object provided in GEO

# ── Inspect top-expressed genes ────────────────────────────────────────────────
C   <- mergedMultiData@assays$Spatial@counts
C@x <- C@x / rep.int(colSums(C), diff(C@p))
most_expressed <- order(Matrix::rowSums(C), decreasing = TRUE)[20:1]
mat <- t(as.matrix(C[most_expressed, ]))

boxplot(mat, cex = 0.1, las = 1, xlab = "% total count per spot",
        col = (scales::hue_pal())(20)[20:1], horizontal = TRUE)

# Remove ALB — constitutively expressed in liver, dominates counts
liver_new <- mergedMultiData[!grepl("^ALB", rownames(mergedMultiData)), ]

# ── Dimensionality reduction ───────────────────────────────────────────────────
liver_new <- RunPCA(liver_new, assay = "SCT", verbose = FALSE)
ElbowPlot(liver_new, ndims = 30)

# Automated PC selection
# https://hbctraining.github.io/scRNA-seq/lessons/sc_exercises_clustering_analysis.html
pct  <- liver_new[["pca"]]@stdev / sum(liver_new[["pca"]]@stdev) * 100
cumu <- cumsum(pct)
co1  <- which(cumu > 90 & pct < 5)[1]
co2  <- sort(which((pct[1:(length(pct) - 1)] - pct[2:length(pct)]) > 0.1), decreasing = TRUE)[1] + 1
pcs  <- min(co1, co2)
message("Selected number of PCs: ", pcs)

plot_df <- data.frame(pct = pct, cumu = cumu, rank = seq_along(pct))
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) +
  geom_text() +
  geom_vline(xintercept = 90, color = "grey") +
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()

liver_new <- FindNeighbors(liver_new, reduction = "pca", dims = 1:pcs)
liver_new <- FindClusters(liver_new,  verbose = FALSE)
liver_new <- RunUMAP(liver_new, reduction = "pca", dims = 1:pcs)

save(file = file.path(output_dir, "01_DownstreamMerged_liver_excluding_HB1364.RData"), liver_new)

# ── Visualization for Figure1 ──────────────────────────────────────────────────────────────
my_colors <- c(
  "Non-relapse" = "#4A90C4",  
  "Relapse"     = "#C0392B", 
  "Resection"   = "#4A8C6A"  
)

png("output_dir.png", 
    res = 300, width = 7, height = 4.5, units = "in")
Idents(liver_new) <- "seurat_clusters"
DimPlot(liver_new, 
        raster = FALSE, 
        label = FALSE,
        cols = c(hcl.colors(palette = "Dynamic", n=10), hcl.colors(palette = "Earth", n=10),hcl.colors(palette = "Dark2", n=11))) + 
  ggtitle("") +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 16),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank())

dev.off()

png("output_dir.png", 
    res = 300, width = 7, height = 4.5, units = "in")
Idents(liver_new) <- "condition"
DimPlot(liver_new, 
        raster = FALSE, 
        label = FALSE, 
        cols = my_colors) + 
  ggtitle("") +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 20),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank())

dev.off()

png("output_dir.png", 
    res = 300, width = 2, height = 2, units = "in")
ggplot(liver_new@meta.data, aes(x = SampleID, fill = seurat_clusters)) +
  geom_bar(position = "fill", width = 0.8) +  # Use 'fill' to normalize by sample
  scale_fill_manual(
    values = c(
      hcl.colors(n = 10, palette = "Dynamic"),
      hcl.colors(n = 10, palette = "Earth"),
      hcl.colors(n = 11, palette = "Dark2")
    )
  ) +
  labs(
    x = NULL,
    y = "Proportion of Spots",
    fill = "Seurat clusters"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 8)
  )

dev.off()

# ── Cluster markers ────────────────────────────────────────────────────────────
Idents(liver_new) <- "seurat_clusters"
cluster_markers <- FindAllMarkers(liver_new, only.pos = TRUE, min.pct = 0.25,
                                  logfc.threshold = 0.25, verbose = FALSE,
                                  recorrect_umi = FALSE)

top100 <- cluster_markers %>% group_by(cluster) %>% top_n(100, avg_log2FC)
write.xlsx(file = file.path(output_dir, "Merged_ClusterMarkers_seuratClusters_top100.xlsx"),
           data.frame(top100))

top20 <- cluster_markers %>% group_by(cluster) %>% top_n(20, avg_log2FC)
liver_scaled <- ScaleData(liver_new, features = rownames(liver_new))
DoHeatmap(subset(liver_scaled, downsample = 1000), features = top20$gene, size = 4) +
  NoLegend() +
  theme(text = element_text(size = 8))
