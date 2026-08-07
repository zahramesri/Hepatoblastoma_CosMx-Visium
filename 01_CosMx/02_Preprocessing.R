###############################################################################
# CosMx SMI - Hepatoblastoma (Sindhi Cohort)
# Script 02: Preprocessing & UMAPs
#
# Description:
#   Assigns patient metadata to merged CosMx object based on FOV ranges,
#   performs QC filtering, SCTransform normalization per patient, clustering,
#   UMAP visualization, and marker gene identification.
#
# Input:  00_AllDataMerge.RData (output of Script 01)
# Output: Per-patient clustered Seurat objects, marker gene tables
#
###############################################################################

# ── Libraries ─────────────────────────────────────────────────────────────────
library(Seurat)
library(dplyr)
library(ggplot2)
library(openxlsx)
library(Polychrome)

# ── Paths (edit to match your environment) ────────────────────────────────────
data_dir   <- "path/to/data"
output_dir <- "path/to/output"

# ── Load merged object ─────────────────────────────────────────────────────────
load(file.path(data_dir, "00_AllDataMerge.RData"))

# ── Assign patient metadata from FOV ranges ────────────────────────────────────
# FOV ranges provided via HB_CosMx_manifest (see script header for ID mapping)
patient_ids      <- c("HB_S1", "HB_S2", "HB_S3", "HB_S4", "HB_S5")
fov_ranges       <- c("1-139,245-277", "35-285", "140-244", "19-173", "1-18,174-215")
tumor_recurrence <- c("Relapse", "Relapse", "non-Relapse", "non-Relapse", "Resection")
cosmx_run        <- c("Hepatoblastoma1", "Hepatoblastoma2", "Hepatoblastoma1",
                      "Hepatoblastoma3", "Hepatoblastoma3")

expand_fov <- function(ranges) {
  unlist(lapply(strsplit(ranges, ",")[[1]], function(x) {
    bounds <- as.numeric(strsplit(x, "-")[[1]])
    seq(bounds[1], bounds[2])
  }))
}

meta_fov <- do.call(rbind, lapply(seq_along(patient_ids), function(i) {
  data.frame(Patient_ID     = patient_ids[i],
             fov            = expand_fov(fov_ranges[i]),
             Tumor_recurrence = tumor_recurrence[i],
             CosMx_run      = cosmx_run[i])
}))

# Extract run prefix from cell barcode and match to FOV manifest
HB_Cosmx_merge@meta.data$sample_short    <- sub("_.*", "", rownames(HB_Cosmx_merge@meta.data))
HB_Cosmx_merge@meta.data$orig_cell_name  <- rownames(HB_Cosmx_merge@meta.data)
HB_Cosmx_merge@meta.data$fov             <- as.numeric(HB_Cosmx_merge@meta.data$fov)
meta_fov$fov                             <- as.numeric(meta_fov$fov)

metaD <- HB_Cosmx_merge@meta.data %>%
  left_join(meta_fov, by = c("fov", "sample_short" = "CosMx_run"))

HB_Cosmx_merge@meta.data              <- metaD
rownames(HB_Cosmx_merge@meta.data)    <- HB_Cosmx_merge@meta.data$orig_cell_name
HB_Cosmx_merge@meta.data$sample_short <- NULL

# ── QC: check unmatched FOVs ───────────────────────────────────────────────────
# FALSE = matched, TRUE = unmatched (FOV not in manifest)
# FALSE   TRUE
# 809654  32127
cat("FOV match summary:\n")
print(table(is.na(HB_Cosmx_merge@meta.data$Patient_ID)))

HB_Cosmx_merge@meta.data %>%
  mutate(match = ifelse(is.na(Patient_ID), "Missing", "Matched")) %>%
  ggplot(aes(factor(fov), fill = match)) +
  geom_bar() +
  facet_wrap(~ orig_cell_name, scales = "free_x") +
  theme_minimal() +
  labs(y = "Cell Count", x = "FOV")

HB_Cosmx_merge_filt <- subset(HB_Cosmx_merge, subset = !is.na(Patient_ID))
# ── QC metrics ────────────────────────────────────────────────────────────────
# percent.NegPrb: ratio of negative probe counts to total RNA (proxy for background)
HB_Cosmx_merge_filt[["percent.NegPrb"]] <-
  (HB_Cosmx_merge_filt$nCount_negprobes / HB_Cosmx_merge_filt$nCount_RNA) * 100

# Combine per-sample count layers into a single counts matrix
counts1 <- HB_Cosmx_merge_filt@assays$Nanostring@layers$`counts.1`
counts2 <- HB_Cosmx_merge_filt@assays$Nanostring@layers$`counts.2`
counts3 <- HB_Cosmx_merge_filt@assays$Nanostring@layers$`counts.3`

colnames(counts1) <- rownames(HB_Cosmx_merge_filt@meta.data)[seq_len(ncol(counts1))]
colnames(counts2) <- rownames(HB_Cosmx_merge_filt@meta.data)[
  (ncol(counts1) + 1):(ncol(counts1) + ncol(counts2))]
colnames(counts3) <- rownames(HB_Cosmx_merge_filt@meta.data)[
  (ncol(counts1) + ncol(counts2) + 1):(ncol(counts1) + ncol(counts2) + ncol(counts3))]

counts_merged <- cbind(counts1, counts2, counts3)
rownames(counts_merged) <- rownames(HB_Cosmx_merge_filt)
HB_Cosmx_merge_filt@assays$Nanostring$counts <- counts_merged

# ── QC plots ───────────────────────────────────────────────────────────────────
ggplot(HB_Cosmx_merge_filt@meta.data, aes(x = factor(fov), y = percent.NegPrb)) +
  geom_violin(fill = "skyblue") +
  geom_jitter(height = 0, width = 0.2, size = 0.1, alpha = 0.3) +
  theme_classic() +
  ylab("% Negative Probes") +
  xlab("FOV")

# ── QC filtering ───────────────────────────────────────────────────────────────
# Thresholds based on distribution summaries:
#   percent.NegPrb: median 0.033, mean 0.092 -> cutoff <= 0.10
#   nCount_RNA:     median 1193, mean 1585   -> keep 500-10000
#   nFeature_RNA:   median 761,  mean 863    -> keep 400-5000
HB_Cosmx_merge_filt <- subset(
  HB_Cosmx_merge_filt,
  subset = percent.NegPrb  <= 0.10 &
           nCount_RNA      >= 500  & nCount_RNA   <= 10000 &
           nFeature_RNA    >= 400  & nFeature_RNA <= 5000
)

VlnPlot(HB_Cosmx_merge_filt,
        features = c("nFeature_RNA", "nCount_RNA", "percent.NegPrb"),
        ncol = 3, pt.size = 0)

# ── Subset to patients of interest ────────────────────────────────────────────
patients_keep <- c("HB_S2", "HB_S4", "HB_S5")

HB_Cosmx_merge_filt_sub <- subset(
  HB_Cosmx_merge_filt,
  subset = Patient_ID %in% patients_keep
)

# ── Normalization (SCTransform per patient) ────────────────────────────────────
cosmx_by_patient <- SplitObject(HB_Cosmx_merge_filt_sub, split.by = "Patient_ID")

cosmx_by_patient <- lapply(cosmx_by_patient, function(obj) {
  SCTransform(obj, assay = "Nanostring", new.assay.name = "SCT",
              variable.features.n = 3000, layer = "counts", verbose = TRUE)
})

# ── Clustering & UMAP ─────────────────────────────────────────────────────────
run_patient_clustering <- function(obj, dims = 1:20,
                                   resolutions = c(0.4, 0.5, 0.6),
                                   npcs = max(dims)) {
  DefaultAssay(obj) <- "SCT"
  obj <- RunPCA(obj, npcs = npcs, verbose = FALSE)
  obj <- RunUMAP(obj, dims = dims, reduction = "pca", verbose = FALSE)
  obj <- FindNeighbors(obj, dims = dims, reduction = "pca", verbose = FALSE)
  for (r in resolutions) obj <- FindClusters(obj, resolution = r, verbose = FALSE)
  obj
}

obj_list <- lapply(cosmx_by_patient, run_patient_clustering)

# Set default clustering resolution (0.5) as active identity
obj_list <- lapply(obj_list, function(obj) {
  obj$seurat_clusters <- obj$SCT_snn_res.0.5
  Idents(obj) <- "seurat_clusters"
  obj
})

cluster_colors <- glasbey(length(0:17))
names(cluster_colors) <- as.character(0:17)

# ── UMAPs ─────────────────────────────────────────────────────────────────────
for (nm in names(obj_list)) {
  print(DimPlot(obj_list[[nm]], reduction = "umap", label = TRUE,
                cols = cluster_colors, raster = FALSE) + ggtitle(nm))
}

# ── Spatial plot (global FOV centroids) ───────────────────────────────────────
for (nm in names(obj_list)) {
  obj <- obj_list[[nm]]
  cell_centroid_df <- data.frame(y    = obj$CenterX_global_px,
                                 x    = obj$CenterY_global_px,
                                 cell = colnames(obj))
  centroid_data <- list("centroids" = CreateCentroids(cell_centroid_df))
  coords        <- CreateFOV(coords = centroid_data, type = "centroids", assay = "RNA")
  obj_list[[nm]][["global"]] <- coords

  print(ImageDimPlot(obj_list[[nm]], fov = "global", axes = TRUE,
                     cols = cluster_colors, size = 0.75, dark.background = FALSE) +
        ggtitle(nm))
}

# ── Marker genes ───────────────────────────────────────────────────────────────
dir.create(file.path(output_dir, "markers"), showWarnings = FALSE)

for (nm in names(obj_list)) {
  Idents(obj_list[[nm]]) <- "seurat_clusters"
  markers_df <- FindAllMarkers(obj_list[[nm]], only.pos = TRUE,
                               min.pct = 0.1, logfc.threshold = 0.25)
  write.xlsx(data.frame(markers_df),
             file = file.path(output_dir, "markers",
                              paste0(nm, "_res0.5_AllMarkers.xlsx")))
}

# ── Dot plots ──────────────────────────────────────────────────────────────────
markers_canonical <- c(
  "EPCAM", "SOX9", "GLUL", "KRT7", "AFP", "ALB", "APOA1", "FABP1", "ALDOB", "HNF4A",
  "PTPRC", "CD2", "CD3D", "CD3E", "CD3G", "CD4", "FOXP3", "CD8A", "CD8B",
  "FCGR3A", "NCAM1", "KLRD1", "KIR2DL3",
  "BLNK", "CD19", "TCL1A", "MS4A1", "CD79A", "CD79B", "CD27",
  "CD38", "SDC1", "JCHAIN",
  "IGHD", "IGHM", "IGHE", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHA1", "IGHA2",
  "LAMP3", "CD14", "CD68", "CD163", "CSF1R", "ITGAM", "ITGAX",
  "PDCD1", "TIGIT", "CTLA4", "LAG3", "HAVCR2", "BTLA", "CD244", "CD160",
  "GZMA", "GZMB", "IGLC2", "IGKC"
)

markers_key <- c(
  "EPCAM", "FGA", "VWF", "CDH5", "LYVE1", "RELN", "SPP1", "GLUL", "CYP2E1", "SDS",
  "ONECUT2", "DLK1", "EFNA1", "KRT7", "KRT8", "KRT18",
  "CD24", "CD44", "CD34", "PROM1", "KIT", "THY1", "ITGB1", "ITGA6",
  "AFP", "CD109", "CDH1",
  "PTPRC", "CD14", "FCGR3A", "ADGRE1", "AXL", "MAFB", "CX3CR1", "C5AR1",
  "CD5L", "VSIG4", "CD3D", "KLRF1", "CD19", "SDC1",
  "GZMA", "NKG7", "CCL4", "KLRB1", "IER2",
  "TBX21", "GZMB", "IL7R", "CXCR6", "EOMES", "GZMK"
)

for (nm in names(obj_list)) {
  obj <- obj_list[[nm]]
  Idents(obj) <- "seurat_clusters"

  # Top 5 markers per cluster dot plot
  top5 <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.1,
                         logfc.threshold = 0.25, verbose = FALSE) %>%
    group_by(cluster) %>%
    top_n(n = 5, wt = avg_log2FC)

  print(
    DotPlot(obj, features = unique(top5$gene), cluster.idents = TRUE, dot.scale = 6) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
      scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
      ggtitle(nm)
  )

  # Canonical marker dot plot
  print(
    DotPlot(obj,
            features = unique(c(markers_key, markers_canonical)),
            cluster.idents = TRUE, dot.scale = 6) +
      theme_classic() +
      coord_flip() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
      scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
      ggtitle(nm)
  )
}

# ── Fig: Manuscript UMAP panel (annotated by cell type) ───────────────────────
# Cell type and state annotation was performed manually based on canonical marker
# gene expression (markers_key, markers_canonical) and top differentially
# expressed genes per cluster. See Methods and supplement table for details.
library(patchwork)

# HB_S2 = Relapse, HB_S4 = non-Relapse, HB_S5 = Resection
Idents(obj_list[["HB_S2"]]) <- "my.annot"
Idents(obj_list[["HB_S4"]]) <- "my.annot"
Idents(obj_list[["HB_S5"]]) <- "my.annot"

# Build a master color palette — consistent color per cell type across all samples
all_celltypes <- unique(c(levels(obj_list[["HB_S2"]]),
                          levels(obj_list[["HB_S4"]]),
                          levels(obj_list[["HB_S5"]])))
all_celltypes <- all_celltypes[!is.na(all_celltypes)]

all_colors <- c(
  "#4E79A7", "#F28E2B", "#76B7B2", "#E15759", "#59A14F",
  "#EDC949", "#AF7AA1", "#FF9DA7", "#9C755F", "#BAB0AC",
  "#86BC3D", "#D37295", "#8B88D8", "#C4AD66", "#6B9AC4",
  "#8B4513", "#A6761D", "#0072B2", "#D55E00", "#009E73",
  "#CC79A7", "#F0E442", "#56B4E9", "#E69F00", "#5A5A5A"
)

master_colors <- setNames(all_colors[seq_along(all_celltypes)], all_celltypes)

# Override specific cell types if needed
master_colors["Cytotoxic NK/T cells"]                    <- "#FF6B6B"
master_colors["inflammatory tumor epithelium (SOX9+)"]   <- "#4ECDC4"

annot_colors_S2 <- master_colors[levels(obj_list[["HB_S2"]])]
annot_colors_S4 <- master_colors[levels(obj_list[["HB_S4"]])]
annot_colors_S5 <- master_colors[levels(obj_list[["HB_S5"]])]

clean_umap_theme <- theme(
  plot.title   = element_text(hjust = 0.5, face = "bold", size = 16),
  axis.text    = element_blank(),
  axis.ticks   = element_blank(),
  axis.title   = element_blank(),
  axis.line    = element_blank(),
  panel.border = element_blank(),
  legend.position = "none",
  plot.margin  = margin(5, 5, 5, 5)
)

p1 <- DimPlot(obj_list[["HB_S2"]], reduction = "umap",
              label = FALSE, raster = FALSE, repel = TRUE, pt.size = 0.3,
              cols = annot_colors_S2) + NoLegend() +
  ggtitle("Relapse") + clean_umap_theme +
  theme(plot.title = element_text(color = "#9B1A1A", face = "bold", size = 16, hjust = 0.5))

p2 <- DimPlot(obj_list[["HB_S4"]], reduction = "umap",
              label = FALSE, raster = FALSE, repel = TRUE, pt.size = 0.3,
              cols = annot_colors_S4) + NoLegend() +
  ggtitle("Non-relapse") + clean_umap_theme +
  theme(plot.title = element_text(color = "#1A4A6B", face = "bold", size = 16, hjust = 0.5))

p3 <- DimPlot(obj_list[["HB_S5"]], reduction = "umap",
              label = FALSE, raster = FALSE, repel = TRUE, pt.size = 0.3,
              cols = annot_colors_S5) + NoLegend() +
  ggtitle("Resection") + clean_umap_theme +
  theme(plot.title = element_text(color = "#2D6B4F", face = "bold", size = 16, hjust = 0.5))

combined <- p1 | p2 | p3

ggsave(file.path(output_dir, "UMAP_panelB.png"),
       combined, width = 18, height = 6, dpi = 300, bg = "white")
