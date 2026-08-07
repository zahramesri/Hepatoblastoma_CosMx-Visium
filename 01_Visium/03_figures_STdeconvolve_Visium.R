###############################################################################
# Spatial Transcriptomics - Hepatoblastoma (Sindhi Cohort)
# Script 03: Manuscript Figures — STdeconvolve
#
# Generates:
#   1. Per-sample STdeconvolve topic proportion plots (vizAllTopics)
#   2. Per-sample UMAP and spatial plots colored by cell-type annotation
#   3. Per-sample dot plots of top STdeconvolve topic markers
#
# Requires ggplot2 version 3.5.1
# Input:
#   - ldas_results_theta&beta.RData  (STdeconvolve results)
#   - 00_STdeconvolve_corpus_count.RData                 (spatial coordinates-refer to script 02)
#   - updated_annotation_data/       (annotated Seurat objects, one per sample) #this is after biological interpretation of the top markers and Topic markers 
#   - TopMarkers_bySample.xlsx       (top markers per STdeconvolve topic)
###############################################################################

# ── Libraries ─────────────────────────────────────────────────────────────────
library(STdeconvolve)
library(Seurat)
library(ggplot2)   # version 3.5.1 required
library(readxl)
library(dplyr)

# ── Paths (edit to match your environment) ────────────────────────────────────
data_dir   <- "path/to/data"    # directory containing STdeconvolve outputs and annotated objects
output_dir <- "path/to/output"  # directory for figure outputs

annot_dir  <- file.path(data_dir, "updated_annotation_data")

# ── Sample metadata ────────────────────────────────────────────────────────────
sample_names <- liver_new$SampleID

# ── Fig: STdeconvolve topic proportions per sample ────────────────────────────
load(file.path(data_dir, "ldas_results_theta&beta.RData")) # results
load(file.path(data_dir, "00_STdeconvolve_corpus_count.RData")) # coords

names(results) <- sample_names

my_colors <- c("#D62728", "#1F77B4", "#2CA02C", "#9467BD", "#FF7F0E")

for (i in seq_along(sample_names)) {
  deconProp <- results[[i]]$theta

  p <- vizAllTopics(deconProp, coords[[i]],
                    r = 4, lwd = 0, showLegend = TRUE) +
    ggplot2::guides(colour = "none") +
    ggplot2::scale_fill_manual(
      name   = "STdec",
      values = my_colors,
      labels = c("Topic.1" = "STdec.1", "Topic.2" = "STdec.2",
                 "Topic.3" = "STdec.3", "Topic.4" = "STdec.4",
                 "Topic.5" = "STdec.5")
    ) +
    ggplot2::ggtitle(sample_names[[i]])

  ggsave(
    filename = file.path(output_dir, paste0("STdec_sample_", i, "_", sample_names[[i]], ".png")),
    plot     = p,
    width    = 8, height = 7, dpi = 300, bg = "white"
  )
}

# ── Load annotated Seurat objects ─────────────────────────────────────────────
annot_files  <- list.files(annot_dir, pattern = "\\.RData$", full.names = TRUE)

for (file in annot_files) {
  object_name <- tools::file_path_sans_ext(basename(file))
  loaded_objs <- load(file)
  assign(object_name, get(loaded_objs[1]))
  rm(list = loaded_objs[1])
  cat("Loaded:", object_name, "\n")
}

# ── Color palette for cell-type annotations ────────────────────────────────────
cluster_colors <- c(
  "Tumor"                      = "#1F78B4",
  "ECM"                        = "#E31A1C",
  "Stroma + InflammatoryMarkers" = "#33A02C",
  "Fibrosis"                   = "#6A3D9A",
  "Stroma"                     = "#B15928",
  "Hepatocyte"                 = "#FDBF2F",
  "Immune + Tumor"             = "#A6CEE3",
  "Hepatocyte + Monocyte"      = "#FB9A99"
)

seurat_objects <- #name of your annotated objects

image_ids <- seurat_objects@images

# ── Fig: UMAP colored by annotation ───────────────────────────────────────────
for (obj_name in seurat_objects) {
  seurat_obj <- get(obj_name)
  Idents(seurat_obj) <- "annotation"

  png(file.path(output_dir, paste0("Region_annotated_UMAP_", obj_name, ".png")),
      res = 300, width = 7, height = 4.5, units = "in")
  print(
    DimPlot(seurat_obj, raster = FALSE, label = FALSE, cols = cluster_colors) +
      ggtitle("") +
      theme(plot.title  = element_text(size = 16, face = "bold", hjust = 0.5),
            axis.text   = element_blank(),
            axis.ticks  = element_blank(),
            axis.title  = element_blank())
  )
  dev.off()
}

# ── Fig: Spatial plot colored by annotation ────────────────────────────────────
for (i in seq_along(seurat_objects)) {
  obj_name   <- seurat_objects[i]
  img_id     <- image_ids[i]
  seurat_obj <- get(obj_name)
  Idents(seurat_obj) <- "annotation"

  png(file.path(output_dir, paste0("SpatialDimPlot_", obj_name, ".png")),
      res = 300, width = 7, height = 4.5, units = "in")
  print(
    SpatialDimPlot(seurat_obj, cols = cluster_colors,
                   images = img_id, pt.size.factor = 3000) +
      ggtitle(paste0("Sample: ", obj_name)) +
      theme(plot.title  = element_text(size = 16, face = "bold", hjust = 0.5),
            axis.text   = element_blank(),
            axis.ticks  = element_blank(),
            axis.title  = element_blank())
  )
  dev.off()
}

# ── Fig: Dot plots of top STdeconvolve topic markers per sample ───────────────
TopMarkers_bySample <- lapply(excel_sheets(file.path(data_dir, "TopMarkers_bySample.xlsx")),
                              function(s) read_excel(file.path(data_dir, "TopMarkers_bySample.xlsx"),
                                                     sheet = s))
names(TopMarkers_bySample) <- excel_sheets(file.path(data_dir, "TopMarkers_bySample.xlsx"))

TopMarkers_bySample <- lapply(TopMarkers_bySample, function(df) {
  colnames(df) <- gsub("CellType_", "STdec.", colnames(df))
  df
})

seurat_list <- list(seurat_objects) #make a list from all the sample objects

seurat_list <- lapply(seurat_list, function(x) { Idents(x) <- "annotation"; x })

marker_list_by_sample <- lapply(TopMarkers_bySample, function(df) {
  markers <- lapply(df, function(x) as.character(na.omit(unique(x))))
  markers[sapply(markers, length) > 0]
})

# Remove genes already used in earlier topics within each sample
marker_list_by_sample_unique <- lapply(marker_list_by_sample, function(sample_markers) {
  seen <- character(0)
  cleaned <- lapply(sample_markers, function(genes) {
    genes <- genes[!genes %in% seen]
    seen  <<- c(seen, genes)
    genes
  })
  cleaned[sapply(cleaned, length) > 0]
})

common_samples <- intersect(names(seurat_list), names(marker_list_by_sample_unique))

dotplot_dir <- file.path(output_dir, "STdeconvolve_DotPlots_bySample")
dir.create(dotplot_dir, showWarnings = FALSE, recursive = TRUE)

for (samp in common_samples) {
  obj     <- seurat_list[[samp]]
  markers <- lapply(marker_list_by_sample_unique[[samp]], head, 3)
  markers <- lapply(markers, function(g) g[g %in% rownames(obj)])
  markers <- markers[sapply(markers, length) > 0]

  p <- DotPlot(object = obj, features = markers) +
    RotatedAxis() +
    ggtitle(samp) +
    theme(plot.title   = element_text(hjust = 0.5, face = "bold"),
          axis.text.x  = element_text(size = 15),
          axis.text.y  = element_text(size = 15))

  ggsave(
    filename = file.path(dotplot_dir, paste0(samp, "_Top3_DotPlot.png")),
    plot     = p,
    width    = 12, height = 4
  )
}
