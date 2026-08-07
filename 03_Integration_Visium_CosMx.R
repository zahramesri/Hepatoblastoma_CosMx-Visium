###############################################################################
# CosMx SMI - Hepatoblastoma (Sindhi Cohort)
# Script 03: Integration — Visium STdeconvolve Topics Projected onto CosMx
#
# Description:
#   Projects Visium STdeconvolve topic marker genes onto CosMx single-cell
#   data using module scoring (AddModuleScore). Topic markers are scored per
#   cell and compared across annotated cell types via violin plots and
#   average score heatmaps.
#
# Input:  Annotated per-patient CosMx Seurat objects (obj_list, from Script 02)
# Output: Module score violin plots, average score heatmaps per patient
#
# STdeconvolve topic-to-biology mapping per sample:
#   HB_S2: CellType 1,2,5 = Tumor  |  CellType 3,4 = Stroma
#   HB_S5: CellType 1,5   = Tumor  |  CellType 2 = Fibrosis  |  CellType 4 = Hepatocyte/Stroma
#   HB_S4: CellType 1,4   = Stroma |  CellType 2,3,5 = Hepatocytes
###############################################################################

# ── Libraries ─────────────────────────────────────────────────────────────────
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(pheatmap)
library(tibble)

# ── Paths (edit to match your environment) ────────────────────────────────────
data_dir   <- "path/to/data"
output_dir <- "path/to/output"

# ── STdeconvolve topic marker lists (derived from Visium, per sample) ─────────
STdeconvolve_markers <- list(

  HB_S2 = list(
    CellType_1 = c("NEFL", "IGKC", "SST", "IGHA1", "IGHG1",
                   "C7", "S100A6", "C1QB", "C1QC", "TAGLN", "ACTA2", "CLU"),
    CellType_2 = c("CYP2E1", "FGG", "FGB", "GLUL", "APOC3",
                   "AMBP", "FABP1", "AHSG", "SCD", "VTN"),
    CellType_3 = c("ASS1", "AGT", "HIST1H1E", "DPEP1", "NOTUM"),
    CellType_4 = c("COL2A1", "ITIH2", "BCAM", "HSPD1", "SCD",
                   "CPS1", "GPC3", "APOA1", "FGG"),
    CellType_5 = c("DEFA5", "MYH7B", "AFP", "CST1", "DLK1",
                   "HSPA1A", "NOTUM")
  ),

  HB_S5 = list(
    CellType_1 = c("GPC3", "CYP2E1", "HPD", "AHSG", "GLUL", "CPS1"),
    CellType_2 = c("IGKC", "COL1A2", "S100A6", "VIM", "C7", "COL1A1",
                   "FOS", "COL3A1", "CCN2", "CCN1", "TIMP1", "EGR1"),
    CellType_3 = c("ELANE", "HBG2", "HBA1", "HBA2"),
    CellType_4 = c("PLA2G2A", "SAA1", "CCL20", "CYP3A7", "APOC2",
                   "APOC3", "CCN1", "APOH", "FABP1", "IGFBP1", "CPB2", "APOB"),
    CellType_5 = c("TF", "HSPD1", "CYP3A5", "HINT1", "FN1",
                   "SPINK1", "COX7C", "ITIH2", "CPS1", "APOE")
  ),

  HB_S4 = list(
    CellType_1 = c("GPC3", "APOB", "DLK1", "GLUL", "PEG10", "ITIH2", "AGT"),
    CellType_2 = c("HRG", "CRP", "PCK1", "SOD2", "FGB", "CEBPD",
                   "SLC39A14", "IGKC", "SERPINA1", "CYP2E1", "APOB", "APOC2"),
    CellType_3 = c("PPIB", "APOA2", "GLUL"),
    CellType_4 = c("AXIN2", "NKD1", "APCDD1", "MYH7B", "SOX13", "CST1"),
    CellType_5 = c("CEBPA", "AMBP", "PABPC1", "AGT")
  )
)

# ── Module scoring and visualization per sample ───────────────────────────────
for (nm in names(STdeconvolve_markers)) {
  obj          <- obj_list[[nm]]
  marker_list  <- STdeconvolve_markers[[nm]]
  score_cols   <- paste0("STdec_", seq_along(marker_list))

  # Filter to genes present in the object
  marker_list_filt <- lapply(marker_list, function(x) intersect(x, rownames(obj)))

  # Spatial feature plot of topic markers
  print(ImageFeaturePlot(obj, features = unlist(marker_list_filt)))

  # Score each topic per cell
  obj <- AddModuleScore(obj, features = marker_list_filt, name = "STdec_")
  obj_list[[nm]] <- obj

  # Violin plot of scores by cell-type annotation
  print(VlnPlot(obj, features = score_cols, group.by = "my.annot", pt.size = 0))

  # Average score per annotation → heatmap
  avg_scores <- obj@meta.data %>%
    group_by(my.annot) %>%
    summarise(across(all_of(score_cols), mean, na.rm = TRUE))

  mat     <- column_to_rownames(as.data.frame(avg_scores), "my.annot") %>% as.matrix()
  max_val <- max(abs(mat))

  pheatmap(mat,
           scale        = "column",
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           color        = colorRampPalette(c("blue", "white", "red"))(100),
           breaks       = seq(-max_val, max_val, length.out = 101),
           main         = nm)
}


