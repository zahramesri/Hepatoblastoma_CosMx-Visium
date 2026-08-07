###############################################################################
# CosMx SMI - Hepatoblastoma (Sindhi Cohort)
# Script 05: Cell-Cell Interaction Analysis using CellChat
#
# Reference:
#   https://rpubs.com/eozcelik/cellchat_spatial
#
# Description:
#   Performs spatial cell-cell interaction analysis on CosMx single-cell data
#   using CellChat. Computes communication probabilities via truncatedMean with
#   distance-weighted interactions (interaction.range = 250 µm, scale.distance
#   = 4.5, contact.range = 20 µm). Summarizes signaling roles by pathway
#   category (Secreted Signaling, ECM-Receptor, Cell-Cell Contact).
#
# Input:  Annotated per-patient CosMx Seurat objects (obj_list, from Script 02)
# Output: CellChat objects, circle plots, heatmaps, signaling role scatter plots
#
# Note: CellChat and presto must be installed before first use:
#   devtools::install_github("jinworks/CellChat")
#   devtools::install_github("immunogenomics/presto")
###############################################################################

# ── Libraries ─────────────────────────────────────────────────────────────────
library(CellChat)
library(Seurat)
library(ggplot2)
library(dplyr)
library(NMF)
library(ggalluvial)
library(patchwork)

# ── Paths (edit to match your environment) ────────────────────────────────────
data_dir   <- "path/to/data"
output_dir <- "path/to/output"

# ── Global settings ───────────────────────────────────────────────────────────
options(future.globals.maxSize = 50 * 1024^3)   # 50 GB; adjust if needed
future::plan("multisession", workers = 4)

# ── CellChatDB ────────────────────────────────────────────────────────────────
CellChatDB     <- CellChatDB.human
CellChatDB.use <- subsetDB(CellChatDB)           # all protein signaling categories

# ── Helper: prepare spatial input from a CosMx Seurat object ─────────────────
prep_spatial_input <- function(obj, sample_name) {
  list(
    expr_mat     = GetAssayData(obj, assay = "SCT", slot = "data"),
    meta         = data.frame(
      labels     = obj$my.annot,
      fov        = obj$fov,
      patient_id = obj$Patient_ID,
      condition  = obj$Tumor_recurrence,
      sample     = sample_name,
      row.names  = colnames(obj)
    ),
    spatial.locs = data.frame(
      x = obj$CenterX_global_px,
      y = obj$CenterY_global_px,
      row.names = colnames(obj)
    )
  )
}

# ── Helper: summarize CellChat results ────────────────────────────────────────
run_cellchat_summary <- function(cellchat_obj,
                                 object_name        = "CellChat",
                                 top_n_interactions = 100,
                                 compute_centrality = TRUE) {

  # Pathway-level summary
  df.netP      <- subsetCommunication(cellchat_obj, slot.name = "netP")
  pathway_df   <- as.data.frame(sort(table(df.netP$pathway_name), decreasing = TRUE))
  colnames(pathway_df) <- c("pathway_name", "count")
  df.netP      <- df.netP[order(-df.netP$prob), ]

  top_interactions <- head(df.netP, top_n_interactions)
  source_counts    <- as.data.frame(table(top_interactions$source))
  target_counts    <- as.data.frame(table(top_interactions$target))
  colnames(source_counts) <- colnames(target_counts) <- c("cell_type", "count")

  p_source <- ggplot(source_counts, aes(x = reorder(cell_type, -count), y = count)) +
    geom_bar(stat = "identity", fill = "skyblue") + theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8)) +
    labs(title = paste0(object_name, ": Top ", top_n_interactions, " Sources"),
         x = "Source Cell Type", y = "Count")

  p_target <- ggplot(target_counts, aes(x = reorder(cell_type, -count), y = count)) +
    geom_bar(stat = "identity", fill = "salmon") + theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8)) +
    labs(title = paste0(object_name, ": Top ", top_n_interactions, " Targets"),
         x = "Target Cell Type", y = "Count")

  # Annotation-based pathway groups
  df.net            <- subsetCommunication(cellchat_obj)
  secreted_pathways <- unique(df.net$pathway_name[df.net$annotation == "Secreted Signaling"])
  ecm_pathways      <- unique(df.net$pathway_name[df.net$annotation == "ECM-Receptor"])
  contact_pathways  <- unique(df.net$pathway_name[df.net$annotation == "Cell-Cell Contact"])

  # Centrality and signaling role scatter plots
  gg_all <- gg_secreted <- gg_ecm <- signaling_role_plot <- NULL
  if (compute_centrality) {
    cellchat_obj <- netAnalysis_computeCentrality(cellchat_obj)
    gg_all       <- netAnalysis_signalingRole_scatter(cellchat_obj,
                      title = paste0(object_name, ": All Signaling"))
    if (length(secreted_pathways) > 0)
      gg_secreted <- netAnalysis_signalingRole_scatter(cellchat_obj,
                       signaling = secreted_pathways,
                       title    = paste0(object_name, ": Secreted Signaling"))
    if (length(ecm_pathways) > 0)
      gg_ecm <- netAnalysis_signalingRole_scatter(cellchat_obj,
                  signaling = ecm_pathways,
                  title     = paste0(object_name, ": ECM-Receptor"))
    plot_list          <- Filter(Negate(is.null), list(gg_all, gg_secreted, gg_ecm))
    if (length(plot_list) > 0) signaling_role_plot <- wrap_plots(plot_list)
  }

  list(
    cellchat_obj_updated        = cellchat_obj,
    netP_table                  = df.netP,
    full_communication_table    = df.net,
    pathway_summary             = pathway_df,
    top_interactions            = top_interactions,
    source_counts               = source_counts,
    target_counts               = target_counts,
    source_target_plot          = p_source + p_target,
    secreted_signaling_pathways = secreted_pathways,
    ecm_receptor_pathways       = ecm_pathways,
    cell_cell_contact_pathways  = contact_pathways,
    signaling_role_all          = gg_all,
    signaling_role_secreted     = gg_secreted,
    signaling_role_ecm          = gg_ecm,
    signaling_role_plot         = signaling_role_plot
  )
}

# ── Spatial conversion factor (CosMx pixels to microns) ──────────────────────
# CosMx coordinates are in pixels; 0.18 µm/pixel for this instrument.
# Skip if your CosMx data already exports coordinates in microns.
conversion.factor <- 0.18

# ── Build CellChat objects per sample ─────────────────────────────────────────
samples       <- c("HB_S2", "HB_S4", "HB_S5")
cellchat_list <- list()

for (nm in samples) {
  obj <- obj_list[[nm]]
  inp <- prep_spatial_input(obj, nm)

  d               <- computeCellDistance(inp$spatial.locs)
  spot.size       <- min(d) * conversion.factor
  spatial.factors <- data.frame(ratio = conversion.factor, tol = spot.size / 2)
  rm(d)

  cc      <- createCellChat(object          = inp$expr_mat,
                            meta            = inp$meta,
                            group.by        = "labels",
                            datatype        = "spatial",
                            coordinates     = inp$spatial.locs,
                            spatial.factors = spatial.factors)
  cc@DB   <- CellChatDB.use
  cc      <- subsetData(cc)
  cc      <- identifyOverExpressedGenes(cc)
  cc      <- identifyOverExpressedInteractions(cc, variable.both = FALSE)

  # trim = 0.1: treat expression as zero if < 10% of cells in a group express
  # the gene; relaxed from the default 25% to suit broad CosMx cell types
  cc <- computeCommunProb(cc,
                          type              = "truncatedMean", trim = 0.1,
                          distance.use      = TRUE, interaction.range = 250,
                          scale.distance    = 4.5,
                          contact.dependent = TRUE, contact.range = 20)

  cc <- filterCommunication(cc, min.cells = 10)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)

  cellchat_list[[nm]] <- cc

  # Interaction circle plots and heatmap
  groupSize <- as.numeric(table(cc@idents))
  par(mfrow = c(1, 2), xpd = TRUE)
  netVisual_circle(cc@net$count,  vertex.weight = groupSize, weight.scale = TRUE,
                   label.edge = FALSE, title.name = "Number of interactions",
                   vertex.label.cex = 0.8)
  netVisual_circle(cc@net$weight, vertex.weight = groupSize, weight.scale = TRUE,
                   label.edge = FALSE, title.name = "Interaction weights/strength",
                   vertex.label.cex = 0.8)
  netVisual_heatmap(cc, measure = "weight", color.heatmap = "Reds")
}

# ── Summary analysis per sample ───────────────────────────────────────────────
cellchat_summary_list <- setNames(
  lapply(names(cellchat_list), function(nm)
    run_cellchat_summary(cellchat_list[[nm]], object_name = nm)),
  names(cellchat_list)
)

# Top 10 interactions per sample
for (nm in names(cellchat_summary_list)) {
  cat("\n--", nm, "top interactions --\n")
  print(cellchat_summary_list[[nm]]$top_interactions[1:10, ])
}

# Source/target frequency plots
for (nm in names(cellchat_summary_list)) {
  print(cellchat_summary_list[[nm]]$source_target_plot)
}

# ── Manuscript figures: pathway chord diagrams per sample ─────────────────────
# Cell colors from the master palette defined in Script 02
cell_colors <- master_colors

# HB_S2 — MK signaling in relapse-associated cell types
cells_HB_S2 <- c("Endothelial cells", "CAF",
                  "AFP+ fetal hepatocyte-like", "Tumor progenitor-like epithelium")
cc_sub <- subsetCellChat(cellchat_list[["HB_S2"]], idents.use = cells_HB_S2)

png(file.path(output_dir, "Figure_MK_chord_HB_S2.png"),
    width = 3000, height = 3000, res = 300)
circlize::circos.clear()
netVisual_aggregate(cc_sub, signaling = "MK", layout = "chord",
                    sources.use = cells_HB_S2, targets.use = cells_HB_S2,
                    color.use = cell_colors,
                    lab.cex = 6, small.gap = 4, big.gap = 10)
dev.off()

# HB_S4 — COLLAGEN signaling in non-relapse-associated cell types
cells_HB_S4 <- c("fibroblast", "antigen presenting TAM",
                  "fibroblast_WNT5A", "Endothelial cells")
cc_sub <- subsetCellChat(cellchat_list[["HB_S4"]], idents.use = cells_HB_S4)

png(file.path(output_dir, "Figure_COLLAGEN_chord_HB_S4.png"),
    width = 3000, height = 3000, res = 300)
circlize::circos.clear()
netVisual_aggregate(cc_sub, signaling = "COLLAGEN", layout = "chord",
                    sources.use = cells_HB_S4, targets.use = cells_HB_S4,
                    color.use = cell_colors,
                    lab.cex = 6, small.gap = 4, big.gap = 10)
dev.off()

# HB_S5 — COLLAGEN signaling in resection-associated cell types
cells_HB_S5 <- c("fibroblast", "antigen presenting TAM",
                  "Tumor progenitor-like epithelium", "Endothelial cells",
                  "inflammatory hepatocyte", "myofibroblast")
cc_sub <- subsetCellChat(cellchat_list[["HB_S5"]], idents.use = cells_HB_S5)

png(file.path(output_dir, "Figure_COLLAGEN_chord_HB_S5.png"),
    width = 3000, height = 3000, res = 300)
circlize::circos.clear()
netVisual_aggregate(cc_sub, signaling = "COLLAGEN", layout = "chord",
                    sources.use = cells_HB_S5, targets.use = cells_HB_S5,
                    color.use = cell_colors,
                    lab.cex = 6, small.gap = 4, big.gap = 10)
dev.off()
