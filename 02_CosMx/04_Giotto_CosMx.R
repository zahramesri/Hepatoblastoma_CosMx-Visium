###############################################################################
# CosMx SMI - Hepatoblastoma (Sindhi Cohort)
# Script 04: Neighborhood Enrichment Analysis using Giotto
#
# Reference:
#   http://giottosuite.com/articles/nanostring_cosmx_lung_cancer.html
#
# Description:
#   Builds a Delaunay spatial network from CosMx cell centroids per patient
#   and computes cell-type proximity enrichment scores (1000 permutations,
#   FDR-adjusted). Results are saved as heatmaps and Excel tables.
#   Cross-FOV edges in the spatial network are reported for QC.
#
# Input:  Annotated per-patient CosMx Seurat objects (obj_list, from Script 02)
# Output: Proximity enrichment heatmaps and result tables per patient
#
# Note: installGiottoEnvironment() must be run once before first use.
###############################################################################

# ── Libraries ─────────────────────────────────────────────────────────────────
library(Giotto)
library(Seurat)
library(openxlsx)

# ── Paths (edit to match your environment) ────────────────────────────────────
data_dir   <- "path/to/data"
output_dir <- "path/to/output"

# ── Giotto instructions ────────────────────────────────────────────────────────
instrs <- createGiottoInstructions(show_plot = FALSE)

# ── Neighborhood enrichment per patient ───────────────────────────────────────
for (nm in names(obj_list)) {
  obj <- obj_list[[nm]]

  # Build spatial locations from global cell centroids
  spatial_locs <- data.frame(
    cell_ID = colnames(obj),
    sdimx   = as.numeric(obj$CenterX_global_px),
    sdimy   = as.numeric(obj$CenterY_global_px),
    stringsAsFactors = FALSE
  )

  # Expression matrix (Nanostring counts)
  expr_mat <- GetAssayData(obj, assay = "Nanostring", slot = "counts")
  expr_mat <- expr_mat[, spatial_locs$cell_ID]

  # Cell metadata
  cell_metadata <- data.frame(
    cell_ID   = colnames(obj),
    cell_type = obj$my.annot,
    fov       = obj$fov,
    patient_id = obj$Patient_ID,
    condition  = obj$Tumor_recurrence,
    stringsAsFactors = FALSE
  )

  # Create Giotto object and spatial network
  gobj <- createGiottoObject(raw_exprs     = expr_mat,
                             spatial_locs  = spatial_locs,
                             cell_metadata = cell_metadata,
                             instructions  = instrs)

  gobj <- createSpatialNetwork(gobject = gobj, method = "Delaunay")

  # Neighborhood enrichment (1000 permutations, FDR-adjusted)
  prox <- cellProximityEnrichment(gobject              = gobj,
                                  cluster_column       = "cell_type",
                                  spatial_network_name = "Delaunay_network",
                                  number_of_simulations = 1000,
                                  adjust_method        = "fdr")

  # Save heatmap
  png(file.path(output_dir, paste0("CosMx_NeighborhoodEnrichment_", nm, ".png")),
      width = 7, height = 7, units = "in", res = 300, bg = "white")
  par(cex.axis = 1.5, cex.lab = 1.5, mar = c(10, 10, 4, 4))
  cellProximityHeatmap(gobject = gobj, CPscore = prox)
  dev.off()

  # Network plot
  cellProximityNetwork(gobject = gobj, CPscore = prox)

  # QC: count cross-FOV edges (should be zero in a valid spatial network)
  spatnet          <- gobj@spatial_network$Delaunay_network
  meta             <- pDataDT(gobj)[, c("cell_ID", "fov")]
  spatnet$fov_from <- meta$fov[match(spatnet$from, meta$cell_ID)]
  spatnet$fov_to   <- meta$fov[match(spatnet$to,   meta$cell_ID)]
  cat(nm, "- cross-FOV edges:", sum(spatnet$fov_from != spatnet$fov_to), "\n")

  # Save enrichment results
  write.xlsx(data.frame(prox$enrichm_res),
             file = file.path(output_dir,
                              paste0(nm, "_GiottoNeighborhoodEnrichment.xlsx")))
}
