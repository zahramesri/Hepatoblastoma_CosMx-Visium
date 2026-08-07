###############################################################################
# Spatial Transcriptomics - Hepatoblastoma (Sindhi Cohort)
# Script 02: STdeconvolve — Unsupervised Cell-Type Deconvolution
#
# Description:
#   Applies STdeconvolve (LDA-based) to deconvolve cell-type proportions
#   from spatial transcriptomics count data across all 11 samples.
#   Extracts per-spot topic proportions (theta) and topic gene expression
#   profiles (beta), exports top marker genes per topic per sample, and
#   compares deconvolution results to transcriptional clusters for one
#   representative sample.
#
# Input:  liver_new — merged Seurat object from Script 01
# Output: LDA models, theta/beta results, LDA pie plots, marker gene tables
#
# Reference:
#   https://github.com/JEFworks-Lab/STdeconvolve/blob/devel/docs/visium_10x.md
###############################################################################

# ── Libraries ─────────────────────────────────────────────────────────────────
library(STdeconvolve)
library(Seurat)
library(ggplot2) #"ggplot2", version = "3.5.1"
library(dplyr)
library(openxlsx)

# ── Paths (edit to match your environment) ────────────────────────────────────
output_dir <- "path/to/output"

# ── Extract per-sample counts and spatial coordinates ─────────────────────────
# Assumes liver_new is already loaded (output of Script 01)
Idents(liver_new)  <- "SampleID"
all_samples        <- levels(liver_new)
images_sample      <- names(liver_new@images)

coords <- list()
count  <- list()

for (i in seq_along(all_samples)) {
  obj       <- subset(liver_new, SampleID == all_samples[[i]])
  count[[i]] <- obj@assays$Spatial@counts
  coords_    <- GetTissueCoordinates(obj, image = images_sample[[i]])
  colnames(coords_) <- c("x", "y")
  coords[[i]] <- coords_
}
rm(obj, coords_)
gc()

names(count)  <- all_samples
names(coords) <- all_samples

# ── Clean counts and build corpus ─────────────────────────────────────────────
counts_cleaned <- list()
corpus         <- list()

for (i in seq_along(all_samples)) {
  count_cleaned   <- cleanCounts(count[[i]], min.lib.size = 100, min.reads = 10)
  corpus_          <- restrictCorpus(count_cleaned, removeAbove = 1.0,
                                     removeBelow = 0.05, nTopOD = 1000)
  counts_cleaned[[i]] <- count_cleaned
  corpus[[i]]         <- corpus_
}

names(counts_cleaned) <- all_samples
names(corpus)         <- all_samples

save(file = file.path(output_dir, "00_STdeconvolve_corpus_count.RData"),
     count, coords, corpus)
save(file = file.path(output_dir, "00_STdeconvolve_countsCleaned.RData"),
     counts_cleaned)

# ── Fit LDA models ─────────────────────────────────────────────────────────────
ldas <- list()

for (i in seq_along(all_samples)) {
  ldas[[i]] <- fitLDA(t(as.matrix(corpus[[i]])), Ks = c(2:10),
                      perc.rare.thresh = 0.1)
}

names(ldas) <- all_samples

save(file = file.path(output_dir, "01_STdeconvolve_ldas.RData"), ldas)

# ── Extract optimal model results and visualize ────────────────────────────────
results <- list()
plt     <- list()

for (i in seq_along(all_samples)) {
  optLDA     <- optimalModel(models = ldas[[i]], opt = 5)
  results[[i]] <- getBetaTheta(optLDA, perc.filt = 0.05, betaScale = 1000)

  deconProp <- results[[i]]$theta
  pos       <- coords[[i]]

  plt[[i]] <- vizAllTopics(
    theta     = deconProp,
    pos       = pos,
    r         = 4,
    lwd       = 0,
    showLegend = TRUE,
    plotTitle  = NA
  ) +
    ggplot2::guides(fill = ggplot2::guide_legend(ncol = 2)) +
    ggplot2::geom_rect(
      data = data.frame(pos),
      ggplot2::aes(xmin = min(x) - 90, xmax = max(x) + 90,
                   ymin = min(y) - 90, ymax = max(y) + 90),
      fill = NA, color = "black", linetype = "solid", linewidth = 0.5
    ) +
    ggplot2::theme(plot.background = ggplot2::element_blank()) +
    ggtitle(all_samples[[i]]) +
    ggplot2::guides(colour = "none")
}

names(results) <- all_samples

save(file = file.path(output_dir, "02_STdeconvolve_results_theta_beta.RData"), results)

pdf(file.path(output_dir, "LDApie_allsamples.pdf"), width = 20, height = 20)
for (i in seq_along(all_samples)) print(plt[[i]])
dev.off()

# ── Export top marker genes per topic per sample ───────────────────────────────
marker_list <- list()

for (sample in all_samples) {
  gexp <- as.matrix(results[[sample]]$beta)
  marker_list[[sample]] <- list()

  for (i in seq_len(nrow(gexp))) {
    highgexp <- names(which(gexp[i, ] > 10))

    if (length(highgexp) > 0 && nrow(gexp) > 1) {
      log2fc     <- log2(gexp[i, highgexp] /
                           colMeans(gexp[-i, highgexp, drop = FALSE]))
      log2fc     <- sort(log2fc, decreasing = TRUE)
      top_log2fc <- head(log2fc, 200)
      markers    <- names(top_log2fc[top_log2fc > 1])
    } else {
      markers <- character(0)
    }

    marker_list[[sample]][[ rownames(gexp)[i] ]] <- markers
  }
}

wb <- createWorkbook()

for (sample in all_samples) {
  celltype_list <- marker_list[[sample]]
  max_len       <- max(sapply(celltype_list, length))

  df            <- data.frame(matrix(ncol = length(celltype_list), nrow = max_len))
  colnames(df)  <- paste0("Topic_", names(celltype_list))

  for (i in seq_along(celltype_list)) {
    df[seq_along(celltype_list[[i]]), i] <- celltype_list[[i]]
  }

  addWorksheet(wb, sheetName = sample)
  writeData(wb, sheet = sample, df)
}

saveWorkbook(wb, file = file.path(output_dir, "TopMarkers_bySample.xlsx"),
             overwrite = TRUE)

# ── Compare to transcriptional clustering (all samples) ───────────────────────
txn_cluster_plots  <- list()
correlation_plots  <- list()
seurat_objs        <- list()
txn_cluster_markers <- list()

for (sample in all_samples) {

  # PCA on cleaned counts
  pcs_info <- stats::prcomp(
    t(log10(as.matrix(counts_cleaned[[sample]]) + 1)),
    center = TRUE
  )
  pcs <- pcs_info$x[, 1:5]

  # t-SNE embedding
  emb <- Rtsne::Rtsne(pcs, is_distance = FALSE, perplexity = 30,
                      num_threads = 1, verbose = FALSE)$Y
  rownames(emb) <- rownames(pcs)
  colnames(emb) <- c("x", "y")

  # Louvain clustering
  com <- MERINGUE::getClusters(pcs, k = 800, weight = TRUE,
                               method = igraph::cluster_louvain)

  # Plot clusters in tissue space
  dat_space <- data.frame(emb1    = coords[[sample]][, "x"],
                          emb2    = coords[[sample]][, "y"],
                          Cluster = com)

  txn_cluster_plots[[sample]][["tissue"]] <-
    ggplot2::ggplot(dat_space, ggplot2::aes(x = emb1, y = emb2, color = Cluster)) +
    ggplot2::geom_point(size = 0.8) +
    ggplot2::scale_color_manual(values = rainbow(length(levels(com)))) +
    ggplot2::labs(title = sample, x = "x", y = "y") +
    ggplot2::theme_classic() +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 2), ncol = 2)) +
    ggplot2::coord_equal()

  # Plot clusters in t-SNE space
  cent_pos <- do.call(rbind, tapply(seq_len(nrow(emb)), com, function(ii) {
    apply(emb[ii, , drop = FALSE], 2, median)
  }))
  cent_pos           <- as.data.frame(cent_pos)
  colnames(cent_pos) <- c("x", "y")
  cent_pos$cluster   <- rownames(cent_pos)
  cent_pos           <- na.omit(cent_pos)

  dat_tsne <- data.frame(emb1 = emb[, 1], emb2 = emb[, 2], Cluster = com)

  txn_cluster_plots[[sample]][["tsne"]] <-
    ggplot2::ggplot(dat_tsne, ggplot2::aes(x = emb1, y = emb2, color = Cluster)) +
    ggplot2::geom_point(size = 0.01) +
    ggplot2::scale_color_manual(values = rainbow(length(levels(com)))) +
    ggplot2::scale_y_continuous(expand = c(0, 0),
                                 limits = c(min(dat_tsne$emb2) - 1, max(dat_tsne$emb2) + 1)) +
    ggplot2::scale_x_continuous(expand = c(0, 0),
                                 limits = c(min(dat_tsne$emb1) - 1, max(dat_tsne$emb1) + 1)) +
    ggplot2::labs(title = sample, x = "t-SNE 1", y = "t-SNE 2") +
    ggplot2::theme_classic() +
    ggplot2::geom_text(data = cent_pos,
                       ggplot2::aes(x = x, y = y, label = cluster),
                       fontface = "bold") +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 2), ncol = 2)) +
    ggplot2::coord_equal()

  # Correlation between transcriptional clusters and STdeconvolve topics
  com_proxyTheta           <- model.matrix(~ 0 + com)
  rownames(com_proxyTheta) <- names(com)
  colnames(com_proxyTheta) <- gsub("^com", "", colnames(com_proxyTheta))
  com_proxyTheta           <- as.data.frame.matrix(com_proxyTheta)

  corMat           <- STdeconvolve::getCorrMtx(m1 = as.matrix(com_proxyTheta),
                                               m2 = as.matrix(results[[sample]]$theta),
                                               type = "t")
  rownames(corMat) <- paste0("com_",   seq_len(nrow(corMat)))
  colnames(corMat) <- paste0("decon_", seq_len(ncol(corMat)))

  pairs  <- STdeconvolve::lsatPairs(corMat)
  m      <- corMat[pairs$rowix, pairs$colsix]

  correlation_plots[[sample]] <-
    STdeconvolve::correlationPlot(mat = m,
                                  colLabs = "Transcriptional clusters",
                                  rowLabs = "STdeconvolve") +
    ggplot2::labs(title = sample) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90))

  # Map transcriptional clusters onto Seurat object for marker gene analysis
  obj_sample           <- subset(liver_new, SampleID == sample)
  com_df               <- data.frame(tempCom = com)
  sample_meta          <- merge(obj_sample@meta.data, com_df, by = "row.names")
  sample_meta          <- column_to_rownames(sample_meta, var = "Row.names")
  obj_sample@meta.data <- sample_meta
  seurat_objs[[sample]] <- obj_sample

  Idents(obj_sample) <- "tempCom"
  txn_cluster_markers[[sample]] <- FindAllMarkers(obj_sample, only.pos = TRUE,
                                                   min.pct = 0.25, logfc.threshold = 0.25,
                                                   verbose = FALSE, recorrect_umi = FALSE)
}

# Save results
save(file = file.path(output_dir, "03_TSNEmodObjects.RData"), seurat_objs)

# Write plots to PDFs
pdf(file.path(output_dir, "TxnClusters_tissue_allsamples.pdf"), width = 10, height = 10)
for (sample in all_samples) print(txn_cluster_plots[[sample]][["tissue"]])
dev.off()

pdf(file.path(output_dir, "TxnClusters_tSNE_allsamples.pdf"), width = 10, height = 10)
for (sample in all_samples) print(txn_cluster_plots[[sample]][["tsne"]])
dev.off()

pdf(file.path(output_dir, "TxnClusters_vs_STdeconvolve_correlation_allsamples.pdf"),
    width = 10, height = 10)
for (sample in all_samples) print(correlation_plots[[sample]])
dev.off()

# Write marker genes to Excel — one sheet per sample
wb_txn <- createWorkbook()
for (sample in all_samples) {
  addWorksheet(wb_txn, sheetName = sample)
  writeData(wb_txn, sheet = sample, data.frame(txn_cluster_markers[[sample]]))
}
saveWorkbook(wb_txn, file = file.path(output_dir, "AllMarkers_tSNEclusters_allsamples.xlsx"),
             overwrite = TRUE)

# ── Cluster annotation ─────────────────────────────────────────────────────────
# tempCom clusters were annotated based on top marker genes per cluster
# (AllMarkers_tSNEclusters_allsamples.xlsx) and STdeconvolve topic markers
# (Supplementary Table S4). Annotated Seurat objects are used as input for
# downstream analysis.

