###############################################################################
# Hepatoblastoma Spatial Transcriptomics (Sindhi Cohort)
# Script 06: Relapse Signature Development and Validation
#
# Description:
#   Relapse-associated signatures are developed from CosMx cell-proximity DEGs
#   (FindAllMarkers on cell types identified as spatially co-localized by Giotto,
#   Script 04) and EnrichR pathway gene sets (colocalized cells). Candidate signatures are scored
#   and compared in the Sindhi Visium cohort (discovery, liver_new). The final
#   composite signature — combining Tumor-progenitor DEGs, pathway genes, and
#   CellChat ligand-receptor genes — is then applied to a public CytAssist FFPE
#   cohort (GSE261958) for independent validation.
#
# Input:
#   - CosMx Seurat objects:  obj_list (from Script 02_Preprocessing&UMAPs_CosMx.R)
#   - Discovery Visium:      liver_new (from Script 01_Preprocessing&UMAPs_Visium.R)
#   - Validation data:       GSE261958 raw files (CytAssist FFPE, GEO)
#
# Output: Cell_CloseProx_fGSEAinput.xlsx, EnrichR_01FDR_Cell_CloseProx.xlsx,
#         per-patient DEG tables, signature comparison barplot, per-patient
#         boxplot, violin plots, SpatialFeaturePlots, Validation.RData
###############################################################################

# ── Libraries ─────────────────────────────────────────────────────────────────
library(Seurat)
library(readxl)
library(tidyverse)
library(ggbeeswarm)
library(ggpubr)
library(rstatix)
library(ggsci)
library(R.utils)
library(enrichR)
library(openxlsx)

# ── Required objects from earlier scripts ─────────────────────────────────────
# obj_list  — CosMx Seurat objects per patient (Script 02_Preprocessing&UMAPs_CosMx.R)
# liver_new — Sindhi Visium Seurat object       (Script 01_Preprocessing&UMAPs_Visium.R)
# These must be present in the environment or loaded here, e.g.:
#   load(file.path(data_dir, "obj_list.RData"))
#   load(file.path(data_dir, "liver_new.RData"))

# ── Paths (edit to match your environment) ────────────────────────────────────
data_dir       <- "path/to/data"
output_dir     <- "path/to/output"
base_path      <- file.path(data_dir, "validation")
gse261958_path <- file.path(base_path, "GSE261958")

# ─────────────────────────────────────────────────────────────────────────────
# 1. CosMx DEA: markers for cell types in spatial proximity (from Script 04)
# ─────────────────────────────────────────────────────────────────────────────
# Cell types identified as spatially co-localized by Giotto neighborhood
# enrichment (Script 04); only these are carried forward for signature building
colocalized_cells <- c("Tumor progenitor-like epithelium", "fibroblast",
                       "TAM", "TAM_SPP1+", "stress responsive tumor cells",
                       "myofibroblast", "Neuronal-like")

condition_prefix <- c("HB_S2" = "R", "HB_S4" = "NR", "HB_S5" = "Tx")

# FindAllMarkers per patient, retain only co-localized cell types
all_markers_list <- setNames(
  lapply(names(condition_prefix), function(nm) {
    obj         <- obj_list[[nm]]
    Idents(obj) <- "my.annot"
    FindAllMarkers(obj, only.pos = TRUE,
                   min.pct = 0.1, logfc.threshold = 0.25) %>%
      filter(cluster %in% colocalized_cells)
  }),
  names(condition_prefix)
)

# Save per-patient DEG tables
for (nm in names(all_markers_list)) {
  write.xlsx(data.frame(all_markers_list[[nm]]),
             file = file.path(output_dir, paste0(nm, "_annotDEGs.xlsx")))
}

# Combine into Cell_CloseProx_fGSEAinput.xlsx — one sheet per condition
Cell_CloseProx <- bind_rows(
  lapply(names(all_markers_list), function(nm)
    all_markers_list[[nm]] %>%
      mutate(sheet = paste0(condition_prefix[[nm]], "_cellCloseProx")))
) %>% filter(p_val_adj <= 0.01, avg_log2FC >= 0.5)

wb <- createWorkbook()
for (nm in names(all_markers_list)) {
  sn <- paste0(condition_prefix[[nm]], "_cellCloseProx")
  addWorksheet(wb, sn)
  writeData(wb, sn,
            data.frame(all_markers_list[[nm]] %>%
                         filter(p_val_adj <= 0.01, avg_log2FC >= 0.5)))
}
saveWorkbook(wb, file = file.path(base_path, "Cell_CloseProx_fGSEAinput.xlsx"),
             overwrite = TRUE)

# ─────────────────────────────────────────────────────────────────────────────
# 2. EnrichR pathway analysis on proximity DEGs
# ─────────────────────────────────────────────────────────────────────────────
setEnrichrSite("Enrichr")
dbs <- c("GO_Biological_Process_2023", "KEGG_2021_Human", "Reactome_2022")

enrichr_res <- Cell_CloseProx %>%
  group_split(sheet, cluster) %>%
  map_df(function(df) {
    genes <- df %>% arrange(desc(avg_log2FC)) %>% pull(gene) %>% unique()
    res   <- enrichr(genes, dbs)
    bind_rows(res, .id = "database") %>%
      mutate(sheet   = unique(df$sheet),
             cluster = unique(df$cluster))
  })

enrichr_plot_df <- enrichr_res %>%
  filter(database == "GO_Biological_Process_2023") %>%
  mutate(
    overlap_num = as.numeric(str_extract(Overlap, "^[0-9]+")),
    overlap_den = as.numeric(str_extract(Overlap, "(?<=/)[0-9]+")),
    gene_ratio  = overlap_num / overlap_den,
    log10_padj  = -log10(Adjusted.P.value),
    sheet_short = str_extract(sheet, "^[^_]+"),
    group       = paste(sheet_short, cluster)
  ) %>%
  filter(Adjusted.P.value < 0.01)

write.xlsx(data.frame(enrichr_plot_df),
           file = file.path(base_path, "EnrichR_01FDR_Cell_CloseProx.xlsx"))

# Top 30 terms per group — dot plot
top_terms <- enrichr_plot_df %>%
  group_by(sheet_short, cluster) %>%
  slice_min(Adjusted.P.value, n = 30, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    Term_short = str_remove(Term, "\\s*\\(GO:\\d+\\)") %>%
                 str_replace_all("\n", " "),
    Term_short = fct_reorder(Term_short, log10_padj)
  )

ggplot(top_terms, aes(x = group, y = Term_short)) +
  geom_point(aes(size = log10_padj, color = Combined.Score)) +
  scale_size(range = c(2, 8)) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title  = element_blank())

# ─────────────────────────────────────────────────────────────────────────────
# 3. Build composite relapse signature
# ─────────────────────────────────────────────────────────────────────────────

# Top 50 relapse-specific DEGs per cluster (not shared with NR or Tx)
other_condition_genes <- Cell_CloseProx %>%
  filter(!grepl("^R_", sheet)) %>%
  pull(gene) %>% unique()

relapse_specific <- Cell_CloseProx %>%
  filter(grepl("^R_", sheet), cluster %in% colocalized_cells) %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 50) %>%
  ungroup() %>%
  filter(!gene %in% other_condition_genes)

relapse_keys <- relapse_specific %>%
  arrange(cluster) %>% group_by(cluster) %>% group_keys() %>% pull(cluster)

relapse_genes_by_cluster <- relapse_specific %>%
  arrange(cluster) %>% group_by(cluster) %>% group_split() %>%
  setNames(relapse_keys) %>%
  lapply(function(x) x$gene)

# EnrichR pathway genes for relapse proximity (selected terms from section 2)
keep_terms <- c(
  "Cytokine-Mediated Signaling Pathway", "T Cell Activation",
  "Regulation Of T Cell Activation", "Response To Type II Interferon",
  "Antigen Processing And Presentation Of Peptide Antigen Via MHC Class II",
  "Regulation Of Interleukin-6 Production", "Regulation Of Interleukin-12 Production",
  "Positive Regulation Of Interleukin-8 Production",
  "B Cell Receptor Signaling Pathway",
  "Regulation Of I-kappaB kinase/NF-kappaB Signaling",
  "Positive Regulation Of MAPK Cascade", "Regulation Of ERK1 And ERK2 Cascade",
  "Transforming Growth Factor Beta Receptor Signaling Pathway",
  "Regulation Of Canonical Wnt Signaling Pathway",
  "DNA Damage Response",
  "Fatty Acid Metabolic Process", "Cholesterol Metabolic Process",
  "Sterol Homeostasis", "Oxidative Phosphorylation", "Glucose Homeostasis",
  "Extracellular Matrix Organization", "Collagen Fibril Organization",
  "Extracellular Structure Organization", "Regulation Of Cell Migration",
  "Cellular Response To Hypoxia", "Response To Unfolded Protein",
  "Negative Regulation Of Programmed Cell Death",
  "Epithelial To Mesenchymal Transition",
  "Negative Regulation Of Cell Differentiation",
  "Regulation Of Angiogenesis", "Axon Guidance", "Phagocytosis"
)

enrichr_gene_list <- enrichr_plot_df %>%
  filter(grepl("^R_", sheet)) %>%
  mutate(Term_short = str_remove(Term, "\\s*\\(GO:\\d+\\)") %>%
           str_replace_all("\n", " ")) %>%
  filter(Term_short %in% keep_terms, sheet == "R_cellCloseProx") %>%
  separate_rows(Genes, sep = ";") %>%
  filter(Genes != "") %>%
  group_by(cluster) %>%
  summarise(genes = list(unique(Genes)), .groups = "drop") %>%
  deframe()

# CellChat ligand-receptor genes (MK and COLLAGEN pathways, from Script 05)
CCC_genes <- unique(c(
  "MDK", "FN1", "COL2A1", "CD99", "AGT", "IGF2",
  "COL9A3", "PPIA", "APP", "MIF", "LAMC1",
  "LAMB1", "COL9A2", "CXCL12", "COL4A5", "MPZL1", "F2", "FGF17",
  "VTN", "DLK1", "FGF2", "LAMA5",
  "CDH1", "GRN", "CADM1", "TGFB2", "LRP1", "SDC2", "SDC1",
  "AGTR1", "IGF2R", "BSG", "ITGA6", "ITGB1", "ITGA5", "ITGA2",
  "CD74", "CXCR4", "ACKR3", "PARD3", "FGFR3", "ITGB4", "ITGAV",
  "NOTCH3", "ITGB5", "CD99L2", "SV2A", "FGFR4",
  "SORT1", "F2RL3", "ACVR1B", "TGFBR2"
))

# Final composite signature: Tumor-progenitor DEGs + pathway genes + LR genes
final_sig <- unique(c(
  relapse_genes_by_cluster$`Tumor progenitor-like epithelium`,
  enrichr_gene_list$`Tumor progenitor-like epithelium`,
  CCC_genes
))

# ─────────────────────────────────────────────────────────────────────────────
# 4. Load GSE261958 metadata and assign instrument type
# ─────────────────────────────────────────────────────────────────────────────
GSE261958_meta <- read_excel(file.path(base_path, "ValidationData_meta.xlsx"),
                             sheet = "GSE261958_meta") %>%
  mutate(
    instrument = case_when(
      `Library name` %in% c("H286", "H287", "H288", "H289",
                             "H386", "H387", "H388", "H389",
                             "H391", "H392") ~ "Visium_FFPE",
      TRUE                                   ~ "CytAssist_FFPE"
    )
  )

# Match GSM accession IDs from filenames
get_gsm_id <- function(library_name, path) {
  matched <- list.files(path, pattern = paste0("_", library_name, "_"),
                        full.names = FALSE)
  if (length(matched) == 0) stop(paste("No files found for:", library_name))
  str_extract(matched[1], "^GSM[0-9]+")
}

GSE261958_meta <- GSE261958_meta %>%
  rowwise() %>%
  mutate(gsm_id = get_gsm_id(`Library name`, gse261958_path)) %>%
  ungroup()

# ─────────────────────────────────────────────────────────────────────────────
# 5. Helper: load one sample from already-decompressed GEO files
# ─────────────────────────────────────────────────────────────────────────────
load_one_sample <- function(library_name, gsm_id, path) {
  cat("Loading:", library_name, "\n")

  tmp         <- file.path(tempdir(), library_name)
  spatial_dir <- file.path(tmp, "spatial")
  dir.create(spatial_dir, recursive = TRUE, showWarnings = FALSE)

  sample_files <- list.files(path,
                             pattern    = paste0(gsm_id, "_", library_name, "_"),
                             full.names = TRUE)
  for (f in sample_files) {
    clean_name <- gsub(paste0("^.*", gsm_id, "_", library_name, "_"), "",
                       basename(f))
    clean_name <- gsub("\\.gz$", "", clean_name)
    R.utils::gunzip(f, destname = file.path(tmp, clean_name),
                    overwrite = TRUE, remove = FALSE)
  }

  spatial_files <- c("tissue_positions.csv",       # Space Ranger >= 2.0
                     "tissue_positions_list.csv",  # Space Ranger < 2.0
                     "scalefactors_json.json",
                     "tissue_hires_image.png", "tissue_lowres_image.png",
                     "aligned_fiducials.jpg", "detected_tissue_image.jpg")
  for (sf in spatial_files) {
    src <- file.path(tmp, sf)
    if (file.exists(src)) file.rename(src, file.path(spatial_dir, sf))
  }

  counts <- ReadMtx(mtx            = file.path(tmp, "matrix.mtx"),
                    cells          = file.path(tmp, "barcodes.tsv"),
                    features       = file.path(tmp, "features.tsv"),
                    feature.column = 2)

  obj              <- CreateSeuratObject(counts = counts, assay = "Spatial",
                                         project = library_name)
  obj$library_name <- library_name
  obj$gsm_id       <- gsm_id

  image <- Read10X_Image(image.dir = spatial_dir, filter.matrix = TRUE)
  image <- image[colnames(obj)]
  DefaultAssay(image) <- "Spatial"
  obj[[paste0("slice_", library_name)]] <- image

  unlink(tmp, recursive = TRUE)
  return(obj)
}

# ─────────────────────────────────────────────────────────────────────────────
# 6. Load samples and subset to chemo-treated
# ─────────────────────────────────────────────────────────────────────────────
visium_meta    <- subset(GSE261958_meta, instrument == "Visium_FFPE")
cytassist_meta <- subset(GSE261958_meta, instrument == "CytAssist_FFPE")

load_platform <- function(meta, path) {
  setNames(
    lapply(seq_len(nrow(meta)), function(i)
      load_one_sample(meta$`Library name`[i], meta$gsm_id[i], path)),
    meta$`Library name`
  )
}

visium_list    <- load_platform(visium_meta,    gse261958_path)
cytassist_list <- load_platform(cytassist_meta, gse261958_path)

visium_chemo    <- visium_meta    %>% filter(grepl("chemo-treated", Condition)) %>%
  pull(`Library name`)
cytassist_chemo <- cytassist_meta %>% filter(grepl("chemo-treated", Condition)) %>%
  pull(`Library name`)

# ─────────────────────────────────────────────────────────────────────────────
# 7. Merge within platform and add sample-level metadata
# ─────────────────────────────────────────────────────────────────────────────
merge_platform <- function(sample_list, platform_meta, instrument_label) {
  merged            <- merge(sample_list[[1]], y = sample_list[-1],
                             add.cell.ids = names(sample_list))
  merged$instrument <- instrument_label

  meta_mod              <- merged@meta.data %>%
    left_join(platform_meta %>% select(-gsm_id, -instrument),
              by = c("orig.ident" = "Library name"))
  rownames(meta_mod)    <- colnames(merged)
  merged@meta.data      <- meta_mod
  merged
}

visium_chemo_combined    <- merge_platform(visium_list[visium_chemo],
                                           visium_meta,    "Visium_FFPE")
cytassist_chemo_combined <- merge_platform(cytassist_list[cytassist_chemo],
                                           cytassist_meta, "CytAssist_FFPE")

# ─────────────────────────────────────────────────────────────────────────────
# 8. QC, filtering, and normalization
# ─────────────────────────────────────────────────────────────────────────────
visium_chemo_combined[["percent.mt"]]    <- PercentageFeatureSet(
  visium_chemo_combined,    pattern = "^MT-")
cytassist_chemo_combined[["percent.mt"]] <- PercentageFeatureSet(
  cytassist_chemo_combined, pattern = "^MT-")

VlnPlot(visium_chemo_combined,
        features = c("nCount_Spatial", "nFeature_Spatial", "percent.mt"),
        group.by = "orig.ident", pt.size = 0, ncol = 3)
VlnPlot(cytassist_chemo_combined,
        features = c("nCount_Spatial", "nFeature_Spatial", "percent.mt"),
        group.by = "orig.ident", pt.size = 0, ncol = 3)

# Visium: no mt filter (FFPE RNA integrity differs from fresh-frozen)
visium_chemo_filtered <- subset(visium_chemo_combined,
                                nCount_Spatial > 500 & nFeature_Spatial > 200)

cytassist_chemo_filtered <- subset(cytassist_chemo_combined,
                                   nCount_Spatial   > 500 &
                                   nFeature_Spatial > 200 &
                                   percent.mt       < 20)

# Derive outcome from survival column
cytassist_chemo_filtered@meta.data <- cytassist_chemo_filtered@meta.data %>%
  mutate(outcome = case_when(
    grepl("^D", `Current status - A = alive, D = dead (age in months)`) ~ "Dead",
    grepl("^A", `Current status - A = alive, D = dead (age in months)`) ~ "Alive",
    TRUE ~ "Unknown"
  ))

visium_chemo_filtered    <- SCTransform(visium_chemo_filtered,    assay = "Spatial",
                                        verbose = TRUE)
cytassist_chemo_filtered <- SCTransform(cytassist_chemo_filtered, assay = "Spatial",
                                        verbose = TRUE)

# ─────────────────────────────────────────────────────────────────────────────
# 9. Score candidate signatures in discovery cohort (liver_new, Sindhi Visium)
# ─────────────────────────────────────────────────────────────────────────────
candidate_signatures <- list(
  "Neuronal_like"     = relapse_genes_by_cluster$`Neuronal-like`,
  "Tumor_progenitor"  = relapse_genes_by_cluster$`Tumor progenitor-like epithelium`,
  "Fibroblast"        = relapse_genes_by_cluster$fibroblast,
  "Stress_responsive" = relapse_genes_by_cluster$`stress responsive tumor cells`,
  "Final_Signature"   = final_sig
)

for (sig in names(candidate_signatures)) {
  genes_use <- intersect(candidate_signatures[[sig]], rownames(liver_new))
  liver_new <- AddModuleScore(liver_new, features = list(genes_use), name = sig)
}

# AddModuleScore appends "1"; rename to clean column names
score_cols <- names(candidate_signatures)
names(liver_new@meta.data)[
  match(paste0(score_cols, "1"), names(liver_new@meta.data))
] <- score_cols

# ─────────────────────────────────────────────────────────────────────────────
# 10. Pseudobulk comparison across candidate signatures (discovery cohort)
# ─────────────────────────────────────────────────────────────────────────────
df_discovery <- liver_new@meta.data %>%
  as_tibble(rownames = "cell_id") %>%
  select(cell_id,
         patient = orig.ident,
         outcome = Tumor.relapse,
         all_of(score_cols))

patient_long <- df_discovery %>%
  group_by(patient, outcome) %>%
  summarise(across(all_of(score_cols), median, na.rm = TRUE), .groups = "drop") %>%
  pivot_longer(cols = all_of(score_cols),
               names_to  = "signature",
               values_to = "median_score")

sig_stats <- patient_long %>%
  group_by(signature) %>%
  rstatix::wilcox_test(median_score ~ outcome) %>%
  rstatix::add_significance() %>%
  left_join(
    patient_long %>%
      group_by(signature) %>%
      rstatix::wilcox_effsize(median_score ~ outcome) %>%
      select(signature, effsize),
    by = "signature"
  ) %>%
  arrange(p)

print(sig_stats)

# Barplot: effect sizes across candidate signatures
p_compare <- sig_stats %>%
  mutate(is_final  = str_detect(signature, "Final"),
         sig_label = paste0(signature, "\n(p=", round(p, 3), ")")) %>%
  ggplot(aes(x = reorder(sig_label, effsize), y = effsize, fill = is_final)) +
  geom_col(width = 0.6, color = "grey30", linewidth = 0.4) +
  geom_text(aes(label = p.signif), hjust = -0.3, size = 5) +
  scale_fill_manual(values = c("FALSE" = "grey75", "TRUE" = "#D6604D"),
                    labels = c("Candidate signatures", "Final signature")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  coord_flip() +
  labs(x = NULL, y = "Effect size", fill = NULL,
       title = "Signature performance — Discovery cohort") +
  theme_classic(base_size = 14) +
  theme(legend.position         = "inside",
        legend.position.inside  = c(0.75, 0.15),
        legend.justification    = c(1, 0),
        axis.text               = element_text(size = 10, color = "black"),
        plot.title              = element_text(size = 12, face = "bold", hjust = 0.5),
        panel.grid.major.x      = element_line(color = "grey92", linewidth = 0.4))

png(file.path(output_dir, "Sig_discovery_comparison.png"),
    res = 300, width = 6, height = 5, units = "in")
print(p_compare)
dev.off()

# Per-patient boxplot for the final signature
final_data <- patient_long %>%
  filter(signature == "Final_Signature") %>%
  mutate(outcome = factor(outcome, levels = c("No", "Yes")))

stat_final <- final_data %>%
  rstatix::wilcox_test(median_score ~ outcome) %>%
  add_significance() %>%
  add_xy_position(x = "outcome")

p_final <- ggplot(final_data, aes(x = outcome, y = median_score, fill = outcome)) +
  geom_boxplot(width = 0.45, outlier.shape = NA,
               alpha = 0.6, color = "grey30", linewidth = 0.5) +
  geom_beeswarm(aes(color = outcome), size = 3, cex = 2.5) +
  stat_pvalue_manual(stat_final, tip.length = 0.02, size = 3.5) +
  scale_fill_manual(values  = c("No" = "#4393C3", "Yes" = "#D6604D")) +
  scale_color_manual(values = c("No" = "#4393C3", "Yes" = "#D6604D")) +
  labs(x = NULL, y = "Relapse signature score\n(per-patient median)",
       title = "Final signature — Discovery cohort") +
  theme_classic(base_size = 12) +
  theme(legend.position    = "none",
        axis.text          = element_text(size = 11, color = "black"),
        plot.title         = element_text(size = 12, face = "bold", hjust = 0.5),
        panel.grid.major.y = element_line(color = "grey92", linewidth = 0.4))

print(p_final)

# ─────────────────────────────────────────────────────────────────────────────
# 11. Score final signature in validation cohort (CytAssist FFPE, GSE261958)
# ─────────────────────────────────────────────────────────────────────────────
liver_new                <- UpdateSeuratObject(liver_new)
liver_new@images         <- lapply(liver_new@images, UpdateSeuratObject)
cytassist_chemo_filtered <- UpdateSeuratObject(cytassist_chemo_filtered)

for (obj_name in c("visium_chemo_filtered", "cytassist_chemo_filtered")) {
  obj <- get(obj_name)
  obj <- AddModuleScore(obj, features = list(final_sig), name = "Relapse_signature")
  assign(obj_name, obj)
}

# ─────────────────────────────────────────────────────────────────────────────
# 12. Visualization: spatial plots and discovery vs. validation violin plots
# ─────────────────────────────────────────────────────────────────────────────
SpatialFeaturePlot(liver_new,                features = "Final_Signature",
                   pt.size.factor = 2000)
SpatialFeaturePlot(visium_chemo_filtered,    features = "Relapse_signature1",
                   pt.size.factor = 3000, ncol = 4)
SpatialFeaturePlot(cytassist_chemo_filtered, features = "Relapse_signature1",
                   pt.size.factor = 4,    ncol = 4)

vln_theme <- theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

V1 <- VlnPlot(liver_new,                features = "Final_Signature",    pt.size = 0) +
  scale_fill_npg() + vln_theme + ggtitle("Discovery")
V2 <- VlnPlot(cytassist_chemo_filtered, features = "Relapse_signature1", pt.size = 0) +
  scale_fill_npg() + vln_theme + ggtitle("Validation")

png(file.path(output_dir, "discovery_validation_relapse_scores_vln.png"),
    res = 300, width = 5, height = 4, units = "in")
print(V1 + V2)
dev.off()

# ── Save ──────────────────────────────────────────────────────────────────────
save(liver_new, cytassist_chemo_filtered,
     file = file.path(base_path, "Validation.RData"))
