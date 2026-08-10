# ============================================================
# Fig. 2b exploration:
# Patient-level fractions of ALL Level 1 major cell types
# Match Fig. 2a disease order and color system
# ============================================================

setwd("data/analysis/atlas_objects")

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(grid)
})

# -------------------------
# 1. Load object
# -------------------------
if (!exists("seurat_obj")) {
  seurat_obj <- readRDS("integrated_pediatric.rds")
}

meta <- as.data.frame(seurat_obj@meta.data)

# -------------------------
# 2. Column settings
# -------------------------
CELL_COL    <- "majorCluster"
DISEASE_COL <- "disease_Level_3"

# 优先按 patient_id 算
# 如果后面你想按 sample 算，可以把 patient_id 换成 sampleID
stopifnot(CELL_COL %in% colnames(meta))
stopifnot(DISEASE_COL %in% colnames(meta))
stopifnot("patient_id" %in% colnames(meta))

# 构建全局唯一 patient ID，避免不同数据集中 patient_id 重名
if ("dataset_id" %in% colnames(meta)) {
  meta$unit_id_global <- paste0(
    trimws(as.character(meta$dataset_id)),
    "__",
    trimws(as.character(meta$patient_id))
  )
} else if ("sampleID" %in% colnames(meta)) {
  meta$unit_id_global <- paste0(
    trimws(as.character(meta$sampleID)),
    "__",
    trimws(as.character(meta$patient_id))
  )
} else {
  meta$unit_id_global <- trimws(as.character(meta$patient_id))
  warning("dataset_id and sampleID were not found. Using raw patient_id only.")
}

UNIT_COL <- "unit_id_global"

# -------------------------
# 3. Disease mapping from Fig. 2a
# -------------------------
disease_short_map <- c(
  "Dengue human infection(DENV-1)" = "Dengue",
  "HFMD(EV-A71)enterovirus" = "HFMD",
  "Influenza A virus (IAV)" = "Flu",
  "Kawasaki disease (KD)" = "KD",
  "Multisystem inflammatory syndrome in children (MIS-C)" = "MIS-C",
  "Pediatric Sepsis" = "Sepsis",
  "bacterial infection" = "Bact.inf",
  "human systemic lupus erythematosus (SLE)" = "SLE",
  "juvenile idopathic arthritis (JIA)" = "JIA",
  "severe acute respiratory syndrome coronavirus 2 (SARS-CoV-2)" = "COVID"
)

disease_order <- c(
  "Dengue", "HFMD", "Flu", "COVID",
  "Bact.inf", "Sepsis", "KD", "MIS-C", "SLE", "JIA"
)

disease_group_map <- c(
  "Dengue"   = "virus",
  "HFMD"     = "virus",
  "Flu"      = "virus",
  "COVID"    = "virus",
  "Bact.inf" = "bacteria",
  "Sepsis"   = "Sepsis",
  "KD"       = "Hyperinflam",
  "MIS-C"    = "Hyperinflam",
  "SLE"      = "Autoimmune",
  "JIA"      = "Autoimmune"
)

group_levels <- c("virus", "bacteria", "Sepsis", "Hyperinflam", "Autoimmune")

disease_cols <- c(
  "Dengue"   = "#2F6FB7",
  "HFMD"     = "#6699CC",
  "Flu"      = "#8BBCE5",
  "COVID"    = "#1F5AA6",
  "Bact.inf" = "#C48EAD",
  "Sepsis"   = "#FF2A2A",
  "KD"       = "#E3A722",
  "MIS-C"    = "#F2C84B",
  "SLE"      = "#36A596",
  "JIA"      = "#168B7A"
)

# -------------------------
# 4. Clean metadata
# -------------------------
meta2 <- meta %>%
  dplyr::mutate(
    disease_raw = trimws(as.character(.data[[DISEASE_COL]])),
    major_raw   = trimws(as.character(.data[[CELL_COL]])),
    unit_id     = trimws(as.character(.data[[UNIT_COL]])),
    disease     = dplyr::recode(disease_raw, !!!disease_short_map, .default = disease_raw),
    disease_group = dplyr::recode(disease, !!!disease_group_map, .default = "Unknown")
  ) %>%
  dplyr::filter(
    !is.na(unit_id), unit_id != "",
    !is.na(disease), disease %in% disease_order,
    !is.na(major_raw), major_raw != ""
  ) %>%
  dplyr::mutate(
    disease = factor(as.character(disease), levels = disease_order),
    disease_group = factor(as.character(disease_group), levels = group_levels),
    major_raw = as.character(major_raw),
    unit_id = as.character(unit_id)
  ) %>%
  tibble::as_tibble()

cat("\nCheck disease counts by cells:\n")
print(table(meta2$disease, useNA = "ifany"))

cat("\nCheck majorCluster names:\n")
print(sort(unique(meta2$major_raw)))

# -------------------------
# 5. All Level 1 cell types from Fig. 2a
# -------------------------
level1_order <- c(
  "B",
  "Mono",
  "pDC",
  "Plasma",
  "DC",
  "T CD4 naive",
  "T CD4 non-naive",
  "T CD8 naive",
  "T CD8 non-naive",
  "UTC",
  "ILC"
)

# 只保留数据里真实存在的
level1_order <- level1_order[level1_order %in% unique(meta2$major_raw)]

label_map <- c(
  "B" = "B cells",
  "Mono" = "Monocytes",
  "pDC" = "pDCs",
  "Plasma" = "Plasma cells",
  "DC" = "DCs",
  "T CD4 naive" = "CD4 T naive",
  "T CD4 non-naive" = "CD4 T non-naive",
  "T CD8 naive" = "CD8 T naive",
  "T CD8 non-naive" = "CD8 T non-naive",
  "UTC" = "UTCs",
  "ILC" = "ILCs"
)

pretty_label <- function(x) {
  ifelse(x %in% names(label_map), unname(label_map[x]), x)
}

selected_order <- pretty_label(level1_order)

cat("\nLevel 1 cell types to plot:\n")
print(data.frame(
  major_raw = level1_order,
  label = selected_order
))

# -------------------------
# 6. Patient-level fractions
# denominator = all PBMCs from each patient
# numerator = each Level 1 cell type count
# -------------------------
total_df <- meta2 %>%
  dplyr::count(unit_id, disease, disease_group, name = "total_cells")

# 检查一个 unit_id 是否对应多个疾病
unit_disease_check <- total_df %>%
  dplyr::count(unit_id, name = "n_disease") %>%
  dplyr::filter(n_disease > 1)

if (nrow(unit_disease_check) > 0) {
  warning("Some unit_id values map to more than one disease. Check patient_id collision.")
  print(head(unit_disease_check, 20))
}

count_df <- meta2 %>%
  dplyr::filter(major_raw %in% level1_order) %>%
  dplyr::count(unit_id, disease, disease_group, major_raw, name = "n")

frac_df <- total_df %>%
  tidyr::expand(
    tidyr::nesting(unit_id, disease, disease_group, total_cells),
    major_raw = level1_order
  ) %>%
  dplyr::left_join(
    count_df,
    by = c("unit_id", "disease", "disease_group", "major_raw")
  ) %>%
  dplyr::mutate(
    n = tidyr::replace_na(n, 0L),
    fraction = n / total_cells * 100,
    celltype_plot = pretty_label(major_raw),
    celltype_plot = factor(celltype_plot, levels = selected_order),
    disease = factor(as.character(disease), levels = disease_order),
    disease_group = factor(as.character(disease_group), levels = group_levels)
  )

cat("\nNumber of patients per disease:\n")
print(
  frac_df %>%
    dplyr::distinct(unit_id, disease) %>%
    dplyr::count(disease)
)

cat("\nCheck fraction summary:\n")
print(summary(frac_df$fraction))

# -------------------------
# 7. Kruskal-Wallis p values for all plotted Level 1 types
# -------------------------
format_p <- function(p) {
  if (is.na(p)) {
    return("P = NA")
  } else if (p < 1e-4) {
    return(paste0("P = ", formatC(p, format = "e", digits = 1)))
  } else {
    return(paste0("P = ", signif(p, 2)))
  }
}

p_df <- frac_df %>%
  dplyr::group_by(major_raw, celltype_plot) %>%
  dplyr::summarise(
    p_value = tryCatch(
      stats::kruskal.test(fraction ~ disease)$p.value,
      error = function(e) NA_real_
    ),
    q_value = NA_real_,
    median_fraction = median(fraction, na.rm = TRUE),
    mean_fraction = mean(fraction, na.rm = TRUE),
    max_fraction = max(fraction, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    q_value = p.adjust(p_value, method = "BH"),
    p_label = vapply(p_value, format_p, character(1)),
    celltype_plot = factor(as.character(celltype_plot), levels = selected_order)
  ) %>%
  dplyr::arrange(p_value)

cat("\nKruskal-Wallis results for all Level 1 cell types:\n")
print(p_df)

# -------------------------
# 8. Plot function
# -------------------------
PSEUDO <- 1e-3

plot_one_celltype <- function(ct, show_y_title = FALSE) {

  df <- frac_df %>%
    dplyr::filter(celltype_plot == ct) %>%
    dplyr::mutate(
      y_plot = ifelse(fraction <= 0, PSEUDO, fraction)
    )

  p_lab <- p_df %>%
    dplyr::filter(celltype_plot == ct) %>%
    dplyr::pull(p_label)

  if (length(p_lab) == 0) {
    p_lab <- "P = NA"
  }

  y_lab <- if (show_y_title) {
    "Fraction of PBMCs (%)"
  } else {
    NULL
  }

  ggplot2::ggplot(
    df,
    ggplot2::aes(x = disease, y = y_plot, color = disease, fill = disease)
  ) +
    ggplot2::geom_boxplot(
      width = 0.55,
      outlier.shape = NA,
      alpha = 0.25,
      linewidth = 0.35
    ) +
    ggplot2::geom_jitter(
      width = 0.14,
      size = 0.65,
      alpha = 0.70,
      stroke = 0
    ) +
    ggplot2::annotate(
      "text",
      x = 1,
      y = 0.0015,
      label = p_lab[1],
      hjust = 0,
      vjust = 0,
      size = 2.6
    ) +
    ggplot2::scale_y_log10(
      limits = c(PSEUDO, 100),
      breaks = c(1e-3, 1e-2, 1e-1, 1, 10, 100),
      labels = c("0.001", "0.01", "0.1", "1", "10", "100")
    ) +
    ggplot2::scale_color_manual(
      values = disease_cols,
      breaks = disease_order,
      drop = FALSE
    ) +
    ggplot2::scale_fill_manual(
      values = disease_cols,
      breaks = disease_order,
      drop = FALSE
    ) +
    ggplot2::labs(
      title = ct,
      x = NULL,
      y = y_lab,
      color = NULL,
      fill = NULL
    ) +
    ggplot2::theme_classic(base_size = 9) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 9.5),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1, size = 6.5),
      axis.text.y = ggplot2::element_text(size = 6.5),
      axis.title.y = ggplot2::element_text(size = 8),
      axis.line = ggplot2::element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = ggplot2::element_line(linewidth = 0.30, colour = "black"),
      legend.position = "bottom",
      legend.text = ggplot2::element_text(size = 8),
      legend.key.size = grid::unit(3.0, "mm"),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.margin = ggplot2::margin(4, 4, 4, 4)
    ) +
    ggplot2::guides(
      fill = "none",
      color = ggplot2::guide_legend(
        nrow = 1,
        override.aes = list(size = 3, alpha = 1)
      )
    )
}

# -------------------------
# 9. Generate all panels
# -------------------------
NCOL <- 4

plots <- lapply(seq_along(selected_order), function(i) {
  plot_one_celltype(
    selected_order[i],
    show_y_title = ((i - 1) %% NCOL == 0)
  )
})

p_combined <- patchwork::wrap_plots(
  plots,
  ncol = NCOL,
  guides = "collect"
) +
  patchwork::plot_annotation(
    title = "Patient-level fractions of all major immune compartments",
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14, hjust = 0.5)
    )
  ) &
  ggplot2::theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    plot.background = ggplot2::element_rect(fill = "white", colour = NA)
  )

# -------------------------
# 10. Save
# -------------------------
outdir <- file.path(getwd(), "fig2_redraw")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

pdf_file <- file.path(outdir, "Fig2b_patient_level_all_Level1_cell_fractions_logY.pdf")
png_file <- file.path(outdir, "Fig2b_patient_level_all_Level1_cell_fractions_logY.png")

ggplot2::ggsave(
  filename = pdf_file,
  plot = p_combined,
  width = 10.5,
  height = 7.2,
  device = grDevices::cairo_pdf
)

ggplot2::ggsave(
  filename = png_file,
  plot = p_combined,
  width = 10.5,
  height = 7.2,
  dpi = 600,
  bg = "white"
)

write.csv(
  frac_df,
  file.path(outdir, "Fig2b_patient_level_all_Level1_cell_fractions_data.csv"),
  row.names = FALSE
)

write.csv(
  p_df,
  file.path(outdir, "Fig2b_patient_level_all_Level1_kruskal_pvalues.csv"),
  row.names = FALSE
)

message("[DONE] Saved all-Level1 Fig. 2b exploration to: ", outdir)








#########fig2补充图############

getwd()


mat_fq_raw <- freq_calculate(meta)


mat_fq_norm<-freq_normalize(mat_fq_raw,normalize="minmax")
dim(mat_fq_norm)
res <- nmf(mat_fq_norm, rank = 2:20, method = "nsNMF", seed = rep(123456, 6), .options = "vp")

plot(res)


pdf("Extended_Data_Fig1c.pdf", width = 15, height = 10)
plot(res)
dev.off()









# ============================================================
# Supplementary analysis based on existing K = 4 NMF model
# Do NOT rerun NMF
# Use the same K4 object as current Fig. 2d/e/f
# ============================================================

suppressPackageStartupMessages({
  library(NMF)
  library(dplyr)
  library(igraph)
  library(CoVarNet)
})

# ------------------------------------------------------------
# 0. Output directories
# ------------------------------------------------------------

# 这里必须和你现在 Fig.2d/e/f 的 outdir 一致
# 你现在 fig2de.R 用的是 fig2_redraw2
outdir <- file.path(getwd(), "fig2_redraw2")

dir.create(file.path(outdir, "NMF_K4_supplement"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "K4_network"), recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 1. Read existing K4 NMF model
# ------------------------------------------------------------

nmf_rds <- file.path(outdir, "NMF_K4_cellular_modules_CM01_CM04.rds")

if (!file.exists(nmf_rds)) {
  stop(
    "Cannot find existing K4 NMF object:\n",
    nmf_rds,
    "\nPlease first run your current Fig.2d/e/f script to generate this file."
  )
}

NMF_K4 <- readRDS(nmf_rds)

cat("[READ] Existing NMF_K4 object:\n")
cat(nmf_rds, "\n")

# Check and standardize names
w <- as.matrix(NMF::basis(NMF_K4))
h0 <- as.matrix(NMF::coef(NMF_K4))

if (ncol(w) == 4) {
  colnames(NMF::basis(NMF_K4)) <- sprintf("CM%02d", 1:4)
}
if (nrow(h0) == 4) {
  rownames(NMF::coef(NMF_K4)) <- sprintf("CM%02d", 1:4)
}

cat("basis dim:\n")
print(dim(NMF::basis(NMF_K4)))

cat("coef dim:\n")
print(dim(NMF::coef(NMF_K4)))

# ------------------------------------------------------------
# 2. Extract sample level CM abundance matrix
# h = CM x sample
# ------------------------------------------------------------

h <- tryCatch(
  NMF::scoef(NMF_K4),
  error = function(e) NMF::coef(NMF_K4)
)

h <- as.matrix(h)

if (nrow(h) == 4) {
  rownames(h) <- sprintf("CM%02d", 1:4)
} else if (ncol(h) == 4) {
  h <- t(h)
  rownames(h) <- sprintf("CM%02d", 1:4)
} else {
  stop("h matrix does not look like CM x sample. Current dim: ", paste(dim(h), collapse = " x "))
}

cat("CM abundance matrix h:\n")
print(dim(h))

cat("Column sum summary:\n")
print(summary(colSums(h)))

# ------------------------------------------------------------
# 3. Build sample level metadata and assign dominant CM
# ------------------------------------------------------------

stopifnot("sampleID" %in% colnames(meta))

meta_sp_sample <- meta %>%
  dplyr::filter(sampleID %in% colnames(h)) %>%
  dplyr::distinct(sampleID, .keep_all = TRUE) %>%
  as.data.frame()

rownames(meta_sp_sample) <- meta_sp_sample$sampleID

common_samples <- intersect(colnames(h), rownames(meta_sp_sample))
h <- h[, common_samples, drop = FALSE]
meta_sp_sample <- meta_sp_sample[common_samples, , drop = FALSE]

stopifnot(all(colnames(h) == rownames(meta_sp_sample)))

max_cm <- apply(h, 2, function(x) rownames(h)[which.max(x)])
meta_sp_sample$CM_dominant <- unname(max_cm[colnames(h)])
meta_sp_sample$CMT <- gsub("^CM", "CMT", meta_sp_sample$CM_dominant)

saveRDS(
  h,
  file.path(outdir, "NMF_K4_supplement", "NMF_h_existing_K4.rds")
)

saveRDS(
  meta_sp_sample,
  file.path(outdir, "NMF_K4_supplement", "meta_sample_CMT_existing_K4.rds")
)

write.csv(
  as.data.frame(h),
  file.path(outdir, "NMF_K4_supplement", "TableS7_sample_level_CM_abundance_existing_K4.csv")
)

write.csv(
  meta_sp_sample,
  file.path(outdir, "NMF_K4_supplement", "TableS7_sample_metadata_with_dominant_CM_existing_K4.csv"),
  row.names = FALSE
)

cat("[DONE] Saved sample level CM abundance and metadata\n")

# ------------------------------------------------------------
# 4. K = 4 consensus matrix from existing NMF object
# ------------------------------------------------------------

pdf(
  file.path(outdir, "NMF_K4_supplement", "FigS2a_consensus_existing_K04.pdf"),
  width = 8,
  height = 7
)

NMF::consensusmap(
  NMF_K4,
  labCol = NA,
  labRow = NA
)

dev.off()

png(
  file.path(outdir, "NMF_K4_supplement", "FigS2a_consensus_existing_K04.png"),
  width = 1600,
  height = 1400,
  res = 200
)

NMF::consensusmap(
  NMF_K4,
  labCol = NA,
  labRow = NA
)

dev.off()

cat("[DONE] Saved K4 consensus matrix from existing NMF object\n")

# ------------------------------------------------------------
# 5. Full Level 2 module weight table
# ------------------------------------------------------------

w <- as.matrix(NMF::basis(NMF_K4))

if (nrow(w) == 4 && ncol(w) > 4) {
  w <- t(w)
}

if (ncol(w) != 4) {
  stop("basis matrix does not look like Level2 x CM. Current dim: ", paste(dim(w), collapse = " x "))
}

colnames(w) <- sprintf("CM%02d", 1:4)

if (is.null(rownames(w)) || any(is.na(rownames(w))) || any(rownames(w) == "")) {
  stop("basis matrix has no Level 2 population rownames.")
}

stopifnot("subCluster" %in% colnames(meta))
stopifnot("majorCluster" %in% colnames(meta))

sub_to_major <- meta %>%
  dplyr::filter(
    !is.na(subCluster), subCluster != "",
    !is.na(majorCluster), majorCluster != ""
  ) %>%
  dplyr::count(subCluster, majorCluster, name = "n") %>%
  dplyr::group_by(subCluster) %>%
  dplyr::slice_max(order_by = n, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::select(subCluster, majorCluster)

module_weight_all <- as.data.frame(w) %>%
  tibble::rownames_to_column("subCluster") %>%
  tidyr::pivot_longer(
    cols = starts_with("CM"),
    names_to = "CM",
    values_to = "weight"
  ) %>%
  dplyr::left_join(sub_to_major, by = "subCluster") %>%
  dplyr::mutate(
    majorCluster = ifelse(is.na(majorCluster), "Other", majorCluster),
    CM = factor(CM, levels = sprintf("CM%02d", 1:4))
  ) %>%
  dplyr::group_by(CM) %>%
  dplyr::arrange(dplyr::desc(weight), .by_group = TRUE) %>%
  dplyr::mutate(rank_within_CM = dplyr::row_number()) %>%
  dplyr::ungroup()

write.csv(
  module_weight_all,
  file.path(outdir, "NMF_K4_supplement", "TableS7_all_Level2_module_weights_existing_K4.csv"),
  row.names = FALSE
)

cat("[DONE] Saved full Level 2 module weight table\n")

# ------------------------------------------------------------
# 6. Built in module weight plots from existing NMF object
# ------------------------------------------------------------

pdf(
  file.path(outdir, "NMF_K4_supplement", "FigS2c_all_Level2_module_weights_existing_K04.pdf"),
  width = 5.5,
  height = 6.5
)

CoVarNet::gr.weight_all(NMF_K4)

dev.off()

pdf(
  file.path(outdir, "NMF_K4_supplement", "FigS2c_top_Level2_module_weights_existing_K04.pdf"),
  width = 11,
  height = 10
)

CoVarNet::gr.weight_top(NMF_K4, num = 15)

dev.off()

cat("[DONE] Saved module weight plots from existing K4 object\n")

# ------------------------------------------------------------
# 7. Read existing raw matrix for network analysis
# Do not recompute NMF
# ------------------------------------------------------------

mat_raw_rds <- file.path(outdir, "CM_input_mat_fq_raw_K4.rds")

if (!file.exists(mat_raw_rds)) {
  stop(
    "Cannot find raw Level 2 frequency matrix:\n",
    mat_raw_rds,
    "\nThis file should be generated by your current Fig.2d/e/f script.\n",
    "Please run fig2de.R first, or provide the correct mat_fq_raw RDS path."
  )
}

mat_fq_raw <- readRDS(mat_raw_rds)

cat("[READ] Existing mat_fq_raw:\n")
cat(mat_raw_rds, "\n")
print(dim(mat_fq_raw))

# ------------------------------------------------------------
# 8. Pairwise Level 2 population co abundance correlation
# ------------------------------------------------------------

cor_pair <- CoVarNet::pair_correlation(
  mat_fq_raw,
  method = "pearson"
)

saveRDS(
  cor_pair,
  file.path(outdir, "K4_network", "cor_pair_Level2_pearson_existing_K04.rds")
)

cat("[DONE] Calculated Level 2 pairwise correlation using existing mat_fq_raw\n")
