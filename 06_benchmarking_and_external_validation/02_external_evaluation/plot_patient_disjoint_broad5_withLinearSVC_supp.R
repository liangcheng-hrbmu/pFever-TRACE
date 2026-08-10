suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(jsonlite)
})

# =========================================================
# Patient-disjoint broad5 supplementary figure, with LinearSVC
#
# a. Level 1 BAS dot plot
# b. Level 1 MCC dot plot
# c. Disease by Level 1 precision heatmap
# d. Disease by Level 1 recall heatmap
#
# Run plot_patient_disjoint_broad5_withLinearSVC_main.R first.
# =========================================================

# -------------------------------
# 0. Paths
# -------------------------------

base_dir <- "data/analysis/09_bulk_deconvolution"
ext_root <- file.path(base_dir, "12_external_validation_multi_cohort_assets_v1")
cohort_id <- "Patient-disjoint set"

eval_dir <- file.path(
  ext_root,
  "results",
  "patient_disjoint_sensitivity",
  "evaluation_withLinearSVC"
)

out_dir <- file.path(
  ext_root,
  "plots",
  "patient_disjoint_broad5_withLinearSVC_supp_full"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

maj_fp <- file.path(eval_dir, "broad5_merged_predictions_with_truth.tsv")
level1_fp <- file.path(eval_dir, "broad5_level1_merged_predictions_with_truth.tsv")
metrics_fp <- file.path(eval_dir, "broad5_external_metrics.tsv")
level1_metrics_fp <- file.path(eval_dir, "broad5_level1_external_metrics.tsv")

need_file <- function(fp) {
  if (!file.exists(fp)) {
    stop(
      "Missing file: ", fp,
      "\nPlease run plot_patient_disjoint_broad5_withLinearSVC_main.R first."
    )
  }
  invisible(fp)
}

need_file(maj_fp)
need_file(level1_fp)
need_file(metrics_fp)
need_file(level1_metrics_fp)

# -------------------------------
# 1. Device and theme
# -------------------------------

base_family <- "sans"

save_gg <- function(plot, stub, width, height) {
  ggplot2::ggsave(
    filename = file.path(out_dir, paste0(stub, ".pdf")),
    plot = plot,
    width = width,
    height = height,
    device = grDevices::cairo_pdf,
    family = base_family
  )

  ggplot2::ggsave(
    filename = file.path(out_dir, paste0(stub, ".png")),
    plot = plot,
    width = width,
    height = height,
    dpi = 300
  )
}

open_pdf <- function(filename, width, height) {
  grDevices::cairo_pdf(
    filename = filename,
    width = width,
    height = height,
    family = base_family
  )
}

paper_theme <- function(base_size = 10) {
  ggplot2::theme_classic(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0, size = base_size + 4),
      plot.subtitle = ggplot2::element_text(size = base_size, color = "black"),
      axis.title = ggplot2::element_text(size = base_size + 1, color = "black"),
      axis.text = ggplot2::element_text(size = base_size, color = "black"),
      axis.line = ggplot2::element_line(linewidth = 0.35, color = "black"),
      axis.ticks = ggplot2::element_line(linewidth = 0.30, color = "black"),
      legend.position = "none",
      plot.margin = ggplot2::margin(4, 4, 4, 4, "pt")
    )
}

# -------------------------------
# 2. Colors and orders
# -------------------------------

celltype_order <- c(
  "Majority vote", "Random", "B", "Plasma", "pDC", "DC", "Mono",
  "T CD4 naive", "T CD4 non-naive",
  "T CD8 naive", "T CD8 non-naive",
  "UTC", "ILC"
)

celltype_order_heatmap <- c(
  "B", "Plasma", "pDC", "DC", "Mono",
  "T CD4 naive", "T CD4 non-naive",
  "T CD8 naive", "T CD8 non-naive",
  "UTC", "ILC"
)

celltype_cols <- c(
  "Majority vote"     = "#535353",
  "Random"            = "#D9D9D9",
  "B"                 = "#7BC6D6",
  "Plasma"            = "#025566",
  "pDC"               = "#A7C957",
  "DC"                = "#6A994E",
  "Mono"              = "#386641",
  "T CD4 naive"       = "#F2E6A2",
  "T CD4 non-naive"   = "#D99B38",
  "T CD8 naive"       = "#B63232",
  "T CD8 non-naive"   = "#6E0000",
  "UTC"               = "#8D6A8C",
  "ILC"               = "#6F2141"
)

celltype_labels <- c(
  "Majority vote"     = "Majority vote",
  "Random"            = "Random",
  "B"                 = "B",
  "Plasma"            = "Plasma",
  "pDC"               = "pDC",
  "DC"                = "DC",
  "Mono"              = "Mono",
  "T CD4 naive"       = "T CD4 naive",
  "T CD4 non-naive"   = "T CD4 non naive",
  "T CD8 naive"       = "T CD8 naive",
  "T CD8 non-naive"   = "T CD8 non naive",
  "UTC"               = "UTC",
  "ILC"               = "ILC"
)

broad_order <- c(
  "Viral infection",
  "Bacterial infection",
  "Malaria",
  "Hyperinflammatory disease",
  "Autoimmune disease",
  "Healthy",
  "Not assigned"
)

broad_group_cols <- c(
  "Viral infection"             = "#6699CC",
  "Bacterial infection"         = "#C48EAD",
  "Malaria"                     = "#6A4C93",
  "Hyperinflammatory disease"   = "#E3A722",
  "Autoimmune disease"          = "#36A596",
  "Healthy"                     = "#8E8E8E",
  "Not assigned"                = "#2F225B",
  "Other"                       = "#BDBDBD"
)

score_breaks <- c(0, 0.20, 0.40, 0.60, 0.80, 1.00)
score_cols <- c("#000004", "#3B0F70", "#8C2981", "#DE4968", "#FE9F6D", "#FCFDBF")
score_col_fun <- circlize::colorRamp2(score_breaks, score_cols)

score_heatmap_legend_param <- function(title) {
  list(
    title = title,
    at = c(0, 0.5, 1),
    labels = c("0", "0.5", "1"),
    direction = "horizontal",
    title_position = "topcenter",
    legend_width = grid::unit(30, "mm"),
    grid_height = grid::unit(2.8, "mm"),
    border = NA,
    labels_gp = grid::gpar(fontsize = 8, fontfamily = base_family),
    title_gp = grid::gpar(fontsize = 9, fontfamily = base_family)
  )
}

complete_named_cols <- function(levels_vec, base_cols, default = "#BDBDBD") {
  missing <- setdiff(levels_vec, names(base_cols))
  out <- c(base_cols, stats::setNames(rep(default, length(missing)), missing))
  out[levels_vec]
}

label_to_group <- function(x) {
  x <- as.character(x)
  out <- ifelse(x %in% names(broad_group_cols), x, "Other")
  factor(out, levels = c(broad_order, "Other"))
}

# -------------------------------
# 3. Metric functions
# -------------------------------

calc_multiclass_metrics <- function(y_true, y_pred, labels = NULL) {
  y_true <- as.character(y_true)
  y_pred <- as.character(y_pred)

  valid <- !is.na(y_true) & !is.na(y_pred)
  y_true <- y_true[valid]
  y_pred <- y_pred[valid]

  if (is.null(labels)) labels <- sort(unique(y_true))
  labels <- labels[labels %in% unique(y_true)]

  per_class <- dplyr::bind_rows(lapply(labels, function(cls) {
    tp <- sum(y_true == cls & y_pred == cls)
    fp <- sum(y_true != cls & y_pred == cls)
    fn <- sum(y_true == cls & y_pred != cls)
    support <- sum(y_true == cls)

    precision <- ifelse(tp + fp == 0, 0, tp / (tp + fp))
    recall <- ifelse(tp + fn == 0, 0, tp / (tp + fn))
    f1 <- ifelse(precision + recall == 0, 0, 2 * precision * recall / (precision + recall))

    data.frame(
      class = cls,
      precision = precision,
      recall = recall,
      F1 = f1,
      support = support,
      stringsAsFactors = FALSE
    )
  }))

  weighted_F1 <- sum(per_class$F1 * per_class$support) / sum(per_class$support)
  macro_F1 <- mean(per_class$F1)
  balanced_accuracy <- mean(per_class$recall)
  accuracy <- mean(y_true == y_pred)

  all_labs <- sort(unique(c(y_true, y_pred)))
  cm <- table(
    factor(y_true, levels = all_labs),
    factor(y_pred, levels = all_labs)
  )

  t_sum <- rowSums(cm)
  p_sum <- colSums(cm)
  s <- sum(cm)
  c_val <- sum(diag(cm))

  numerator <- c_val * s - sum(t_sum * p_sum)
  denominator <- sqrt((s^2 - sum(p_sum^2)) * (s^2 - sum(t_sum^2)))
  mcc <- ifelse(denominator == 0, NA_real_, numerator / denominator)

  list(
    metrics = data.frame(
      accuracy = accuracy,
      weighted_F1 = weighted_F1,
      macro_F1 = macro_F1,
      balanced_accuracy = balanced_accuracy,
      MCC = mcc,
      n_samples = length(y_true),
      n_truth_classes = length(unique(y_true)),
      n_pred_classes = length(unique(y_pred)),
      stringsAsFactors = FALSE
    ),
    per_class = per_class
  )
}

calc_binary_class_metrics <- function(y_true, y_pred, positive_class) {
  y_true <- as.character(y_true)
  y_pred <- as.character(y_pred)

  valid <- !is.na(y_true) & !is.na(y_pred)
  y_true <- y_true[valid]
  y_pred <- y_pred[valid]

  yt <- y_true == positive_class
  yp <- y_pred == positive_class

  tp <- sum(yt & yp)
  fp <- sum(!yt & yp)
  fn <- sum(yt & !yp)
  tn <- sum(!yt & !yp)

  precision <- ifelse(tp + fp == 0, 0, tp / (tp + fp))
  recall <- ifelse(tp + fn == 0, 0, tp / (tp + fn))
  specificity <- ifelse(tn + fp == 0, 0, tn / (tn + fp))
  F1 <- ifelse(precision + recall == 0, 0, 2 * precision * recall / (precision + recall))
  BAS <- (recall + specificity) / 2

  denom <- sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
  MCC <- ifelse(denom == 0, NA_real_, ((tp * tn) - (fp * fn)) / denom)

  data.frame(
    class = positive_class,
    precision = precision,
    recall = recall,
    F1 = F1,
    specificity = specificity,
    balanced_accuracy_binary = BAS,
    MCC_binary = MCC,
    support = sum(yt),
    tp = tp,
    fp = fp,
    fn = fn,
    tn = tn,
    stringsAsFactors = FALSE
  )
}

# -------------------------------
# 4. Read data
# -------------------------------

maj_df <- data.table::fread(maj_fp, sep = "\t", data.table = FALSE)
pred_df <- data.table::fread(level1_fp, sep = "\t", data.table = FALSE)
metrics_df <- data.table::fread(metrics_fp, sep = "\t", data.table = FALSE)
level1_metrics_df <- data.table::fread(level1_metrics_fp, sep = "\t", data.table = FALSE)

required_maj <- c("sample_id", "true_short", "pred_short")
required_l1 <- c("sample_id", "Level1_plot", "true_short", "pred_short")
required_metrics <- c("weighted_F1", "balanced_accuracy", "MCC")
required_l1_metrics <- c("Level1", "weighted_F1", "balanced_accuracy", "MCC")

miss <- setdiff(required_maj, colnames(maj_df))
if (length(miss) > 0) stop("maj_df missing columns: ", paste(miss, collapse = ", "))

miss <- setdiff(required_l1, colnames(pred_df))
if (length(miss) > 0) stop("pred_df missing columns: ", paste(miss, collapse = ", "))

miss <- setdiff(required_metrics, colnames(metrics_df))
if (length(miss) > 0) stop("metrics_df missing columns: ", paste(miss, collapse = ", "))

miss <- setdiff(required_l1_metrics, colnames(level1_metrics_df))
if (length(miss) > 0) stop("level1_metrics_df missing columns: ", paste(miss, collapse = ", "))

disease_order <- broad_order[broad_order %in% unique(maj_df$true_short)]
disease_order <- setdiff(disease_order, "Not assigned")

celltype_order_existing <- celltype_order[
  celltype_order %in% c("Majority vote", "Random", unique(level1_metrics_df$Level1), unique(pred_df$Level1_plot))
]

celltype_heatmap_existing <- celltype_order_heatmap[
  celltype_order_heatmap %in% unique(pred_df$Level1_plot)
]

# -------------------------------
# 5. Random baseline for panel a and b
# -------------------------------

set.seed(1234)
n_random <- 100

random_metrics <- dplyr::bind_rows(lapply(seq_len(n_random), function(i) {
  y_random <- sample(maj_df$true_short, length(maj_df$true_short), replace = FALSE)
  calc_multiclass_metrics(
    y_true = maj_df$true_short,
    y_pred = y_random,
    labels = disease_order
  )$metrics %>%
    dplyr::mutate(iteration = i, .before = 1)
}))

supp_metric_source <- dplyr::bind_rows(
  data.frame(
    Celltype = "Majority vote",
    Metric = "Balanced accuracy score",
    Score = metrics_df$balanced_accuracy[1],
    ymin = NA_real_,
    ymax = NA_real_
  ),
  data.frame(
    Celltype = "Majority vote",
    Metric = "Matthew correlation coefficient",
    Score = metrics_df$MCC[1],
    ymin = NA_real_,
    ymax = NA_real_
  ),
  data.frame(
    Celltype = "Random",
    Metric = "Balanced accuracy score",
    Score = mean(random_metrics$balanced_accuracy, na.rm = TRUE),
    ymin = mean(random_metrics$balanced_accuracy, na.rm = TRUE) - stats::sd(random_metrics$balanced_accuracy, na.rm = TRUE),
    ymax = mean(random_metrics$balanced_accuracy, na.rm = TRUE) + stats::sd(random_metrics$balanced_accuracy, na.rm = TRUE)
  ),
  data.frame(
    Celltype = "Random",
    Metric = "Matthew correlation coefficient",
    Score = mean(random_metrics$MCC, na.rm = TRUE),
    ymin = mean(random_metrics$MCC, na.rm = TRUE) - stats::sd(random_metrics$MCC, na.rm = TRUE),
    ymax = mean(random_metrics$MCC, na.rm = TRUE) + stats::sd(random_metrics$MCC, na.rm = TRUE)
  ),
  level1_metrics_df %>%
    dplyr::transmute(
      Celltype = as.character(.data$Level1),
      Metric = "Balanced accuracy score",
      Score = as.numeric(.data$balanced_accuracy),
      ymin = NA_real_,
      ymax = NA_real_
    ),
  level1_metrics_df %>%
    dplyr::transmute(
      Celltype = as.character(.data$Level1),
      Metric = "Matthew correlation coefficient",
      Score = as.numeric(.data$MCC),
      ymin = NA_real_,
      ymax = NA_real_
    )
) %>%
  dplyr::mutate(
    Celltype = factor(Celltype, levels = celltype_order_existing)
  )

data.table::fwrite(
  supp_metric_source,
  file.path(out_dir, "SourceData_PatientDisjointSuppFig_ab_BAS_MCC_random_level1.tsv"),
  sep = "\t"
)

plot_one_metric <- function(metric_name, panel_label, y_label, outfile_prefix) {
  plot_df <- supp_metric_source %>%
    dplyr::filter(Metric == metric_name) %>%
    dplyr::mutate(
      Celltype = factor(as.character(Celltype), levels = celltype_order_existing)
    )

  y_min <- min(plot_df$Score, plot_df$ymin, na.rm = TRUE)
  y_max <- max(plot_df$Score, plot_df$ymax, na.rm = TRUE)

  if (grepl("MCC|Matthew", metric_name)) {
    y_lower <- min(-0.1, floor(y_min * 10) / 10)
    y_upper <- min(1, ceiling(y_max * 10) / 10 + 0.05)
    y_breaks <- seq(-0.2, 1, 0.2)
  } else {
    y_lower <- max(0, floor(y_min * 10) / 10 - 0.02)
    y_upper <- min(1, ceiling(y_max * 10) / 10 + 0.05)
    y_breaks <- seq(0, 1, 0.1)
  }

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = Celltype, y = Score, color = Celltype)) +
    ggplot2::geom_errorbar(
      data = dplyr::filter(plot_df, as.character(Celltype) == "Random"),
      ggplot2::aes(ymin = ymin, ymax = ymax),
      width = 0.18,
      linewidth = 0.45,
      color = "#BDBDBD"
    ) +
    ggplot2::geom_vline(xintercept = 2.5, color = "black", linewidth = 0.35) +
    ggplot2::geom_point(size = 2.8) +
    ggplot2::scale_color_manual(values = complete_named_cols(levels(plot_df$Celltype), celltype_cols), drop = FALSE) +
    ggplot2::scale_x_discrete(labels = celltype_labels) +
    ggplot2::coord_cartesian(ylim = c(y_lower, y_upper), clip = "off") +
    ggplot2::scale_y_continuous(
      breaks = y_breaks,
      expand = ggplot2::expansion(mult = c(0.02, 0.04))
    ) +
    paper_theme(base_size = 10) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)
    ) +
    ggplot2::labs(
      title = panel_label,
      subtitle = paste0(cohort_id, ", n = ", nrow(maj_df)),
      x = NULL,
      y = y_label
    )

  save_gg(p, outfile_prefix, width = 4.6, height = 3.3)
  p
}

p_a <- plot_one_metric(
  metric_name = "Balanced accuracy score",
  panel_label = "a",
  y_label = "Balanced accuracy score",
  outfile_prefix = "PatientDisjointSuppFig_a_BAS_random_level1_withLinearSVC"
)

p_b <- plot_one_metric(
  metric_name = "Matthew correlation coefficient",
  panel_label = "b",
  y_label = "Matthew correlation coefficient",
  outfile_prefix = "PatientDisjointSuppFig_b_MCC_random_level1_withLinearSVC"
)

# -------------------------------
# 6. Disease by Level 1 class metrics
# -------------------------------

heat_rows <- list()

for (ct in celltype_heatmap_existing) {
  sub_ct <- pred_df %>%
    dplyr::filter(as.character(.data$Level1_plot) == ct)

  for (d in disease_order) {
    heat_rows[[length(heat_rows) + 1]] <- calc_binary_class_metrics(
      y_true = sub_ct$true_short,
      y_pred = sub_ct$pred_short,
      positive_class = d
    ) %>%
      dplyr::mutate(
        Disease = d,
        Celltype = ct,
        .before = 1
      )
  }
}

heat_long <- dplyr::bind_rows(heat_rows) %>%
  dplyr::mutate(
    Disease = factor(Disease, levels = disease_order),
    Celltype = factor(Celltype, levels = celltype_heatmap_existing)
  )

data.table::fwrite(
  heat_long,
  file.path(out_dir, "SourceData_PatientDisjointSuppFig_cd_disease_celltype_precision_recall.tsv"),
  sep = "\t"
)

plot_heatmap_metric <- function(metric_col, legend_title, panel_label, file_stub) {
  mat_df <- heat_long %>%
    dplyr::select(Disease, Celltype, value = .data[[metric_col]]) %>%
    tidyr::pivot_wider(names_from = Celltype, values_from = value) %>%
    as.data.frame()

  rownames(mat_df) <- as.character(mat_df$Disease)
  mat_df$Disease <- NULL

  mat <- as.matrix(mat_df)
  mat[is.na(mat)] <- 0

  row_group <- droplevels(label_to_group(rownames(mat)))

  ha_bottom <- ComplexHeatmap::HeatmapAnnotation(
    `Cell types` = factor(colnames(mat), levels = celltype_heatmap_existing),
    col = list(`Cell types` = complete_named_cols(colnames(mat), celltype_cols)),
    show_annotation_name = FALSE,
    show_legend = FALSE,
    which = "column",
    simple_anno_size = grid::unit(3.2, "mm"),
    border = FALSE
  )

  ra_left <- ComplexHeatmap::rowAnnotation(
    `Disease group` = row_group,
    col = list(`Disease group` = broad_group_cols[levels(row_group)]),
    show_annotation_name = FALSE,
    show_legend = FALSE,
    simple_anno_size = grid::unit(3.2, "mm"),
    border = FALSE
  )

  ht <- ComplexHeatmap::Heatmap(
    mat,
    name = legend_title,
    col = score_col_fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    left_annotation = ra_left,
    bottom_annotation = ha_bottom,
    row_split = row_group,
    row_gap = grid::unit(1.1, "mm"),
    column_labels = unname(celltype_labels[colnames(mat)]),
    column_names_rot = 90,
    column_names_gp = grid::gpar(fontsize = 8, fontfamily = base_family),
    row_names_gp = grid::gpar(fontsize = 8.5, fontfamily = base_family),
    rect_gp = grid::gpar(col = NA),
    na_col = "#000004",
    heatmap_legend_param = score_heatmap_legend_param(legend_title),
    column_title = NULL,
    row_title = NULL,
    width = grid::unit(ncol(mat) * 4.8, "mm"),
    height = grid::unit(nrow(mat) * 6.0, "mm")
  )

  w <- max(5.0, ncol(mat) * 0.25 + 2.4)
  h <- max(3.2, nrow(mat) * 0.40 + 1.8)

  pdf_fp <- file.path(out_dir, paste0(file_stub, ".pdf"))
  png_fp <- file.path(out_dir, paste0(file_stub, ".png"))

  open_pdf(pdf_fp, width = w, height = h)
  grid::grid.text(
    panel_label,
    x = grid::unit(0.02, "npc"),
    y = grid::unit(0.98, "npc"),
    just = c("left", "top"),
    gp = grid::gpar(fontsize = 14, fontface = "bold", fontfamily = base_family)
  )
  ComplexHeatmap::draw(
    ht,
    heatmap_legend_side = "top",
    annotation_legend_side = "bottom",
    newpage = FALSE
  )
  dev.off()

  png(png_fp, width = w, height = h, units = "in", res = 300)
  grid::grid.text(
    panel_label,
    x = grid::unit(0.02, "npc"),
    y = grid::unit(0.98, "npc"),
    just = c("left", "top"),
    gp = grid::gpar(fontsize = 14, fontface = "bold", fontfamily = base_family)
  )
  ComplexHeatmap::draw(
    ht,
    heatmap_legend_side = "top",
    annotation_legend_side = "bottom",
    newpage = FALSE
  )
  dev.off()

  invisible(ht)
}

ht_c <- plot_heatmap_metric(
  metric_col = "precision",
  legend_title = "Precision",
  panel_label = "c",
  file_stub = "PatientDisjointSuppFig_c_disease_celltype_precision_heatmap_withLinearSVC"
)

ht_d <- plot_heatmap_metric(
  metric_col = "recall",
  legend_title = "Recall",
  panel_label = "d",
  file_stub = "PatientDisjointSuppFig_d_disease_celltype_recall_heatmap_withLinearSVC"
)

# -------------------------------
# 7. Optional combined a and b
# -------------------------------

if (requireNamespace("patchwork", quietly = TRUE)) {
  combined_ab <- p_a | p_b

  ggplot2::ggsave(
    filename = file.path(out_dir, "PatientDisjointSuppFig_ab_BAS_MCC_combined_withLinearSVC.pdf"),
    plot = combined_ab,
    width = 9.2,
    height = 3.4,
    device = grDevices::cairo_pdf,
    family = base_family
  )

  ggplot2::ggsave(
    filename = file.path(out_dir, "PatientDisjointSuppFig_ab_BAS_MCC_combined_withLinearSVC.png"),
    plot = combined_ab,
    width = 9.2,
    height = 3.4,
    dpi = 300
  )
}

summary_list <- list(
  cohort_id = cohort_id,
  eval_dir = eval_dir,
  out_dir = out_dir,
  n_majority_samples = nrow(maj_df),
  n_level1_rows = nrow(pred_df),
  disease_order = disease_order,
  celltype_heatmap_existing = celltype_heatmap_existing,
  panels = c(
    "a: BAS dot plot",
    "b: MCC dot plot",
    "c: disease by Level 1 precision heatmap",
    "d: disease by Level 1 recall heatmap"
  )
)

jsonlite::write_json(
  summary_list,
  file.path(out_dir, "patient_disjoint_supp_withLinearSVC_summary.json"),
  auto_unbox = TRUE,
  pretty = TRUE
)

cat("\n[DONE] Patient-disjoint supplementary figure saved to:\n")
cat(out_dir, "\n")
cat("\nGenerated files:\n")
cat("  - PatientDisjointSuppFig_a_BAS_random_level1_withLinearSVC.pdf/png\n")
cat("  - PatientDisjointSuppFig_b_MCC_random_level1_withLinearSVC.pdf/png\n")
cat("  - PatientDisjointSuppFig_c_disease_celltype_precision_heatmap_withLinearSVC.pdf/png\n")
cat("  - PatientDisjointSuppFig_d_disease_celltype_recall_heatmap_withLinearSVC.pdf/png\n")
cat("  - SourceData_PatientDisjointSuppFig_cd_disease_celltype_precision_recall.tsv\n")