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
# Bulk classifier main figure c to e, fine 13 disease task
#
# c: WF1 score distribution for majority vote and each Level 1 classifier
# d: F1 score for each fine disease and cell type combination after pooling all folds
# e: Column normalized confusion matrix after pooling all folds
#
# This script keeps the analysis logic from the original bulk scripts, but
# restructures the plotting into independent, clearly annotated blocks.
# =========================================================

# -------------------------------
# 0. Input and output
# -------------------------------
fine_dir <- "data/analysis/09_bulk_deconvolution/08_bulk_classifier/03_patient_classifier_fine_bulkZ_L1x11_final_withLinearSVC"

out_dir <- "data/analysis/09_bulk_deconvolution/08_bulk_classifier/05_final_bulk_plots_R1/Fig_fine13_cde_paper_mimic_separate"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

split_summary_tsv  <- file.path(fine_dir, "scenario1_split_summary.tsv")
overall_json       <- file.path(fine_dir, "scenario1_overall_metrics.json")
level1_metrics_tsv <- file.path(fine_dir, "scenario1_level1_metrics_all.tsv")
majority_pred_tsv  <- file.path(fine_dir, "scenario1_majority_predictions_all.tsv")
level1_pred_tsv    <- file.path(fine_dir, "query_level1_predictions_all.tsv")

# -------------------------------
# 1. Safe PDF device
# -------------------------------
# Use cairo_pdf and sans to avoid font family errors on Linux servers.
base_family <- "sans"

save_gg_pdf <- function(filename, plot, width, height) {
  ggplot2::ggsave(
    filename = filename,
    plot = plot,
    width = width,
    height = height,
    device = grDevices::cairo_pdf,
    family = base_family
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

# -------------------------------
# 2. Shared palettes
# -------------------------------
# Continuous score palette: inferno like, matching the current paper mimic main figure style.
score_breaks <- c(0, 0.20, 0.40, 0.60, 0.80, 1.00)
score_cols   <- c("#000004", "#3B0F70", "#8C2981", "#DE4968", "#FE9F6D", "#FCFDBF")
score_ticks  <- c(0, 0.5, 1)
score_col_fun <- circlize::colorRamp2(score_breaks, score_cols)

score_heatmap_legend_param <- function(title) {
  list(
    title = title,
    at = score_ticks,
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

celltype_cols <- c(
  "Majority vote"     = "#535353",
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

celltype_order <- c(
  "Majority vote", "B", "Plasma", "pDC", "DC", "Mono",
  "T CD4 naive", "T CD4 non-naive",
  "T CD8 naive", "T CD8 non-naive",
  "UTC", "ILC"
)

celltype_labels <- c(
  "Majority vote"     = "Majority\nvote",
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

bulk_group_cols <- c(
  "Viral infection"             = "#6699CC",
  "Bacterial infection"         = "#C48EAD",
  "Malaria"                     = "#6A4C93",
  "Hyperinflammatory disease"   = "#E3A722",
  "Autoimmune disease"          = "#36A596",
  "Other"                       = "#BDBDBD",
  "Not assigned"                = "#2F225B"
)





bulk_group_order <- c(
  "Viral infection",
  "Bacterial infection",
  "Malaria",
  "Hyperinflammatory disease",
  "Autoimmune disease",
  "Other",
  "Not assigned"
)

fine_order <- c(
  "IAV", "RSV", "Adenovirus", "Rhinovirus",
  "Meningococcus", "Pneumococcus", "S. aureus", "GAS", "E. coli", "Tuberculosis",
  "Malaria", "KD", "JIA"
)

fine_group_map <- c(
  "IAV"           = "Viral infection",
  "RSV"           = "Viral infection",
  "Adenovirus"    = "Viral infection",
  "Rhinovirus"    = "Viral infection",
  "Meningococcus" = "Bacterial infection",
  "Pneumococcus"  = "Bacterial infection",
  "S. aureus"     = "Bacterial infection",
  "GAS"           = "Bacterial infection",
  "E. coli"       = "Bacterial infection",
  "Tuberculosis"  = "Bacterial infection",
  "Malaria"       = "Malaria",
  "KD"            = "Hyperinflammatory disease",
  "JIA"           = "Autoimmune disease",
  "Not assigned"  = "Not assigned"
)

paper_theme <- function(base_size = 10) {
  ggplot2::theme_classic(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0, size = base_size + 4),
      axis.title = ggplot2::element_text(size = base_size + 1, color = "black"),
      axis.text = ggplot2::element_text(size = base_size, color = "black"),
      axis.line = ggplot2::element_line(linewidth = 0.35, color = "black"),
      axis.ticks = ggplot2::element_line(linewidth = 0.30, color = "black"),
      legend.position = "none",
      plot.margin = ggplot2::margin(4, 4, 4, 4, "pt")
    )
}

# -------------------------------
# 3. Helper functions
# -------------------------------
pick_col <- function(df, candidates, required = TRUE) {
  hit <- intersect(candidates, colnames(df))
  if (length(hit) == 0) {
    if (required) {
      stop("Cannot find any of columns: ", paste(candidates, collapse = ", "))
    }
    return(NA_character_)
  }
  hit[1]
}

filter_primary <- function(df) {
  if ("benchmark_status" %in% colnames(df)) {
    df <- dplyr::filter(df, as.character(.data$benchmark_status) == "primary_benchmark")
  }
  df
}

short_celltype <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  map <- c(
    "Monocyte" = "Mono",
    "Mono" = "Mono",
    "T CD4 non naive" = "T CD4 non-naive",
    "T CD8 non naive" = "T CD8 non-naive",
    "T CD4 non-naive" = "T CD4 non-naive",
    "T CD8 non-naive" = "T CD8 non-naive"
  )
  out <- ifelse(x %in% names(map), map[x], x)
  unname(out)
}

short_fine <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[is.na(x) | x == "" | x %in% c("NA", "NaN", "nan", "NULL", "None")] <- "Not assigned"

  map <- c(
    "Influenza A virus (IAV)" = "IAV",
    "IAV" = "IAV",
    "Respiratory Syncytial Virus Infection" = "RSV",
    "RSV" = "RSV",
    "adenovirus infection" = "Adenovirus",
    "Adenovirus" = "Adenovirus",
    "rhinovirus infection" = "Rhinovirus",
    "Rhinovirus" = "Rhinovirus",
    "meningococcal infection" = "Meningococcus",
    "Meningococcus" = "Meningococcus",
    "pneumococcal infection" = "Pneumococcus",
    "Pneumococcus" = "Pneumococcus",
    "Staphylococcus aureus infection" = "S. aureus",
    "S. aureus" = "S. aureus",
    "group A streptococcal infection" = "GAS",
    "GAS" = "GAS",
    "escherichia coli infection" = "E. coli",
    "E. coli" = "E. coli",
    "Tuberculosis" = "Tuberculosis",
    "malaria" = "Malaria",
    "Malaria" = "Malaria",
    "Kawasaki disease (KD)" = "KD",
    "KD" = "KD",
    "juvenile idiopathic arthritis (JIA)" = "JIA",
    "juvenile idopathic arthritis (JIA)" = "JIA",
    "JIA" = "JIA",
    "TIE" = "Not assigned",
    "Not assigned" = "Not assigned"
  )

  out <- ifelse(x %in% names(map), map[x], x)
  unname(out)
}

label_to_group <- function(x) {
  x <- as.character(x)
  out <- fine_group_map[x]
  out[is.na(out) & x == "Not assigned"] <- "Not assigned"
  out[is.na(out)] <- "Other"
  factor(out, levels = bulk_group_order)
}

complete_named_cols <- function(levels_vec, base_cols, default = "#BDBDBD") {
  missing <- setdiff(levels_vec, names(base_cols))
  out <- c(base_cols, stats::setNames(rep(default, length(missing)), missing))
  out[levels_vec]
}

compute_f1_bin <- function(y_true, y_pred, positive_class) {
  y_true <- as.character(y_true)
  y_pred <- as.character(y_pred)
  valid <- !is.na(y_true) & !is.na(y_pred)
  y_true <- y_true[valid]
  y_pred <- y_pred[valid]

  yt <- ifelse(y_true == positive_class, 1, 0)
  yp <- ifelse(y_pred == positive_class, 1, 0)

  tp <- sum(yt == 1 & yp == 1)
  fp <- sum(yt == 0 & yp == 1)
  fn <- sum(yt == 1 & yp == 0)

  if ((2 * tp + fp + fn) == 0) return(NA_real_)
  2 * tp / (2 * tp + fp + fn)
}

read_level1_predictions <- function(cls_dir) {
  merged_fp <- file.path(cls_dir, "query_level1_predictions_all.tsv")
  pred_dir <- file.path(cls_dir, "predictions")

  if (file.exists(merged_fp)) {
    message("[READ] ", merged_fp)
    return(data.table::fread(merged_fp, sep = "\t", data.table = FALSE))
  }

  fps <- list.files(
    pred_dir,
    pattern = "^query_level1_predictions_SPLIT_.*\\.tsv$",
    full.names = TRUE
  )

  if (length(fps) == 0) {
    stop("Cannot find query_level1_predictions_all.tsv or split prediction files in: ", cls_dir)
  }

  message("[READ] merging ", length(fps), " split level prediction files")
  df <- dplyr::bind_rows(lapply(fps, function(fp) {
    x <- data.table::fread(fp, sep = "\t", data.table = FALSE)
    x$.source_file <- basename(fp)
    x
  }))

  data.table::fwrite(df, merged_fp, sep = "\t")
  message("[WRITE] ", merged_fp)
  df
}

# -------------------------------
# 4. Read and standardize data
# -------------------------------
if (!file.exists(split_summary_tsv))  stop("Missing split summary file: ", split_summary_tsv)
if (!file.exists(level1_metrics_tsv)) stop("Missing Level 1 metrics file: ", level1_metrics_tsv)
if (!file.exists(majority_pred_tsv))  stop("Missing majority prediction file: ", majority_pred_tsv)

split_df  <- data.table::fread(split_summary_tsv, sep = "\t", data.table = FALSE)
level1_df <- data.table::fread(level1_metrics_tsv, sep = "\t", data.table = FALSE)
maj_df    <- data.table::fread(majority_pred_tsv, sep = "\t", data.table = FALSE)
pred_df   <- read_level1_predictions(fine_dir)

maj_df <- filter_primary(maj_df)
pred_df <- filter_primary(pred_df)

true_col <- pick_col(maj_df, c("true_disease", "disease_clean", "disease", "disease_Level_3", "true_label"))
pred_col <- pick_col(maj_df, c("pred_disease_majority", "pred_majority", "majority_pred", "pred_disease", "prediction", "pred_label"))

level1_col <- pick_col(level1_df, c("Level1", "classifier_ct", "celltype", "Celltype", "cell_type", "Cell_type"))
pred_level1_col <- pick_col(pred_df, c("Level1", "classifier_ct", "celltype", "Celltype", "cell_type", "Cell_type"))
pred_true_col <- pick_col(pred_df, c("true_disease", "disease_clean", "disease", "disease_Level_3", "true_label"))
pred_pred_col <- pick_col(pred_df, c("pred_disease", "pred_label", "prediction", "pred_disease_majority", "pred_majority", "majority_pred"))

level1_df$Level1_plot <- short_celltype(level1_df[[level1_col]])
pred_df$Level1_plot <- short_celltype(pred_df[[pred_level1_col]])

maj_df$true_short <- short_fine(maj_df[[true_col]])
maj_df$pred_short <- short_fine(maj_df[[pred_col]])
pred_df$true_short <- short_fine(pred_df[[pred_true_col]])
pred_df$pred_short <- short_fine(pred_df[[pred_pred_col]])

truth_labels <- setdiff(unique(maj_df$true_short), "Not assigned")
disease_order <- c(fine_order[fine_order %in% truth_labels], setdiff(sort(truth_labels), fine_order))

pred_labels <- unique(maj_df$pred_short)
pred_row_order <- unique(c(disease_order, setdiff(sort(pred_labels), c(disease_order, "Not assigned")), "Not assigned"))

celltype_order_existing <- celltype_order[celltype_order %in% c("Majority vote", unique(level1_df$Level1_plot), unique(pred_df$Level1_plot))]
celltype_cols_existing <- complete_named_cols(celltype_order_existing, celltype_cols)

# =========================================================
# 5. Panel c: WF1 boxplot
# =========================================================
# Each box contains the five split level WF1 scores. Whiskers use the default
# 1.5 x IQR rule; outliers are shown as individual points.
plot_panel_c <- function() {
  majority_wf1_col <- pick_col(split_df, c("WF1_majority", "weighted_F1_majority", "weighted_f1_majority", "WeightedF1_majority", "WF1"))
  level1_wf1_col <- pick_col(level1_df, c("WF1", "weighted_F1", "weighted_f1", "WeightedF1", "F1"))

  box_df <- dplyr::bind_rows(
    dplyr::transmute(split_df, split = as.character(.data$split), Celltype = "Majority vote", WF1 = as.numeric(.data[[majority_wf1_col]])),
    dplyr::transmute(level1_df, split = as.character(.data$split), Celltype = .data$Level1_plot, WF1 = as.numeric(.data[[level1_wf1_col]]))
  ) %>%
    dplyr::filter(!is.na(WF1)) %>%
    dplyr::mutate(Celltype = factor(Celltype, levels = celltype_order_existing))

  data.table::fwrite(box_df, file.path(out_dir, "SourceData_BulkFig_c_WF1_by_split.tsv"), sep = "\t")

  y_min <- max(0, floor(min(box_df$WF1, na.rm = TRUE) * 10) / 10)
  y_max <- min(1, ceiling(max(box_df$WF1, na.rm = TRUE) * 20) / 20)
  if (y_max < 0.98) y_max <- 1

  p_c <- ggplot2::ggplot(box_df, ggplot2::aes(x = Celltype, y = WF1, fill = Celltype)) +
    ggplot2::stat_boxplot(geom = "errorbar", width = 0.34, linewidth = 0.35, color = "black") +
    ggplot2::geom_boxplot(
      width = 0.65,
      linewidth = 0.35,
      color = "black",
      outlier.shape = 21,
      outlier.size = 1.15,
      outlier.stroke = 0.25,
      outlier.fill = "white",
      outlier.colour = "black"
    ) +
    ggplot2::stat_summary(
      fun.data = function(y) {
        m <- stats::median(y, na.rm = TRUE)
        data.frame(y = m, ymin = m, ymax = m)
      },
      geom = "crossbar",
      width = 0.50,
      fatten = 0,
      linewidth = 0.38,
      color = "black",
      fill = NA,
      show.legend = FALSE
    ) +
    ggplot2::geom_vline(xintercept = 1.5, color = "black", linewidth = 0.35) +
    ggplot2::scale_fill_manual(values = celltype_cols_existing, drop = FALSE) +
    ggplot2::scale_x_discrete(labels = celltype_labels) +
    ggplot2::coord_cartesian(ylim = c(y_min, y_max)) +
    ggplot2::scale_y_continuous(breaks = seq(0, 1, 0.05), expand = ggplot2::expansion(mult = c(0.02, 0.02))) +
    paper_theme(base_size = 10) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    ggplot2::labs(title = "c", x = NULL, y = "WF1 score")

  save_gg_pdf(file.path(out_dir, "BulkFig_c_fine13_WF1_boxplot.pdf"), p_c, width = 4.2, height = 3.2)
  ggplot2::ggsave(file.path(out_dir, "BulkFig_c_fine13_WF1_boxplot.png"), p_c, width = 4.2, height = 3.2, dpi = 300)

  invisible(p_c)
}

# =========================================================
# 6. Panel d: fine disease x cell type F1 heatmap
# =========================================================
# F1 is computed after pooling predictions from all five folds.
plot_panel_d <- function() {
  rows <- list()

  for (d in disease_order) {
    rows[[length(rows) + 1]] <- data.frame(
      Disease = d,
      Celltype = "Majority vote",
      F1 = compute_f1_bin(maj_df$true_short, maj_df$pred_short, d),
      stringsAsFactors = FALSE
    )
  }

  for (ct in setdiff(celltype_order_existing, "Majority vote")) {
    sub_ct <- dplyr::filter(pred_df, as.character(.data$Level1_plot) == ct)
    for (d in disease_order) {
      rows[[length(rows) + 1]] <- data.frame(
        Disease = d,
        Celltype = ct,
        F1 = compute_f1_bin(sub_ct$true_short, sub_ct$pred_short, d),
        stringsAsFactors = FALSE
      )
    }
  }

  d_long <- dplyr::bind_rows(rows) %>%
    dplyr::mutate(
      Disease = factor(Disease, levels = disease_order),
      Celltype = factor(Celltype, levels = celltype_order_existing)
    )

  data.table::fwrite(d_long, file.path(out_dir, "SourceData_BulkFig_d_fine13_disease_celltype_F1.tsv"), sep = "\t")

  d_mat <- d_long %>%
    dplyr::select(Disease, Celltype, F1) %>%
    tidyr::pivot_wider(names_from = Celltype, values_from = F1) %>%
    as.data.frame()

  rownames(d_mat) <- as.character(d_mat$Disease)
  d_mat$Disease <- NULL
  d_mat <- as.matrix(d_mat)
  d_mat[is.na(d_mat)] <- 0

  row_group <- droplevels(label_to_group(rownames(d_mat)))

  ha_bottom <- ComplexHeatmap::HeatmapAnnotation(
    `Cell types` = factor(colnames(d_mat), levels = celltype_order_existing),
    col = list(`Cell types` = celltype_cols_existing),
    show_annotation_name = FALSE,
    show_legend = FALSE,
    which = "column",
    simple_anno_size = grid::unit(3.2, "mm"),
    border = FALSE
  )

  ra_left <- ComplexHeatmap::rowAnnotation(
    `Disease group` = row_group,
    col = list(`Disease group` = bulk_group_cols[levels(row_group)]),
    show_annotation_name = FALSE,
    show_legend = FALSE,
    simple_anno_size = grid::unit(3.2, "mm"),
    border = FALSE
  )

  ht_d <- ComplexHeatmap::Heatmap(
    d_mat,
    name = "F1 score",
    col = score_col_fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    left_annotation = ra_left,
    bottom_annotation = ha_bottom,
    row_split = row_group,
    row_gap = grid::unit(1.1, "mm"),
    column_labels = unname(celltype_labels[colnames(d_mat)]),
    column_names_rot = 90,
    column_names_gp = grid::gpar(fontsize = 8, fontfamily = base_family),
    row_names_gp = grid::gpar(fontsize = 8, fontfamily = base_family),
    rect_gp = grid::gpar(col = NA),
    na_col = "#000004",
    heatmap_legend_param = score_heatmap_legend_param("F1 score"),
    column_title = NULL,
    row_title = NULL,
    width = grid::unit(ncol(d_mat) * 4.8, "mm"),
    height = grid::unit(nrow(d_mat) * 4.8, "mm")
  )

  d_w <- max(4.2, ncol(d_mat) * 0.24 + 2.3)
  d_h <- max(4.4, nrow(d_mat) * 0.25 + 1.9)

  open_pdf(file.path(out_dir, "BulkFig_d_fine13_disease_celltype_F1_heatmap.pdf"), width = d_w, height = d_h)
  ComplexHeatmap::draw(ht_d, heatmap_legend_side = "top", annotation_legend_side = "bottom", newpage = TRUE)
  dev.off()

  png(file.path(out_dir, "BulkFig_d_fine13_disease_celltype_F1_heatmap.png"), width = d_w, height = d_h, units = "in", res = 300)
  ComplexHeatmap::draw(ht_d, heatmap_legend_side = "top", annotation_legend_side = "bottom", newpage = TRUE)
  dev.off()

  invisible(ht_d)
}

# =========================================================
# 7. Panel e: fine13 column normalized confusion matrix
# =========================================================
# Columns are true diseases; rows are predicted diseases. Values are normalized
# within each true disease column, so diagonal values correspond to recall.
plot_panel_e <- function() {
  cm <- table(
    factor(maj_df$pred_short, levels = pred_row_order),
    factor(maj_df$true_short, levels = disease_order)
  )
  cm <- as.matrix(cm)
  cm_col <- sweep(cm, 2, colSums(cm), "/")
  cm_col[is.na(cm_col)] <- 0

  cm_count_df <- as.data.frame.matrix(cm)
  cm_count_df <- cbind(Predicted_condition = rownames(cm_count_df), cm_count_df)
  data.table::fwrite(cm_count_df, file.path(out_dir, "SourceData_BulkFig_e_confusion_counts.tsv"), sep = "\t")

  cm_col_df <- as.data.frame.matrix(cm_col)
  cm_col_df <- cbind(Predicted_condition = rownames(cm_col_df), cm_col_df)
  data.table::fwrite(cm_col_df, file.path(out_dir, "SourceData_BulkFig_e_confusion_colnorm.tsv"), sep = "\t")

  top_group <- droplevels(label_to_group(colnames(cm_col)))
  right_group <- droplevels(label_to_group(rownames(cm_col)))

  all_group_levels <- unique(c(levels(top_group), levels(right_group)))
  annotation_group_cols <- bulk_group_cols[all_group_levels]
  if (any(is.na(annotation_group_cols))) {
    missing <- names(annotation_group_cols)[is.na(annotation_group_cols)]
    stop("Missing disease group colors for: ", paste(missing, collapse = ", "))
  }

  ha_top <- ComplexHeatmap::HeatmapAnnotation(
    `True condition group` = top_group,
    col = list(`True condition group` = annotation_group_cols),
    show_annotation_name = FALSE,
    show_legend = FALSE,
    simple_anno_size = grid::unit(3.2, "mm"),
    border = FALSE
  )

  ra_right <- ComplexHeatmap::rowAnnotation(
    `Predicted group` = right_group,
    col = list(`Predicted group` = annotation_group_cols),
    show_annotation_name = FALSE,
    show_legend = FALSE,
    simple_anno_size = grid::unit(3.2, "mm"),
    border = FALSE
  )

  ht_e <- ComplexHeatmap::Heatmap(
    cm_col,
    name = "True condition proportion",
    col = score_col_fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    top_annotation = ha_top,
    right_annotation = ra_right,
    row_split = right_group,
    column_split = top_group,
    row_gap = grid::unit(1.1, "mm"),
    column_gap = grid::unit(1.1, "mm"),
    column_names_rot = 90,
    column_names_gp = grid::gpar(fontsize = 8, fontfamily = base_family),
    row_names_gp = grid::gpar(fontsize = 8, fontfamily = base_family),
    rect_gp = grid::gpar(col = "#1B1B1B", lwd = 0.20),
    na_col = "#000004",
    heatmap_legend_param = score_heatmap_legend_param("True condition proportion"),
    column_title = NULL,
    row_title = NULL,
    width = grid::unit(ncol(cm_col) * 4.9, "mm"),
    height = grid::unit(nrow(cm_col) * 4.9, "mm")
  )

  e_w <- max(4.6, ncol(cm_col) * 0.26 + 2.4)
  e_h <- max(4.8, nrow(cm_col) * 0.27 + 1.9)

  open_pdf(file.path(out_dir, "BulkFig_e_fine13_confusion_matrix.pdf"), width = e_w, height = e_h)
  ComplexHeatmap::draw(ht_e, heatmap_legend_side = "top", annotation_legend_side = "bottom", merge_legends = TRUE, newpage = TRUE)
  dev.off()

  png(file.path(out_dir, "BulkFig_e_fine13_confusion_matrix.png"), width = e_w, height = e_h, units = "in", res = 300)
  ComplexHeatmap::draw(ht_e, heatmap_legend_side = "top", annotation_legend_side = "bottom", merge_legends = TRUE, newpage = TRUE)
  dev.off()

  invisible(ht_e)
}

# -------------------------------
# 8. Run panels separately
# -------------------------------
panel_c <- plot_panel_c()
panel_d <- plot_panel_d()
panel_e <- plot_panel_e()

# -------------------------------
# 9. Summary
# -------------------------------
summary_list <- list(
  fine_dir = fine_dir,
  out_dir = out_dir,
  n_majority_predictions = nrow(maj_df),
  diseases = disease_order,
  predicted_rows = pred_row_order,
  cell_types = celltype_order_existing,
  note = "Bulk fine13 c-e plotted as separated paper mimic panels. Panel e is column normalized by true disease."
)

if (file.exists(overall_json)) {
  summary_list$overall <- jsonlite::fromJSON(overall_json)
}

jsonlite::write_json(
  summary_list,
  file.path(out_dir, "bulk_fine13_cde_separate_summary.json"),
  auto_unbox = TRUE,
  pretty = TRUE
)

cat("[DONE] Bulk fine13 c-e separated paper mimic plots saved to:\n")
cat(out_dir, "\n")
cat("Generated files:\n")
cat("  - BulkFig_c_fine13_WF1_boxplot.pdf/png\n")
cat("  - BulkFig_d_fine13_disease_celltype_F1_heatmap.pdf/png\n")
cat("  - BulkFig_e_fine13_confusion_matrix.pdf/png\n")
cat("  - SourceData_BulkFig_c/d/e TSV files\n")
