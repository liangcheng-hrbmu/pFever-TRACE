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
# Bulk supplementary figures S6 and S7
#
# S6: Fine 13 disease classifier detailed performance
#   a. Per disease majority vote metrics across five folds
#   b. Majority vote and Level 1 classifier BAS and MCC across five folds
#   c. Pooled recall and precision heatmaps by disease and cell type
#
# S7: Broad 5 category classifier detailed performance
#   a. Per category majority vote metrics across five folds
#   b. Majority vote and Level 1 classifier BAS and MCC across five folds
#   c. Pooled recall and precision heatmaps by category and cell type
#
# This script follows the current paper mimic style used for the bulk main figure.
# =========================================================

# -------------------------------
# 0. Paths
# -------------------------------
fine_dir  <- "data/analysis/09_bulk_deconvolution/08_bulk_classifier/03_patient_classifier_fine_bulkZ_L1x11_final"
broad_dir <- "data/analysis/09_bulk_deconvolution/08_bulk_classifier/03_patient_classifier_broad5_bulkZ_L1x11_final"

out_dir <- "data/analysis/09_bulk_deconvolution/08_bulk_classifier/06_bulk_supplementary_S6_S7"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------------
# 1. Devices and theme
# -------------------------------
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

paper_theme <- function(base_size = 10) {
  ggplot2::theme_classic(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0, size = base_size + 4),
      axis.title = ggplot2::element_text(size = base_size + 1, color = "black"),
      axis.text = ggplot2::element_text(size = base_size, color = "black"),
      axis.line = ggplot2::element_line(linewidth = 0.35, color = "black"),
      axis.ticks = ggplot2::element_line(linewidth = 0.30, color = "black"),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(size = base_size + 1, color = "black"),
      legend.position = "none",
      plot.margin = ggplot2::margin(4, 4, 4, 4, "pt")
    )
}

# -------------------------------
# 2. Shared palettes
# -------------------------------
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
    legend_width = grid::unit(26, "mm"),
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

fine_order <- c(
  "IAV", "RSV", "Adenovirus", "Rhinovirus",
  "Meningococcus", "Pneumococcus", "S. aureus", "GAS", "E. coli", "Tuberculosis",
  "Malaria", "KD", "JIA"
)

broad_order <- c(
  "Viral infection",
  "Bacterial infection",
  "Malaria",
  "Hyperinflammatory disease",
  "Autoimmune disease"
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

broad_group_map <- c(
  "Viral infection"           = "Viral infection",
  "Bacterial infection"       = "Bacterial infection",
  "Malaria"                   = "Malaria",
  "Hyperinflammatory disease" = "Hyperinflammatory disease",
  "Autoimmune disease"        = "Autoimmune disease",
  "Not assigned"              = "Not assigned"
)

group_cols <- c(
  "Viral infection"           = "#6699CC",
  "Bacterial infection"       = "#C48EAD",
  "Malaria"                   = "#6A4C93",
  "Hyperinflammatory disease" = "#E3A722",
  "Autoimmune disease"        = "#36A596",
  "Other"                     = "#BDBDBD",
  "Not assigned"              = "#2F225B"
)

group_order <- c(
  "Viral infection",
  "Bacterial infection",
  "Malaria",
  "Hyperinflammatory disease",
  "Autoimmune disease",
  "Other",
  "Not assigned"
)

fine_box_cols <- c(
  "IAV"           = "#D1BA5A",
  "RSV"           = "#6AA7A2",
  "Adenovirus"    = "#84B97A",
  "Rhinovirus"    = "#5E8E6A",
  "Meningococcus" = "#C98E7A",
  "Pneumococcus"  = "#E7A789",
  "S. aureus"     = "#C86A6A",
  "GAS"           = "#B85B63",
  "E. coli"       = "#D97A50",
  "Tuberculosis"  = "#A97C50",
  "Malaria"       = "#8F6FB2",
  "KD"            = "#C9A33A",
  "JIA"           = "#4B6B8F"
)

broad_box_cols <- c(
  "Viral infection"           = "#6699CC",
  "Bacterial infection"       = "#C48EAD",
  "Malaria"                   = "#6A4C93",
  "Hyperinflammatory disease" = "#E3A722",
  "Autoimmune disease"        = "#36A596"
)

# -------------------------------
# 3. Helper functions
# -------------------------------
pick_col <- function(df, candidates, required = TRUE) {
  hit <- intersect(candidates, colnames(df))
  if (length(hit) == 0) {
    if (required) stop("Cannot find any of columns: ", paste(candidates, collapse = ", "))
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

short_broad <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[is.na(x) | x == "" | x %in% c("NA", "NaN", "nan", "NULL", "None")] <- "Not assigned"
  map <- c(
    "virus" = "Viral infection",
    "Virus" = "Viral infection",
    "Viral infection" = "Viral infection",
    "bacteria" = "Bacterial infection",
    "Bacteria" = "Bacterial infection",
    "Bacterial infection" = "Bacterial infection",
    "malaria" = "Malaria",
    "Malaria" = "Malaria",
    "Hyperinflam" = "Hyperinflammatory disease",
    "Hyperinflammatory disease" = "Hyperinflammatory disease",
    "Autoimmune" = "Autoimmune disease",
    "Autoimmune disease" = "Autoimmune disease",
    "TIE" = "Not assigned",
    "Not assigned" = "Not assigned"
  )
  out <- ifelse(x %in% names(map), map[x], x)
  unname(out)
}

complete_named_cols <- function(levels_vec, base_cols, default = "#BDBDBD") {
  missing <- setdiff(levels_vec, names(base_cols))
  out <- c(base_cols, stats::setNames(rep(default, length(missing)), missing))
  out[levels_vec]
}

label_to_group <- function(x, task = c("fine", "broad")) {
  task <- match.arg(task)
  x <- as.character(x)
  if (task == "fine") {
    out <- fine_group_map[x]
  } else {
    out <- broad_group_map[x]
  }
  out[is.na(out) & x == "Not assigned"] <- "Not assigned"
  out[is.na(out)] <- "Other"
  factor(out, levels = group_order)
}

compute_one_vs_rest <- function(y_true, y_pred, positive_class) {
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
  tn <- sum(yt == 0 & yp == 0)

  precision <- ifelse((tp + fp) == 0, 0, tp / (tp + fp))
  recall    <- ifelse((tp + fn) == 0, 0, tp / (tp + fn))
  f1        <- ifelse((2 * tp + fp + fn) == 0, 0, 2 * tp / (2 * tp + fp + fn))
  specificity <- ifelse((tn + fp) == 0, 0, tn / (tn + fp))
  bas <- (recall + specificity) / 2

  data.frame(
    support = sum(yt == 1),
    precision = precision,
    recall = recall,
    F1 = f1,
    BAS = bas,
    stringsAsFactors = FALSE
  )
}

read_level1_predictions <- function(cls_dir) {
  merged_fp <- file.path(cls_dir, "query_level1_predictions_all.tsv")
  pred_dir <- file.path(cls_dir, "predictions")

  if (file.exists(merged_fp)) {
    message("[READ] ", merged_fp)
    return(data.table::fread(merged_fp, sep = "\t", data.table = FALSE))
  }

  fps <- list.files(pred_dir, pattern = "^query_level1_predictions_SPLIT_.*\\.tsv$", full.names = TRUE)
  if (length(fps) == 0) stop("Cannot find merged or split Level 1 prediction files in: ", cls_dir)

  message("[READ] merging ", length(fps), " split Level 1 prediction files")
  df <- dplyr::bind_rows(lapply(fps, function(fp) {
    x <- data.table::fread(fp, sep = "\t", data.table = FALSE)
    x$.source_file <- basename(fp)
    x
  }))
  data.table::fwrite(df, merged_fp, sep = "\t")
  message("[WRITE] ", merged_fp)
  df
}

make_heat_matrix <- function(df, value_col) {
  mat_df <- df %>%
    dplyr::select(Disease, Celltype, value = dplyr::all_of(value_col)) %>%
    tidyr::pivot_wider(names_from = Celltype, values_from = value) %>%
    as.data.frame()
  rownames(mat_df) <- as.character(mat_df$Disease)
  mat_df$Disease <- NULL
  mat <- as.matrix(mat_df)
  mat[is.na(mat)] <- 0
  mat
}

# -------------------------------
# 4. Core plotting function
# -------------------------------
plot_bulk_supplementary <- function(task = c("fine", "broad")) {
  task <- match.arg(task)

  if (task == "fine") {
    cls_dir <- fine_dir
    prefix <- "SuppFigS6"
    short_fun <- short_fine
    disease_order0 <- fine_order
    disease_fill_cols <- fine_box_cols
    task_label <- "fine"
  } else {
    cls_dir <- broad_dir
    prefix <- "SuppFigS7"
    short_fun <- short_broad
    disease_order0 <- broad_order
    disease_fill_cols <- broad_box_cols
    task_label <- "broad"
  }

  split_df  <- data.table::fread(file.path(cls_dir, "scenario1_split_summary.tsv"), sep = "\t", data.table = FALSE)
  level1_df <- data.table::fread(file.path(cls_dir, "scenario1_level1_metrics_all.tsv"), sep = "\t", data.table = FALSE)
  maj_df    <- data.table::fread(file.path(cls_dir, "scenario1_majority_predictions_all.tsv"), sep = "\t", data.table = FALSE)
  pred_df   <- read_level1_predictions(cls_dir)

  maj_df <- filter_primary(maj_df)
  pred_df <- filter_primary(pred_df)

  true_col <- pick_col(maj_df, c("true_disease", "disease_broad", "disease_clean", "disease", "disease_Level_3", "true_label"))
  pred_col <- pick_col(maj_df, c("pred_disease_majority", "pred_majority", "majority_pred", "pred_disease", "prediction", "pred_label"))
  level1_col <- pick_col(level1_df, c("Level1", "classifier_ct", "celltype", "Celltype", "cell_type", "Cell_type"))
  pred_level1_col <- pick_col(pred_df, c("Level1", "classifier_ct", "celltype", "Celltype", "cell_type", "Cell_type"))
  pred_true_col <- pick_col(pred_df, c("disease_broad", "true_disease", "disease_clean", "disease", "disease_Level_3", "true_label"))
  pred_pred_col <- pick_col(pred_df, c("pred_disease", "pred_label", "prediction", "pred_disease_majority", "pred_majority", "majority_pred"))

  split_col_maj <- pick_col(maj_df, c("split", "fold", "cv_split"))
  split_col_splitdf <- pick_col(split_df, c("split", "fold", "cv_split"))
  split_col_level1df <- pick_col(level1_df, c("split", "fold", "cv_split"))

  maj_df$split_id <- as.character(maj_df[[split_col_maj]])
  level1_df$split_id <- as.character(level1_df[[split_col_level1df]])
  split_df$split_id <- as.character(split_df[[split_col_splitdf]])

  level1_df$Level1_plot <- short_celltype(level1_df[[level1_col]])
  pred_df$Level1_plot <- short_celltype(pred_df[[pred_level1_col]])
  maj_df$true_short <- short_fun(maj_df[[true_col]])
  maj_df$pred_short <- short_fun(maj_df[[pred_col]])
  pred_df$true_short <- short_fun(pred_df[[pred_true_col]])
  pred_df$pred_short <- short_fun(pred_df[[pred_pred_col]])

  truth_labels <- setdiff(unique(maj_df$true_short), "Not assigned")
  disease_order <- c(disease_order0[disease_order0 %in% truth_labels], setdiff(sort(truth_labels), disease_order0))
  celltype_order_existing <- celltype_order[celltype_order %in% c("Majority vote", unique(level1_df$Level1_plot), unique(pred_df$Level1_plot))]
  celltype_cols_existing <- complete_named_cols(celltype_order_existing, celltype_cols)
  disease_fill_cols <- complete_named_cols(disease_order, disease_fill_cols)

  # ======================================
  # Panel a: per disease majority vote metrics
  # ======================================
  metric_rows <- list()
  for (sp in sort(unique(maj_df$split_id))) {
    sub_sp <- dplyr::filter(maj_df, .data$split_id == sp)
    for (d in disease_order) {
      met <- compute_one_vs_rest(sub_sp$true_short, sub_sp$pred_short, d)
      metric_rows[[length(metric_rows) + 1]] <- data.frame(
        split = sp,
        disease = d,
        support = met$support,
        BalancedAccuracy = met$BAS,
        F1 = met$F1,
        Precision = met$precision,
        Recall = met$recall,
        stringsAsFactors = FALSE
      )
    }
  }

  metric_df_a <- dplyr::bind_rows(metric_rows) %>%
    dplyr::mutate(disease = factor(disease, levels = disease_order))

  data.table::fwrite(metric_df_a, file.path(out_dir, paste0(prefix, "a_source_data_majority_per_disease_metrics.tsv")), sep = "\t")

  avg_support_df <- metric_df_a %>%
    dplyr::group_by(disease) %>%
    dplyr::summarise(avg_support = mean(support, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(
      disease = factor(disease, levels = disease_order),
      avg_support_lab = sprintf("%.1f", avg_support)
    )

  metric_long_a <- metric_df_a %>%
    dplyr::select(split, disease, support, BalancedAccuracy, F1, Precision, Recall) %>%
    tidyr::pivot_longer(
      cols = c(BalancedAccuracy, F1, Precision, Recall),
      names_to = "Metric",
      values_to = "Score"
    ) %>%
    dplyr::mutate(
      Metric = factor(Metric,
                      levels = c("BalancedAccuracy", "F1", "Precision", "Recall"),
                      labels = c("Balanced accuracy", "F1", "Precision", "Recall")),
      disease = factor(disease, levels = disease_order)
    )

  label_df_a <- tidyr::crossing(
    disease = factor(disease_order, levels = disease_order),
    Metric = levels(metric_long_a$Metric)
  ) %>%
    dplyr::left_join(avg_support_df, by = "disease") %>%
    dplyr::mutate(y = 1.05)

  p_a <- ggplot2::ggplot(metric_long_a, ggplot2::aes(x = disease, y = Score, fill = disease)) +
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
    ggplot2::geom_text(
      data = label_df_a,
      ggplot2::aes(x = disease, y = y, label = avg_support_lab),
      inherit.aes = FALSE,
      size = 2.7,
      family = base_family,
      vjust = 0
    ) +
    ggplot2::annotate("text", x = 0.25, y = 1.05, label = "average #\nof samples", hjust = 0,
                      vjust = 0, size = 3.0, family = base_family) +
    ggplot2::facet_wrap(~Metric, nrow = 2, scales = "fixed") +
    ggplot2::scale_fill_manual(values = disease_fill_cols[disease_order], guide = "none") +
    ggplot2::coord_cartesian(ylim = c(-0.05, 1.08), clip = "off") +
    ggplot2::scale_y_continuous(breaks = seq(0, 1, 0.2), expand = ggplot2::expansion(mult = c(0.02, 0.02))) +
    paper_theme(base_size = 9.5) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1),
      strip.text = ggplot2::element_text(size = 10),
      panel.spacing = grid::unit(4, "mm")
    ) +
    ggplot2::labs(title = "a", x = NULL, y = "score")

  width_a <- if (task == "fine") 8.8 else 6.8
  height_a <- if (task == "fine") 7.6 else 6.8
  save_gg_pdf(file.path(out_dir, paste0(prefix, "a_majority_per_disease_metrics_boxplots.pdf")), p_a, width = width_a, height = height_a)
  ggplot2::ggsave(file.path(out_dir, paste0(prefix, "a_majority_per_disease_metrics_boxplots.png")), p_a, width = width_a, height = height_a, dpi = 300)

  # ======================================
  # Panel b: majority vote and Level 1 classifier BAS and MCC
  # ======================================
  bas_majority_col <- pick_col(split_df, c("BAS_majority", "balanced_accuracy_majority", "BalancedAccuracy_majority", "BAS"))
  mcc_majority_col <- pick_col(split_df, c("MCC_majority", "MCC"))
  bas_level1_col <- pick_col(level1_df, c("BAS", "balanced_accuracy", "BalancedAccuracy"))
  mcc_level1_col <- pick_col(level1_df, c("MCC"))

  metric_df_b <- dplyr::bind_rows(
    split_df %>%
      dplyr::transmute(
        split = .data$split_id,
        Celltype = "Majority vote",
        BAS = as.numeric(.data[[bas_majority_col]]),
        MCC = as.numeric(.data[[mcc_majority_col]])
      ),
    level1_df %>%
      dplyr::transmute(
        split = .data$split_id,
        Celltype = .data$Level1_plot,
        BAS = as.numeric(.data[[bas_level1_col]]),
        MCC = as.numeric(.data[[mcc_level1_col]])
      )
  ) %>%
    dplyr::mutate(Celltype = factor(Celltype, levels = celltype_order_existing)) %>%
    tidyr::pivot_longer(cols = c(BAS, MCC), names_to = "Metric", values_to = "Score") %>%
    dplyr::filter(!is.na(Score)) %>%
    dplyr::mutate(
      Metric = factor(Metric,
                      levels = c("BAS", "MCC"),
                      labels = c("Balanced accuracy score", "Matthew correlation coefficient"))
    )

  data.table::fwrite(metric_df_b, file.path(out_dir, paste0(prefix, "b_source_data_classifier_metrics.tsv")), sep = "\t")

  p_b <- ggplot2::ggplot(metric_df_b, ggplot2::aes(x = Celltype, y = Score, fill = Celltype)) +
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
    ggplot2::facet_wrap(~Metric, ncol = 1) +
    ggplot2::scale_fill_manual(values = celltype_cols_existing, drop = FALSE) +
    ggplot2::scale_x_discrete(labels = celltype_labels) +
    ggplot2::coord_cartesian(ylim = c(min(metric_df_b$Score, na.rm = TRUE) - 0.03, 1.0)) +
    ggplot2::scale_y_continuous(breaks = seq(0, 1, 0.1), expand = ggplot2::expansion(mult = c(0.02, 0.02))) +
    paper_theme(base_size = 9.5) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1),
      strip.text = ggplot2::element_text(size = 10),
      panel.spacing = grid::unit(4, "mm")
    ) +
    ggplot2::labs(title = "b", x = NULL, y = NULL)

  save_gg_pdf(file.path(out_dir, paste0(prefix, "b_classifier_metrics_boxplots.pdf")), p_b, width = 3.3, height = 4.9)
  ggplot2::ggsave(file.path(out_dir, paste0(prefix, "b_classifier_metrics_boxplots.png")), p_b, width = 3.3, height = 4.9, dpi = 300)

  # ======================================
  # Panel c: pooled recall and precision heatmaps
  # ======================================
  heat_rows <- list()
  for (d in disease_order) {
    met <- compute_one_vs_rest(maj_df$true_short, maj_df$pred_short, d)
    heat_rows[[length(heat_rows) + 1]] <- data.frame(
      Disease = d,
      Celltype = "Majority vote",
      Recall = met$recall,
      Precision = met$precision,
      stringsAsFactors = FALSE
    )
  }

  for (ct in setdiff(celltype_order_existing, "Majority vote")) {
    sub_ct <- dplyr::filter(pred_df, as.character(.data$Level1_plot) == ct)
    for (d in disease_order) {
      met <- compute_one_vs_rest(sub_ct$true_short, sub_ct$pred_short, d)
      heat_rows[[length(heat_rows) + 1]] <- data.frame(
        Disease = d,
        Celltype = ct,
        Recall = met$recall,
        Precision = met$precision,
        stringsAsFactors = FALSE
      )
    }
  }

  heat_long <- dplyr::bind_rows(heat_rows) %>%
    dplyr::mutate(
      Disease = factor(Disease, levels = disease_order),
      Celltype = factor(Celltype, levels = celltype_order_existing)
    )

  data.table::fwrite(heat_long, file.path(out_dir, paste0(prefix, "c_source_data_recall_precision.tsv")), sep = "\t")

  recall_mat <- make_heat_matrix(heat_long, "Recall")
  precision_mat <- make_heat_matrix(heat_long, "Precision")

  ha_bottom <- ComplexHeatmap::HeatmapAnnotation(
    `Cell types` = factor(colnames(recall_mat), levels = celltype_order_existing),
    col = list(`Cell types` = celltype_cols_existing),
    show_annotation_name = FALSE,
    show_legend = FALSE,
    which = "column",
    simple_anno_size = grid::unit(3.2, "mm"),
    border = FALSE
  )

  row_group <- droplevels(label_to_group(rownames(recall_mat), task = task_label))
  row_group_cols <- group_cols[levels(row_group)]

  ra_left <- ComplexHeatmap::rowAnnotation(
    `Disease group` = row_group,
    col = list(`Disease group` = row_group_cols),
    show_annotation_name = FALSE,
    show_legend = FALSE,
    simple_anno_size = grid::unit(3.2, "mm"),
    border = FALSE
  )

  ht_recall <- ComplexHeatmap::Heatmap(
    recall_mat,
    name = "Recall score",
    col = score_col_fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    left_annotation = ra_left,
    bottom_annotation = ha_bottom,
    row_split = row_group,
    row_gap = grid::unit(1.1, "mm"),
    column_labels = unname(celltype_labels[colnames(recall_mat)]),
    column_names_rot = 90,
    column_names_gp = grid::gpar(fontsize = 8, fontfamily = base_family),
    row_names_gp = grid::gpar(fontsize = 8.5, fontfamily = base_family),
    rect_gp = grid::gpar(col = NA),
    heatmap_legend_param = score_heatmap_legend_param("Recall score"),
    column_title = "Recall score",
    column_title_gp = grid::gpar(fontsize = 11, fontfamily = base_family),
    row_title = NULL,
    width = grid::unit(ncol(recall_mat) * 4.6, "mm"),
    height = grid::unit(nrow(recall_mat) * ifelse(task == "fine", 4.6, 6.0), "mm")
  )

  ht_precision <- ComplexHeatmap::Heatmap(
    precision_mat,
    name = "Precision score",
    col = score_col_fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    left_annotation = ra_left,
    bottom_annotation = ha_bottom,
    row_split = row_group,
    row_gap = grid::unit(1.1, "mm"),
    column_labels = unname(celltype_labels[colnames(precision_mat)]),
    column_names_rot = 90,
    column_names_gp = grid::gpar(fontsize = 8, fontfamily = base_family),
    row_names_gp = grid::gpar(fontsize = 8.5, fontfamily = base_family),
    rect_gp = grid::gpar(col = NA),
    heatmap_legend_param = score_heatmap_legend_param("Precision score"),
    column_title = "Precision score",
    column_title_gp = grid::gpar(fontsize = 11, fontfamily = base_family),
    row_title = NULL,
    width = grid::unit(ncol(precision_mat) * 4.6, "mm"),
    height = grid::unit(nrow(precision_mat) * ifelse(task == "fine", 4.6, 6.0), "mm")
  )

  w_heat <- if (task == "fine") max(6.6, ncol(recall_mat) * 0.23 + 4.2) else max(6.2, ncol(recall_mat) * 0.24 + 3.8)
  h_heat <- if (task == "fine") max(4.8, nrow(recall_mat) * 0.21 + 1.7) else max(4.0, nrow(recall_mat) * 0.35 + 1.4)

  open_pdf(file.path(out_dir, paste0(prefix, "c_recall_heatmap.pdf")), width = w_heat, height = h_heat)
  ComplexHeatmap::draw(ht_recall, heatmap_legend_side = "top", annotation_legend_side = "bottom", newpage = TRUE)
  dev.off()
  png(file.path(out_dir, paste0(prefix, "c_recall_heatmap.png")), width = w_heat, height = h_heat, units = "in", res = 300)
  ComplexHeatmap::draw(ht_recall, heatmap_legend_side = "top", annotation_legend_side = "bottom", newpage = TRUE)
  dev.off()

  open_pdf(file.path(out_dir, paste0(prefix, "c_precision_heatmap.pdf")), width = w_heat, height = h_heat)
  ComplexHeatmap::draw(ht_precision, heatmap_legend_side = "top", annotation_legend_side = "bottom", newpage = TRUE)
  dev.off()
  png(file.path(out_dir, paste0(prefix, "c_precision_heatmap.png")), width = w_heat, height = h_heat, units = "in", res = 300)
  ComplexHeatmap::draw(ht_precision, heatmap_legend_side = "top", annotation_legend_side = "bottom", newpage = TRUE)
  dev.off()

  # combined quick check
  combined_w <- if (task == "fine") 10.8 else 10.2
  open_pdf(file.path(out_dir, paste0(prefix, "c_recall_precision_combined.pdf")), width = combined_w, height = h_heat)
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow = 1, ncol = 2, widths = grid::unit(c(1, 1), "null"))))
  grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  ComplexHeatmap::draw(ht_recall, heatmap_legend_side = "top", annotation_legend_side = "bottom", newpage = FALSE)
  grid::upViewport()
  grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
  ComplexHeatmap::draw(ht_precision, heatmap_legend_side = "top", annotation_legend_side = "bottom", newpage = FALSE)
  grid::upViewport(2)
  dev.off()

  png(file.path(out_dir, paste0(prefix, "c_recall_precision_combined.png")), width = combined_w, height = h_heat, units = "in", res = 300)
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow = 1, ncol = 2, widths = grid::unit(c(1, 1), "null"))))
  grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  ComplexHeatmap::draw(ht_recall, heatmap_legend_side = "top", annotation_legend_side = "bottom", newpage = FALSE)
  grid::upViewport()
  grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
  ComplexHeatmap::draw(ht_precision, heatmap_legend_side = "top", annotation_legend_side = "bottom", newpage = FALSE)
  grid::upViewport(2)
  dev.off()

  summary_list <- list(
    task = task,
    cls_dir = cls_dir,
    diseases = disease_order,
    cell_types = celltype_order_existing,
    n_majority_predictions = nrow(maj_df),
    files = c(
      paste0(prefix, "a_majority_per_disease_metrics_boxplots.pdf/png"),
      paste0(prefix, "b_classifier_metrics_boxplots.pdf/png"),
      paste0(prefix, "c_recall_heatmap.pdf/png"),
      paste0(prefix, "c_precision_heatmap.pdf/png"),
      paste0(prefix, "c_recall_precision_combined.pdf/png")
    )
  )

  jsonlite::write_json(
    summary_list,
    file.path(out_dir, paste0(prefix, "_summary.json")),
    auto_unbox = TRUE,
    pretty = TRUE
  )
}

# -------------------------------
# 5. Run
# -------------------------------
plot_bulk_supplementary("fine")
plot_bulk_supplementary("broad")

cat("[DONE] Bulk supplementary figures S6 and S7 saved to:\n")
cat(out_dir, "\n")
cat("Generated figure groups:\n")
cat("  - SuppFigS6a/b/c for fine 13 disease task\n")
cat("  - SuppFigS7a/b/c for broad 5 category task\n")
