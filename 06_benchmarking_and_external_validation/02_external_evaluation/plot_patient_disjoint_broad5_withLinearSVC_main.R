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
# Patient-disjoint broad5 evaluation, with LinearSVC
# 222 disease samples remaining after exclusion of participants
# linked to the model-development cohort
#
# a. WF1 for pFever majority vote, published Med139 signature,
#    random baseline, and each Level 1 classifier
# b. Disease by Level 1 F1 heatmap
# c. True-condition-normalized confusion matrix
#
# Metrics are recomputed from the patient-disjoint prediction outputs.
# =========================================================

# -------------------------------
# 0. Input and output
# -------------------------------

base_dir <- "data/analysis/09_bulk_deconvolution"
ext_root <- file.path(base_dir, "12_external_validation_multi_cohort_assets_v1")

cohort_id <- "Patient-disjoint evaluation"
task_name <- "broad5"

run_dir <- file.path(
  ext_root,
  "results",
  "patient_disjoint_sensitivity",
  "pfever_trace_output_withLinearSVC"
)

include_healthy <- FALSE

metadata_candidates <- c(
  file.path(
    base_dir,
    "external_validation",
    "EUCLIDS_PERFORM_external_v2",
    "metadata",
    "patient_disjoint_sensitivity_disease_metadata.tsv"
  ),
  file.path(
    ext_root,
    "metadata",
    "patient_disjoint_sensitivity_disease_metadata.tsv"
  )
)

majority_candidates <- c(
  file.path(run_dir, "pFever_Immune_predictions.csv"),
  file.path(run_dir, "prediction_primary_broad5", "pFever_Immune_predictions.csv"),
  file.path(run_dir, "prediction_primary_broad5", "broad5_majority_predictions.tsv"),
  file.path(run_dir, "broad5_majority_predictions.tsv")
)

level1_candidates <- c(
  file.path(run_dir, "prediction_primary_broad5", "broad5_level1_predictions.tsv"),
  file.path(run_dir, "broad5_level1_predictions.tsv"),
  file.path(run_dir, "prediction_primary_broad5", "level1_predictions.tsv")
)

med139_metrics_candidates <- c(
  file.path(
    ext_root,
    "results",
    "patient_disjoint_sensitivity",
    "Med139_evaluation",
    "broad5_external_metrics.tsv"
  )
)

out_dir <- file.path(
  ext_root,
  "plots",
  "patient_disjoint_broad5_withLinearSVC_abc"
)

eval_out_dir <- file.path(
  ext_root,
  "results",
  "patient_disjoint_sensitivity",
  "evaluation_withLinearSVC"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(eval_out_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------------
# 1. Device and theme
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
# 2. Palettes and orders
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
    legend_width = grid::unit(30, "mm"),
    grid_height = grid::unit(2.8, "mm"),
    border = NA,
    labels_gp = grid::gpar(fontsize = 8, fontfamily = base_family),
    title_gp = grid::gpar(fontsize = 9, fontfamily = base_family)
  )
}

celltype_order <- c(
  "Majority vote", "Published signature", "Random",
  "B", "Plasma", "pDC", "DC", "Mono",
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
  "Majority vote"       = "#535353",
  "Published signature" = "#1F78B4",
  "Random"              = "#D9D9D9",
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
  "Majority vote"       = "pFever TRACE",
  "Published signature" = "Published signature",
  "Random"              = "Random",
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

# -------------------------------
# 3. Helper functions
# -------------------------------

pick_existing <- function(paths, label) {
  hit <- unique(paths)[file.exists(unique(paths))]
  if (length(hit) == 0) {
    stop(
      "Missing ", label, ". Tried:\n",
      paste0("  - ", unique(paths), collapse = "\n"),
      call. = FALSE
    )
  }
  hit[1]
}

read_any <- function(fp) {
  if (grepl("\\.csv$", fp, ignore.case = TRUE)) {
    data.table::fread(fp, data.table = FALSE)
  } else {
    data.table::fread(fp, sep = "\t", data.table = FALSE)
  }
}

pick_col <- function(df, candidates, required = TRUE, label = "column") {
  cols <- colnames(df)
  lower_map <- stats::setNames(cols, tolower(cols))
  for (cc in candidates) {
    if (cc %in% cols) return(cc)
    if (tolower(cc) %in% names(lower_map)) return(lower_map[[tolower(cc)]])
  }
  if (required) {
    stop(
      "Cannot find ", label, ". Tried: ",
      paste(candidates, collapse = ", "),
      "\nCurrent columns: ",
      paste(cols, collapse = " | "),
      call. = FALSE
    )
  }
  NA_character_
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

short_broad <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  xl <- tolower(x)

  out <- rep(NA_character_, length(x))

  out[xl %in% c("viral infection", "viral", "virus", "viral disease", "viral infections")] <- "Viral infection"
  out[xl %in% c("bacterial infection", "bacterial", "bacteria", "bacterial disease", "bacterial infections")] <- "Bacterial infection"
  out[xl %in% c("malaria", "plasmodium")] <- "Malaria"
  out[xl %in% c("hyperinflammatory disease", "hyperinflammatory", "hyperinflam", "hyperinflammation")] <- "Hyperinflammatory disease"
  out[xl %in% c("autoimmune disease", "autoimmune", "autoimmunity")] <- "Autoimmune disease"
  out[xl %in% c("healthy", "normal", "control", "pediatric healthy")] <- "Healthy"
  out[xl %in% c("tie", "not assigned", "unassigned", "unknown", "", "na", "nan", "none", "null")] <- "Not assigned"

  out[is.na(out) & grepl("influenza|iav|respiratory syncytial|rsv|adeno|adenovirus|rhinovirus|rhino|dengue|denv|sars|covid|hfmd|enterovirus", xl)] <- "Viral infection"
  out[is.na(out) & grepl("mening|pneumo|staph|aureus|streptococcal|streptococcus|\\bgas\\b|escherichia|e\\.coli|e coli|ecoli|tuberculosis|\\btb\\b|sepsis", xl)] <- "Bacterial infection"
  out[is.na(out) & grepl("kawasaki|\\bkd\\b|mis-c|misc|multisystem inflammatory", xl)] <- "Hyperinflammatory disease"
  out[is.na(out) & grepl("juvenile idiopathic|juvenile idopathic|\\bjia\\b|arthritis|systemic lupus|\\bsle\\b", xl)] <- "Autoimmune disease"
  out[is.na(out) & grepl("malaria|plasmodium", xl)] <- "Malaria"

  out[is.na(out)] <- x[is.na(out)]
  out
}

label_to_group <- function(x) {
  x <- short_broad(x)
  out <- ifelse(x %in% names(broad_group_cols), x, "Other")
  factor(out, levels = c(broad_order, "Other"))
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

calc_multiclass_metrics <- function(y_true, y_pred, labels = NULL) {
  y_true <- as.character(y_true)
  y_pred <- as.character(y_pred)

  valid <- !is.na(y_true) & !is.na(y_pred)
  y_true <- y_true[valid]
  y_pred <- y_pred[valid]

  if (is.null(labels)) {
    labels <- sort(unique(y_true))
  }
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

complete_named_cols <- function(levels_vec, base_cols, default = "#BDBDBD") {
  missing <- setdiff(levels_vec, names(base_cols))
  out <- c(base_cols, stats::setNames(rep(default, length(missing)), missing))
  out[levels_vec]
}

# -------------------------------
# 4. Read LinearSVC predictions and metadata
# -------------------------------

metadata_fp <- pick_existing(metadata_candidates, "metadata")
majority_fp <- pick_existing(majority_candidates, "LinearSVC majority prediction")
level1_fp <- pick_existing(level1_candidates, "LinearSVC Level 1 prediction")
med139_metrics_fp <- pick_existing(
  med139_metrics_candidates,
  "Med139 patient-disjoint metrics"
)

cat("[READ] metadata:\n", metadata_fp, "\n")
cat("[READ] majority predictions:\n", majority_fp, "\n")
cat("[READ] level1 predictions:\n", level1_fp, "\n")
cat("[READ] Med139 metrics:\n", med139_metrics_fp, "\n")

meta_raw <- read_any(metadata_fp)
maj_raw <- read_any(majority_fp)
l1_raw <- read_any(level1_fp)
med139_metrics_raw <- read_any(med139_metrics_fp)

if (!"weighted_F1" %in% colnames(med139_metrics_raw)) {
  stop("Med139 metrics file does not contain weighted_F1.", call. = FALSE)
}
med139_wf1 <- as.numeric(med139_metrics_raw$weighted_F1[1])

meta_sample_col <- pick_col(meta_raw, c("sample_id", "Sample", "sample", "patient_id", "patient"), label = "metadata sample column")
meta_truth_col <- pick_col(
  meta_raw,
  c("disease_broad_external", "disease_broad", "broad5_truth", "disease_group", "disease", "condition", "label"),
  label = "metadata broad truth column"
)

maj_sample_col <- pick_col(maj_raw, c("sample_id", "Sample", "sample", "patient_id", "patient"), label = "majority sample column")
maj_pred_col <- pick_col(
  maj_raw,
  c("broad5_prediction", "prediction", "predicted_category", "pred_label", "pred_disease"),
  label = "majority prediction column"
)

l1_sample_col <- pick_col(l1_raw, c("sample_id", "Sample", "sample", "patient_id", "patient"), label = "Level1 sample column")
l1_cell_col <- pick_col(l1_raw, c("Level1", "classifier_ct", "celltype", "Celltype", "cell_type", "Cell_type"), label = "Level1 cell column")
l1_pred_col <- pick_col(l1_raw, c("pred_label", "pred_disease", "prediction", "predicted_category"), label = "Level1 prediction column")

meta_df <- meta_raw %>%
  dplyr::transmute(
    sample_id = trimws(as.character(.data[[meta_sample_col]])),
    true_short = short_broad(.data[[meta_truth_col]])
  ) %>%
  dplyr::distinct(sample_id, .keep_all = TRUE)

maj_df <- maj_raw %>%
  dplyr::transmute(
    sample_id = trimws(as.character(.data[[maj_sample_col]])),
    pred_short = short_broad(.data[[maj_pred_col]])
  ) %>%
  dplyr::left_join(meta_df, by = "sample_id") %>%
  dplyr::filter(!is.na(true_short), !is.na(pred_short), true_short != "Not assigned")

if (!include_healthy) {
  maj_df <- maj_df %>%
    dplyr::filter(true_short != "Healthy")
}

pred_df <- l1_raw %>%
  dplyr::transmute(
    sample_id = trimws(as.character(.data[[l1_sample_col]])),
    Level1_plot = short_celltype(.data[[l1_cell_col]]),
    pred_short = short_broad(.data[[l1_pred_col]])
  ) %>%
  dplyr::left_join(meta_df, by = "sample_id") %>%
  dplyr::filter(!is.na(true_short), !is.na(pred_short), true_short != "Not assigned") %>%
  dplyr::filter(sample_id %in% maj_df$sample_id)

if (!include_healthy) {
  pred_df <- pred_df %>%
    dplyr::filter(true_short != "Healthy")
}

truth_labels <- setdiff(unique(maj_df$true_short), "Not assigned")
disease_order <- c(
  broad_order[broad_order %in% truth_labels],
  setdiff(sort(truth_labels), broad_order)
)

cat("[INFO] n majority samples used: ", nrow(maj_df), "\n")
cat("[INFO] truth labels: ", paste(disease_order, collapse = " | "), "\n")
cat("[INFO] predicted labels: ", paste(sort(unique(maj_df$pred_short)), collapse = " | "), "\n")

# -------------------------------
# 5. Recompute and save evaluation files
# -------------------------------

maj_metrics <- calc_multiclass_metrics(
  y_true = maj_df$true_short,
  y_pred = maj_df$pred_short,
  labels = disease_order
)

metrics_df <- maj_metrics$metrics %>%
  dplyr::mutate(
    cohort_id = cohort_id,
    task = task_name,
    model = "pFever_Immune_withLinearSVC_majority",
    .before = 1
  )

per_class_majority <- maj_metrics$per_class %>%
  dplyr::mutate(
    cohort_id = cohort_id,
    task = task_name,
    model = "pFever_Immune_withLinearSVC_majority",
    .before = 1
  )

level1_metrics_df <- dplyr::bind_rows(lapply(sort(unique(pred_df$Level1_plot)), function(ct) {
  sub <- pred_df %>% dplyr::filter(Level1_plot == ct)
  met <- calc_multiclass_metrics(
    y_true = sub$true_short,
    y_pred = sub$pred_short,
    labels = disease_order
  )$metrics

  met %>%
    dplyr::mutate(
      cohort_id = cohort_id,
      task = task_name,
      Level1 = ct,
      .before = 1
    )
}))

data.table::fwrite(metrics_df, file.path(eval_out_dir, "broad5_external_metrics.tsv"), sep = "\t")
data.table::fwrite(level1_metrics_df, file.path(eval_out_dir, "broad5_level1_external_metrics.tsv"), sep = "\t")
data.table::fwrite(per_class_majority, file.path(eval_out_dir, "broad5_external_per_class_metrics.tsv"), sep = "\t")
data.table::fwrite(maj_df, file.path(eval_out_dir, "broad5_merged_predictions_with_truth.tsv"), sep = "\t")
data.table::fwrite(pred_df, file.path(eval_out_dir, "broad5_level1_merged_predictions_with_truth.tsv"), sep = "\t")

# -------------------------------
# 6. Panel a: WF1 dot plot with random baseline
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

random_df <- data.frame(
  Celltype = "Random",
  WF1 = mean(random_metrics$weighted_F1, na.rm = TRUE),
  ymin = mean(random_metrics$weighted_F1, na.rm = TRUE) - stats::sd(random_metrics$weighted_F1, na.rm = TRUE),
  ymax = mean(random_metrics$weighted_F1, na.rm = TRUE) + stats::sd(random_metrics$weighted_F1, na.rm = TRUE),
  stringsAsFactors = FALSE
)

majority_df <- data.frame(
  Celltype = "Majority vote",
  WF1 = as.numeric(metrics_df$weighted_F1[1]),
  ymin = NA_real_,
  ymax = NA_real_,
  stringsAsFactors = FALSE
)

med139_df <- data.frame(
  Celltype = "Published signature",
  WF1 = med139_wf1,
  ymin = NA_real_,
  ymax = NA_real_,
  stringsAsFactors = FALSE
)

level1_df_a <- level1_metrics_df %>%
  dplyr::transmute(
    Celltype = short_celltype(.data$Level1),
    WF1 = as.numeric(.data$weighted_F1),
    ymin = NA_real_,
    ymax = NA_real_
  )

celltype_order_existing <- celltype_order[
  celltype_order %in% c(
    "Majority vote", "Published signature", "Random",
    unique(level1_df_a$Celltype), unique(pred_df$Level1_plot)
  )
]

a_df <- dplyr::bind_rows(majority_df, med139_df, random_df, level1_df_a) %>%
  dplyr::filter(!is.na(WF1)) %>%
  dplyr::mutate(
    Celltype = factor(Celltype, levels = celltype_order_existing)
  )

data.table::fwrite(
  a_df,
  file.path(out_dir, "SourceData_PatientDisjointBroad5Fig_a_WF1_random_level1.tsv"),
  sep = "\t"
)

data.table::fwrite(
  random_metrics,
  file.path(out_dir, "SourceData_PatientDisjointBroad5Fig_a_random_100_iterations.tsv"),
  sep = "\t"
)

# Source data for updated supplementary figure
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
      Celltype = short_celltype(.data$Level1),
      Metric = "Balanced accuracy score",
      Score = as.numeric(.data$balanced_accuracy),
      ymin = NA_real_,
      ymax = NA_real_
    ),
  level1_metrics_df %>%
    dplyr::transmute(
      Celltype = short_celltype(.data$Level1),
      Metric = "Matthew correlation coefficient",
      Score = as.numeric(.data$MCC),
      ymin = NA_real_,
      ymax = NA_real_
    )
) %>%
  dplyr::mutate(Celltype = factor(Celltype, levels = celltype_order_existing))

data.table::fwrite(
  supp_metric_source,
  file.path(out_dir, "SourceData_PatientDisjointSuppFig_a_BAS_MCC_random_level1.tsv"),
  sep = "\t"
)

p_a <- ggplot2::ggplot(a_df, ggplot2::aes(x = Celltype, y = WF1, color = Celltype)) +
  ggplot2::geom_errorbar(
    data = dplyr::filter(a_df, as.character(Celltype) == "Random"),
    ggplot2::aes(ymin = ymin, ymax = ymax),
    width = 0.18,
    linewidth = 0.45,
    color = "#BDBDBD"
  ) +
  ggplot2::geom_vline(xintercept = 3.5, color = "black", linewidth = 0.35) +
  ggplot2::geom_point(size = 2.6) +
  ggplot2::scale_color_manual(values = complete_named_cols(levels(a_df$Celltype), celltype_cols), drop = FALSE) +
  ggplot2::scale_x_discrete(labels = celltype_labels) +
  ggplot2::coord_cartesian(
    ylim = c(max(0, min(a_df$WF1, na.rm = TRUE) - 0.12), min(1, max(a_df$WF1, na.rm = TRUE) + 0.12))
  ) +
  ggplot2::scale_y_continuous(
    breaks = seq(0, 1, 0.1),
    expand = ggplot2::expansion(mult = c(0.02, 0.04))
  ) +
  paper_theme(base_size = 10) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)
  ) +
  ggplot2::labs(
    title = "a",
    subtitle = paste0("Patient-disjoint set, n = ", nrow(maj_df)),
    x = NULL,
    y = "WF1 score"
  )

save_gg_pdf(
  file.path(out_dir, "PatientDisjointBroad5Fig_a_WF1_random_level1_withLinearSVC.pdf"),
  p_a,
  width = 4.6,
  height = 3.2
)

ggplot2::ggsave(
  file.path(out_dir, "PatientDisjointBroad5Fig_a_WF1_random_level1_withLinearSVC.png"),
  p_a,
  width = 4.6,
  height = 3.2,
  dpi = 300
)

# -------------------------------
# 7. Panel b: disease by Level 1 F1 heatmap
# -------------------------------

celltype_order_heatmap_existing <- celltype_order_heatmap[
  celltype_order_heatmap %in% unique(pred_df$Level1_plot)
]

b_rows <- list()

for (ct in celltype_order_heatmap_existing) {
  sub_ct <- dplyr::filter(pred_df, as.character(.data$Level1_plot) == ct)
  for (d in disease_order) {
    b_rows[[length(b_rows) + 1]] <- data.frame(
      Disease = d,
      Celltype = ct,
      F1 = compute_f1_bin(sub_ct$true_short, sub_ct$pred_short, d),
      stringsAsFactors = FALSE
    )
  }
}

b_long <- dplyr::bind_rows(b_rows) %>%
  dplyr::mutate(
    Disease = factor(Disease, levels = disease_order),
    Celltype = factor(Celltype, levels = celltype_order_heatmap_existing)
  )

data.table::fwrite(
  b_long,
  file.path(out_dir, "SourceData_PatientDisjointBroad5Fig_b_disease_celltype_F1.tsv"),
  sep = "\t"
)

b_mat <- b_long %>%
  dplyr::select(Disease, Celltype, F1) %>%
  tidyr::pivot_wider(names_from = Celltype, values_from = F1) %>%
  as.data.frame()

rownames(b_mat) <- as.character(b_mat$Disease)
b_mat$Disease <- NULL
b_mat <- as.matrix(b_mat)
b_mat[is.na(b_mat)] <- 0

ha_bottom <- ComplexHeatmap::HeatmapAnnotation(
  `Cell types` = factor(colnames(b_mat), levels = celltype_order_heatmap_existing),
  col = list(`Cell types` = complete_named_cols(colnames(b_mat), celltype_cols)),
  show_annotation_name = FALSE,
  show_legend = FALSE,
  which = "column",
  simple_anno_size = grid::unit(3.2, "mm"),
  border = FALSE
)

row_group <- droplevels(label_to_group(rownames(b_mat)))

ra_left <- ComplexHeatmap::rowAnnotation(
  `Disease group` = row_group,
  col = list(`Disease group` = broad_group_cols[levels(row_group)]),
  show_annotation_name = FALSE,
  show_legend = FALSE,
  simple_anno_size = grid::unit(3.2, "mm"),
  border = FALSE
)

ht_b <- ComplexHeatmap::Heatmap(
  b_mat,
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
  column_labels = unname(celltype_labels[colnames(b_mat)]),
  column_names_rot = 90,
  column_names_gp = grid::gpar(fontsize = 8, fontfamily = base_family),
  row_names_gp = grid::gpar(fontsize = 8.5, fontfamily = base_family),
  rect_gp = grid::gpar(col = NA),
  na_col = "#000004",
  heatmap_legend_param = score_heatmap_legend_param("F1 score"),
  column_title = NULL,
  row_title = NULL,
  width = grid::unit(ncol(b_mat) * 4.8, "mm"),
  height = grid::unit(nrow(b_mat) * 6.0, "mm")
)

b_w <- max(5.0, ncol(b_mat) * 0.25 + 2.4)
b_h <- max(3.2, nrow(b_mat) * 0.40 + 1.8)

open_pdf(file.path(out_dir, "PatientDisjointBroad5Fig_b_disease_celltype_F1_heatmap_withLinearSVC.pdf"), width = b_w, height = b_h)
ComplexHeatmap::draw(ht_b, heatmap_legend_side = "top", annotation_legend_side = "bottom", newpage = TRUE)
dev.off()

png(file.path(out_dir, "PatientDisjointBroad5Fig_b_disease_celltype_F1_heatmap_withLinearSVC.png"), width = b_w, height = b_h, units = "in", res = 300)
ComplexHeatmap::draw(ht_b, heatmap_legend_side = "top", annotation_legend_side = "bottom", newpage = TRUE)
dev.off()

# -------------------------------
# 8. Panel c: column normalized confusion matrix
# -------------------------------

true_order <- disease_order

pred_order <- c(
  broad_order[broad_order %in% unique(maj_df$pred_short)],
  setdiff(sort(unique(maj_df$pred_short)), broad_order)
)

pred_order <- pred_order[pred_order %in% unique(maj_df$pred_short)]

cm <- table(
  factor(maj_df$pred_short, levels = pred_order),
  factor(maj_df$true_short, levels = true_order)
)

cm <- as.matrix(cm)
cm <- cm[rowSums(cm) > 0, colSums(cm) > 0, drop = FALSE]

cm_col <- sweep(cm, 2, colSums(cm), "/")
cm_col[is.na(cm_col)] <- 0

cm_count_df <- as.data.frame.matrix(cm)
cm_count_df <- cbind(Predicted_condition = rownames(cm_count_df), cm_count_df)
data.table::fwrite(
  cm_count_df,
  file.path(out_dir, "SourceData_PatientDisjointBroad5Fig_c_confusion_counts.tsv"),
  sep = "\t"
)

cm_col_df <- as.data.frame.matrix(cm_col)
cm_col_df <- cbind(Predicted_condition = rownames(cm_col_df), cm_col_df)
data.table::fwrite(
  cm_col_df,
  file.path(out_dir, "SourceData_PatientDisjointBroad5Fig_c_confusion_colnorm.tsv"),
  sep = "\t"
)

top_group <- droplevels(label_to_group(colnames(cm_col)))
right_group <- droplevels(label_to_group(rownames(cm_col)))

annotation_group_cols <- broad_group_cols[unique(c(levels(top_group), levels(right_group)))]

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

ht_c <- ComplexHeatmap::Heatmap(
  cm_col,
  name = "True condition\nproportion",
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
  column_names_rot = 45,
  column_names_gp = grid::gpar(fontsize = 8.5, fontfamily = base_family),
  row_names_gp = grid::gpar(fontsize = 8.5, fontfamily = base_family),
  rect_gp = grid::gpar(col = "#1B1B1B", lwd = 0.20),
  na_col = "#000004",
  heatmap_legend_param = score_heatmap_legend_param("True condition proportion"),
  column_title = NULL,
  row_title = NULL,
  width = grid::unit(ncol(cm_col) * 8.0, "mm"),
  height = grid::unit(nrow(cm_col) * 8.0, "mm")
)

c_w <- max(4.2, ncol(cm_col) * 0.42 + 2.3)
c_h <- max(4.1, nrow(cm_col) * 0.42 + 1.8)

open_pdf(file.path(out_dir, "PatientDisjointBroad5Fig_c_confusion_matrix_withLinearSVC.pdf"), width = c_w, height = c_h)
ComplexHeatmap::draw(ht_c, heatmap_legend_side = "top", annotation_legend_side = "bottom", merge_legends = TRUE, newpage = TRUE)
dev.off()

png(file.path(out_dir, "PatientDisjointBroad5Fig_c_confusion_matrix_withLinearSVC.png"), width = c_w, height = c_h, units = "in", res = 300)
ComplexHeatmap::draw(ht_c, heatmap_legend_side = "top", annotation_legend_side = "bottom", merge_legends = TRUE, newpage = TRUE)
dev.off()

# -------------------------------
# 9. Summary
# -------------------------------

summary_list <- list(
  cohort_id = cohort_id,
  task_name = task_name,
  run_dir = run_dir,
  metadata_fp = metadata_fp,
  majority_fp = majority_fp,
  level1_fp = level1_fp,
  med139_metrics_fp = med139_metrics_fp,
  med139_weighted_F1 = med139_wf1,
  include_healthy = include_healthy,
  eval_out_dir = eval_out_dir,
  out_dir = out_dir,
  n_majority_predictions = nrow(maj_df),
  diseases = disease_order,
  predicted_rows = rownames(cm_col),
  cell_types = celltype_order_heatmap_existing,
  metrics = as.list(metrics_df[1, ]),
  random_baseline = list(
    n_iterations = n_random,
    mean_WF1 = mean(random_metrics$weighted_F1, na.rm = TRUE),
    sd_WF1 = stats::sd(random_metrics$weighted_F1, na.rm = TRUE),
    mean_BAS = mean(random_metrics$balanced_accuracy, na.rm = TRUE),
    sd_BAS = stats::sd(random_metrics$balanced_accuracy, na.rm = TRUE),
    mean_MCC = mean(random_metrics$MCC, na.rm = TRUE),
    sd_MCC = stats::sd(random_metrics$MCC, na.rm = TRUE)
  ),
  note = "Patient-disjoint broad5 figures recomputed from frozen pFever TRACE predictions. Med139 is evaluated in the same 222 samples. Panel c is normalized within each true category."
)

jsonlite::write_json(
  summary_list,
  file.path(out_dir, "patient_disjoint_broad5_withLinearSVC_summary.json"),
  auto_unbox = TRUE,
  pretty = TRUE
)

cat("\n[DONE] Patient-disjoint broad5 main plots saved to:\n")
cat(out_dir, "\n")
cat("\n[UPDATED EVALUATION FILES]\n")
cat(eval_out_dir, "\n")
cat("\n[OVERALL METRICS]\n")
print(metrics_df)