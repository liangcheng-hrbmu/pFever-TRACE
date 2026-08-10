suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(grid)
})

# =========================================================
# Patient-disjoint benchmark figure
# pFever Immune with LinearSVC versus Med139 ridge baseline
#
# This script:
#   1. Reads current pFever withLinearSVC external predictions
#   2. Reads current Med139 ridge multinomial baseline predictions
#   3. Uses Med139 imbalance only lambda.min as the reference prediction
#   4. Recomputes paired external metrics
#   5. Recomputes bootstrap delta and McNemar test
#   6. Draws final external comparison plots
# =========================================================

# -------------------------------
# 0. Paths
# -------------------------------

base_dir <- "data/analysis/09_bulk_deconvolution"
ext_root <- file.path(base_dir, "12_external_validation_multi_cohort_assets_v1")

# Evaluation set used in the final manuscript figure
evaluation_id <- "patient_disjoint_sensitivity"
expected_n <- 222

# The Med139 predictions were originally generated for the complete
# E-MTAB-12793 matrix and are restricted below to the same 222 samples
# through an inner join with the pFever patient-disjoint predictions.
med139_source_id <- "E-MTAB-12793"

# pFever evaluation generated for the patient-disjoint set
pfever_eval_dir <- file.path(
  ext_root,
  "results",
  evaluation_id,
  "evaluation_withLinearSVC"
)

pfever_pred_fp <- file.path(
  pfever_eval_dir,
  "broad5_merged_predictions_with_truth.tsv"
)

# New Med139 ridge multinomial baseline directory
med139_dir <- file.path(
  ext_root,
  "results",
  med139_source_id,
  "med139_ridge_multinomial_baseline"
)

# Use the all weight schemes file and explicitly select imbalance only lambda.min
med139_prediction_fp_manual <- file.path(
  med139_dir,
  "med139_external_predictions_all_weight_schemes.tsv"
)

med139_prediction_col_manual <- "med139_imbalance_only_lambda_min_prediction"

# New output directories
cmp_dir <- file.path(
  ext_root,
  "results",
  evaluation_id,
  "comparison_pfever_withLinearSVC_vs_med139_ridge"
)

out_dir <- file.path(
  ext_root,
  "plots",
  "patient_disjoint_222_pfever_withLinearSVC_vs_med139_ridge_final"
)

dir.create(cmp_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(pfever_pred_fp)) {
  stop(
    "Missing pFever withLinearSVC prediction file:\n",
    pfever_pred_fp,
    "\nGenerate the Med139 external-evaluation metrics before running this comparison."
  )
}

if (!dir.exists(med139_dir)) {
  stop("Missing Med139 baseline directory:\n", med139_dir)
}

if (!file.exists(med139_prediction_fp_manual)) {
  stop(
    "Missing manually specified Med139 prediction file:\n",
    med139_prediction_fp_manual
  )
}

# -------------------------------
# 1. Helpers
# -------------------------------

base_family <- "sans"

read_any <- function(path) {
  if (grepl("\\.csv$", path, ignore.case = TRUE)) {
    data.table::fread(path, data.table = FALSE)
  } else {
    data.table::fread(path, sep = "\t", data.table = FALSE)
  }
}

read_tsv <- function(path) {
  data.table::fread(path, sep = "\t", data.table = FALSE)
}

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

paper_theme <- function(base_size = 9.5) {
  ggplot2::theme_classic(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0, size = base_size + 5),
      plot.subtitle = ggplot2::element_text(size = base_size + 1, color = "black"),
      axis.title = ggplot2::element_text(size = base_size + 1, color = "black"),
      axis.text = ggplot2::element_text(size = base_size, color = "black"),
      axis.line = ggplot2::element_line(linewidth = 0.35, color = "black"),
      axis.ticks = ggplot2::element_line(linewidth = 0.30, color = "black"),
      legend.title = ggplot2::element_text(size = base_size),
      legend.text = ggplot2::element_text(size = base_size - 1),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(size = base_size + 1, color = "black"),
      panel.spacing = grid::unit(4, "mm"),
      plot.margin = ggplot2::margin(4, 4, 4, 4, "pt")
    )
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
      "\nCurrent columns are:\n",
      paste(cols, collapse = " | "),
      call. = FALSE
    )
  }

  NA_character_
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

  out[is.na(out) & grepl("influenza|\\biav\\b|respiratory syncytial|\\brsv\\b|adeno|adenovirus|rhinovirus|rhino|dengue|denv|sars|covid|hfmd|enterovirus", xl)] <- "Viral infection"
  out[is.na(out) & grepl("mening|pneumo|staph|aureus|streptococcal|streptococcus|\\bgas\\b|escherichia|e\\.coli|e coli|ecoli|tuberculosis|\\btb\\b|sepsis", xl)] <- "Bacterial infection"
  out[is.na(out) & grepl("kawasaki|\\bkd\\b|mis-c|misc|multisystem inflammatory", xl)] <- "Hyperinflammatory disease"
  out[is.na(out) & grepl("juvenile idiopathic|juvenile idopathic|\\bjia\\b|arthritis|systemic lupus|\\bsle\\b", xl)] <- "Autoimmune disease"
  out[is.na(out) & grepl("malaria|plasmodium", xl)] <- "Malaria"

  out[is.na(out)] <- x[is.na(out)]
  out
}

calc_metrics <- function(truth, pred, labels = NULL) {
  truth <- as.character(truth)
  pred <- as.character(pred)

  valid <- !is.na(truth) & !is.na(pred)
  truth <- truth[valid]
  pred <- pred[valid]

  if (is.null(labels)) {
    labels <- sort(unique(truth))
  }

  labels <- labels[labels %in% unique(truth)]

  per_class <- dplyr::bind_rows(lapply(labels, function(cls) {
    tp <- sum(truth == cls & pred == cls)
    fp <- sum(truth != cls & pred == cls)
    fn <- sum(truth == cls & pred != cls)
    support <- sum(truth == cls)

    precision <- ifelse(tp + fp == 0, 0, tp / (tp + fp))
    recall <- ifelse(tp + fn == 0, 0, tp / (tp + fn))
    f1 <- ifelse(
      precision + recall == 0,
      0,
      2 * precision * recall / (precision + recall)
    )

    data.frame(
      class = cls,
      precision = precision,
      recall = recall,
      F1 = f1,
      support = support,
      stringsAsFactors = FALSE
    )
  }))

  accuracy <- mean(truth == pred)
  weighted_F1 <- sum(per_class$F1 * per_class$support) / sum(per_class$support)
  macro_F1 <- mean(per_class$F1)
  balanced_accuracy <- mean(per_class$recall)

  all_labs <- sort(unique(c(truth, pred)))

  cm <- table(
    factor(truth, levels = all_labs),
    factor(pred, levels = all_labs)
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
      n_samples = length(truth),
      n_truth_classes = length(unique(truth)),
      n_pred_classes = length(unique(pred)),
      stringsAsFactors = FALSE
    ),
    per_class = per_class
  )
}

mcnemar_exact <- function(med139_correct, pfever_correct) {
  both_correct <- sum(med139_correct & pfever_correct)
  both_wrong <- sum(!med139_correct & !pfever_correct)
  pfever_correct_med139_wrong <- sum(!med139_correct & pfever_correct)
  pfever_wrong_med139_correct <- sum(med139_correct & !pfever_correct)

  discordant <- pfever_correct_med139_wrong + pfever_wrong_med139_correct

  if (discordant == 0) {
    p <- 1
  } else {
    p <- stats::binom.test(
      x = min(pfever_correct_med139_wrong, pfever_wrong_med139_correct),
      n = discordant,
      p = 0.5,
      alternative = "two.sided"
    )$p.value
  }

  data.frame(
    both_correct = both_correct,
    pFever_correct_med139_wrong = pfever_correct_med139_wrong,
    pFever_wrong_med139_correct = pfever_wrong_med139_correct,
    both_wrong = both_wrong,
    discordant_total = discordant,
    mcnemar_exact_p = p,
    stringsAsFactors = FALSE
  )
}

# -------------------------------
# 2. Colors and labels
# -------------------------------

method_cols <- c(
  "Med139 ridge baseline" = "#B784A7",
  "pFever Immune" = "#E3A722"
)

method_short_cols <- c(
  "Med139" = "#B784A7",
  "pFever Immune" = "#E3A722"
)

group_cols <- c(
  "Viral infection" = "#6699CC",
  "Bacterial infection" = "#C48EAD",
  "Malaria" = "#6A4C93",
  "Hyperinflammatory disease" = "#E3A722",
  "Autoimmune disease" = "#36A596",
  "Healthy" = "#8E8E8E",
  "Other" = "#BDBDBD"
)

class_order_all <- c(
  "Viral infection",
  "Bacterial infection",
  "Malaria",
  "Hyperinflammatory disease",
  "Autoimmune disease",
  "Healthy"
)

metric_order <- c(
  "Weighted F1",
  "Balanced accuracy",
  "MCC",
  "Macro F1",
  "Accuracy"
)

metric_map <- c(
  "weighted_F1" = "Weighted F1",
  "balanced_accuracy" = "Balanced accuracy",
  "MCC" = "MCC",
  "macro_F1" = "Macro F1",
  "accuracy" = "Accuracy"
)

# -------------------------------
# 3. Read pFever withLinearSVC predictions
# -------------------------------

pfever_raw <- read_tsv(pfever_pred_fp)

pfever_sample_col <- pick_col(
  pfever_raw,
  c("sample_id", "Sample", "sample", "patient_id", "patient"),
  label = "pFever sample column"
)

pfever_truth_col <- pick_col(
  pfever_raw,
  c("true_short", "truth", "true_disease", "disease_broad", "disease_broad_external", "disease", "label"),
  label = "pFever truth column"
)

pfever_pred_col <- pick_col(
  pfever_raw,
  c("pred_short", "broad5_prediction", "prediction", "predicted_category", "pred_label", "pred_disease"),
  label = "pFever prediction column"
)

pfever_df <- pfever_raw %>%
  dplyr::transmute(
    sample_id = trimws(as.character(.data[[pfever_sample_col]])),
    truth = short_broad(.data[[pfever_truth_col]]),
    pfever_prediction = short_broad(.data[[pfever_pred_col]])
  ) %>%
  dplyr::filter(
    !is.na(sample_id),
    !is.na(truth),
    !is.na(pfever_prediction),
    truth != "Not assigned",
    truth != "Healthy"
  ) %>%
  dplyr::distinct(sample_id, .keep_all = TRUE)

# -------------------------------
# 4. Read Med139 ridge baseline predictions
# -------------------------------

med139_prediction_fp <- med139_prediction_fp_manual

cat("[READ] pFever withLinearSVC predictions:\n", pfever_pred_fp, "\n")
cat("[READ] Med139 ridge baseline predictions:\n", med139_prediction_fp, "\n")
cat("[READ] Med139 prediction column:\n", med139_prediction_col_manual, "\n")

med139_raw <- read_any(med139_prediction_fp)

med139_sample_col <- pick_col(
  med139_raw,
  c("sample_id", "Sample", "sample", "patient_id", "patient"),
  label = "Med139 sample column"
)

if (!med139_prediction_col_manual %in% colnames(med139_raw)) {
  stop(
    "Cannot find manually selected Med139 prediction column:\n",
    med139_prediction_col_manual,
    "\nCurrent columns are:\n",
    paste(colnames(med139_raw), collapse = " | ")
  )
}

med139_pred_col <- med139_prediction_col_manual

med139_truth_col <- pick_col(
  med139_raw,
  c(
    "truth",
    "true_short",
    "true_disease",
    "disease_broad",
    "disease_broad_external",
    "disease",
    "label"
  ),
  required = FALSE,
  label = "Med139 truth column"
)

med139_df <- med139_raw %>%
  dplyr::transmute(
    sample_id = trimws(as.character(.data[[med139_sample_col]])),
    med139_prediction = short_broad(.data[[med139_pred_col]]),
    med139_truth = if (!is.na(med139_truth_col)) short_broad(.data[[med139_truth_col]]) else NA_character_
  ) %>%
  dplyr::filter(
    !is.na(sample_id),
    !is.na(med139_prediction)
  ) %>%
  dplyr::distinct(sample_id, .keep_all = TRUE)

# -------------------------------
# 5. Pair samples
# -------------------------------

paired <- pfever_df %>%
  dplyr::inner_join(med139_df, by = "sample_id") %>%
  dplyr::filter(
    !is.na(truth),
    !is.na(pfever_prediction),
    !is.na(med139_prediction),
    truth != "Healthy",
    truth != "Not assigned"
  )

if (nrow(paired) == 0) {
  stop("No matched samples between pFever withLinearSVC and Med139 ridge predictions.")
}

if (nrow(paired) != expected_n) {
  stop(
    "Expected ", expected_n,
    " paired patient-disjoint disease samples, but observed ",
    nrow(paired),
    ". Check the pFever evaluation file and sample IDs before plotting."
  )
}

if (!all(is.na(paired$med139_truth))) {
  truth_disagree <- paired %>%
    dplyr::filter(!is.na(med139_truth), med139_truth != truth)

  if (nrow(truth_disagree) > 0) {
    warning(
      "There are ", nrow(truth_disagree),
      " samples with different truth labels between pFever file and Med139 file. ",
      "Using truth from pFever withLinearSVC evaluation file."
    )
  }
}

labels <- class_order_all[class_order_all %in% unique(paired$truth)]

cat("\n[INFO] matched paired samples: ", nrow(paired), "\n")
cat("[INFO] truth labels: ", paste(labels, collapse = " | "), "\n")
cat("[INFO] pFever prediction labels: ", paste(sort(unique(paired$pfever_prediction)), collapse = " | "), "\n")
cat("[INFO] Med139 prediction labels: ", paste(sort(unique(paired$med139_prediction)), collapse = " | "), "\n")

data.table::fwrite(
  paired,
  file.path(cmp_dir, "external_paired_predictions_pfever_withLinearSVC_vs_med139.tsv"),
  sep = "\t"
)

# -------------------------------
# 6. Recompute metrics
# -------------------------------

med139_met <- calc_metrics(
  truth = paired$truth,
  pred = paired$med139_prediction,
  labels = labels
)

pfever_met <- calc_metrics(
  truth = paired$truth,
  pred = paired$pfever_prediction,
  labels = labels
)

overall <- dplyr::bind_rows(
  med139_met$metrics %>%
    dplyr::mutate(model = "Med139 ridge baseline", .before = 1),
  pfever_met$metrics %>%
    dplyr::mutate(model = "pFever Immune", .before = 1)
)

per_class <- dplyr::bind_rows(
  med139_met$per_class %>%
    dplyr::mutate(model = "Med139 ridge baseline", .before = 1),
  pfever_met$per_class %>%
    dplyr::mutate(model = "pFever Immune", .before = 1)
)

mc <- mcnemar_exact(
  med139_correct = paired$med139_prediction == paired$truth,
  pfever_correct = paired$pfever_prediction == paired$truth
)

data.table::fwrite(
  overall,
  file.path(cmp_dir, "external_metrics_pfever_withLinearSVC_vs_med139.tsv"),
  sep = "\t"
)

data.table::fwrite(
  per_class,
  file.path(cmp_dir, "per_class_metrics_pfever_withLinearSVC_vs_med139.tsv"),
  sep = "\t"
)

data.table::fwrite(
  mc,
  file.path(cmp_dir, "mcnemar_pfever_withLinearSVC_vs_med139.tsv"),
  sep = "\t"
)

# -------------------------------
# 7. Bootstrap delta
# -------------------------------

set.seed(2026)
n_boot <- 2000

metrics_to_boot <- c("weighted_F1", "balanced_accuracy", "MCC", "macro_F1", "accuracy")

boot_rows <- vector("list", n_boot)

for (i in seq_len(n_boot)) {
  idx <- sample(seq_len(nrow(paired)), size = nrow(paired), replace = TRUE)
  sub <- paired[idx, , drop = FALSE]

  med139_i <- calc_metrics(sub$truth, sub$med139_prediction, labels)$metrics
  pfever_i <- calc_metrics(sub$truth, sub$pfever_prediction, labels)$metrics

  med139_vec <- as.numeric(unlist(med139_i[1, metrics_to_boot], use.names = FALSE))
  pfever_vec <- as.numeric(unlist(pfever_i[1, metrics_to_boot], use.names = FALSE))

  boot_rows[[i]] <- data.frame(
    iteration = i,
    metric = metrics_to_boot,
    delta = pfever_vec - med139_vec,
    stringsAsFactors = FALSE
  )
}

boot_all <- dplyr::bind_rows(boot_rows)

boot <- boot_all %>%
  dplyr::group_by(metric) %>%
  dplyr::summarise(
    delta_median = stats::median(delta, na.rm = TRUE),
    delta_mean = mean(delta, na.rm = TRUE),
    delta_ci_low = stats::quantile(delta, 0.025, na.rm = TRUE),
    delta_ci_high = stats::quantile(delta, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

data.table::fwrite(
  boot_all,
  file.path(cmp_dir, "bootstrap_delta_all_pfever_withLinearSVC_minus_med139.tsv"),
  sep = "\t"
)

data.table::fwrite(
  boot,
  file.path(cmp_dir, "bootstrap_delta_summary_pfever_withLinearSVC_minus_med139.tsv"),
  sep = "\t"
)

# -------------------------------
# 8. Panel a: overall metrics
# -------------------------------

overall_plot <- overall %>%
  dplyr::filter(model %in% c("Med139 ridge baseline", "pFever Immune")) %>%
  dplyr::mutate(
    Method = dplyr::case_when(
      model == "Med139 ridge baseline" ~ "Med139",
      model == "pFever Immune" ~ "pFever Immune",
      TRUE ~ model
    )
  ) %>%
  dplyr::select(
    Method,
    `Weighted F1` = weighted_F1,
    `Balanced accuracy` = balanced_accuracy,
    MCC = MCC,
    `Macro F1` = macro_F1,
    Accuracy = accuracy
  ) %>%
  tidyr::pivot_longer(
    cols = c(`Weighted F1`, `Balanced accuracy`, MCC, `Macro F1`, Accuracy),
    names_to = "Metric",
    values_to = "Score"
  ) %>%
  dplyr::mutate(
    Method = factor(Method, levels = c("Med139", "pFever Immune")),
    Metric = factor(Metric, levels = metric_order)
  )

data.table::fwrite(
  overall_plot,
  file.path(out_dir, "SourceData_a_overall_metrics.tsv"),
  sep = "\t"
)

p_a <- ggplot2::ggplot(overall_plot, ggplot2::aes(x = Metric, y = Score, fill = Method)) +
  ggplot2::geom_col(
    position = ggplot2::position_dodge(width = 0.68),
    width = 0.58,
    color = "black",
    linewidth = 0.22
  ) +
  ggplot2::scale_fill_manual(values = method_short_cols, drop = FALSE) +
  ggplot2::coord_cartesian(ylim = c(0, 1)) +
  ggplot2::scale_y_continuous(
    breaks = seq(0, 1, 0.2),
    expand = ggplot2::expansion(mult = c(0, 0.04))
  ) +
  paper_theme(base_size = 9.4) +
  ggplot2::theme(
    legend.position = "top",
    legend.title = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(angle = 35, hjust = 1)
  ) +
  ggplot2::labs(
    title = "a",
    subtitle = paste0("Patient-disjoint evaluation, n = ", nrow(paired)),
    x = NULL,
    y = "Performance"
  )

save_gg(p_a, "PatientDisjoint_Fig_a_overall_metrics", width = 5.6, height = 4.2)

# -------------------------------
# 9. Per disease point plots
# -------------------------------

plot_per_class_points <- function(metric_col, panel_letter, stub, x_title) {
  class_order <- labels

  d0 <- per_class %>%
    dplyr::filter(model %in% c("Med139 ridge baseline", "pFever Immune")) %>%
    dplyr::mutate(
      Method = dplyr::case_when(
        model == "Med139 ridge baseline" ~ "Med139",
        model == "pFever Immune" ~ "pFever Immune",
        TRUE ~ model
      ),
      class = factor(class, levels = rev(class_order))
    ) %>%
    dplyr::select(class, Method, value = .data[[metric_col]], support)

  support_df <- d0 %>%
    dplyr::group_by(class) %>%
    dplyr::summarise(
      support = suppressWarnings(max(as.numeric(support), na.rm = TRUE)),
      .groups = "drop"
    )

  d <- d0 %>%
    dplyr::select(class, Method, value) %>%
    tidyr::pivot_wider(names_from = Method, values_from = value) %>%
    dplyr::left_join(support_df, by = "class") %>%
    dplyr::mutate(
      class_chr = as.character(class),
      label = paste0(class_chr, "  n=", as.integer(support)),
      group = class_chr
    )

  label_map <- d$label
  names(label_map) <- as.character(d$class)

  data.table::fwrite(
    d,
    file.path(out_dir, paste0("SourceData_", stub, ".tsv")),
    sep = "\t"
  )

  p <- ggplot2::ggplot(d, ggplot2::aes(y = class)) +
    ggplot2::geom_segment(
      ggplot2::aes(x = Med139, xend = `pFever Immune`, yend = class),
      color = "grey70",
      linewidth = 1.1
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x = Med139),
      shape = 21,
      size = 4.4,
      stroke = 0.85,
      fill = "white",
      color = "#4D4D4D"
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x = `pFever Immune`, fill = group),
      shape = 21,
      size = 4.5,
      stroke = 0.35,
      color = "black"
    ) +
    ggplot2::scale_fill_manual(values = group_cols, drop = FALSE) +
    ggplot2::scale_y_discrete(labels = label_map) +
    ggplot2::coord_cartesian(xlim = c(0, 1.02)) +
    ggplot2::scale_x_continuous(
      breaks = seq(0, 1, 0.2),
      expand = ggplot2::expansion(mult = c(0.01, 0.03))
    ) +
    paper_theme(base_size = 9.3) +
    ggplot2::theme(
      legend.position = "none",
      panel.grid.major.x = ggplot2::element_line(color = "#D9D9D9", linewidth = 0.35),
      panel.grid.minor.x = ggplot2::element_blank()
    ) +
    ggplot2::labs(
      title = panel_letter,
      subtitle = x_title,
      x = x_title,
      y = NULL
    )

  save_gg(p, stub, width = 5.4, height = 3.8)
  p
}

p_b <- plot_per_class_points(
  metric_col = "F1",
  panel_letter = "b",
  stub = "PatientDisjoint_Fig_b_per_disease_F1_points",
  x_title = "F1"
)

p_c <- plot_per_class_points(
  metric_col = "recall",
  panel_letter = "c",
  stub = "PatientDisjoint_Fig_c_per_disease_recall_points",
  x_title = "Recall"
)

p_c2 <- plot_per_class_points(
  metric_col = "precision",
  panel_letter = "c2",
  stub = "PatientDisjoint_Supp_per_disease_precision_points",
  x_title = "Precision"
)

# -------------------------------
# 10. Bootstrap delta forest plot
# -------------------------------

boot_plot <- boot %>%
  dplyr::filter(metric %in% c("weighted_F1", "balanced_accuracy", "MCC", "macro_F1", "accuracy")) %>%
  dplyr::mutate(
    Metric = dplyr::recode(metric, !!!metric_map),
    Metric = factor(Metric, levels = rev(metric_order))
  )

data.table::fwrite(
  boot_plot,
  file.path(out_dir, "SourceData_d_bootstrap_delta.tsv"),
  sep = "\t"
)

x_min <- min(boot_plot$delta_ci_low, 0, na.rm = TRUE)
x_max <- max(boot_plot$delta_ci_high, 0, na.rm = TRUE)
x_pad <- max(0.02, (x_max - x_min) * 0.10)

p_d <- ggplot2::ggplot(boot_plot, ggplot2::aes(y = Metric)) +
  ggplot2::annotate(
    "rect",
    xmin = 0,
    xmax = Inf,
    ymin = -Inf,
    ymax = Inf,
    fill = "#FFF5D6",
    alpha = 0.6
  ) +
  ggplot2::geom_vline(
    xintercept = 0,
    color = "grey35",
    linewidth = 0.35,
    linetype = "dashed"
  ) +
  ggplot2::geom_segment(
    ggplot2::aes(
      x = delta_ci_low,
      xend = delta_ci_high,
      yend = Metric
    ),
    linewidth = 0.65,
    color = "black"
  ) +
  ggplot2::geom_point(
    ggplot2::aes(x = delta_median),
    shape = 21,
    fill = "#E3A722",
    color = "black",
    size = 3.8,
    stroke = 0.35
  ) +
  ggplot2::coord_cartesian(
    xlim = c(x_min - x_pad, x_max + x_pad),
    clip = "off"
  ) +
  paper_theme(base_size = 9.4) +
  ggplot2::theme(
    panel.grid.major.x = ggplot2::element_line(color = "#D9D9D9", linewidth = 0.35),
    panel.grid.minor.x = ggplot2::element_blank()
  ) +
  ggplot2::labs(
    title = "d",
    subtitle = "Bootstrap difference",
    x = "Performance difference\npFever Immune minus Med139",
    y = NULL
  )

save_gg(p_d, "PatientDisjoint_Fig_d_bootstrap_delta_forest", width = 5.3, height = 3.8)

# -------------------------------
# 11. Paired correctness
# -------------------------------

mc_long <- mc %>%
  tidyr::pivot_longer(
    cols = c(
      both_correct,
      pFever_correct_med139_wrong,
      pFever_wrong_med139_correct,
      both_wrong
    ),
    names_to = "Category",
    values_to = "n"
  ) %>%
  dplyr::mutate(
    Category = dplyr::recode(
      Category,
      "both_correct" = "Both correct",
      "pFever_correct_med139_wrong" = "pFever correct only",
      "pFever_wrong_med139_correct" = "Med139 correct only",
      "both_wrong" = "Both wrong"
    ),
    Category = factor(
      Category,
      levels = c(
        "Both correct",
        "pFever correct only",
        "Med139 correct only",
        "Both wrong"
      )
    )
  )

data.table::fwrite(
  mc_long,
  file.path(out_dir, "SourceData_e_paired_correctness.tsv"),
  sep = "\t"
)

p_e <- ggplot2::ggplot(mc_long, ggplot2::aes(x = Category, y = n, fill = Category)) +
  ggplot2::geom_col(color = "black", linewidth = 0.25, width = 0.62) +
  ggplot2::geom_text(
    ggplot2::aes(label = n),
    vjust = -0.35,
    size = 3.2,
    family = base_family
  ) +
  ggplot2::scale_fill_manual(
    values = c(
      "Both correct" = "#BDBDBD",
      "pFever correct only" = "#E3A722",
      "Med139 correct only" = "#B784A7",
      "Both wrong" = "#6E6E6E"
    ),
    drop = FALSE
  ) +
  ggplot2::coord_cartesian(
    ylim = c(0, max(mc_long$n, na.rm = TRUE) * 1.15)
  ) +
  paper_theme(base_size = 9.3) +
  ggplot2::theme(
    legend.position = "none",
    axis.text.x = ggplot2::element_text(angle = 35, hjust = 1)
  ) +
  ggplot2::labs(
    title = "e",
    subtitle = paste0("McNemar p = ", signif(mc$mcnemar_exact_p[1], 3)),
    x = NULL,
    y = "Number of samples"
  )

save_gg(p_e, "PatientDisjoint_Fig_e_paired_correctness", width = 4.8, height = 3.8)

# -------------------------------
# 12. Combined figure
# -------------------------------

if (requireNamespace("patchwork", quietly = TRUE)) {
  combined <- (p_a | p_d) / (p_b | p_c)

  ggplot2::ggsave(
    filename = file.path(out_dir, "PatientDisjoint_pfever_withLinearSVC_vs_med139_ridge_combined.pdf"),
    plot = combined,
    width = 10.8,
    height = 8.2,
    device = grDevices::cairo_pdf,
    family = base_family
  )

  ggplot2::ggsave(
    filename = file.path(out_dir, "PatientDisjoint_pfever_withLinearSVC_vs_med139_ridge_combined.png"),
    plot = combined,
    width = 10.8,
    height = 8.2,
    dpi = 300
  )
}

# -------------------------------
# 13. Summary
# -------------------------------

cat("\n[DONE] Patient-disjoint comparison pFever withLinearSVC vs Med139 ridge plots saved to:\n")
cat(out_dir, "\n")

cat("\n[UPDATED COMPARISON TABLES]\n")
cat(cmp_dir, "\n")

cat("\n[EVALUATION SET]\n")
cat("evaluation_id: ", evaluation_id, "\n")
cat("expected paired n: ", expected_n, "\n")
cat("observed paired n: ", nrow(paired), "\n")

cat("\n[INPUT FILES]\n")
cat("pFever: ", pfever_pred_fp, "\n")
cat("Med139: ", med139_prediction_fp, "\n")
cat("Med139 column: ", med139_prediction_col_manual, "\n")

cat("\n[Overall metrics]\n")
print(overall)

cat("\n[McNemar]\n")
print(mc)

cat("\n[Files]\n")
cat("  external_paired_predictions_pfever_withLinearSVC_vs_med139.tsv\n")
cat("  external_metrics_pfever_withLinearSVC_vs_med139.tsv\n")
cat("  per_class_metrics_pfever_withLinearSVC_vs_med139.tsv\n")
cat("  bootstrap_delta_summary_pfever_withLinearSVC_minus_med139.tsv\n")
cat("  mcnemar_pfever_withLinearSVC_vs_med139.tsv\n")
cat("  External_pfever_withLinearSVC_vs_med139_ridge_combined.pdf/png\n")
