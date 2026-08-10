suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(grid)
})

# =========================================================
# Internal matched CV benchmark supplementary figure
# Broad 5 + Fine 13
#
# a. Overall metric grouped bars
# b. Per class raw F1 point comparison
# c. Per class raw precision point comparison
# d. Per class raw recall point comparison
# =========================================================

# -------------------------------
# 0. Paths
# -------------------------------

base_dir <- "data/analysis/09_bulk_deconvolution"

broad_dir <- file.path(
  base_dir,
  "08_bulk_classifier/04_internal_matched_cv_med139_vs_pfever"
)

fine_dir <- file.path(
  base_dir,
  "08_bulk_classifier/04_internal_matched_cv_med139_vs_pfever_fine13"
)

out_dir <- file.path(
  base_dir,
  "08_bulk_classifier/04_internal_matched_cv_med139_vs_pfever_combined_plots_supp"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------------
# 1. Helpers
# -------------------------------

base_family <- "sans"

read_tsv <- function(path) {
  data.table::fread(path, sep = "\t", data.table = FALSE)
}

require_cols <- function(df, cols, name) {
  miss <- setdiff(cols, colnames(df))
  if (length(miss) > 0) {
    stop(
      name, " is missing required columns: ",
      paste(miss, collapse = ", "),
      "\nCurrent columns are: ",
      paste(colnames(df), collapse = " | "),
      call. = FALSE
    )
  }
  invisible(TRUE)
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

paper_theme <- function(base_size = 9.3) {
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

model_label <- function(x) {
  x0 <- as.character(x)
  xl <- tolower(trimws(x0))
  dplyr::case_when(
    grepl("pfever", xl) ~ "pFever Immune",
    grepl("med139", xl) | grepl("med", xl) | grepl("signature", xl) ~ "Med139 signature",
    TRUE ~ x0
  )
}

model_short <- function(x) {
  x <- model_label(x)
  dplyr::case_when(
    x == "Med139 signature" ~ "Med139",
    x == "pFever Immune" ~ "pFever Immune",
    TRUE ~ x
  )
}

class_label <- function(x) {
  x <- trimws(as.character(x))
  map <- c(
    "Viral infection" = "Viral infection",
    "Bacterial infection" = "Bacterial infection",
    "Malaria" = "Malaria",
    "Hyperinflammatory disease" = "Hyperinflammatory disease",
    "Autoimmune disease" = "Autoimmune disease",
    "IAV" = "IAV",
    "RSV" = "RSV",
    "Adenovirus" = "Adenovirus",
    "Rhinovirus" = "Rhinovirus",
    "Meningococcus" = "Meningococcus",
    "Pneumococcus" = "Pneumococcus",
    "S. aureus" = "S. aureus",
    "GAS" = "GAS",
    "E. coli" = "E. coli",
    "Tuberculosis" = "Tuberculosis",
    "KD" = "KD",
    "JIA" = "JIA"
  )
  out <- ifelse(x %in% names(map), map[x], x)
  unname(out)
}

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
  "Other" = "#BDBDBD",
  "Not assigned" = "#2F225B"
)

broad_order <- c(
  "Viral infection",
  "Bacterial infection",
  "Malaria",
  "Hyperinflammatory disease",
  "Autoimmune disease"
)

fine_order <- c(
  "IAV", "RSV", "Adenovirus", "Rhinovirus",
  "Meningococcus", "Pneumococcus", "S. aureus", "GAS", "E. coli", "Tuberculosis",
  "Malaria", "KD", "JIA"
)

group_map <- c(
  "Viral infection" = "Viral infection",
  "Bacterial infection" = "Bacterial infection",
  "Malaria" = "Malaria",
  "Hyperinflammatory disease" = "Hyperinflammatory disease",
  "Autoimmune disease" = "Autoimmune disease",
  "IAV" = "Viral infection",
  "RSV" = "Viral infection",
  "Adenovirus" = "Viral infection",
  "Rhinovirus" = "Viral infection",
  "Meningococcus" = "Bacterial infection",
  "Pneumococcus" = "Bacterial infection",
  "S. aureus" = "Bacterial infection",
  "GAS" = "Bacterial infection",
  "E. coli" = "Bacterial infection",
  "Tuberculosis" = "Bacterial infection",
  "KD" = "Hyperinflammatory disease",
  "JIA" = "Autoimmune disease",
  "Not assigned" = "Not assigned"
)

get_group <- function(x) {
  x <- as.character(x)
  out <- ifelse(x %in% names(group_map), group_map[x], "Other")
  factor(out, levels = names(group_cols))
}

detect_support_col <- function(df) {
  candidates <- c("support", "n", "n_samples", "class_n", "N", "count")
  hit <- candidates[candidates %in% colnames(df)]
  if (length(hit) == 0) return(NULL)
  hit[1]
}

read_task_data <- function(task_name, task_dir) {
  overall_fp <- file.path(task_dir, "internal_med139_vs_pfever_overall_metrics.tsv")
  per_class_fp <- file.path(task_dir, "internal_med139_vs_pfever_per_class_metrics.tsv")

  if (!file.exists(overall_fp)) stop("Missing overall metrics for ", task_name, ": ", overall_fp)
  if (!file.exists(per_class_fp)) stop("Missing per class metrics for ", task_name, ": ", per_class_fp)

  overall_df <- read_tsv(overall_fp)
  per_class_df <- read_tsv(per_class_fp)

  require_cols(
    overall_df,
    c("model", "weighted_F1", "balanced_accuracy", "MCC", "macro_F1", "accuracy"),
    paste0(task_name, " overall_df")
  )

  require_cols(
    per_class_df,
    c("model", "class", "F1", "precision", "recall"),
    paste0(task_name, " per_class_df")
  )

  support_col <- detect_support_col(per_class_df)

  overall_df <- overall_df %>%
    dplyr::mutate(
      Task = task_name,
      Method = factor(model_label(model), levels = c("Med139 signature", "pFever Immune")),
      Method_short = factor(model_short(model), levels = c("Med139", "pFever Immune"))
    ) %>%
    dplyr::filter(!is.na(Method))

  per_class_df <- per_class_df %>%
    dplyr::mutate(
      Task = task_name,
      Method = factor(model_label(model), levels = c("Med139 signature", "pFever Immune")),
      Method_short = factor(model_short(model), levels = c("Med139", "pFever Immune")),
      Class_display = class_label(class),
      Group = get_group(Class_display)
    ) %>%
    dplyr::filter(!is.na(Method))

  if (!is.null(support_col)) {
    per_class_df$n_support <- suppressWarnings(as.numeric(per_class_df[[support_col]]))
  } else {
    per_class_df$n_support <- NA_real_
  }

  list(overall = overall_df, per_class = per_class_df)
}

# -------------------------------
# 2. Read data
# -------------------------------

broad <- read_task_data("Broad 5", broad_dir)
fine <- read_task_data("Fine 13", fine_dir)

overall_df <- dplyr::bind_rows(broad$overall, fine$overall)
per_class_df <- dplyr::bind_rows(broad$per_class, fine$per_class)

# =========================================================
# a. Overall metric grouped bars
# =========================================================

overall_long <- overall_df %>%
  dplyr::select(
    Task,
    Method_short,
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
    Task = factor(Task, levels = c("Broad 5", "Fine 13")),
    Metric = factor(
      Metric,
      levels = c("Weighted F1", "Balanced accuracy", "MCC", "Macro F1", "Accuracy")
    )
  )

data.table::fwrite(overall_long, file.path(out_dir, "SourceData_a_overall_metrics.tsv"), sep = "\t")

p_a <- ggplot2::ggplot(overall_long, ggplot2::aes(x = Metric, y = Score, fill = Method_short)) +
  ggplot2::geom_col(
    position = ggplot2::position_dodge(width = 0.68),
    width = 0.58,
    color = "black",
    linewidth = 0.22
  ) +
  ggplot2::facet_grid(. ~ Task) +
  ggplot2::scale_fill_manual(values = method_short_cols, drop = FALSE) +
  ggplot2::coord_cartesian(ylim = c(0, 1.05)) +
  ggplot2::scale_y_continuous(
    breaks = seq(0, 1, 0.2),
    expand = ggplot2::expansion(mult = c(0, 0.04))
  ) +
  paper_theme(base_size = 9.3) +
  ggplot2::theme(
    legend.position = "top",
    legend.title = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(angle = 35, hjust = 1)
  ) +
  ggplot2::labs(
    title = "a",
    x = NULL,
    y = "Internal matched CV performance"
  )

save_gg(p_a, "Internal_Supp_a_overall_metric_grouped_bars", width = 7.2, height = 3.8)

# =========================================================
# b, c, d. Absolute class level point comparison
# =========================================================

plot_absolute_points <- function(metric_col, panel_letter, file_stub, x_title) {
  d0 <- per_class_df %>%
    dplyr::select(
      Task,
      Class_display,
      Group,
      Method_short,
      n_support,
      value = .data[[metric_col]]
    )

  support_df <- d0 %>%
    dplyr::group_by(Task, Class_display, Group) %>%
    dplyr::summarise(
      n_support = suppressWarnings(max(as.numeric(n_support), na.rm = TRUE)),
      .groups = "drop"
    )

  support_df$n_support[!is.finite(support_df$n_support)] <- NA_real_

  d <- d0 %>%
    dplyr::select(Task, Class_display, Group, Method_short, value) %>%
    dplyr::group_by(Task, Class_display, Group, Method_short) %>%
    dplyr::summarise(value = mean(value, na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = Method_short, values_from = value)

  require_cols(d, c("Med139", "pFever Immune"), paste0(metric_col, " point plot data"))

  d <- d %>%
    dplyr::left_join(support_df, by = c("Task", "Class_display", "Group")) %>%
    dplyr::mutate(
      Task = factor(Task, levels = c("Broad 5", "Fine 13")),
      Class_key = paste(as.character(Task), Class_display, sep = " | ")
    )

  class_key_order <- unique(c(
    paste("Broad 5", broad_order[broad_order %in% d$Class_display], sep = " | "),
    paste("Fine 13", fine_order[fine_order %in% d$Class_display], sep = " | ")
  ))

  d <- d %>%
    dplyr::mutate(
      Class_key = factor(Class_key, levels = rev(class_key_order)),
      Label = dplyr::if_else(
        is.na(n_support),
        as.character(Class_display),
        paste0(as.character(Class_display), "  n=", as.integer(n_support))
      )
    )

  label_map <- d$Label
  names(label_map) <- as.character(d$Class_key)

  data.table::fwrite(d, file.path(out_dir, paste0("SourceData_", file_stub, ".tsv")), sep = "\t")

  p <- ggplot2::ggplot(d, ggplot2::aes(y = Class_key)) +
    ggplot2::geom_segment(
      ggplot2::aes(x = Med139, xend = `pFever Immune`, yend = Class_key),
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
      ggplot2::aes(x = `pFever Immune`, fill = Group),
      shape = 21,
      size = 4.5,
      stroke = 0.35,
      color = "black"
    ) +
    ggplot2::facet_grid(Task ~ ., scales = "free_y", space = "free_y") +
    ggplot2::scale_fill_manual(values = group_cols, drop = FALSE) +
    ggplot2::scale_y_discrete(labels = label_map) +
    ggplot2::coord_cartesian(xlim = c(0, 1.02)) +
    ggplot2::scale_x_continuous(
      breaks = seq(0, 1, 0.2),
      expand = ggplot2::expansion(mult = c(0.01, 0.03))
    ) +
    paper_theme(base_size = 9.2) +
    ggplot2::theme(
      legend.position = "top",
      legend.title = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_line(color = "#D9D9D9", linewidth = 0.35),
      panel.grid.minor.x = ggplot2::element_blank()
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 3.8))
    ) +
    ggplot2::labs(
      title = panel_letter,
      x = x_title,
      y = NULL
    )

  save_gg(p, file_stub, width = 6.8, height = 6.3)
  p
}

p_b <- plot_absolute_points(
  metric_col = "F1",
  panel_letter = "b",
  file_stub = "Internal_Supp_b_absolute_F1_points",
  x_title = "F1"
)

p_c <- plot_absolute_points(
  metric_col = "precision",
  panel_letter = "c",
  file_stub = "Internal_Supp_c_absolute_precision_points",
  x_title = "Precision"
)

p_d <- plot_absolute_points(
  metric_col = "recall",
  panel_letter = "d",
  file_stub = "Internal_Supp_d_absolute_recall_points",
  x_title = "Recall"
)

if (requireNamespace("patchwork", quietly = TRUE)) {
  combined <- p_a / (p_b | p_c | p_d) +
    patchwork::plot_layout(heights = c(0.8, 1.35))

  ggplot2::ggsave(
    filename = file.path(out_dir, "Internal_Supp_combined_absolute_points.pdf"),
    plot = combined,
    width = 17,
    height = 9.5,
    device = grDevices::cairo_pdf,
    family = base_family
  )

  ggplot2::ggsave(
    filename = file.path(out_dir, "Internal_Supp_combined_absolute_points.png"),
    plot = combined,
    width = 17,
    height = 9.5,
    dpi = 300
  )
}

cat("\n[DONE] Combined Broad 5 and Fine 13 internal benchmark supplementary figure saved to:\n")
cat(out_dir, "\n")