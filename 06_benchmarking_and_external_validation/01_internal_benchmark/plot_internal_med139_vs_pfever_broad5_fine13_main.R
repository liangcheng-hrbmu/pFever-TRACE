suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(grid)
})

# =========================================================
# Internal matched CV benchmark main figure
# Broad 5 + Fine 13
#
# a. Broad 5 BAS and WF1 trajectory
# b. Fine 13 BAS and WF1 trajectory
# c. Per class F1 difference
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
  "08_bulk_classifier/04_internal_matched_cv_med139_vs_pfever_combined_plots"
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

model_label <- function(x) {
  x0 <- as.character(x)
  xl <- tolower(trimws(x0))
  dplyr::case_when(
    grepl("pfever", xl) ~ "pFever Immune",
    grepl("med139", xl) | grepl("med", xl) | grepl("signature", xl) ~ "Med139 signature",
    TRUE ~ x0
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

method_cols <- c(
  "Med139 signature" = "#B784A7",
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

read_task_data <- function(task_name, task_dir) {
  split_fp <- file.path(task_dir, "internal_med139_vs_pfever_split_metrics.tsv")
  per_class_fp <- file.path(task_dir, "internal_med139_vs_pfever_per_class_metrics.tsv")

  if (!file.exists(split_fp)) stop("Missing split metrics for ", task_name, ": ", split_fp)
  if (!file.exists(per_class_fp)) stop("Missing per class metrics for ", task_name, ": ", per_class_fp)

  split_df <- read_tsv(split_fp)
  per_class_df <- read_tsv(per_class_fp)

  require_cols(split_df, c("model", "split", "weighted_F1", "balanced_accuracy"), paste0(task_name, " split_df"))
  require_cols(per_class_df, c("model", "class", "F1"), paste0(task_name, " per_class_df"))

  split_df <- split_df %>%
    dplyr::mutate(
      Task = task_name,
      Method = factor(model_label(model), levels = c("Med139 signature", "pFever Immune")),
      split = as.character(split)
    ) %>%
    dplyr::filter(!is.na(Method)) %>%
    dplyr::group_by(Task, split, Method) %>%
    dplyr::summarise(
      weighted_F1 = mean(weighted_F1, na.rm = TRUE),
      balanced_accuracy = mean(balanced_accuracy, na.rm = TRUE),
      .groups = "drop"
    )

  per_class_df <- per_class_df %>%
    dplyr::mutate(
      Task = task_name,
      Method = factor(model_label(model), levels = c("Med139 signature", "pFever Immune")),
      Class_display = class_label(class),
      Group = get_group(Class_display)
    ) %>%
    dplyr::filter(!is.na(Method)) %>%
    dplyr::group_by(Task, Method, Class_display, Group) %>%
    dplyr::summarise(
      F1 = mean(F1, na.rm = TRUE),
      .groups = "drop"
    )

  list(split = split_df, per_class = per_class_df)
}

# -------------------------------
# 2. Read data
# -------------------------------

broad <- read_task_data("Broad 5", broad_dir)
fine <- read_task_data("Fine 13", fine_dir)

split_df <- dplyr::bind_rows(broad$split, fine$split)
per_class_df <- dplyr::bind_rows(broad$per_class, fine$per_class)

# =========================================================
# a and b. Trajectory plots
# =========================================================

make_trajectory_data <- function(task_name) {
  d <- split_df %>%
    dplyr::filter(Task == task_name) %>%
    dplyr::select(Task, split, Method, balanced_accuracy, weighted_F1) %>%
    tidyr::pivot_wider(
      names_from = Method,
      values_from = c(balanced_accuracy, weighted_F1),
      names_sep = "__"
    )

  required <- c(
    "balanced_accuracy__Med139 signature",
    "balanced_accuracy__pFever Immune",
    "weighted_F1__Med139 signature",
    "weighted_F1__pFever Immune"
  )

  require_cols(d, required, paste0("trajectory data for ", task_name))

  d %>%
    dplyr::rename(
      BAS_ref = .data[["balanced_accuracy__Med139 signature"]],
      BAS_pfever = .data[["balanced_accuracy__pFever Immune"]],
      WF1_ref = .data[["weighted_F1__Med139 signature"]],
      WF1_pfever = .data[["weighted_F1__pFever Immune"]]
    )
}

make_summary_data <- function(task_name) {
  split_df %>%
    dplyr::filter(Task == task_name) %>%
    dplyr::group_by(Method) %>%
    dplyr::summarise(
      mean_BAS = mean(balanced_accuracy, na.rm = TRUE),
      sd_BAS = stats::sd(balanced_accuracy, na.rm = TRUE),
      mean_WF1 = mean(weighted_F1, na.rm = TRUE),
      sd_WF1 = stats::sd(weighted_F1, na.rm = TRUE),
      .groups = "drop"
    )
}

make_limits <- function(traj, summary_df) {
  xs <- c(
    traj$BAS_ref,
    traj$BAS_pfever,
    summary_df$mean_BAS - summary_df$sd_BAS,
    summary_df$mean_BAS + summary_df$sd_BAS
  )

  ys <- c(
    traj$WF1_ref,
    traj$WF1_pfever,
    summary_df$mean_WF1 - summary_df$sd_WF1,
    summary_df$mean_WF1 + summary_df$sd_WF1
  )

  x_range <- range(xs, na.rm = TRUE)
  y_range <- range(ys, na.rm = TRUE)
  x_pad <- max(0.04, diff(x_range) * 0.20)
  y_pad <- max(0.04, diff(y_range) * 0.20)

  list(
    xlim = c(max(0, x_range[1] - x_pad), min(1, x_range[2] + x_pad)),
    ylim = c(max(0, y_range[1] - y_pad), min(1, y_range[2] + y_pad))
  )
}

plot_trajectory <- function(task_name, panel_letter) {
  traj <- make_trajectory_data(task_name)
  summary_df <- make_summary_data(task_name)
  lims <- make_limits(traj, summary_df)

  data.table::fwrite(
    traj,
    file.path(out_dir, paste0("SourceData_", panel_letter, "_", gsub(" ", "_", task_name), "_trajectory.tsv")),
    sep = "\t"
  )

  data.table::fwrite(
    summary_df,
    file.path(out_dir, paste0("SourceData_", panel_letter, "_", gsub(" ", "_", task_name), "_summary.tsv")),
    sep = "\t"
  )

  ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = traj,
      ggplot2::aes(x = BAS_ref, y = WF1_ref, xend = BAS_pfever, yend = WF1_pfever),
      color = "grey72",
      linewidth = 0.55,
      arrow = grid::arrow(length = grid::unit(2.2, "mm"), type = "closed")
    ) +
    ggplot2::geom_point(
      data = traj,
      ggplot2::aes(x = BAS_ref, y = WF1_ref),
      shape = 21,
      size = 2.4,
      stroke = 0.70,
      fill = "white",
      color = method_cols[["Med139 signature"]]
    ) +
    ggplot2::geom_point(
      data = traj,
      ggplot2::aes(x = BAS_pfever, y = WF1_pfever),
      shape = 21,
      size = 2.6,
      stroke = 0.45,
      fill = method_cols[["pFever Immune"]],
      color = method_cols[["pFever Immune"]]
    ) +
    ggplot2::geom_segment(
      data = summary_df,
      ggplot2::aes(
        x = pmax(0, mean_BAS - sd_BAS),
        xend = pmin(1, mean_BAS + sd_BAS),
        y = mean_WF1,
        yend = mean_WF1,
        color = Method
      ),
      linewidth = 0.60
    ) +
    ggplot2::geom_segment(
      data = summary_df,
      ggplot2::aes(
        x = mean_BAS,
        xend = mean_BAS,
        y = pmax(0, mean_WF1 - sd_WF1),
        yend = pmin(1, mean_WF1 + sd_WF1),
        color = Method
      ),
      linewidth = 0.60
    ) +
    ggplot2::geom_point(
      data = summary_df,
      ggplot2::aes(x = mean_BAS, y = mean_WF1, fill = Method, color = Method),
      shape = 21,
      size = 3.4,
      stroke = 0.65
    ) +
    ggplot2::scale_color_manual(values = method_cols, drop = FALSE) +
    ggplot2::scale_fill_manual(values = method_cols, drop = FALSE) +
    ggplot2::coord_cartesian(xlim = lims$xlim, ylim = lims$ylim) +
    ggplot2::scale_x_continuous(breaks = pretty(lims$xlim, n = 4)) +
    ggplot2::scale_y_continuous(breaks = pretty(lims$ylim, n = 4)) +
    paper_theme(base_size = 9.5) +
    ggplot2::theme(
      legend.position = "top",
      legend.title = ggplot2::element_blank(),
      legend.margin = ggplot2::margin(0, 0, 0, 0)
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(override.aes = list(size = 2.8)),
      fill = "none"
    ) +
    ggplot2::labs(
      title = panel_letter,
      subtitle = task_name,
      x = "BAS",
      y = "WF1 score"
    )
}

p_a <- plot_trajectory("Broad 5", "a")
save_gg(p_a, "Internal_Fig_a_Broad5_BAS_WF1_trajectory", width = 4.5, height = 4.5)

p_b <- plot_trajectory("Fine 13", "b")
save_gg(p_b, "Internal_Fig_b_Fine13_BAS_WF1_trajectory", width = 4.5, height = 4.5)

# =========================================================
# c. Per class F1 difference, lollipop style like example
# =========================================================

delta_df <- per_class_df %>%
  dplyr::filter(Method %in% c("Med139 signature", "pFever Immune")) %>%
  dplyr::select(Task, Class_display, Group, Method, F1) %>%
  tidyr::pivot_wider(names_from = Method, values_from = F1) %>%
  dplyr::mutate(
    delta_F1 = .data[["pFever Immune"]] - .data[["Med139 signature"]],
    Task = factor(Task, levels = c("Broad 5", "Fine 13"))
  )

# 你想要的显示顺序
broad_show_order <- c(
  "Viral infection",
  "Bacterial infection",
  "Malaria",
  "Hyperinflammatory disease",
  "Autoimmune disease"
)
broad_show_order <- broad_show_order[broad_show_order %in% delta_df$Class_display]

fine_show_order <- c(
  "Malaria",
  "IAV",
  "RSV",
  "Adenovirus",
  "Rhinovirus",
  "Meningococcus",
  "Pneumococcus",
  "S. aureus",
  "GAS",
  "E. coli",
  "Tuberculosis",
  "KD",
  "JIA"
)
fine_show_order <- fine_show_order[fine_show_order %in% delta_df$Class_display]

# 为了 facet 后还能保持各自顺序，拼一个 Task + Class 的联合因子
delta_df <- delta_df %>%
  dplyr::mutate(
    class_facet = dplyr::case_when(
      Task == "Broad 5" ~ paste0("Broad 5__", Class_display),
      Task == "Fine 13" ~ paste0("Fine 13__", Class_display),
      TRUE ~ Class_display
    )
  )

class_levels <- c(
  paste0("Broad 5__", rev(broad_show_order)),
  paste0("Fine 13__", rev(fine_show_order))
)

delta_df$class_facet <- factor(delta_df$class_facet, levels = class_levels)

# x 轴范围
x_min <- min(delta_df$delta_F1, na.rm = TRUE)
x_max <- max(delta_df$delta_F1, na.rm = TRUE)

left_lim  <- min(-0.05, floor((x_min - 0.01) * 20) / 20)
right_lim <- max(0.20, ceiling((x_max + 0.01) * 20) / 20)

group_breaks <- c(
  "Autoimmune disease",
  "Bacterial infection",
  "Hyperinflammatory disease",
  "Malaria",
  "Viral infection"
)
group_breaks <- group_breaks[group_breaks %in% unique(delta_df$Group)]

p_c <- ggplot2::ggplot(
  delta_df,
  ggplot2::aes(x = delta_F1, y = class_facet, color = Group)
) +
  ggplot2::geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.5,
    color = "grey45"
  ) +
  ggplot2::geom_segment(
    ggplot2::aes(x = 0, xend = delta_F1, yend = class_facet),
    linewidth = 0.9,
    lineend = "round",
    alpha = 0.95
  ) +
  ggplot2::geom_point(size = 3.2) +
  ggplot2::facet_grid(
    rows = ggplot2::vars(Task),
    scales = "free_y",
    space = "free_y"
  ) +
  ggplot2::scale_y_discrete(
    labels = function(x) sub("^.*__", "", x)
  ) +
  ggplot2::scale_color_manual(
    values = group_cols,
    breaks = group_breaks,
    drop = FALSE
  ) +
  ggplot2::coord_cartesian(
    xlim = c(left_lim, right_lim),
    clip = "off"
  ) +
  ggplot2::scale_x_continuous(
    breaks = pretty(c(left_lim, right_lim), n = 4)
  ) +
  ggplot2::labs(
    title = "c",
    x = "F1 difference, pFever TRACE minus reference signature",
    y = NULL,
    color = NULL
  ) +
  ggplot2::theme_classic(base_size = 12) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      face = "bold",
      size = 20,
      hjust = 0.5
    ),

    axis.title.x = ggplot2::element_text(size = 14, color = "black"),
    axis.text.x  = ggplot2::element_text(size = 12, color = "black"),
    axis.text.y  = ggplot2::element_text(size = 12, color = "black"),
    axis.title.y = ggplot2::element_blank(),

    legend.position = "top",
    legend.direction = "horizontal",
    legend.justification = "center",
    legend.title = ggplot2::element_blank(),
    legend.text = ggplot2::element_text(size = 11),
    legend.key.width = grid::unit(1.4, "lines"),
    legend.key.height = grid::unit(0.9, "lines"),
    legend.margin = ggplot2::margin(0, 0, 0, 0),
    legend.box.margin = ggplot2::margin(0, 0, 4, 0),

    panel.spacing.y = grid::unit(1.2, "lines"),

    strip.background = ggplot2::element_blank(),
    strip.text.y.right = ggplot2::element_text(
      angle = -90,
      size = 16,
      color = "black",
      face = "plain"
    ),

    plot.margin = ggplot2::margin(8, 12, 8, 8)
  ) +
  ggplot2::guides(
    color = ggplot2::guide_legend(
      nrow = 1,
      byrow = TRUE,
      override.aes = list(size = 4, linewidth = 1)
    )
  )

ggplot2::ggsave(
  filename = file.path(out_dir, "Internal_Fig_c_per_class_F1_difference_lollipop.pdf"),
  plot = p_c,
  width = 8.4,
  height = 9.2,
  device = grDevices::cairo_pdf,
  family = base_family
)

ggplot2::ggsave(
  filename = file.path(out_dir, "Internal_Fig_c_per_class_F1_difference_lollipop.png"),
  plot = p_c,
  width = 8.4,
  height = 9.2,
  dpi = 300
)