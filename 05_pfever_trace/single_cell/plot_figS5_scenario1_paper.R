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
# Supplementary Fig. S5, Scenario 1
# Five fold cross validation performance details
#
# S5a: Majority vote per disease metrics across five splits
#      Balanced accuracy, F1, precision, recall
# S5b: Classifier performance across five splits
#      Balanced accuracy and MCC for majority vote and each Level 1 classifier
# S5c: Disease x cell type heatmaps after pooling all folds
#      Recall and precision
#
# Style target: use the same paper mimic style as Fig 5c-e.
# =========================================================

# -------------------------------
# 0. Input and output
# -------------------------------
base_dir <- "data/analysis/08_patient_classifier/03_patient_classifier_scenario1_v3"

split_summary_tsv  <- file.path(base_dir, "scenario1_split_summary.tsv")
level1_metrics_tsv <- file.path(base_dir, "scenario1_level1_metrics_all.tsv")
majority_pred_tsv  <- file.path(base_dir, "scenario1_majority_predictions_all.tsv")
level1_pred_tsv    <- file.path(base_dir, "query_level1_predictions_all.tsv")

out_dir <- "data/analysis/08_patient_classifier/05_supplementary_figS5_paper_mimic"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------------
# 1. Safe PDF device
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

# -------------------------------
# 2. Shared palettes and labels
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
  "Majority vote"   = "Majority\nvote",
  "B"               = "B",
  "Plasma"          = "Plasma",
  "pDC"             = "pDC",
  "DC"              = "DC",
  "Mono"            = "Mono",
  "T CD4 naive"     = "T CD4 naive",
  "T CD4 non-naive" = "T CD4 non naive",
  "T CD8 naive"     = "T CD8 naive",
  "T CD8 non-naive" = "T CD8 non naive",
  "UTC"             = "UTC",
  "ILC"             = "ILC"
)

preferred_disease_order <- c(
  "Healthy",
  "Bact.inf",
  "Dengue",
  "HFMD",
  "Flu",
  "COVID",
  "Sepsis",
  "KD",
  "MIS-C",
  "SLE",
  "JIA"
)

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
# 3. Helper functions
# -------------------------------
short_celltype <- function(x) {
  x <- as.character(x)
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

short_disease <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[is.na(x) | x == ""] <- "Not assigned"

  map <- c(
    "Pediatric healthy" = "Healthy",
    "Pediatric Sepsis" = "Sepsis",
    "severe acute respiratory syndrome coronavirus 2 (SARS-CoV-2)" = "COVID",
    "SARS-CoV-2" = "COVID",
    "COVID-19" = "COVID",
    "Multisystem inflammatory syndrome in children (MIS-C)" = "MIS-C",
    "Kawasaki disease (KD)" = "KD",
    "human systemic lupus erythematosus (SLE)" = "SLE",
    "systemic lupus erythematosus" = "SLE",
    "juvenile idopathic arthritis (JIA)" = "JIA",
    "juvenile idiopathic arthritis (JIA)" = "JIA",
    "HFMD(EV-A71)enterovirus" = "HFMD",
    "HFMD" = "HFMD",
    "Influenza A virus (IAV)" = "Flu",
    "IAV" = "Flu",
    "bacterial infection" = "Bact.inf",
    "Bacterial infection" = "Bact.inf",
    "Dengue human infection(DENV-1)" = "Dengue",
    "DENV" = "Dengue",
    "TIE" = "Not assigned",
    "NA" = "Not assigned",
    "nan" = "Not assigned",
    "Not assigned" = "Not assigned"
  )

  out <- ifelse(x %in% names(map), map[x], x)
  unname(out)
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

# -------------------------------
# 4. Read data
# -------------------------------
split_df  <- data.table::fread(split_summary_tsv, data.table = FALSE)
level1_df <- data.table::fread(level1_metrics_tsv, data.table = FALSE)
maj_df    <- data.table::fread(majority_pred_tsv, data.table = FALSE)

if (!file.exists(level1_pred_tsv)) {
  stop("Cannot find query_level1_predictions_all.tsv. Please merge split level predictions first.")
}
pred_df <- data.table::fread(level1_pred_tsv, data.table = FALSE)

if ("benchmark_status" %in% colnames(maj_df)) {
  maj_df <- dplyr::filter(maj_df, as.character(benchmark_status) == "primary_benchmark")
}
if ("benchmark_status" %in% colnames(pred_df)) {
  pred_df <- dplyr::filter(pred_df, as.character(benchmark_status) == "primary_benchmark")
}

level1_df$Level1 <- short_celltype(level1_df$Level1)
pred_df$Level1 <- short_celltype(pred_df$Level1)
maj_df$true_short <- short_disease(maj_df$true_disease)
maj_df$pred_short <- short_disease(maj_df$pred_disease_majority)
pred_df$true_short <- short_disease(pred_df$disease_Level_3)
pred_df$pred_short <- short_disease(pred_df$pred_disease)

truth_diseases <- setdiff(unique(maj_df$true_short), "Not assigned")
disease_order <- c(
  preferred_disease_order[preferred_disease_order %in% truth_diseases],
  setdiff(sort(truth_diseases), preferred_disease_order)
)

celltype_order_existing <- celltype_order[celltype_order %in% c("Majority vote", unique(level1_df$Level1), unique(pred_df$Level1))]

# -------------------------------
# 5. S5a Majority vote per disease across splits
# -------------------------------
if (!"split" %in% colnames(maj_df)) {
  stop("scenario1_majority_predictions_all.tsv must contain a split column for S5a.")
}

maj_split_metric_rows <- list()
for (sp in sort(unique(maj_df$split))) {
  sub_sp <- dplyr::filter(maj_df, split == sp)
  for (d in disease_order) {
    met <- compute_one_vs_rest(sub_sp$true_short, sub_sp$pred_short, d)
    maj_split_metric_rows[[length(maj_split_metric_rows) + 1]] <- data.frame(
      split = as.character(sp),
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

maj_split_metric <- dplyr::bind_rows(maj_split_metric_rows) %>%
  dplyr::mutate(
    disease = factor(disease, levels = disease_order)
  )

data.table::fwrite(maj_split_metric, file.path(out_dir, "SourceData_S5a_majority_vote_per_disease_by_split.tsv"), sep = "\t")

avg_support_df <- maj_split_metric %>%
  dplyr::group_by(disease) %>%
  dplyr::summarise(avg_support = mean(support, na.rm = TRUE), .groups = "drop") %>%
  dplyr::mutate(
    disease = factor(disease, levels = disease_order),
    avg_support_lab = sprintf("%.1f", avg_support)
  )

metric_long_a <- maj_split_metric %>%
  dplyr::select(split, disease, support, BalancedAccuracy, F1, Precision, Recall) %>%
  tidyr::pivot_longer(
    cols = c(BalancedAccuracy, F1, Precision, Recall),
    names_to = "Metric",
    values_to = "Score"
  ) %>%
  dplyr::mutate(
    Metric = factor(Metric, levels = c("BalancedAccuracy", "F1", "Precision", "Recall"),
                    labels = c("Balanced accuracy", "F1", "Precision", "Recall")),
    disease = factor(disease, levels = disease_order)
  )

# top labels for average sample count per disease
label_df_a <- tidyr::crossing(
  disease = factor(disease_order, levels = disease_order),
  Metric = levels(metric_long_a$Metric)
) %>%
  dplyr::left_join(avg_support_df, by = "disease") %>%
  dplyr::mutate(y = 1.05)

# custom disease fills for a
# reassign a stable palette without changing the main figure style
supp_disease_cols <- c(
  "Healthy" = "#B4B4B4",
  "Bact.inf" = "#E0866A",
  "Dengue" = "#94B37A",
  "HFMD" = "#6F9A72",
  "Flu" = "#6AA7A2",
  "COVID" = "#D1BA5A",
  "Sepsis" = "#DD4A61",
  "KD" = "#9D6EA8",
  "MIS-C" = "#A6C6DD",
  "SLE" = "#4B6B8F",
  "JIA" = "#C39A3C"
)
supp_disease_cols <- c(supp_disease_cols,
                       stats::setNames(rep("#BDBDBD", length(setdiff(disease_order, names(supp_disease_cols)))),
                                       setdiff(disease_order, names(supp_disease_cols))))

p_s5a <- ggplot(metric_long_a, aes(x = disease, y = Score, fill = disease)) +
  stat_boxplot(geom = "errorbar", width = 0.34, linewidth = 0.35, color = "black") +
  geom_boxplot(
    width = 0.65,
    linewidth = 0.35,
    color = "black",
    outlier.shape = 21,
    outlier.size = 1.15,
    outlier.stroke = 0.25,
    outlier.fill = "white",
    outlier.colour = "black"
  ) +
  stat_summary(
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
  geom_text(
    data = label_df_a,
    aes(x = disease, y = y, label = avg_support_lab),
    inherit.aes = FALSE,
    size = 2.7,
    family = base_family,
    vjust = 0
  ) +
  annotate("text", x = 0.25, y = 1.05, label = "average #\nof samples", hjust = 0,
           vjust = 0, size = 3.0, family = base_family) +
  facet_wrap(~Metric, nrow = 2, scales = "fixed") +
  scale_fill_manual(values = supp_disease_cols[disease_order], guide = "none") +
  coord_cartesian(ylim = c(-0.05, 1.08), clip = "off") +
  scale_y_continuous(breaks = seq(0, 1, 0.2), expand = expansion(mult = c(0.02, 0.02))) +
  paper_theme(base_size = 9.5) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    strip.text = element_text(size = 10),
    panel.spacing = unit(4, "mm")
  ) +
  labs(title = "a", x = NULL, y = "score")

save_gg_pdf(file.path(out_dir, "SuppFigS5a_majority_per_disease_metrics_boxplots.pdf"), p_s5a, width = 8.6, height = 7.6)
ggsave(file.path(out_dir, "SuppFigS5a_majority_per_disease_metrics_boxplots.png"), p_s5a, width = 8.6, height = 7.6, dpi = 300)

# -------------------------------
# 6. S5b Majority vote and Level 1 classifier metrics
# -------------------------------
metric_df_b <- bind_rows(
  split_df %>%
    transmute(
      split = as.character(split),
      Celltype = "Majority vote",
      BAS = as.numeric(BAS_majority),
      MCC = as.numeric(MCC_majority)
    ),
  level1_df %>%
    transmute(
      split = as.character(split),
      Celltype = short_celltype(Level1),
      BAS = as.numeric(BAS),
      MCC = as.numeric(MCC)
    )
) %>%
  mutate(
    Celltype = factor(Celltype, levels = celltype_order_existing)
  ) %>%
  pivot_longer(
    cols = c(BAS, MCC),
    names_to = "Metric",
    values_to = "Score"
  ) %>%
  filter(!is.na(Score)) %>%
  mutate(
    Metric = factor(Metric, levels = c("BAS", "MCC"),
                    labels = c("Balanced accuracy score", "Matthew correlation coefficient"))
  )

data.table::fwrite(metric_df_b, file.path(out_dir, "SourceData_S5b_classifier_metrics_long.tsv"), sep = "\t")

p_s5b <- ggplot(metric_df_b, aes(x = Celltype, y = Score, fill = Celltype)) +
  stat_boxplot(geom = "errorbar", width = 0.34, linewidth = 0.35, color = "black") +
  geom_boxplot(
    width = 0.65,
    linewidth = 0.35,
    color = "black",
    outlier.shape = 21,
    outlier.size = 1.15,
    outlier.stroke = 0.25,
    outlier.fill = "white",
    outlier.colour = "black"
  ) +
  stat_summary(
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
  facet_wrap(~Metric, ncol = 1) +
  scale_fill_manual(values = celltype_cols, drop = FALSE) +
  scale_x_discrete(labels = celltype_labels) +
  coord_cartesian(ylim = c(min(metric_df_b$Score, na.rm = TRUE) - 0.03, 1.0)) +
  scale_y_continuous(breaks = seq(0, 1, 0.1), expand = expansion(mult = c(0.02, 0.02))) +
  paper_theme(base_size = 9.5) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    strip.text = element_text(size = 10),
    panel.spacing = unit(4, "mm")
  ) +
  labs(title = "b", x = NULL, y = NULL)

save_gg_pdf(file.path(out_dir, "SuppFigS5b_classifier_metrics_boxplots.pdf"), p_s5b, width = 3.3, height = 4.9)
ggsave(file.path(out_dir, "SuppFigS5b_classifier_metrics_boxplots.png"), p_s5b, width = 3.3, height = 4.9, dpi = 300)

# -------------------------------
# 7. S5c pooled recall and precision heatmaps
# -------------------------------
heat_rows <- list()

# Majority vote
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

# Per Level 1 classifier
for (ct in setdiff(celltype_order_existing, "Majority vote")) {
  sub_ct <- dplyr::filter(pred_df, as.character(Level1) == ct)
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

heat_long <- bind_rows(heat_rows) %>%
  mutate(
    Disease = factor(Disease, levels = disease_order),
    Celltype = factor(Celltype, levels = celltype_order_existing)
  )

data.table::fwrite(heat_long, file.path(out_dir, "SourceData_S5c_recall_precision_heatmaps.tsv"), sep = "\t")

make_heat_matrix <- function(df, value_col) {
  mat_df <- df %>%
    dplyr::select(Disease, Celltype, !!rlang::sym(value_col)) %>%
    tidyr::pivot_wider(names_from = Celltype, values_from = !!rlang::sym(value_col)) %>%
    as.data.frame()
  rownames(mat_df) <- as.character(mat_df$Disease)
  mat_df$Disease <- NULL
  mat <- as.matrix(mat_df)
  mat[is.na(mat)] <- 0
  mat
}

recall_mat <- make_heat_matrix(heat_long, "Recall")
precision_mat <- make_heat_matrix(heat_long, "Precision")

ha_bottom <- ComplexHeatmap::HeatmapAnnotation(
  `Cell types` = factor(colnames(recall_mat), levels = celltype_order_existing),
  col = list(`Cell types` = celltype_cols),
  show_annotation_name = FALSE,
  show_legend = FALSE,
  which = "column",
  simple_anno_size = unit(3.2, "mm"),
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
  bottom_annotation = ha_bottom,
  column_labels = unname(celltype_labels[colnames(recall_mat)]),
  column_names_rot = 90,
  column_names_gp = gpar(fontsize = 8, fontfamily = base_family),
  row_names_gp = gpar(fontsize = 8.5, fontfamily = base_family),
  rect_gp = gpar(col = NA),
  heatmap_legend_param = score_heatmap_legend_param("Recall score"),
  column_title = "Recall score",
  column_title_gp = gpar(fontsize = 11, fontfamily = base_family),
  row_title = NULL,
  width = unit(ncol(recall_mat) * 4.6, "mm"),
  height = unit(nrow(recall_mat) * 4.6, "mm")
)

ht_precision <- ComplexHeatmap::Heatmap(
  precision_mat,
  name = "Precision score",
  col = score_col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  bottom_annotation = ha_bottom,
  column_labels = unname(celltype_labels[colnames(precision_mat)]),
  column_names_rot = 90,
  column_names_gp = gpar(fontsize = 8, fontfamily = base_family),
  row_names_gp = gpar(fontsize = 8.5, fontfamily = base_family),
  rect_gp = gpar(col = NA),
  heatmap_legend_param = score_heatmap_legend_param("Precision score"),
  column_title = "Precision score",
  column_title_gp = gpar(fontsize = 11, fontfamily = base_family),
  row_title = NULL,
  width = unit(ncol(precision_mat) * 4.6, "mm"),
  height = unit(nrow(precision_mat) * 4.6, "mm")
)

w_heat <- max(6.4, ncol(recall_mat) * 0.23 + 4.0)
h_heat <- max(4.6, nrow(recall_mat) * 0.21 + 1.6)

open_pdf(file.path(out_dir, "SuppFigS5c_recall_heatmap.pdf"), width = w_heat, height = h_heat)
draw(ht_recall, heatmap_legend_side = "top", annotation_legend_side = "bottom", newpage = TRUE)
dev.off()

png(file.path(out_dir, "SuppFigS5c_recall_heatmap.png"), width = w_heat, height = h_heat, units = "in", res = 300)
draw(ht_recall, heatmap_legend_side = "top", annotation_legend_side = "bottom", newpage = TRUE)
dev.off()

open_pdf(file.path(out_dir, "SuppFigS5c_precision_heatmap.pdf"), width = w_heat, height = h_heat)
draw(ht_precision, heatmap_legend_side = "top", annotation_legend_side = "bottom", newpage = TRUE)
dev.off()

png(file.path(out_dir, "SuppFigS5c_precision_heatmap.png"), width = w_heat, height = h_heat, units = "in", res = 300)
draw(ht_precision, heatmap_legend_side = "top", annotation_legend_side = "bottom", newpage = TRUE)
dev.off()

# optional combined panel c for quick checking
open_pdf(file.path(out_dir, "SuppFigS5c_recall_precision_combined.pdf"), width = 10.6, height = h_heat)
grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow = 1, ncol = 2, widths = unit(c(1, 1), "null"))))
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
draw(ht_recall, heatmap_legend_side = "top", annotation_legend_side = "bottom", newpage = FALSE)
upViewport()
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
draw(ht_precision, heatmap_legend_side = "top", annotation_legend_side = "bottom", newpage = FALSE)
upViewport(2)
dev.off()

png(file.path(out_dir, "SuppFigS5c_recall_precision_combined.png"), width = 10.6, height = h_heat, units = "in", res = 300)
grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow = 1, ncol = 2, widths = unit(c(1, 1), "null"))))
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
draw(ht_recall, heatmap_legend_side = "top", annotation_legend_side = "bottom", newpage = FALSE)
upViewport()
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
draw(ht_precision, heatmap_legend_side = "top", annotation_legend_side = "bottom", newpage = FALSE)
upViewport(2)
dev.off()

# -------------------------------
# 8. Summary
# -------------------------------
summary_list <- list(
  n_cv_splits = dplyr::n_distinct(split_df$split),
  diseases = disease_order,
  cell_types = celltype_order_existing,
  generated = c(
    "SuppFigS5a_majority_per_disease_metrics_boxplots.pdf/png",
    "SuppFigS5b_classifier_metrics_boxplots.pdf/png",
    "SuppFigS5c_recall_heatmap.pdf/png",
    "SuppFigS5c_precision_heatmap.pdf/png",
    "SuppFigS5c_recall_precision_combined.pdf/png"
  )
)

jsonlite::write_json(
  summary_list,
  file.path(out_dir, "suppfigS5_paper_mimic_summary.json"),
  auto_unbox = TRUE,
  pretty = TRUE
)

cat("[DONE] Supplementary Fig. S5 paper mimic plots saved to:\n")
cat(out_dir, "\n")
cat("Generated files:\n")
cat("  - SuppFigS5a_majority_per_disease_metrics_boxplots.pdf/png\n")
cat("  - SuppFigS5b_classifier_metrics_boxplots.pdf/png\n")
cat("  - SuppFigS5c_recall_heatmap.pdf/png\n")
cat("  - SuppFigS5c_precision_heatmap.pdf/png\n")
cat("  - SuppFigS5c_recall_precision_combined.pdf/png\n")
cat("  - SourceData tables\n")
