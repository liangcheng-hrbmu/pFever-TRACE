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
# Fig 5c, Fig 5d and Fig 5e
# Scenario 1: five fold cross validation performance
#
# 本脚本只输出三张独立图，不再拼合 c d e 总图。
#
# Fig 5c:
#   每个五折拆分中的 WF1 分数分布。
#   箱体表示 IQR，中位数为箱体内横线，胡须为 1.5 x IQR，离群值为单点。
#
# Fig 5d:
#   汇总所有五折预测后，计算每种疾病和每个 Level 1 细胞类型分类器的 one versus rest F1。
#
# Fig 5e:
#   汇总所有五折预测后，绘制列归一化混淆矩阵。
#   每一列是真实疾病，每一格表示该真实疾病被预测为某一类别的比例。
#   对角线值可理解为各疾病的 recall。
# =========================================================


# =========================================================
# 0. 用户需要修改的路径
# =========================================================
base_dir <- "data/analysis/08_patient_classifier/03_patient_classifier_scenario1_v3"

split_summary_tsv  <- file.path(base_dir, "scenario1_split_summary.tsv")
overall_json       <- file.path(base_dir, "scenario1_overall_metrics.json")
level1_metrics_tsv <- file.path(base_dir, "scenario1_level1_metrics_all.tsv")
majority_pred_tsv  <- file.path(base_dir, "scenario1_majority_predictions_all.tsv")
level1_pred_tsv    <- file.path(base_dir, "query_level1_predictions_all.tsv")

out_dir <- "data/analysis/08_patient_classifier/04_scenario1_plots_R_paper_mimic_separate"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


# =========================================================
# 1. 输出设备
# =========================================================
# 说明：
# 服务器上 Arial 经常会导致 PDF 字体错误。
# 这里统一使用 cairo_pdf + sans，PDF 文字通常仍可在 AI 中编辑。
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


# =========================================================
# 2. 颜色和显示顺序
# =========================================================
# 2.1 连续分数色板
# 参考原文性能热图风格：低值接近黑紫，高值接近浅黄。
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

# 2.2 Level 1 细胞类型颜色
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

# 2.3 疾病大类颜色，用于 Fig 5e 顶部和右侧注释条
disease_group_cols_full <- c(
  "Healthy"           = "#808080",
  "Bacteria"          = "#C48EAD",
  "Virus"             = "#6699CC",
  "Sepsis"            = "#FF2A2A",
  "Hyperinflammation" = "#E3A722",
  "Autoimmune"        = "#36A596",
  "Other"             = "#BDBDBD",
  "Not assigned"      = "#2F225B"
)

disease_group_levels <- c(
  "Healthy", "Bacteria", "Virus", "Sepsis",
  "Hyperinflammation", "Autoimmune", "Other", "Not assigned"
)

# 2.4 疾病显示顺序
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


# =========================================================
# 3. 标签清洗和统计函数
# =========================================================
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

get_disease_group <- function(x) {
  x <- as.character(x)

  map <- c(
    "Healthy" = "Healthy",
    "Bact.inf" = "Bacteria",
    "Dengue" = "Virus",
    "HFMD" = "Virus",
    "Flu" = "Virus",
    "COVID" = "Virus",
    "Sepsis" = "Sepsis",
    "KD" = "Hyperinflammation",
    "MIS-C" = "Hyperinflammation",
    "SLE" = "Autoimmune",
    "JIA" = "Autoimmune",
    "Not assigned" = "Not assigned"
  )

  out <- ifelse(x %in% names(map), map[x], "Other")
  factor(out, levels = disease_group_levels)
}

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

get_annotation_cols <- function(group_factor) {
  used_levels <- levels(droplevels(group_factor))
  missing_levels <- setdiff(used_levels, names(disease_group_cols_full))

  if (length(missing_levels) > 0) {
    extra_cols <- stats::setNames(rep("#BDBDBD", length(missing_levels)), missing_levels)
    cols <- c(disease_group_cols_full, extra_cols)
  } else {
    cols <- disease_group_cols_full
  }

  cols[used_levels]
}


# =========================================================
# 4. 数据读入和统一预处理
# =========================================================
split_df  <- data.table::fread(split_summary_tsv, sep = "\t", data.table = FALSE)
level1_df <- data.table::fread(level1_metrics_tsv, sep = "\t", data.table = FALSE)
maj_df    <- data.table::fread(majority_pred_tsv, sep = "\t", data.table = FALSE)

overall <- NULL
if (file.exists(overall_json)) {
  overall <- jsonlite::fromJSON(overall_json)
}

if (!file.exists(level1_pred_tsv)) {
  stop("Cannot find query_level1_predictions_all.tsv. Please merge split level predictions first.")
}
pred_df <- data.table::fread(level1_pred_tsv, sep = "\t", data.table = FALSE)

# 如果存在 benchmark_status，只保留 primary benchmark，避免混入其他附加预测。
if ("benchmark_status" %in% colnames(maj_df)) {
  maj_df <- dplyr::filter(maj_df, as.character(benchmark_status) == "primary_benchmark")
}
if ("benchmark_status" %in% colnames(pred_df)) {
  pred_df <- dplyr::filter(pred_df, as.character(benchmark_status) == "primary_benchmark")
}

# 统一细胞类型和疾病名称。
level1_df$Level1 <- short_celltype(level1_df$Level1)
pred_df$Level1 <- short_celltype(pred_df$Level1)

maj_df$true_short <- short_disease(maj_df$true_disease)
maj_df$pred_short <- short_disease(maj_df$pred_disease_majority)

pred_df$true_short <- short_disease(pred_df$disease_Level_3)
pred_df$pred_short <- short_disease(pred_df$pred_disease)

# 按真实疾病决定 Fig 5d 和 Fig 5e 的列顺序。
truth_diseases <- setdiff(unique(maj_df$true_short), "Not assigned")
disease_order <- c(
  preferred_disease_order[preferred_disease_order %in% truth_diseases],
  setdiff(sort(truth_diseases), preferred_disease_order)
)

# 预测类别中可能存在 Not assigned 或额外类别，因此行顺序单独确定。
pred_present <- unique(maj_df$pred_short)
pred_row_order <- unique(c(
  disease_order,
  setdiff(sort(pred_present), c(disease_order, "Not assigned")),
  "Not assigned"
))

# 如果某次运行缺失某个 Level 1 细胞类型，只保留实际存在的列。
celltype_order_existing <- celltype_order[
  celltype_order %in% c("Majority vote", unique(level1_df$Level1), unique(pred_df$Level1))
]


# =========================================================
# 5. Fig 5c: WF1 箱线图
# =========================================================
make_fig5c <- function() {
  box_df <- dplyr::bind_rows(
    dplyr::transmute(
      split_df,
      split = as.character(split),
      Celltype = "Majority vote",
      WF1 = as.numeric(WF1_majority)
    ),
    dplyr::transmute(
      level1_df,
      split = as.character(split),
      Celltype = short_celltype(Level1),
      WF1 = as.numeric(WF1)
    )
  ) %>%
    dplyr::filter(!is.na(WF1)) %>%
    dplyr::mutate(Celltype = factor(Celltype, levels = celltype_order_existing))

  data.table::fwrite(
    box_df,
    file.path(out_dir, "SourceData_Fig5c_WF1_by_split.tsv"),
    sep = "\t"
  )

  y_min <- max(0, floor(min(box_df$WF1, na.rm = TRUE) * 10) / 10)
  y_max <- min(1, ceiling(max(box_df$WF1, na.rm = TRUE) * 20) / 20)
  if (y_max < 0.98) y_max <- 1

  p_c <- ggplot2::ggplot(box_df, ggplot2::aes(x = Celltype, y = WF1, fill = Celltype)) +
    # 胡须末端横杠。geom_boxplot 默认不显示横杠，所以这里单独加。
    ggplot2::stat_boxplot(
      geom = "errorbar",
      width = 0.34,
      linewidth = 0.35,
      color = "black"
    ) +
    # 箱体。默认统计就是 IQR，中位数和 1.5 x IQR 胡须。
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
    # 显式加中位数横线，防止颜色较深时看不清。
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
    # 原文风格中 Majority vote 和各细胞类型之间有一条分隔线。
    ggplot2::geom_vline(xintercept = 1.5, color = "black", linewidth = 0.35) +
    ggplot2::scale_fill_manual(values = celltype_cols, drop = FALSE) +
    ggplot2::scale_x_discrete(labels = celltype_labels) +
    ggplot2::coord_cartesian(ylim = c(y_min, y_max)) +
    ggplot2::scale_y_continuous(
      breaks = seq(0, 1, 0.05),
      expand = ggplot2::expansion(mult = c(0.02, 0.02))
    ) +
    paper_theme(base_size = 10) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1),
      legend.position = "none"
    ) +
    ggplot2::labs(
      title = "c",
      x = NULL,
      y = "WF1 score"
    )

  save_gg_pdf(
    file.path(out_dir, "Fig5c_scenario1_WF1_boxplot.pdf"),
    p_c,
    width = 4.2,
    height = 3.2
  )
  ggplot2::ggsave(
    file.path(out_dir, "Fig5c_scenario1_WF1_boxplot.png"),
    p_c,
    width = 4.2,
    height = 3.2,
    dpi = 300
  )

  invisible(list(plot = p_c, data = box_df))
}


# =========================================================
# 6. Fig 5d: 疾病 x 细胞类型 F1 热图
# =========================================================
make_fig5d <- function() {
  d_rows <- list()

  # Majority vote 作为第一列。
  for (d in disease_order) {
    d_rows[[length(d_rows) + 1]] <- data.frame(
      Disease = d,
      Celltype = "Majority vote",
      F1 = compute_f1_bin(maj_df$true_short, maj_df$pred_short, d),
      stringsAsFactors = FALSE
    )
  }

  # 每个 Level 1 细胞类型分类器分别计算 one versus rest F1。
  for (ct in setdiff(celltype_order_existing, "Majority vote")) {
    sub_ct <- dplyr::filter(pred_df, as.character(Level1) == ct)

    for (d in disease_order) {
      d_rows[[length(d_rows) + 1]] <- data.frame(
        Disease = d,
        Celltype = ct,
        F1 = compute_f1_bin(sub_ct$true_short, sub_ct$pred_short, d),
        stringsAsFactors = FALSE
      )
    }
  }

  d_long <- dplyr::bind_rows(d_rows) %>%
    dplyr::mutate(
      Disease = factor(Disease, levels = disease_order),
      Celltype = factor(Celltype, levels = celltype_order_existing)
    )

  data.table::fwrite(
    d_long,
    file.path(out_dir, "SourceData_Fig5d_F1_by_disease_celltype.tsv"),
    sep = "\t"
  )

  d_mat <- d_long %>%
    tidyr::pivot_wider(names_from = Celltype, values_from = F1) %>%
    as.data.frame()

  rownames(d_mat) <- as.character(d_mat$Disease)
  d_mat$Disease <- NULL
  d_mat <- as.matrix(d_mat)
  d_mat[is.na(d_mat)] <- 0

  ha_bottom <- ComplexHeatmap::HeatmapAnnotation(
    `Cell types` = factor(colnames(d_mat), levels = celltype_order_existing),
    col = list(`Cell types` = celltype_cols),
    show_annotation_name = FALSE,
    show_legend = FALSE,
    which = "column",
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
    bottom_annotation = ha_bottom,
    column_labels = unname(celltype_labels[colnames(d_mat)]),
    column_names_rot = 90,
    column_names_gp = grid::gpar(fontsize = 8, fontfamily = base_family),
    row_names_gp = grid::gpar(fontsize = 8.5, fontfamily = base_family),
    rect_gp = grid::gpar(col = NA),
    na_col = "#000004",
    heatmap_legend_param = score_heatmap_legend_param("F1 score"),
    column_title = NULL,
    row_title = NULL,
    width = grid::unit(ncol(d_mat) * 4.6, "mm"),
    height = grid::unit(nrow(d_mat) * 4.6, "mm")
  )

  d_w <- max(3.2, ncol(d_mat) * 0.22 + 2.2)
  d_h <- max(3.0, nrow(d_mat) * 0.22 + 1.6)

  open_pdf(file.path(out_dir, "Fig5d_scenario1_F1_heatmap.pdf"), width = d_w, height = d_h)
  ComplexHeatmap::draw(
    ht_d,
    heatmap_legend_side = "top",
    annotation_legend_side = "bottom",
    newpage = TRUE
  )
  dev.off()

  png(
    file.path(out_dir, "Fig5d_scenario1_F1_heatmap.png"),
    width = d_w,
    height = d_h,
    units = "in",
    res = 300
  )
  ComplexHeatmap::draw(
    ht_d,
    heatmap_legend_side = "top",
    annotation_legend_side = "bottom",
    newpage = TRUE
  )
  dev.off()

  invisible(list(heatmap = ht_d, data_long = d_long, matrix = d_mat))
}


# =========================================================
# 7. Fig 5e: 列归一化混淆矩阵
# =========================================================
make_fig5e <- function() {
  # 行为预测类别，列为真实疾病。
  cm <- table(
    factor(maj_df$pred_short, levels = pred_row_order),
    factor(maj_df$true_short, levels = disease_order)
  )
  cm <- as.matrix(cm)

  # 列归一化：每个真实疾病内部，预测到各类别的比例。
  cm_col <- sweep(cm, 2, colSums(cm), "/")
  cm_col[is.na(cm_col)] <- 0

  cm_out <- as.data.frame(cm_col)
  cm_out$Predicted_condition <- rownames(cm_out)
  cm_out <- cm_out[, c("Predicted_condition", setdiff(colnames(cm_out), "Predicted_condition")), drop = FALSE]

  data.table::fwrite(
    cm_out,
    file.path(out_dir, "SourceData_Fig5e_confusion_colnorm.tsv"),
    sep = "\t"
  )

  top_group <- droplevels(get_disease_group(colnames(cm_col)))
  right_group <- droplevels(get_disease_group(rownames(cm_col)))

  # ComplexHeatmap 要求注释条颜色名和 factor level 完全一致。
  annotation_group_cols <- get_annotation_cols(c(top_group, right_group))

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

  e_w <- max(3.4, ncol(cm_col) * 0.25 + 2.2)
  e_h <- max(4.0, nrow(cm_col) * 0.25 + 1.7)

  open_pdf(file.path(out_dir, "Fig5e_scenario1_confusion_matrix.pdf"), width = e_w, height = e_h)
  ComplexHeatmap::draw(
    ht_e,
    heatmap_legend_side = "top",
    annotation_legend_side = "bottom",
    merge_legends = TRUE,
    newpage = TRUE
  )
  dev.off()

  png(
    file.path(out_dir, "Fig5e_scenario1_confusion_matrix.png"),
    width = e_w,
    height = e_h,
    units = "in",
    res = 300
  )
  ComplexHeatmap::draw(
    ht_e,
    heatmap_legend_side = "top",
    annotation_legend_side = "bottom",
    merge_legends = TRUE,
    newpage = TRUE
  )
  dev.off()

  invisible(list(heatmap = ht_e, matrix = cm_col, raw_table = cm))
}


# =========================================================
# 8. 分别运行三张图
# =========================================================
res_c <- make_fig5c()
res_d <- make_fig5d()
res_e <- make_fig5e()


# =========================================================
# 9. 运行总结
# =========================================================
summary_list <- list(
  n_pooled_majority_predictions = nrow(maj_df),
  n_cv_splits_in_metric_tables = dplyr::n_distinct(res_c$data$split),
  diseases = disease_order,
  predicted_rows = pred_row_order,
  cell_types = celltype_order_existing,
  note = "Fig 5c, 5d and 5e are exported as separate files. Fig 5c uses IQR boxplots with 1.5 x IQR whiskers. Fig 5d uses pooled one versus rest F1. Fig 5e uses column normalized confusion proportions."
)

if (!is.null(overall)) {
  summary_list$overall <- overall
}

jsonlite::write_json(
  summary_list,
  file.path(out_dir, "scenario1_cde_separate_summary.json"),
  auto_unbox = TRUE,
  pretty = TRUE
)

cat("[DONE] Separate Scenario 1 plots saved to:\n")
cat(out_dir, "\n")
cat("Generated files:\n")
cat("  - Fig5c_scenario1_WF1_boxplot.pdf/png\n")
cat("  - Fig5d_scenario1_F1_heatmap.pdf/png\n")
cat("  - Fig5e_scenario1_confusion_matrix.pdf/png\n")
cat("  - SourceData_Fig5c_WF1_by_split.tsv\n")
cat("  - SourceData_Fig5d_F1_by_disease_celltype.tsv\n")
cat("  - SourceData_Fig5e_confusion_colnorm.tsv\n")
cat("  - scenario1_cde_separate_summary.json\n")
