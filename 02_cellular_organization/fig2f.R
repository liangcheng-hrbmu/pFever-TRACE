# ============================================================
# Fig. 2f: Top Level 2 populations defining each cellular module
# K = 4
# Full label version with original majorCluster colors
# Legend on the right, fully shown
# ============================================================

setwd("data/analysis/atlas_objects")

suppressPackageStartupMessages({
  library(Seurat)
  library(NMF)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(grid)
  library(stringr)
})

# ============================================================
# 1. Load Seurat object and NMF_K4
# ============================================================

outdir <- file.path(getwd(), "fig2_redraw")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

seurat_obj <- readRDS("integrated_pediatric.rds")
meta <- as.data.frame(seurat_obj@meta.data)

# 兼容列名
if ("sample_id" %in% colnames(meta) && !"sampleID" %in% colnames(meta)) {
  colnames(meta)[colnames(meta) == "sample_id"] <- "sampleID"
}
if ("aifi_L1_paper" %in% colnames(meta) && !"majorCluster" %in% colnames(meta)) {
  colnames(meta)[colnames(meta) == "aifi_L1_paper"] <- "majorCluster"
}
if ("aifi_L2_paper" %in% colnames(meta) && !"subCluster" %in% colnames(meta)) {
  colnames(meta)[colnames(meta) == "aifi_L2_paper"] <- "subCluster"
}

stopifnot("majorCluster" %in% colnames(meta))
stopifnot("subCluster" %in% colnames(meta))

meta$majorCluster <- trimws(as.character(meta$majorCluster))
meta$subCluster   <- trimws(as.character(meta$subCluster))

nmf_rds <- file.path(outdir, "NMF_K4_cellular_modules_CM01_CM04.rds")

if (!file.exists(nmf_rds)) {
  stop(
    "Cannot find NMF_K4 file:\n",
    nmf_rds,
    "\nPlease run the K=4 NMF script first."
  )
}

NMF_K4 <- readRDS(nmf_rds)

# ============================================================
# 2. Extract NMF basis matrix
# rows = Level 2 populations
# columns = CM01-CM04
# ============================================================

w <- as.matrix(NMF::basis(NMF_K4))

# 如果方向反了，自动转置
if (nrow(w) == 4 && ncol(w) > 4) {
  w <- t(w)
}

if (ncol(w) != 4) {
  stop(
    "basis matrix does not look like Level2 x 4 CM matrix. Current dim: ",
    paste(dim(w), collapse = " x ")
  )
}

colnames(w) <- sprintf("CM%02d", 1:4)

if (is.null(rownames(w)) || any(is.na(rownames(w))) || any(rownames(w) == "")) {
  stop("basis matrix has no Level 2 population rownames. Please check NMF input matrix.")
}

cat("basis matrix dim:\n")
print(dim(w))

cat("basis colnames:\n")
print(colnames(w))

cat("basis row examples:\n")
print(head(rownames(w), 20))

# ============================================================
# 3. Build Level 2 -> majorCluster mapping
# ============================================================

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

cat("Level 2 to Level 1 mapping examples:\n")
print(head(sub_to_major, 20))

# ============================================================
# 4. MajorCluster color system
# Use the same celltype_colors as Fig. 1 / Fig. 2a
# ============================================================

celltype_col <- c(
  "B"               = "#1F77B4",
  "Plasma"          = "#AEC7E8",
  "pDC"             = "#FF7F0E",
  "DC"              = "#FFBB78",
  "Mono"            = "#2CA02C",
  "T CD4 naive"     = "#98DF8A",
  "T CD4 non-naive" = "#D62728",
  "T CD8 naive"     = "#FF9896",
  "T CD8 non-naive" = "#9467BD",
  "UTC"             = "#C5B0D5",
  "ILC"             = "#8C564B",
  "Other"           = "#BDBDBD"
)

cm_title_cols <- c(
  "CM01" = "#F26B6C",
  "CM02" = "#F4A261",
  "CM03" = "#4CB7A5",
  "CM04" = "#A78BDA"
)

# ============================================================
# 5. Prepare top Level 2 populations per CM
# ============================================================

TOP_N <- 10

top_df <- as.data.frame(w) %>%
  tibble::rownames_to_column("subCluster") %>%
  tidyr::pivot_longer(
    cols = starts_with("CM"),
    names_to = "CM",
    values_to = "weight"
  ) %>%
  dplyr::group_by(CM) %>%
  dplyr::slice_max(order_by = weight, n = TOP_N, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(sub_to_major, by = "subCluster") %>%
  dplyr::mutate(
    majorCluster = ifelse(is.na(majorCluster), "Other", majorCluster),
    majorCluster = factor(
      majorCluster,
      levels = c(
        "B", "Plasma", "pDC", "DC", "Mono",
        "T CD4 naive", "T CD4 non-naive",
        "T CD8 naive", "T CD8 non-naive",
        "UTC", "ILC", "Other"
      )
    ),
    CM = factor(CM, levels = sprintf("CM%02d", 1:4))
  )

# 如果出现未设定颜色的 majorCluster，自动补色
missing_major <- setdiff(
  as.character(unique(top_df$majorCluster)),
  names(celltype_col)
)

missing_major <- missing_major[!is.na(missing_major)]

if (length(missing_major) > 0) {
  extra_cols <- grDevices::hcl.colors(length(missing_major), "Set 2")
  names(extra_cols) <- missing_major
  celltype_col <- c(celltype_col, extra_cols)
  message("Added colors for missing majorCluster: ", paste(missing_major, collapse = ", "))
}

# 每个 CM 内按 weight 排序
top_df <- top_df %>%
  dplyr::group_by(CM) %>%
  dplyr::arrange(weight, .by_group = TRUE) %>%
  dplyr::mutate(
    feature_id = paste(CM, subCluster, sep = "___"),
    feature_id = factor(feature_id, levels = unique(feature_id))
  ) %>%
  dplyr::ungroup()

cat("Top Level 2 populations per CM:\n")
print(top_df %>% dplyr::select(CM, subCluster, majorCluster, weight))

write.csv(
  top_df,
  file.path(outdir, "Fig2f_top_Level2_populations_defining_CM_K4_fullLabels_data.csv"),
  row.names = FALSE
)

# ============================================================
# 6. Plot each CM separately
# ============================================================

label_fun <- function(x) {
  x <- sub("^.*___", "", x)
  stringr::str_wrap(x, width = 30)
}

# 统一 x 轴范围
global_xmax <- max(top_df$weight, na.rm = TRUE) * 1.15

make_one_cm_plot <- function(cm_id) {
  df <- top_df %>%
    dplyr::filter(CM == cm_id)

  ggplot(df, aes(x = weight, y = feature_id, fill = majorCluster)) +
    geom_col(width = 0.72, color = "white", linewidth = 0.25) +
    scale_y_discrete(labels = label_fun) +
    scale_x_continuous(
      limits = c(0, global_xmax),
      expand = expansion(mult = c(0, 0.02))
    ) +
    scale_fill_manual(
      values = celltype_col,
      drop = FALSE,
      name = "Parent Level 1\ncompartment"
    ) +
    labs(
      title = cm_id,
      x = "Module weight",
      y = NULL
    ) +
    theme_classic(base_size = 9) +
    theme(
      plot.title = element_text(
        face = "bold",
        size = 13,
        hjust = 0.5,
        color = cm_title_cols[cm_id]
      ),

      axis.text.y = element_text(
        size = 7.5,
        color = "black",
        lineheight = 0.9
      ),
      axis.text.x = element_text(size = 7, color = "black"),
      axis.title.x = element_text(size = 8),

      axis.line = element_line(linewidth = 0.35),
      axis.ticks = element_line(linewidth = 0.3),

      legend.position = "right",
      legend.title = element_text(size = 8, face = "bold"),
      legend.text = element_text(size = 7.5),
      legend.key.size = unit(3.5, "mm"),
      legend.spacing.y = unit(1.5, "mm"),

      plot.margin = margin(t = 4, r = 10, b = 4, l = 14)
    ) +
    guides(
      fill = guide_legend(
        ncol = 1,
        byrow = TRUE,
        override.aes = list(alpha = 1)
      )
    )
}

p_list <- lapply(sprintf("CM%02d", 1:4), make_one_cm_plot)

# ============================================================
# 7. Combine plots
# Legend collected on the right
# ============================================================

p_final <- patchwork::wrap_plots(
  p_list,
  nrow = 1,
  guides = "collect",
  widths = c(1.25, 1.25, 1.25, 1.25)
) +
  patchwork::plot_annotation(
    title = "Top Level 2 populations defining each cellular module",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0)
    )
  ) &
  theme(
    legend.position = "right",
    legend.box = "vertical",
    legend.margin = margin(0, 0, 0, 4),
    plot.background = element_rect(fill = "white", color = NA)
  )

# ============================================================
# 8. Save
# ============================================================

pdf_file <- file.path(outdir, "Fig2f_top_Level2_populations_defining_CM_K4_rightLegend_fullLabels.pdf")
png_file <- file.path(outdir, "Fig2f_top_Level2_populations_defining_CM_K4_rightLegend_fullLabels.png")

ggsave(
  filename = pdf_file,
  plot = p_final,
  width = 18.5,
  height = 5.6,
  device = cairo_pdf,
  limitsize = FALSE
)

ggsave(
  filename = png_file,
  plot = p_final,
  width = 18.5,
  height = 5.6,
  dpi = 600,
  limitsize = FALSE
)

message("[DONE] Saved:")
message(pdf_file)
message(png_file)
