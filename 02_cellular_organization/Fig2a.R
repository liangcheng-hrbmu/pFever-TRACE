

library(dplyr)
library(tidyr)

meta <- seurat_obj@meta.data






##########################################
### Fig1e-style: Disease prevalence heatmap (Ro/e)
### 使用新的疾病简称 + 新的疾病分组注释体系
##########################################


library(ComplexHeatmap)
library(circlize)
library(grid)
library(dplyr)

## =========================================================
## 0) 输入对象检查
## =========================================================
# 你需要提前保证环境中已有:
# meta : data.frame / metadata
# 且至少包含以下两列:
#   1) disease_Level_3
#   2) majorCluster
#
# 如果你有 celltype_colors，也会自动用于列顺序
#
# 例如：
# meta <- seurat_obj@meta.data

stopifnot(exists("meta"))
stopifnot("disease_Level_3" %in% colnames(meta))
stopifnot("majorCluster" %in% colnames(meta))

## =========================================================
## 1) 参数
## =========================================================
row_var_raw <- "disease_Level_3"   # 原始长疾病名列（保留）
col_var     <- "majorCluster"      # 列变量：主细胞群

# 小格子阈值：某 disease x celltype 计数低于此值则设为 NA
#min_cell_per_celltype_in_group <- 5

# 是否把健康组纳入主图
include_healthy <- FALSE

# 输出文件名
out_pdf <- "Fig1e_style_RoE_DiseaseShort_x_majorCluster_newGroup.pdf"

## =========================================================
## 2) 疾病长名 -> 短名映射
## =========================================================
disease_short_map <- c(
  "Dengue human infection(DENV-1)" = "Dengue",
  "HFMD(EV-A71)enterovirus" = "HFMD",
  "Influenza A virus (IAV)" = "Flu",
  "Kawasaki disease (KD)" = "KD",
  "Multisystem inflammatory syndrome in children (MIS-C)" = "MIS-C",
  "Pediatric Sepsis" = "Sepsis",
  "bacterial infection" = "Bact.inf",
  "human systemic lupus erythematosus (SLE)" = "SLE",
  "juvenile idopathic arthritis (JIA)" = "JIA",
  "severe acute respiratory syndrome coronavirus 2 (SARS-CoV-2)" = "COVID"
)

# 如果你后面想纳入健康组，在这里加
# disease_short_map <- c(
#   disease_short_map,
#   "Healthy children" = "Healthy"
# )

## 固定疾病顺序
disease_order <- c(
  "Dengue", "HFMD", "Flu", "COVID",
  "Bact.inf", "Sepsis", "KD", "MIS-C", "SLE", "JIA"
)

# 若 include_healthy = TRUE，可用：
# disease_order <- c(
#   "Healthy",
#   "Dengue", "HFMD", "Flu", "COVID",
#   "Bact.inf", "Sepsis", "KD", "MIS-C", "SLE", "JIA"
# )

## =========================================================
## 3) 短名 -> 疾病大类映射
## =========================================================
disease_group_map <- c(
  "Dengue"   = "virus",
  "HFMD"     = "virus",
  "Flu"      = "virus",
  "COVID"    = "virus",
  "Bact.inf" = "bacteria",
  "Sepsis"   = "Sepsis",
  "KD"       = "Hyperinflam",
  "MIS-C"    = "Hyperinflam",
  "SLE"      = "Autoimmune",
  "JIA"      = "Autoimmune"
)

# 如果纳入健康组：
# disease_group_map <- c(
#   "Healthy"  = "Healthy",
#   disease_group_map
# )

## 疾病大类显示顺序
group_levels <- c("virus", "bacteria", "Sepsis", "Hyperinflam", "Autoimmune")

# 如果纳入健康组：
# group_levels <- c("Healthy", "virus", "bacteria", "Sepsis", "Hyperinflam", "Autoimmune")

## 颜色：尽量贴近你现在新图的配色
group_col <- c(
  "virus"       = "#6699CC",
  "bacteria"    = "#C48EAD",
  "Sepsis"      = "#FF2A2A",
  "Hyperinflam" = "#E3A722",
  "Autoimmune"  = "#36A596"
)

# 如果纳入健康组：
# group_col <- c(
#   "Healthy"     = "#BDBDBD",
#   "virus"       = "#6699CC",
#   "bacteria"    = "#C48EAD",
#   "Sepsis"      = "#FF2A2A",
#   "Hyperinflam" = "#E3A722",
#   "Autoimmune"  = "#36A596"
# )

## =========================================================
## 4) 新建 disease_short 和 disease_group（不覆盖原始列）
## =========================================================
meta <- meta %>%
  mutate(
    disease_short = dplyr::recode(.data[[row_var_raw]], !!!disease_short_map),
    disease_short = ifelse(is.na(disease_short), .data[[row_var_raw]], disease_short),
    disease_group = dplyr::recode(disease_short, !!!disease_group_map),
    disease_group = ifelse(is.na(disease_group), "Unknown", disease_group)
  )


meta <- meta %>%
  mutate(
    disease_raw = as.character(.data[[row_var_raw]])
  ) %>%
  mutate(
    disease_short = dplyr::recode(disease_raw, !!!disease_short_map, .default = disease_raw),
    disease_group = dplyr::recode(disease_short, !!!disease_group_map, .default = "Unknown")
  )

## 因子化，固定顺序
meta$disease_short <- factor(meta$disease_short, levels = disease_order)
meta$disease_group <- factor(meta$disease_group, levels = group_levels)

## =========================================================
## 5) 映射检查
## =========================================================
cat("========================================\n")
cat("Check 1: 原始长名中未被 disease_short_map 覆盖的项\n")
print(base::setdiff(unique(meta[[row_var_raw]]), names(disease_short_map)))

cat("========================================\n")
cat("Check 2: disease_short 中未被 disease_group_map 覆盖的项\n")
print(base::setdiff(unique(as.character(meta$disease_short)), names(disease_group_map)))

cat("========================================\n")
cat("Check 3: disease_short x disease_group 交叉表\n")
print(table(meta$disease_short, meta$disease_group, useNA = "ifany"))

## =========================================================
## 6) 是否纳入健康组
## =========================================================
if (!include_healthy) {
  meta_use <- meta %>%
    dplyr::filter(!is.na(disease_short)) %>%
    dplyr::filter(as.character(disease_short) %in% disease_order)
} else {
  meta_use <- meta %>%
    dplyr::filter(!is.na(disease_short))
}

## 再次确保顺序正确
meta_use$disease_short <- factor(as.character(meta_use$disease_short), levels = disease_order)
meta_use$disease_group <- factor(as.character(meta_use$disease_group), levels = group_levels)

## =========================================================
## 7) 计数矩阵 + Ro/e
## =========================================================
row_var  <- "disease_short"
row_anno <- "disease_group"

mat <- base::table(meta_use[[row_var]], meta_use[[col_var]])

# 去掉全空行/列
mat <- mat[rowSums(mat) > 0, colSums(mat) > 0, drop = FALSE]

# 小格子设为 NA
mat2 <- mat
#mat2[mat2 < min_cell_per_celltype_in_group] <- NA

# Row over expected
roe <- base::t(
  base::t(mat2 / base::rowSums(mat2, na.rm = TRUE)) /
    base::colSums(mat2, na.rm = TRUE) * base::sum(mat2, na.rm = TRUE)
)

# 非有限值设为 NA
roe[!is.finite(roe)] <- NA

## =========================================================
## 8) 行注释：用 disease_group，并固定行顺序
## =========================================================
map_df <- data.frame(
  disease_short = meta_use[[row_var]],
  disease_group = meta_use[[row_anno]],
  stringsAsFactors = FALSE
) %>%
  dplyr::count(disease_short, disease_group, name = "n") %>%
  dplyr::group_by(disease_short) %>%
  dplyr::slice_max(order_by = n, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup()

row_group <- map_df$disease_group
names(row_group) <- map_df$disease_short
row_group <- row_group[rownames(roe)]
row_group <- factor(row_group, levels = group_levels)

# 行顺序严格按 disease_order
target_row_order <- disease_order[disease_order %in% rownames(roe)]
roe <- roe[target_row_order, , drop = FALSE]
row_group <- row_group[match(rownames(roe), names(row_group))]

## 右侧注释条
row_ha <- ComplexHeatmap::rowAnnotation(
  `Disease Group` = row_group,
  col = list(`Disease Group` = group_col),
  border = FALSE,
  simple_anno_size = grid::unit(3, "mm"),
  annotation_name_gp = grid::gpar(fontsize = 9),
  annotation_legend_param = list(
    title = "Disease Group",
    at = group_levels[group_levels %in% levels(row_group)],
    title_gp = grid::gpar(fontsize = 9),
    labels_gp = grid::gpar(fontsize = 8),
    grid_height = grid::unit(3, "mm"),
    grid_width  = grid::unit(3, "mm")
  )
)

## =========================================================
## 9) 列顺序：若存在 celltype_colors，则按其顺序
## =========================================================
if (exists("celltype_colors")) {
  common_ct <- base::intersect(colnames(roe), names(celltype_colors))
  common_ct <- common_ct[match(common_ct, names(celltype_colors))]
  common_ct <- common_ct[!is.na(common_ct)]
  roe <- roe[, common_ct, drop = FALSE]
}

## =========================================================
## 10) 热图颜色
## =========================================================
roe_col_fun <- circlize::colorRamp2(
  c(0, 1, 2, 3),
  c("#FFFFFF", "#F7C6C6", "#E07A7A", "#B2182B")
)

## =========================================================
## 11) 画热图
## =========================================================



cell_size_cm <- 0.45

ht <- ComplexHeatmap::Heatmap(
  matrix = roe,
  name = "Ro/e",
  col = roe_col_fun,
  na_col = "#FFFFFF",
  border = TRUE,

  width  = grid::unit(ncol(roe) * cell_size_cm, "cm"),
  height = grid::unit(nrow(roe) * cell_size_cm, "cm"),

  right_annotation = row_ha,

  cluster_rows = FALSE,
  cluster_columns = FALSE,

  row_names_gp = grid::gpar(fontsize = 8),
  column_names_gp = grid::gpar(fontsize = 8),
  column_names_rot = 45,

  column_title = "Disease prevalence of major cell types",
  column_title_gp = grid::gpar(fontsize = 10, fontface = "plain"),

  cell_fun = function(j, i, x, y, width, height, fill) {
    v <- roe[i, j]
    if (!is.na(v)) {
      if (v > 3) {
        grid::grid.text("+++", x, y, gp = grid::gpar(fontsize = 7))
      } else if (v > 2) {
        grid::grid.text("++", x, y, gp = grid::gpar(fontsize = 7))
      } else if (v > 1) {
        grid::grid.text("+", x, y, gp = grid::gpar(fontsize = 7))
      }
    }
  },

  heatmap_legend_param = list(
    title = expression(R[oe]),
    at = c(0, 1, 2, 3),
    labels = c("0", "1 (+)", "2 (++)", "3 (+++)"),
    title_gp = grid::gpar(fontsize = 9),
    labels_gp = grid::gpar(fontsize = 8),
    grid_width = grid::unit(4, "mm"),
    legend_height = grid::unit(25, "mm")
  )
)

pdf_w <- max(7, ncol(roe) * cell_size_cm + 3)
pdf_h <- max(5, nrow(roe) * cell_size_cm + 2)

grDevices::pdf(out_pdf, width = pdf_w, height = pdf_h)

ComplexHeatmap::draw(
  ht,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)

grDevices::dev.off()







suppressPackageStartupMessages({
  library(dplyr)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

plot_roe_heatmap_one_major <- function(
    meta,
    major_name,
    row_var_raw = "disease_Level_3",
    col_state   = "subCluster",

    disease_short_map = c(
      "Dengue human infection(DENV-1)" = "Dengue",
      "HFMD(EV-A71)enterovirus" = "HFMD",
      "Influenza A virus (IAV)" = "Flu",
      "Kawasaki disease (KD)" = "KD",
      "Multisystem inflammatory syndrome in children (MIS-C)" = "MIS-C",
      "Pediatric Sepsis" = "Sepsis",
      "bacterial infection" = "Bact.inf",
      "human systemic lupus erythematosus (SLE)" = "SLE",
      "juvenile idopathic arthritis (JIA)" = "JIA",
      "severe acute respiratory syndrome coronavirus 2 (SARS-CoV-2)" = "COVID"
    ),
    disease_order = c(
      "Dengue", "HFMD", "Flu", "COVID",
      "Bact.inf", "Sepsis", "KD", "MIS-C", "SLE", "JIA"
    ),
    disease_group_map = c(
      "Dengue"   = "virus",
      "HFMD"     = "virus",
      "Flu"      = "virus",
      "COVID"    = "virus",
      "Bact.inf" = "bacteria",
      "Sepsis"   = "Sepsis",
      "KD"       = "Hyperinflam",
      "MIS-C"    = "Hyperinflam",
      "SLE"      = "Autoimmune",
      "JIA"      = "Autoimmune"
    ),
    group_levels = c("virus", "bacteria", "Sepsis", "Hyperinflam", "Autoimmune"),
    group_col = c(
      "virus"       = "#6699CC",
      "bacteria"    = "#C48EAD",
      "Sepsis"      = "#FF2A2A",
      "Hyperinflam" = "#E3A722",
      "Autoimmune"  = "#36A596"
    ),

    include_healthy = FALSE,

    min_cells_per_row = 30,
    min_cells_per_state_total = 30,
    min_cell_per_celltype_in_group = 3,

    cap_for_clust = 3,

    out_pdf = NULL,
    title = NULL,

    column_names_rot = 60,
    row_font = 7,
    col_font = 7,
    plus_font = 7,
    title_font = 10,
    cell_size_cm = 0.42,
    pdf_margin_w = 3.2,
    pdf_margin_h = 2.2,
    verbose = TRUE
) {

  stopifnot(row_var_raw %in% base::colnames(meta))
  stopifnot("majorCluster" %in% base::colnames(meta))
  stopifnot(col_state %in% base::colnames(meta))

  meta2 <- meta %>%
    dplyr::mutate(
      disease_short = dplyr::recode(base::as.character(.data[[row_var_raw]]), !!!disease_short_map),
      disease_short = base::ifelse(base::is.na(disease_short), base::as.character(.data[[row_var_raw]]), disease_short),
      disease_group = dplyr::recode(base::as.character(disease_short), !!!disease_group_map),
      disease_group = base::ifelse(base::is.na(disease_group), "Unknown", disease_group)
    )

  if (verbose) {
    base::cat("========================================\n")
    base::cat("[", major_name, "] Unmapped long disease names:\n", sep = "")
    base::print(
      base::setdiff(
        base::unique(base::as.character(meta[[row_var_raw]])),
        base::names(disease_short_map)
      )
    )

    base::cat("[", major_name, "] Unmapped disease_short to disease_group:\n", sep = "")
    base::print(
      base::setdiff(
        base::unique(base::as.character(meta2$disease_short)),
        base::names(disease_group_map)
      )
    )
  }

  if (!include_healthy) {
    meta2 <- meta2 %>%
      dplyr::filter(base::as.character(disease_short) %in% disease_order)
  }

  meta2$disease_short <- base::factor(meta2$disease_short, levels = disease_order)
  meta2$disease_group <- base::factor(meta2$disease_group, levels = group_levels)

  df <- meta2 %>%
    dplyr::mutate(
      .row  = base::trimws(base::as.character(.data[["disease_short"]])),
      .anno = base::trimws(base::as.character(.data[["disease_group"]])),
      .maj  = base::trimws(base::as.character(.data[["majorCluster"]])),
      .col  = base::trimws(base::as.character(.data[[col_state]]))
    ) %>%
    dplyr::filter(
      !base::is.na(.row), .row != "",
      !base::is.na(.anno), .anno != "",
      !base::is.na(.col), .col != "",
      .maj == major_name
    )

  if (base::nrow(df) == 0) {
    base::message("[SKIP] No cells for majorCluster = ", major_name)
    return(base::invisible(NULL))
  }

  state_total <- df %>%
    dplyr::count(.col, name = "n") %>%
    dplyr::arrange(dplyr::desc(n))

  keep_states <- state_total %>%
    dplyr::filter(n >= min_cells_per_state_total) %>%
    dplyr::pull(.col)

  df <- df %>% dplyr::filter(.col %in% keep_states)

  if (base::nrow(df) == 0) {
    base::message("[SKIP] After state filtering, no data left for ", major_name)
    return(base::invisible(NULL))
  }

  row_total <- df %>%
    dplyr::count(.row, name = "n") %>%
    dplyr::arrange(dplyr::desc(n))

  keep_rows <- row_total %>%
    dplyr::filter(n >= min_cells_per_row) %>%
    dplyr::pull(.row)

  df <- df %>% dplyr::filter(.row %in% keep_rows)

  if (base::nrow(df) == 0) {
    base::message("[SKIP] After row filtering, no data left for ", major_name)
    return(base::invisible(NULL))
  }

  keep_row_order <- disease_order[disease_order %in% base::unique(df$.row)]
  df$.row <- base::factor(df$.row, levels = keep_row_order)

  mat <- base::table(df$.row, df$.col)
  mat <- mat[base::rowSums(mat) > 0, base::colSums(mat) > 0, drop = FALSE]

  mat2 <- mat
  mat2[mat2 < min_cell_per_celltype_in_group] <- NA

  roe <- base::t(
    base::t(mat2 / base::rowSums(mat2, na.rm = TRUE)) /
      base::colSums(mat2, na.rm = TRUE) * base::sum(mat2, na.rm = TRUE)
  )

  roe[!base::is.finite(roe)] <- NA

  keep_c <- base::colSums(!base::is.na(roe)) > 0
  keep_r <- base::rowSums(!base::is.na(roe)) > 0
  roe <- roe[keep_r, keep_c, drop = FALSE]

  if (base::nrow(roe) < 1 || base::ncol(roe) < 1) {
    base::message("[SKIP] roe matrix empty after NA trimming for ", major_name)
    return(base::invisible(NULL))
  }

  map_df <- df %>%
    dplyr::count(.row, .anno, name = "n") %>%
    dplyr::group_by(.row) %>%
    dplyr::slice_max(order_by = n, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup()

  group_vec <- map_df$.anno
  base::names(group_vec) <- map_df$.row
  row_group <- group_vec[base::rownames(roe)]

  target_row_order <- disease_order[disease_order %in% base::rownames(roe)]
  roe <- roe[target_row_order, , drop = FALSE]
  row_group <- row_group[base::match(base::rownames(roe), base::names(row_group))]
  row_group <- base::factor(row_group, levels = group_levels)

  row_ha <- ComplexHeatmap::rowAnnotation(
    `Disease Group` = row_group,
    col = base::list(`Disease Group` = group_col),
    border = FALSE,
    simple_anno_size = grid::unit(3, "mm"),
    annotation_name_gp = grid::gpar(fontsize = 9),
    annotation_legend_param = base::list(
      title = "Disease Group",
      at = group_levels[group_levels %in% base::as.character(base::unique(row_group))],
      title_gp  = grid::gpar(fontsize = 9),
      labels_gp = grid::gpar(fontsize = 8),
      grid_height = grid::unit(3, "mm"),
      grid_width  = grid::unit(3, "mm")
    )
  )

  do_col_clust <- (base::ncol(roe) >= 2)
  col_hc <- NULL

  if (do_col_clust) {
    roe_for_clust <- roe
    roe_for_clust[base::is.na(roe_for_clust)] <- 0
    roe_for_clust <- base::pmin(roe_for_clust, cap_for_clust)

    col_hc <- stats::hclust(
      stats::dist(base::t(roe_for_clust)),
      method = "ward.D2"
    )
  }

  roe_col_fun <- circlize::colorRamp2(
    c(0, 1, 2, 3),
    c("#FFFFFF", "#F7C6C6", "#E07A7A", "#B2182B")
  )

  if (base::is.null(title)) {
    title <- base::paste0("Disease prevalence of ", major_name, " states")
  }
  if (base::is.null(out_pdf)) {
    out_pdf <- base::paste0("Fig1e_style_", base::gsub("[ /]", "_", major_name), "_states_newGroup.pdf")
  }

  ht <- ComplexHeatmap::Heatmap(
    matrix = roe,
    col = roe_col_fun,
    na_col = "#FFFFFF",
    border = TRUE,

    width  = grid::unit(base::ncol(roe) * cell_size_cm, "cm"),
    height = grid::unit(base::nrow(roe) * cell_size_cm, "cm"),

    right_annotation = row_ha,
    cluster_rows = FALSE,
    cluster_columns = if (do_col_clust) col_hc else FALSE,
    column_dend_height = if (do_col_clust) grid::unit(6, "mm") else grid::unit(0, "mm"),
    row_names_gp = grid::gpar(fontsize = row_font),
    column_names_gp = grid::gpar(fontsize = col_font),
    column_names_rot = column_names_rot,
    column_title = title,
    column_title_gp = grid::gpar(fontsize = title_font),
    cell_fun = function(j, i, x, y, width, height, fill) {
      v <- roe[i, j]
      if (!base::is.na(v)) {
        if (v > 3) grid::grid.text("+++", x, y, gp = grid::gpar(fontsize = plus_font))
        else if (v > 2) grid::grid.text("++", x, y, gp = grid::gpar(fontsize = plus_font))
        else if (v > 1) grid::grid.text("+", x, y, gp = grid::gpar(fontsize = plus_font))
      }
    },
    heatmap_legend_param = base::list(
      title = expression(R[oe]),
      at = c(0, 1, 2, 3),
      labels = c("0", "1 (+)", "2 (++)", "3 (+++)"),
      title_gp = grid::gpar(fontsize = 9),
      labels_gp = grid::gpar(fontsize = 8),
      grid_width = grid::unit(4, "mm"),
      legend_height = grid::unit(25, "mm")
    )
  )

  pdf_w <- base::ncol(roe) * cell_size_cm + pdf_margin_w
  pdf_h <- base::nrow(roe) * cell_size_cm + pdf_margin_h

  grDevices::pdf(out_pdf, width = pdf_w, height = pdf_h)
  ComplexHeatmap::draw(
    ht,
    heatmap_legend_side = "right",
    annotation_legend_side = "right",
    padding = grid::unit(c(6, 32, 6, 6), "mm")
  )
  grDevices::dev.off()
  if (verbose) {
    base::cat("[DONE] ", major_name, " -> ", out_pdf, "\n", sep = "")
  }

  return(base::invisible(base::list(
    roe = roe,
    out_pdf = out_pdf,
    row_group = row_group,
    state_total = state_total,
    row_total = row_total
  )))
}


meta <- seurat_obj@meta.data

major_df <- meta %>%
  dplyr::mutate(majorCluster = base::trimws(base::as.character(majorCluster))) %>%
  dplyr::filter(!base::is.na(majorCluster), majorCluster != "") %>%
  dplyr::count(majorCluster, name = "n") %>%
  dplyr::arrange(dplyr::desc(n))

major_list <- major_df %>%
  dplyr::filter(n >= 50) %>%
  dplyr::pull(majorCluster)

base::cat("Major clusters to plot:\n")
base::print(major_df)
base::print(major_list)

res_list <- base::lapply(major_list, function(mj) {
  plot_roe_heatmap_one_major(
    meta = meta,
    major_name = mj,
    row_var_raw = "disease_Level_3",
    col_state = "subCluster",
    min_cells_per_row = 30,
    min_cells_per_state_total = 30,
    min_cell_per_celltype_in_group = 3,
    include_healthy = FALSE,
    out_pdf = base::paste0("Fig1e_style_", base::gsub("[ /]", "_", mj), "_states_newGroup.pdf"),
    title = base::paste0("Disease prevalence of ", mj, " states"),
    verbose = TRUE
  )
})

base::names(res_list) <- major_list
