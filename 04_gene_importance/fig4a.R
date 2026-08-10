suppressPackageStartupMessages({
  library(data.table)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

# =========================================================
# 1. 输入输出
# =========================================================
base_dir <- "data/analysis/08_gene_importance"
out_dir  <- file.path(base_dir, "48_Fig3a_pretty_R_separate")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

candidate_dirs <- c(
  file.path(base_dir, "44_fig3a_confusion_v4"),
  file.path(base_dir, "44_fig3a_confusion_v2"),
  file.path(base_dir, "44_fig3a_confusion")
)

find_first_existing <- function(fname, dirs) {
  for (d in dirs) {
    fp <- file.path(d, fname)
    if (file.exists(fp)) return(fp)
  }
  # 如果上面没找到，就全局搜一下
  hits <- Sys.glob(file.path(base_dir, "**", fname))
  if (length(hits) > 0) return(hits[1])
  stop("Cannot find file: ", fname)
}

corrected_tsv <- find_first_existing("confusion_matrix_corrected_colnorm.tsv", candidate_dirs)
raw_tsv       <- find_first_existing("confusion_matrix_unintegrated_colnorm.tsv", candidate_dirs)

cat("[INFO] corrected_tsv   =", corrected_tsv, "\n")
cat("[INFO] raw_tsv         =", raw_tsv, "\n")

# =========================================================
# 2. 疾病短名 / 顺序 / 分组
# =========================================================
disease_short_map <- c(
  "Pediatric healthy" = "Healthy",
  "Dengue human infection(DENV-1)" = "Dengue",
  "HFMD(EV-A71)enterovirus" = "HFMD",
  "Influenza A virus (IAV)" = "Flu",
  "severe acute respiratory syndrome coronavirus 2 (SARS-CoV-2)" = "COVID",
  "bacterial infection" = "Bact.inf",
  "Pediatric Sepsis" = "Sepsis",
  "Kawasaki disease (KD)" = "KD",
  "Multisystem inflammatory syndrome in children (MIS-C)" = "MIS-C",
  "human systemic lupus erythematosus (SLE)" = "SLE",
  "juvenile idopathic arthritis (JIA)" = "JIA"
)

disease_order <- c(
  "Healthy",
  "Dengue", "HFMD", "Flu", "COVID",
  "Bact.inf", "Sepsis",
  "KD", "MIS-C",
  "SLE", "JIA"
)

disease_group_map <- c(
  "Healthy"  = "Healthy",
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

group_order <- c("Healthy", "Autoimmune", "Sepsis", "bacteria", "virus", "Hyperinflam")

disease_group_cols <- c(
  "virus"       = "#9DB2E5",
  "bacteria"    = "#C98FB1",
  "Sepsis"      = "#FF1F2D",
  "Hyperinflam" = "#E9A826",
  "Autoimmune"  = "#2F9E8F",
  "Healthy"     = "#808080"
)

# =========================================================
# 3. 图形参数
# =========================================================
pdf_width  <- 5.8
pdf_height <- 4.9

row_name_size   <- 10
col_name_size   <- 9
title_size      <- 12
legend_title_sz <- 10
legend_label_sz <- 9

col_fun <- circlize::colorRamp2(
  c(0, 0.2, 0.5, 0.8, 1.0),
  c("#3B4CC0", "#A7C7E7", "#F2F2E8", "#F4A261", "#B40426")
)

# =========================================================
# 4. 工具函数
# =========================================================
clean_names <- function(x) {
  trimws(as.character(x))
}

read_cm <- function(fp) {
  x <- fread(fp, data.table = FALSE)
  rn <- x[[1]]
  x[[1]] <- NULL
  mat <- as.matrix(x)
  rownames(mat) <- clean_names(rn)
  colnames(mat) <- clean_names(colnames(mat))
  storage.mode(mat) <- "numeric"
  mat
}

pretty_mat <- function(mat) {
  rownames(mat) <- ifelse(
    rownames(mat) %in% names(disease_short_map),
    disease_short_map[rownames(mat)],
    rownames(mat)
  )
  colnames(mat) <- ifelse(
    colnames(mat) %in% names(disease_short_map),
    disease_short_map[colnames(mat)],
    colnames(mat)
  )

  full <- matrix(
    NA_real_,
    nrow = length(disease_order),
    ncol = length(disease_order),
    dimnames = list(disease_order, disease_order)
  )

  rr <- base::intersect(rownames(mat), disease_order)
  cc <- base::intersect(colnames(mat), disease_order)
  full[rr, cc] <- mat[rr, cc, drop = FALSE]
  full
}

make_group_bar <- function(row_names) {
  grp <- factor(disease_group_map[row_names], levels = group_order)

  rowAnnotation(
    DiseaseGroup = grp,
    col = list(DiseaseGroup = disease_group_cols),
    show_annotation_name = FALSE,
    width = unit(3.2, "mm"),
    border = FALSE,
    annotation_legend_param = list(
      title = "Disease group:",
      title_gp = gpar(fontsize = legend_title_sz, fontface = "bold"),
      labels_gp = gpar(fontsize = legend_label_sz)
    )
  )
}

make_heatmap <- function(mat, title_text) {
  row_grp <- factor(disease_group_map[rownames(mat)], levels = group_order)
  col_grp <- factor(disease_group_map[colnames(mat)], levels = group_order)

  ht <- Heatmap(
    mat,
    name = "Proportion\n(col.norm.)",
    col = col_fun,
    na_col = "#ECECEC",

    cluster_rows = FALSE,
    cluster_columns = FALSE,

    row_split = row_grp,
    column_split = col_grp,

    row_gap = unit(1.2, "mm"),
    column_gap = unit(1.2, "mm"),

    rect_gp = gpar(col = "#D9D9D9", lwd = 0.6),
    border = FALSE,

    row_names_side = "left",
    row_names_gp = gpar(fontsize = row_name_size),
    column_names_rot = 90,
    column_names_gp = gpar(fontsize = col_name_size),

    row_title = "Predicted conditions",
    column_title = title_text,
    column_title_gp = gpar(fontsize = title_size, fontface = "bold"),
    row_title_gp = gpar(fontsize = 11),

    heatmap_legend_param = list(
      at = c(0, 0.5, 1.0),
      labels = c("0", "0.5", "1.0"),
      title_gp = gpar(fontsize = legend_title_sz, fontface = "bold"),
      labels_gp = gpar(fontsize = legend_label_sz),
      grid_width = unit(4.5, "mm"),
      grid_height = unit(4.5, "mm")
    )
  )

  make_group_bar(rownames(mat)) + ht
}

save_one_panel <- function(ht_obj, out_pdf, out_png, letter = "a") {
  pdf(out_pdf, width = pdf_width, height = pdf_height)
  ComplexHeatmap::draw(
    ht_obj,
    heatmap_legend_side = "right",
    annotation_legend_side = "right",
    merge_legends = FALSE,
    padding = unit(c(3, 4, 3, 4), "mm")
  )
  grid.text(
    letter,
    x = unit(2, "mm"),
    y = unit(1, "npc") - unit(2, "mm"),
    just = c("left", "top"),
    gp = gpar(fontsize = 16, fontface = "bold")
  )
  grid.text(
    "True conditions",
    x = unit(0.48, "npc"),
    y = unit(2.2, "mm"),
    just = c("center", "bottom"),
    gp = gpar(fontsize = 11)
  )
  dev.off()

  png(out_png, width = pdf_width, height = pdf_height, units = "in", res = 500)
  ComplexHeatmap::draw(
    ht_obj,
    heatmap_legend_side = "right",
    annotation_legend_side = "right",
    merge_legends = FALSE,
    padding = unit(c(3, 4, 3, 4), "mm")
  )
  grid.text(
    letter,
    x = unit(2, "mm"),
    y = unit(1, "npc") - unit(2, "mm"),
    just = c("left", "top"),
    gp = gpar(fontsize = 16, fontface = "bold")
  )
  grid.text(
    "True conditions",
    x = unit(0.48, "npc"),
    y = unit(2.2, "mm"),
    just = c("center", "bottom"),
    gp = gpar(fontsize = 11)
  )
  dev.off()
}

# =========================================================
# 5. 读取矩阵
# =========================================================
mat_corr <- pretty_mat(read_cm(corrected_tsv))
mat_raw  <- pretty_mat(read_cm(raw_tsv))

# =========================================================
# 6. 构建两张图
# =========================================================
ht_corr <- make_heatmap(
  mat_corr,
  "XGBoost trained on\nscANVI-corrected expression"
)

ht_raw <- make_heatmap(
  mat_raw,
  "XGBoost trained on\nnon-corrected expression"
)

# =========================================================
# 7. 分开保存
# =========================================================
save_one_panel(
  ht_corr,
  file.path(out_dir, "Fig3a_corrected_pretty_R.pdf"),
  file.path(out_dir, "Fig3a_corrected_pretty_R.png"),
  letter = "a"
)

save_one_panel(
  ht_raw,
  file.path(out_dir, "Fig3a_uncorrected_pretty_R.pdf"),
  file.path(out_dir, "Fig3a_uncorrected_pretty_R.png"),
  letter = "a"
)

write.table(
  mat_corr,
  file = file.path(out_dir, "Fig3a_corrected_matrix_pretty.tsv"),
  sep = "\t", quote = FALSE, col.names = NA
)

write.table(
  mat_raw,
  file = file.path(out_dir, "Fig3a_uncorrected_matrix_pretty.tsv"),
  sep = "\t", quote = FALSE, col.names = NA
)

cat("[DONE]\n")
cat("Output dir:", out_dir, "\n")
