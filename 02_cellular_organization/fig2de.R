library(CoVarNet)
library(sceasy)
library(Matrix)
library(SummarizedExperiment)
library(data.table)
library(dplyr)
library(igraph)
library(Seurat)
library(tidyverse)
library(stringr)
library(magick)
library(Metrics)
library(circlize)
library(filesstrings)



setwd("data/analysis/atlas_objects")

library(SeuratDisk)

#Convert("data/analysis/scANVI/08_downstream_object/mainObject_unionGenes_dropBadCT.h5ad", dest = "h5seurat", overwrite = TRUE)
#跑个1.py

seurat_obj <- LoadH5Seurat("data/analysis/scANVI/08_downstream_object/mainObject_unionGenes_dropBadCT.h5seurat", assays = "RNA")

cat("Assays:\n")
print(Assays(seurat_obj))

cat("\nDefaultAssay:\n")
print(DefaultAssay(seurat_obj))

rna_counts <- GetAssayData(seurat_obj, assay = "RNA", slot = "counts")
rna_data   <- GetAssayData(seurat_obj, assay = "RNA", slot = "data")

cat("\ncounts [1:5,1:5]:\n")
print(as.matrix(rna_counts[1:5, 1:5]))

cat("\ndata [1:5,1:5]:\n")
print(as.matrix(rna_data[1:5, 1:5]))

vals <- as.vector(rna_counts[1:min(100, nrow(rna_counts)), 1:min(20, ncol(rna_counts))])
vals <- vals[!is.na(vals)]

cat("\ncounts integer-like:\n")
print(all(abs(vals - round(vals)) < 1e-8))

cat("\ncounts summary:\n")
print(summary(vals))

# 假设 seurat_obj 是你的 Seurat 对象
# 修改列名
colnames(seurat_obj@meta.data) <- gsub("sample_id", "sampleID", colnames(seurat_obj@meta.data))
colnames(seurat_obj@meta.data) <- gsub("aifi_L1_paper", "majorCluster", colnames(seurat_obj@meta.data))
colnames(seurat_obj@meta.data) <- gsub("aifi_L2_paper", "subCluster", colnames(seurat_obj@meta.data))



# 确保 subCluster 和 majorCluster 一一对应，且匹配正确
meta <- seurat_obj@meta.data


# 打印对应关系，检查是否正确匹配
table(meta$subCluster, meta$majorCluster)



table(meta$subCluster)  # 检查 subCluster 的类别数量

# ============================================================
# Step 1. Rebuild and save NMF_K4 for cellular modules
# ============================================================

setwd("data/analysis/atlas_objects")

suppressPackageStartupMessages({
  library(Seurat)
  library(CoVarNet)
  library(NMF)
  library(dplyr)
})

# 读对象
seurat_obj <- readRDS("integrated_pediatric.rds")
meta <- as.data.frame(seurat_obj@meta.data)

# 检查列名
print(colnames(meta))
stopifnot("sampleID" %in% colnames(meta))
stopifnot("majorCluster" %in% colnames(meta))
stopifnot("subCluster" %in% colnames(meta))
stopifnot("disease_Level_3" %in% colnames(meta))

# 清洗空格
meta$sampleID        <- trimws(as.character(meta$sampleID))
meta$majorCluster    <- trimws(as.character(meta$majorCluster))
meta$subCluster      <- trimws(as.character(meta$subCluster))
meta$disease_Level_3 <- trimws(as.character(meta$disease_Level_3))

# 计算 sample × Level 2 cell population abundance
mat_fq_raw  <- CoVarNet::freq_calculate(meta)
mat_fq_norm <- CoVarNet::freq_normalize(mat_fq_raw, normalize = "minmax")

cat("mat_fq_raw dim:\n")
print(dim(mat_fq_raw))

cat("mat_fq_norm dim:\n")
print(dim(mat_fq_norm))

# 跑 NMF，K = 4
set.seed(77)
K <- 4

NMF_K12 <- NMF::nmf(
  mat_fq_norm,
  rank = K,
  method = "nsNMF",
  seed = rep(77, 6),
  nrun = 30
)

# 统一命名为 CM01-CM04
colnames(NMF::basis(NMF_K12)) <- sprintf("CM%02d", 1:K)
rownames(NMF::coef(NMF_K12))  <- sprintf("CM%02d", 1:K)

# 保存，后面直接 load
saveRDS(NMF_K12, "NMF_K4_cellular_modules_CM01_CM04.rds")
saveRDS(mat_fq_raw, "CM_input_mat_fq_raw.rds")
saveRDS(mat_fq_norm, "CM_input_mat_fq_norm_minmax.rds")

cat("[DONE] Saved:\n")
cat("NMF_K4_cellular_modules_CM01_CM04.rds\n")
cat("CM_input_mat_fq_raw.rds\n")
cat("CM_input_mat_fq_norm_minmax.rds\n")







# ============================================================
# Fig. 2d/e/f for cellular modules, K = 4
# d: sample-level CM abundance heatmap
# e: disease-level CM Ro/e heatmap
# f: top Level 2 populations defining each CM
# ============================================================

setwd("data/analysis/atlas_objects")

suppressPackageStartupMessages({
  library(Seurat)
  library(CoVarNet)
  library(NMF)
  library(dplyr)
  library(tidyr)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

# ============================================================
# 0. Load object
# ============================================================

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

stopifnot("sampleID" %in% colnames(meta))
stopifnot("disease_Level_3" %in% colnames(meta))
stopifnot("majorCluster" %in% colnames(meta))
stopifnot("subCluster" %in% colnames(meta))

meta$sampleID        <- trimws(as.character(meta$sampleID))
meta$disease_Level_3 <- trimws(as.character(meta$disease_Level_3))
meta$majorCluster    <- trimws(as.character(meta$majorCluster))
meta$subCluster      <- trimws(as.character(meta$subCluster))

outdir <- file.path(getwd(), "fig2_redraw2")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 1. Disease and color settings
# ============================================================

INCLUDE_HEALTHY <- FALSE
USE_CLUSTER_FOR_E <- FALSE

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
  "severe acute respiratory syndrome coronavirus 2 (SARS-CoV-2)" = "COVID",
  "Healthy children" = "Pediatric healthy",
  "Pediatric healthy" = "Pediatric healthy"
)

if (INCLUDE_HEALTHY) {
  disease_order <- c(
    "Dengue", "HFMD", "Flu", "COVID",
    "Bact.inf", "Sepsis", "KD", "MIS-C",
    "Pediatric healthy", "SLE", "JIA"
  )
} else {
  disease_order <- c(
    "Dengue", "HFMD", "Flu", "COVID",
    "Bact.inf", "Sepsis", "KD", "MIS-C", "SLE", "JIA"
  )
}

disease_group_map <- c(
  "Dengue" = "virus",
  "HFMD" = "virus",
  "Flu" = "virus",
  "COVID" = "virus",
  "Bact.inf" = "bacteria",
  "Sepsis" = "Sepsis",
  "KD" = "Hyperinflam",
  "MIS-C" = "Hyperinflam",
  "Pediatric healthy" = "Healthy",
  "SLE" = "Autoimmune",
  "JIA" = "Autoimmune"
)

if (INCLUDE_HEALTHY) {
  group_levels <- c("virus", "bacteria", "Sepsis", "Hyperinflam", "Healthy", "Autoimmune")
} else {
  group_levels <- c("virus", "bacteria", "Sepsis", "Hyperinflam", "Autoimmune")
}

group_col <- c(
  "virus"       = "#6699CC",
  "bacteria"    = "#C48EAD",
  "Sepsis"      = "#FF2A2A",
  "Hyperinflam" = "#E3A722",
  "Healthy"     = "#BDBDBD",
  "Autoimmune"  = "#36A596"
)

cm_col <- c(
  "CM01" = "#F26B6C",
  "CM02" = "#F4A261",
  "CM03" = "#4CB7A5",
  "CM04" = "#A78BDA"
)

celltype_col <- c(
  "B"                = "#4C78A8",
  "Plasma"           = "#9ECAE1",
  "Mono"             = "#59A14F",
  "DC"               = "#F28E2B",
  "pDC"              = "#E15759",
  "T CD4 naive"      = "#A1D99B",
  "T CD4 non-naive"  = "#E15759",
  "T CD8 naive"      = "#F4A261",
  "T CD8 non-naive"  = "#8E63CE",
  "UTC"              = "#76B7B2",
  "ILC"              = "#C77CFF",
  "Other"            = "#BDBDBD"
)

# ============================================================
# 2. Prepare metadata for K=4 NMF
# ============================================================

meta_mod <- meta %>%
  dplyr::mutate(
    disease_raw = as.character(disease_Level_3),
    disease = dplyr::recode(disease_raw, !!!disease_short_map, .default = disease_raw),
    disease_group = dplyr::recode(disease, !!!disease_group_map, .default = "Unknown")
  ) %>%
  dplyr::filter(
    !is.na(sampleID), sampleID != "",
    !is.na(subCluster), subCluster != "",
    !is.na(majorCluster), majorCluster != "",
    disease %in% disease_order
  )

# 和你原来模块分析一致，去掉非主分析类型
drop_ct <- c("Erythrocyte", "Platelet")
meta_mod <- meta_mod %>%
  dplyr::filter(!(majorCluster %in% drop_ct))

cat("Samples used:\n")
print(length(unique(meta_mod$sampleID)))
cat("Disease counts:\n")
print(table(meta_mod$disease))
cat("Major cluster counts:\n")
print(table(meta_mod$majorCluster))

is_fac <- vapply(meta_mod, is.factor, logical(1))
meta_mod[is_fac] <- lapply(meta_mod[is_fac], as.character)

# ============================================================
# 3. Load or compute K=4 NMF
# ============================================================

nmf_rds <- file.path(outdir, "NMF_K4_cellular_modules_CM01_CM04.rds")
mat_raw_rds <- file.path(outdir, "CM_input_mat_fq_raw_K4.rds")
mat_norm_rds <- file.path(outdir, "CM_input_mat_fq_norm_minmax_K4.rds")

if (file.exists(nmf_rds)) {
  message("[READ] ", nmf_rds)
  NMF_K4 <- readRDS(nmf_rds)
} else {
  message("[RUN] CoVarNet freq_calculate + freq_normalize + NMF K=4")

  mat_fq_raw <- CoVarNet::freq_calculate(meta_mod)
  mat_fq_norm <- CoVarNet::freq_normalize(mat_fq_raw, normalize = "minmax")

  cat("mat_fq_raw dim:\n")
  print(dim(mat_fq_raw))
  cat("mat_fq_norm dim:\n")
  print(dim(mat_fq_norm))

  set.seed(77)
  NMF_K4 <- NMF::nmf(
    mat_fq_norm,
    rank = 4,
    method = "nsNMF",
    seed = rep(77, 6),
    nrun = 30
  )

  colnames(NMF::basis(NMF_K4)) <- sprintf("CM%02d", 1:4)
  rownames(NMF::coef(NMF_K4))  <- sprintf("CM%02d", 1:4)

  saveRDS(NMF_K4, nmf_rds)
  saveRDS(mat_fq_raw, mat_raw_rds)
  saveRDS(mat_fq_norm, mat_norm_rds)
}

# ============================================================
# 4. Extract h matrix: CM x sample
# ============================================================

h <- tryCatch(
  NMF::scoef(NMF_K4),
  error = function(e) NMF::coef(NMF_K4)
)

h <- as.matrix(h)
rownames(h) <- sprintf("CM%02d", 1:4)
h <- h[rownames(h), , drop = FALSE]

meta_sp <- meta_mod %>%
  dplyr::select(sampleID, disease, disease_group) %>%
  dplyr::distinct(sampleID, .keep_all = TRUE)

rownames(meta_sp) <- meta_sp$sampleID

common_samples <- intersect(colnames(h), rownames(meta_sp))
h <- h[, common_samples, drop = FALSE]
meta_sp <- meta_sp[common_samples, , drop = FALSE]

stopifnot(all(colnames(h) == rownames(meta_sp)))

meta_sp$disease <- factor(meta_sp$disease, levels = disease_order)
meta_sp$disease_group <- factor(meta_sp$disease_group, levels = group_levels)

cat("h dim:\n")
print(dim(h))
cat("Column sum summary:\n")
print(summary(colSums(h)))

# ============================================================
# Fig. 2d: sample-level CM abundance heatmap
# ============================================================

dominant_cm <- apply(h, 2, function(x) rownames(h)[which.max(x)])
meta_sp$dominant_CM <- dominant_cm[rownames(meta_sp)]
meta_sp$dominant_CM <- factor(meta_sp$dominant_CM, levels = rownames(h))

# 排序：先疾病大类，再疾病，再 dominant CM，再 CM01 abundance
sample_order <- meta_sp %>%
  dplyr::mutate(
    disease_group = factor(as.character(disease_group), levels = group_levels),
    disease = factor(as.character(disease), levels = disease_order)
  ) %>%
  dplyr::arrange(disease_group, disease, dominant_CM, dplyr::desc(h["CM01", rownames(.)])) %>%
  rownames()

h_order <- h[, sample_order, drop = FALSE]
meta_order <- meta_sp[sample_order, , drop = FALSE]

# 行内 z score，展示样本间模块偏向
h_z <- t(scale(t(h_order)))
h_z[!is.finite(h_z)] <- 0
h_z <- pmax(pmin(h_z, 2), -2)

ha_d <- HeatmapAnnotation(
  `Disease group` = meta_order$disease_group,
  `Dominant CM` = meta_order$dominant_CM,
  col = list(
    `Disease group` = group_col,
    `Dominant CM` = cm_col
  ),
  border = FALSE,
  simple_anno_size = unit(3, "mm"),
  annotation_name_gp = gpar(fontsize = 8),
  annotation_legend_param = list(
    `Disease group` = list(
      title = "Disease group",
      title_gp = gpar(fontsize = 8),
      labels_gp = gpar(fontsize = 7),
      grid_width = unit(3, "mm"),
      grid_height = unit(3, "mm")
    ),
    `Dominant CM` = list(
      title = "Dominant CM",
      title_gp = gpar(fontsize = 8),
      labels_gp = gpar(fontsize = 7),
      grid_width = unit(3, "mm"),
      grid_height = unit(3, "mm")
    )
  )
)

z_col_fun <- colorRamp2(
  c(-2, 0, 2),
  c("#4A78C2", "#FFFFFF", "#D73027")
)

row_gp_d <- gpar(
  fontsize = 9,
  fontface = "bold",
  col = cm_col[rownames(h_z)]
)

ht_d <- Heatmap(
  h_z,
  name = "CM abundance\nz score",
  col = z_col_fun,
  na_col = "#FFFFFF",
  border = TRUE,
  rect_gp = gpar(col = "white", lwd = 0.3),
  top_annotation = ha_d,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_column_names = FALSE,
  row_names_gp = row_gp_d,
  column_title = "Cellular module abundance across samples",
  column_title_gp = gpar(fontsize = 10, fontface = "bold"),
  width = unit(7.4, "cm"),
  height = unit(2.4, "cm"),
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 8),
    labels_gp = gpar(fontsize = 7),
    at = c(-2, -1, 0, 1, 2),
    legend_height = unit(22, "mm")
  )
)

pdf_file_d <- file.path(outdir, "Fig2d_sample_level_CM_abundance_heatmap_K4.pdf")
png_file_d <- file.path(outdir, "Fig2d_sample_level_CM_abundance_heatmap_K4.png")

pdf(pdf_file_d, width = 7.0, height = 3.2)
draw(ht_d, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

png(png_file_d, width = 7.0, height = 3.2, units = "in", res = 600)
draw(ht_d, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

# ============================================================
# Fig. 2e: disease-level CM Ro/e heatmap
# ============================================================

cm_names <- rownames(h)

sample_cm_df <- as.data.frame(t(h)) %>%
  tibble::rownames_to_column("sampleID") %>%
  dplyr::left_join(
    meta_sp %>%
      dplyr::mutate(sampleID = rownames(meta_sp)) %>%
      dplyr::select(sampleID, disease, disease_group),
    by = "sampleID"
  ) %>%
  dplyr::filter(!is.na(disease), disease %in% disease_order)

disease_cm <- sample_cm_df %>%
  dplyr::group_by(disease) %>%
  dplyr::summarise(
    dplyr::across(dplyr::all_of(cm_names), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  as.data.frame()

rownames(disease_cm) <- disease_cm$disease
disease_cm$disease <- NULL

disease_cm <- disease_cm[
  disease_order[disease_order %in% rownames(disease_cm)],
  cm_names,
  drop = FALSE
]

# Fig. 2a 同款 Ro/e 逻辑
mat <- as.matrix(disease_cm)
mat_prop <- mat / rowSums(mat, na.rm = TRUE)

roe_disease_by_cm <- t(
  t(mat_prop / rowSums(mat_prop, na.rm = TRUE)) /
    colSums(mat_prop, na.rm = TRUE) * sum(mat_prop, na.rm = TRUE)
)

roe_disease_by_cm[!is.finite(roe_disease_by_cm)] <- NA
roe_plot <- t(roe_disease_by_cm)

roe_plot <- roe_plot[
  cm_names,
  disease_order[disease_order %in% colnames(roe_plot)],
  drop = FALSE
]

roe_col_fun <- colorRamp2(
  c(0, 1, 2, 3),
  c("#FFFFFF", "#F7C6C6", "#E07A7A", "#B2182B")
)

col_group <- disease_group_map[colnames(roe_plot)]
col_group <- factor(col_group, levels = group_levels)

ha_e <- HeatmapAnnotation(
  `Disease group` = col_group,
  col = list(`Disease group` = group_col),
  border = FALSE,
  simple_anno_size = unit(3, "mm"),
  annotation_name_gp = gpar(fontsize = 8),
  annotation_legend_param = list(
    title = "Disease group",
    title_gp = gpar(fontsize = 8),
    labels_gp = gpar(fontsize = 7),
    grid_width = unit(3, "mm"),
    grid_height = unit(3, "mm")
  )
)

row_gp_e <- gpar(
  fontsize = 9,
  fontface = "bold",
  col = cm_col[rownames(roe_plot)]
)

ht_e <- Heatmap(
  roe_plot,
  name = "Ro/e",
  col = roe_col_fun,
  na_col = "#FFFFFF",
  border = TRUE,
  rect_gp = gpar(col = "white", lwd = 0.8),
  top_annotation = ha_e,
  cluster_rows = USE_CLUSTER_FOR_E,
  cluster_columns = USE_CLUSTER_FOR_E,
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  column_dend_height = if (USE_CLUSTER_FOR_E) unit(5, "mm") else unit(0, "mm"),
  row_dend_width = if (USE_CLUSTER_FOR_E) unit(5, "mm") else unit(0, "mm"),
  row_names_gp = row_gp_e,
  column_names_gp = gpar(fontsize = 8),
  column_names_rot = 45,
  column_title = "Disease prevalence of cellular modules",
  column_title_gp = gpar(fontsize = 10, fontface = "bold"),
  width = unit(ncol(roe_plot) * 0.55, "cm"),
  height = unit(nrow(roe_plot) * 0.55, "cm"),
  cell_fun = function(j, i, x, y, width, height, fill) {
    v <- roe_plot[i, j]
    if (!is.na(v)) {
      if (v > 3) {
        grid.text("+++", x, y, gp = gpar(fontsize = 7))
      } else if (v > 2) {
        grid.text("++", x, y, gp = gpar(fontsize = 7))
      } else if (v > 1) {
        grid.text("+", x, y, gp = gpar(fontsize = 7))
      }
    }
  },
  heatmap_legend_param = list(
    title = expression(R[oe]),
    at = c(0, 1, 2, 3),
    labels = c("0", "1 (+)", "2 (++)", "3 (+++)"),
    title_gp = gpar(fontsize = 8),
    labels_gp = gpar(fontsize = 7),
    grid_width = unit(4, "mm"),
    legend_height = unit(22, "mm")
  )
)

pdf_file_e <- file.path(outdir, "Fig2e_disease_prevalence_CM_RoE_K4.pdf")
png_file_e <- file.path(outdir, "Fig2e_disease_prevalence_CM_RoE_K4.png")

pdf(pdf_file_e, width = 6.7, height = 3.5)
draw(ht_e, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

png(png_file_e, width = 6.7, height = 3.5, units = "in", res = 600)
draw(ht_e, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

# ============================================================
# Fig. 2f: top Level 2 populations defining each CM
# ============================================================

w <- NMF::basis(NMF_K4)
w <- as.matrix(w)
colnames(w) <- sprintf("CM%02d", 1:4)

# 如果 basis 的列是 CM，行是 Level 2 population
if (ncol(w) != 4 && nrow(w) == 4) {
  w <- t(w)
  colnames(w) <- sprintf("CM%02d", 1:4)
}

cat("basis matrix dim:\n")
print(dim(w))
cat("basis row examples:\n")
print(rownames(w)[1:min(10, nrow(w))])

# subCluster -> majorCluster majority map
sub_to_major <- meta_mod %>%
  dplyr::count(subCluster, majorCluster, name = "n") %>%
  dplyr::group_by(subCluster) %>%
  dplyr::slice_max(order_by = n, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::select(subCluster, majorCluster)

top_n <- 10

top_df <- as.data.frame(w) %>%
  tibble::rownames_to_column("subCluster") %>%
  tidyr::pivot_longer(
    cols = starts_with("CM"),
    names_to = "CM",
    values_to = "weight"
  ) %>%
  dplyr::group_by(CM) %>%
  dplyr::slice_max(order_by = weight, n = top_n, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(sub_to_major, by = "subCluster") %>%
  dplyr::mutate(
    majorCluster = ifelse(is.na(majorCluster), "Other", majorCluster),
    CM = factor(CM, levels = sprintf("CM%02d", 1:4))
  ) %>%
  dplyr::group_by(CM) %>%
  dplyr::arrange(weight, .by_group = TRUE) %>%
  dplyr::mutate(
    feature_id = paste(CM, subCluster, sep = "___"),
    feature_id = factor(feature_id, levels = unique(feature_id))
  ) %>%
  dplyr::ungroup()

# 防止有些 majorCluster 没颜色
missing_major <- setdiff(unique(top_df$majorCluster), names(celltype_col))
if (length(missing_major) > 0) {
  extra_cols <- grDevices::hcl.colors(length(missing_major), "Set 2")
  names(extra_cols) <- missing_major
  celltype_col <- c(celltype_col, extra_cols)
}

p_f <- ggplot(top_df, aes(x = weight, y = feature_id, fill = majorCluster)) +
  geom_col(width = 0.72, color = "white", linewidth = 0.25) +
  facet_wrap(~ CM, scales = "free_y", nrow = 1) +
  scale_y_discrete(labels = function(x) sub("^.*___", "", x)) +
  scale_fill_manual(values = celltype_col, name = "Parent Level 1 compartment") +
  labs(
    title = "Top Level 2 populations defining each cellular module",
    x = "Module weight",
    y = NULL
  ) +
  theme_classic(base_size = 9) +
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 10),
    axis.text.y = element_text(size = 6.5),
    axis.text.x = element_text(size = 7),
    axis.title.x = element_text(size = 9),
    axis.line = element_line(linewidth = 0.35),
    axis.ticks = element_line(linewidth = 0.3),
    legend.position = "bottom",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    legend.key.size = unit(3.2, "mm"),
    panel.spacing.x = unit(5, "mm")
  ) +
  guides(fill = guide_legend(nrow = 1))

pdf_file_f <- file.path(outdir, "Fig2f_top_Level2_populations_defining_CM_K4.pdf")
png_file_f <- file.path(outdir, "Fig2f_top_Level2_populations_defining_CM_K4.png")

ggsave(pdf_file_f, p_f, width = 8.6, height = 3.4, device = cairo_pdf)
ggsave(png_file_f, p_f, width = 8.6, height = 3.4, dpi = 600)

# ============================================================
# Save data
# ============================================================

write.csv(
  as.data.frame(h_order),
  file.path(outdir, "Fig2d_sample_level_CM_abundance_matrix_ordered.csv")
)

write.csv(
  roe_plot,
  file.path(outdir, "Fig2e_CM_RoE_matrix.csv")
)

write.csv(
  top_df,
  file.path(outdir, "Fig2f_top_Level2_populations_defining_CM.csv"),
  row.names = FALSE
)

message("[DONE] Saved Fig. 2d/e/f to: ", outdir)
