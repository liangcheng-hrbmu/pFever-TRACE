





setwd("data/analysis/inflammation_signatures/14_LMEM_Corr_L1")



suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(RColorBrewer)
})

# =========================================================
# 1. 输入输出
# =========================================================
all_path <- "data/analysis/inflammation_signatures/14_LMEM_Corr_L1_v2/LMEM_L1_all_results.tsv"

out_dir <- "data/analysis/inflammation_signatures/14_LMEM_Corr_L1_v2/LMEM_L1_reference_style_heatmap_narrow_asym_final"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# =========================================================
# 2. 参数
# =========================================================
clip_quantile <- 0.90

# 收窄“宽度”
pdf_width  <- 10.8
pdf_height <- 6

sig_point_fontsize <- 8
row_name_size <- 9.0
col_name_size <- 5.0

# 用原始 p_value 判显著
p_cutoff <- 0.05

# 只保留“至少一个疾病中显著”的因子
keep_only_any_sig <- TRUE

# =========================================================
# 3. 疾病短名与顺序
# =========================================================
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

disease_order <- c(
  "Dengue", "HFMD", "Flu", "COVID",
  "Bact.inf", "Sepsis", "KD", "MIS-C", "SLE", "JIA"
)

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

row_group_order <- c("virus", "bacteria", "Sepsis", "Hyperinflam", "Autoimmune")

row_group_cols <- c(
  "virus" = "#4F81BD",
  "bacteria" = "#2C7FB8",
  "Sepsis" = "#1F4E79",
  "Hyperinflam" = "#F4A261",
  "Autoimmune" = "#C1121F"
)

# =========================================================
# 4. cell type 顺序与颜色
# =========================================================
celltype_order <- c(
  "B", "Plasma",
  "pDC", "DC", "Mono",
  "T CD4 naive", "T CD4 non-naive",
  "T CD8 naive", "T CD8 non-naive",
  "UTC", "ILC"
)

celltype_cols <- c(
  "B" = "#7bc6d6",
  "Plasma" = "#025566",
  "pDC" = "#a7c957",
  "DC" = "#6a994e",
  "Mono" = "#386641",
  "T CD4 naive" = "#e9dd8f",
  "T CD4 non-naive" = "#d99b38",
  "T CD8 naive" = "#b63232",
  "T CD8 non-naive" = "#6e0000",
  "UTC" = "#8d6a8c",
  "ILC" = "#6f2141"
)

# =========================================================
# 5. 因子名规范化
# =========================================================
clean_factor_feature <- function(x) {
  x <- as.character(x)
  x <- gsub("^global_all_", "", x)
  x <- gsub("^B_", "", x)
  x <- gsub("^Plasma_", "", x)
  x <- gsub("^DC_", "", x)
  x <- gsub("^pDC_", "", x)
  x <- gsub("^Mono_", "", x)
  x <- gsub("^ILC_", "", x)
  x <- gsub("^T_CD4_Naive_", "", x)
  x <- gsub("^T_CD4_non_naive_", "", x)
  x <- gsub("^T_CD8_Naive_", "", x)
  x <- gsub("^T_CD8_non_naive_", "", x)
  x <- gsub("^UTC_", "", x)
  x <- gsub("cytokine_and receptors_ antiinflammatory", "cytokine_andreceptors_antiinflammatory", x, fixed = TRUE)
  x <- gsub("  +", " ", x)
  trimws(x)
}

short_factor_feature <- function(x) {
  x <- clean_factor_feature(x)

  x <- gsub("adhesion_molecules", "Adhesion\nmolecules", x, fixed = TRUE)
  x <- gsub("antigen_presentation_molecules", "Ag. presenting\nmolecules", x, fixed = TRUE)
  x <- gsub("cytokine_and_receptors_proinflammatory", "Pro-inflamm\ncytokines and R", x, fixed = TRUE)
  x <- gsub("cytokine_andreceptors_antiinflammatory", "Anti-inflamm\ncytokines and R", x, fixed = TRUE)
  x <- gsub("chemokines", "Chemokines", x, fixed = TRUE)
  x <- gsub("chemokine_receptors", "Chemokine\nreceptors", x, fixed = TRUE)
  x <- gsub("IFN_Type_1_2_Lambda", "IFN type\n1 and 2", x, fixed = TRUE)
  x <- gsub("IFN_response", "IFN-induced", x, fixed = TRUE)
  x <- gsub("TNF_receptors_ligands", "TNF\nreceptor/ligand", x, fixed = TRUE)
  x <- gsub("CD4T_TFH_UP", "Tfh", x, fixed = TRUE)
  x <- gsub("CD4T_TH1_UP", "Th1", x, fixed = TRUE)
  x <- gsub("CD4T_TH2_UP", "Th2", x, fixed = TRUE)
  x <- gsub("CD4T_TH17_UP", "Th17", x, fixed = TRUE)
  x <- gsub("Tregs_FoxP3_stabilization", "FOXP3", x, fixed = TRUE)
  x <- gsub("CD8T_exhaustion", "Exhaustion", x, fixed = TRUE)
  x <- gsub("CD8T_tcr_activation", "TCR activation", x, fixed = TRUE)
  x <- gsub("effector", "Effector", x, fixed = TRUE)
  x <- gsub("IFNG_response", "IFNG response", x, fixed = TRUE)
  x <- gsub("IL4-IL13_response", "IL4/13 response", x, fixed = TRUE)
  x <- gsub("antigen-crosspresentation", "Ag crosspresentation", x, fixed = TRUE)
  x <- gsub("TNF-via-NFkB_signaling", "TNF via NFκB sig.", x, fixed = TRUE)
  x <- gsub("IL6-JAK-STAT3_signaling", "IL6-JAK-STAT3", x, fixed = TRUE)
  x <- gsub("JAK-STAT_signaling", "JAK-STAT", x, fixed = TRUE)
  x <- gsub("MHC-II-I presentation", "MHC-I/II", x, fixed = TRUE)
  x <- gsub("complement_production", "Complement", x, fixed = TRUE)

  x <- gsub("_", " ", x, fixed = TRUE)
  x <- gsub("\\s+", " ", x)
  trimws(x)
}

infer_function_group <- function(x) {
  x <- clean_factor_feature(x)

  dplyr::case_when(
    grepl("adhesion_molecules$", x) ~ "Adhesion\nmolecules",
    grepl("antigen_presentation_molecules$", x) ~ "Ag. presenting\nmolecules",
    grepl("cytokine_and_receptors_proinflammatory$", x) ~ "Pro-inflamm\ncytokines and R",
    grepl("cytokine_andreceptors_antiinflammatory$", x) ~ "Anti-inflamm\ncytokines and R",
    grepl("chemokines$", x) ~ "Chemokines",
    grepl("chemokine_receptors$", x) ~ "Chemokine\nreceptors",
    grepl("IFN_Type_1_2_Lambda$", x) ~ "IFN type\n1 and 2",
    grepl("IFN_response$", x) ~ "IFN-induced",
    grepl("TNF_receptors_ligands$", x) ~ "TNF\nreceptor/ligand",
    grepl("CD4T_TFH_UP$", x) ~ "Tfh",
    grepl("CD4T_TH1_UP$", x) ~ "Th1",
    grepl("CD4T_TH2_UP$", x) ~ "Th2",
    grepl("CD4T_TH17_UP$", x) ~ "Th17",
    grepl("Tregs_FoxP3_stabilization$", x) ~ "FOXP3",
    grepl("CD8T_exhaustion$", x) ~ "Exhaustion",
    grepl("CD8T_tcr_activation$", x) ~ "TCR activation",
    grepl("effector$", x) ~ "Effector",
    grepl("IFNG_response$", x) ~ "IFNG response",
    grepl("IL4-IL13_response$", x) ~ "IL4/13 response",
    grepl("antigen-crosspresentation$", x) ~ "Ag crosspresentation",
    grepl("TNF-via-NFkB_signaling$", x) ~ "TNF via NFκB sig.",
    grepl("IL6-JAK-STAT3_signaling$", x) ~ "IFN-induced",
    grepl("JAK-STAT_signaling$", x) ~ "IFN-induced",
    TRUE ~ NA_character_
  )
}

is_bad_factor_label <- function(x) {
  x <- trimws(as.character(x))
  is.na(x) | x == "" | x == "NA" |
    grepl("^[0-9]+$", x) |
    grepl("^[0-9]+\\.[0-9]+$", x)
}

# =========================================================
# 6. 读取数据
# =========================================================
all_df <- fread(all_path, data.table = FALSE)
colnames(all_df) <- trimws(colnames(all_df))
all_df <- all_df %>% dplyr::filter(group0 == "Pediatric healthy")

# =========================================================
# 7. 解析 factor
# =========================================================
parse_factor_info <- function(df) {
  sp <- strsplit(df$factor, "-X-", fixed = TRUE)
  df$factor_id <- vapply(sp, function(x) ifelse(length(x) >= 1, x[1], NA_character_), character(1))
  df$factor_src <- vapply(sp, function(x) ifelse(length(x) >= 2, x[2], NA_character_), character(1))
  df$factor_feature <- vapply(sp, function(x) ifelse(length(x) >= 3, x[3], NA_character_), character(1))
  df
}
all_df <- parse_factor_info(all_df)

# =========================================================
# 8. 构造长表 + 清理无效因子
# =========================================================
all_df <- all_df %>%
  dplyr::mutate(
    disease = group1,
    disease_short = disease_short_map[disease],
    RowID = paste0(celltype, "-X-", factor),
    factor_feature_clean = clean_factor_feature(factor_feature),
    function_group = infer_function_group(factor_feature),
    factor_short = short_factor_feature(factor_feature)
  ) %>%
  dplyr::filter(!is.na(disease_short)) %>%
  dplyr::filter(!is_bad_factor_label(factor_feature_clean)) %>%
  dplyr::filter(!is_bad_factor_label(factor_short)) %>%
  dplyr::filter(!is.na(function_group))

label_df <- all_df %>%
  dplyr::select(RowID, celltype, factor_id, factor_short, function_group) %>%
  dplyr::distinct() %>%
  dplyr::arrange(function_group, factor_short, factor_id) %>%
  dplyr::group_by(celltype, factor_short) %>%
  dplyr::mutate(n_same = dplyr::n()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    factor_label = ifelse(n_same > 1, paste0(factor_short, " (F", factor_id, ")"), factor_short)
  )

all_df <- all_df %>%
  dplyr::left_join(label_df %>% dplyr::select(RowID, factor_label), by = "RowID") %>%
  dplyr::filter(!is_bad_factor_label(factor_label))

# 用原始 p_value 判显著，标记用 *
sig_key <- all_df %>%
  dplyr::select(RowID, disease_short, p_value) %>%
  dplyr::group_by(RowID, disease_short) %>%
  dplyr::summarise(p_value = min(p_value, na.rm = TRUE), .groups = "drop") %>%
  dplyr::mutate(sig_mark = ifelse(is.na(p_value), "", ifelse(p_value < p_cutoff, "*", "")))

plot_df <- all_df %>%
  dplyr::select(RowID, disease_short, beta) %>%
  dplyr::distinct() %>%
  dplyr::left_join(sig_key %>% dplyr::select(RowID, disease_short, sig_mark),
                   by = c("RowID", "disease_short"))

col_meta <- all_df %>%
  dplyr::select(RowID, celltype, function_group, factor_label) %>%
  dplyr::distinct()

# =========================================================
# 9. 只保留“至少一个疾病中显著”的因子
# =========================================================
if (keep_only_any_sig) {
  sig_rowid_df <- all_df %>%
    dplyr::select(RowID, disease_short, p_value) %>%
    dplyr::distinct() %>%
    dplyr::group_by(RowID) %>%
    dplyr::summarise(
      min_p_value = min(p_value, na.rm = TRUE),
      any_sig = any(p_value < p_cutoff, na.rm = TRUE),
      .groups = "drop"
    )

  keep_rowid <- sig_rowid_df %>%
    dplyr::filter(any_sig) %>%
    dplyr::pull(RowID)

  plot_df <- plot_df %>% dplyr::filter(RowID %in% keep_rowid)
  col_meta <- col_meta %>% dplyr::filter(RowID %in% keep_rowid)
}

# =========================================================
# 10. 按模块列数多少自动排序：多的放前面
# =========================================================
group_count_df <- col_meta %>%
  dplyr::count(function_group, name = "n_factor") %>%
  dplyr::arrange(dplyr::desc(n_factor), function_group)

function_order_dynamic <- group_count_df$function_group

# =========================================================
# 11. 矩阵
# =========================================================
mat_df <- plot_df %>%
  dplyr::select(RowID, disease_short, beta) %>%
  tidyr::pivot_wider(names_from = disease_short, values_from = beta)

sig_df2 <- plot_df %>%
  dplyr::select(RowID, disease_short, sig_mark) %>%
  tidyr::pivot_wider(names_from = disease_short, values_from = sig_mark)

wanted_diseases <- disease_order[disease_order %in% colnames(mat_df)]

for (d in wanted_diseases) {
  if (!(d %in% colnames(mat_df))) mat_df[[d]] <- NA_real_
  if (!(d %in% colnames(sig_df2))) sig_df2[[d]] <- ""
}

mat <- as.matrix(as.data.frame(mat_df[, wanted_diseases, drop = FALSE]))
rownames(mat) <- mat_df$RowID

sig_mat <- as.matrix(as.data.frame(sig_df2[, wanted_diseases, drop = FALSE]))
rownames(sig_mat) <- sig_df2$RowID

mat <- t(mat)
sig_mat <- t(sig_mat)

# =========================================================
# 12. 列顺序：模块按数量多少排序，模块内再按 celltype / label
# =========================================================
col_meta <- col_meta %>%
  dplyr::mutate(
    function_group = factor(function_group, levels = function_order_dynamic),
    celltype = factor(celltype, levels = celltype_order)
  ) %>%
  dplyr::arrange(function_group, celltype, factor_label)

keep_cols <- col_meta$RowID[col_meta$RowID %in% colnames(mat)]
col_meta <- col_meta %>% dplyr::filter(RowID %in% keep_cols)
mat <- mat[, col_meta$RowID, drop = FALSE]
sig_mat <- sig_mat[, col_meta$RowID, drop = FALSE]

colnames(mat) <- col_meta$factor_label
colnames(sig_mat) <- col_meta$factor_label

# =========================================================
# 13. 行顺序与分组
# =========================================================
row_meta <- data.frame(
  disease_short = rownames(mat),
  disease_group = disease_group_map[rownames(mat)],
  stringsAsFactors = FALSE
)

row_meta$disease_group <- factor(row_meta$disease_group, levels = row_group_order)
row_meta <- row_meta %>% dplyr::arrange(disease_group, disease_short)

mat <- mat[row_meta$disease_short, , drop = FALSE]
sig_mat <- sig_mat[row_meta$disease_short, , drop = FALSE]

# =========================================================
# 14. 颜色增强（不对称，贴近原图）
# =========================================================
abs_vals <- abs(mat[is.finite(mat)])
if (length(abs_vals) > 0) {
  max_abs <- as.numeric(stats::quantile(abs_vals, probs = clip_quantile, na.rm = TRUE))
} else {
  max_abs <- 1
}
if (!is.finite(max_abs) || max_abs <= 0) max_abs <- 1

mat[mat >  max_abs] <-  max_abs
mat[mat < -max_abs] <- -max_abs

neg_lim <- -0.55 * max_abs
neg_mid <- -0.16 * max_abs
pos_mid1 <- 0.18 * max_abs
pos_mid2 <- 0.42 * max_abs
pos_hi1  <- 0.72 * max_abs
pos_hi2  <- max_abs

col_fun <- circlize::colorRamp2(
  c(neg_lim, neg_mid, 0, pos_mid1, pos_mid2, pos_hi1, pos_hi2),
  c("#0B4FA2", "#9CC3E5", "#F3F0DE", "#F7D97A", "#F4A259", "#E85D2A", "#B3001B")
)





col_fun <- circlize::colorRamp2(
  c(neg_lim, neg_mid, 0, pos_mid, pos_hi1, pos_hi2),
  c("#0B4FA2", "#8CB9DF", "#F3F0DE", "#F3C55A", "#F07C3E", "#B3001B")
)

# =========================================================
# 15. 注释
# =========================================================
bottom_ha <- HeatmapAnnotation(
  CellType = col_meta$celltype,
  col = list(CellType = celltype_cols),
  annotation_name_side = "left",
  annotation_name_gp = gpar(fontsize = 8.8, fontface = "bold"),
  simple_anno_size = unit(3.2, "mm"),
  show_annotation_name = TRUE,
  show_legend = TRUE,
  border = TRUE,
  gp = gpar(col = "black", lwd = 0.9)
)

left_ha <- rowAnnotation(
  DiseaseGroup = row_meta$disease_group,
  col = list(DiseaseGroup = row_group_cols),
  show_annotation_name = FALSE,
  width = unit(2.8, "mm"),
  border = TRUE,
  gp = gpar(col = "black", lwd = 0.9)
)

col_split <- factor(as.character(col_meta$function_group), levels = function_order_dynamic)
col_split <- droplevels(col_split)

row_split <- factor(as.character(row_meta$disease_group), levels = row_group_order)
row_split <- droplevels(row_split)

# =========================================================
# 16. 热图
# =========================================================
ht <- left_ha + Heatmap(
  mat,
  name = "Corrected activity score\n(ULM coefficient)",
  col = col_fun,
  na_col = "#ECECEC",

  cluster_rows = FALSE,
  cluster_columns = FALSE,

  row_split = row_split,
  column_split = col_split,


  row_gap = unit(1.6, "mm"),
  column_gap = unit(0.7, "mm"),

  show_row_dend = FALSE,
  show_column_dend = FALSE,

  row_names_side = "right",
  row_names_gp = gpar(fontsize = row_name_size, col = "#2B2B2B"),
  column_names_side = "top",
  column_names_rot = 90,
  column_names_gp = gpar(fontsize = col_name_size, col = "#2B2B2B"),

  bottom_annotation = bottom_ha,

  rect_gp = gpar(col = "#6F6F6F", lwd = 0.4),
  border = FALSE,

  row_title = NULL,
  column_title = NULL,

  heatmap_legend_param = list(
    at = c(neg_lim, 0, pos_hi2),
    labels = c(paste0("<", round(neg_lim, 1)), "0", paste0(">", round(pos_hi2, 1))),
    legend_direction = "horizontal",
    title_gp = gpar(fontsize = 8.8, fontface = "bold"),
    labels_gp = gpar(fontsize = 8),
    border = "#666666",
    grid_height = unit(3.2, "mm"),
    grid_width = unit(14, "mm")
  ),

  cell_fun = function(j, i, x, y, w, h, fill) {
    if (!is.na(sig_mat[i, j]) && sig_mat[i, j] != "") {
      grid.text(
        sig_mat[i, j],
        x = x, y = y,
        gp = gpar(fontsize = sig_point_fontsize, col = "#111111", fontface = "bold")
      )
    }
  }
)

# =========================================================
# 17. 输出
# =========================================================
pdf(file.path(out_dir, "LMEM_L1_reference_style_heatmap_narrow_asym_final.pdf"), width = pdf_width, height = pdf_height)
draw(
  ht,
  heatmap_legend_side = "bottom",
  annotation_legend_side = "bottom",
  merge_legends = FALSE,
  padding = unit(c(0.8, 0.6, 0.8, 0.6), "mm")
)
dev.off()

png(file.path(out_dir, "LMEM_L1_reference_style_heatmap_narrow_asym_final.png"),
    width = pdf_width, height = pdf_height, units = "in", res = 500)
draw(
  ht,
  heatmap_legend_side = "bottom",
  annotation_legend_side = "bottom",
  merge_legends = FALSE,
  padding = unit(c(0.8, 0.6, 0.8, 0.6), "mm")
)
dev.off()

write.table(
  mat,
  file = file.path(out_dir, "LMEM_L1_reference_style_matrix_narrow_asym_final.tsv"),
  sep = "\t", quote = FALSE, col.names = NA
)

write.table(
  sig_mat,
  file = file.path(out_dir, "LMEM_L1_reference_style_sig_narrow_asym_final.tsv"),
  sep = "\t", quote = FALSE, col.names = NA
)

write.table(
  col_meta,
  file = file.path(out_dir, "LMEM_L1_reference_style_column_metadata_narrow_asym_final.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(
  row_meta,
  file = file.path(out_dir, "LMEM_L1_reference_style_row_metadata_narrow_asym_final.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

cat("\n[DONE]\n")
cat("Output dir: ", out_dir, "\n")
cat("PDF: ", file.path(out_dir, "LMEM_L1_reference_style_heatmap_narrow_asym_final.pdf"), "\n")
cat("PNG: ", file.path(out_dir, "LMEM_L1_reference_style_heatmap_narrow_asym_final.png"), "\n")
cat("Kept columns: ", ncol(mat), "\n")
cat("Kept rows: ", nrow(mat), "\n")
