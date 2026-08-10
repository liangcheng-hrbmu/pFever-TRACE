suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
})

# =========================================================
# 1. 输入输出
# =========================================================
infile <- "data/analysis/inflammation_signatures/14_LMEM_Corr_L1/LMEM_L1_all_results.tsv"
outdir <- "data/analysis/inflammation_signatures/14_LMEM_Corr_L1/L1_disease_clustering_onlytree"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# =========================================================
# 2. 参数
# =========================================================
use_sig_only <- TRUE      # 只保留至少一个疾病中 p<0.05 的特征
run_zscore   <- TRUE      # 是否同时跑 z-score 版本

# 想看的聚类组合
cluster_configs <- list(
  list(name = "euclidean_complete",   dist = "euclidean",  hclust = "complete"),
  list(name = "euclidean_average",    dist = "euclidean",  hclust = "average"),
  list(name = "euclidean_wardD2",     dist = "euclidean",  hclust = "ward.D2"),
  list(name = "manhattan_complete",   dist = "manhattan",  hclust = "complete"),
  list(name = "manhattan_average",    dist = "manhattan",  hclust = "average"),
  list(name = "correlation_complete", dist = "correlation", hclust = "complete"),
  list(name = "correlation_average",  dist = "correlation", hclust = "average")
)

# =========================================================
# 3. 读取数据
# =========================================================
df <- fread(infile, data.table = FALSE)
colnames(df) <- trimws(colnames(df))

# 只保留 disease vs healthy
df <- df %>%
  dplyr::filter(group0 == "Pediatric healthy")

# =========================================================
# 4. 疾病短名
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

# =========================================================
# 5. 构造 feature ID
# =========================================================
df <- df %>%
  dplyr::mutate(
    disease_short = disease_short_map[group1],
    FeatureID = paste(celltype, factor, sep = " | ")
  ) %>%
  dplyr::filter(!is.na(disease_short))

# =========================================================
# 6. 是否只用显著特征
# =========================================================
if (use_sig_only) {
  keep_features <- df %>%
    dplyr::group_by(FeatureID) %>%
    dplyr::summarise(any_sig = any(p_value < 0.05, na.rm = TRUE), .groups = "drop") %>%
    dplyr::filter(any_sig) %>%
    dplyr::pull(FeatureID)

  df <- df %>% dplyr::filter(FeatureID %in% keep_features)
}

# =========================================================
# 7. 宽表：行=疾病，列=feature，值=beta
# =========================================================
mat_df <- df %>%
  dplyr::select(disease_short, FeatureID, beta) %>%
  dplyr::distinct() %>%
  tidyr::pivot_wider(names_from = FeatureID, values_from = beta)

mat <- as.data.frame(mat_df)
rownames(mat) <- mat$disease_short
mat$disease_short <- NULL
mat <- as.matrix(mat)

# 缺失填 0
mat[is.na(mat)] <- 0

# 去掉全0列
keep_cols <- colSums(abs(mat)) > 0
mat <- mat[, keep_cols, drop = FALSE]

# =========================================================
# 8. 构造 raw / zscore 两套矩阵
# =========================================================
mat_list <- list(raw = mat)

if (run_zscore) {
  mat_z <- scale(mat, center = TRUE, scale = TRUE)
  keep_cols_z <- colSums(is.na(mat_z)) == 0
  mat_z <- mat_z[, keep_cols_z, drop = FALSE]
  mat_list$zscore <- mat_z
}

# =========================================================
# 9. 距离函数
# =========================================================
get_dist <- function(x, method) {
  if (method %in% c("euclidean", "manhattan", "maximum", "canberra", "binary", "minkowski")) {
    return(dist(x, method = method))
  }

  if (method == "correlation") {
    cmat <- suppressWarnings(cor(t(x), method = "pearson", use = "pairwise.complete.obs"))
    cmat[is.na(cmat)] <- 0
    return(as.dist(1 - cmat))
  }

  stop("Unsupported distance method: ", method)
}

# =========================================================
# 10. 画树函数（只画聚类树）
# =========================================================
plot_dendrogram_only <- function(hc_obj,
                                 main_letter = "b",
                                 cex_lab = 1.1,
                                 line_width = 1.0) {
  par(mar = c(4.5, 2.8, 1.2, 0.5), xpd = NA, lwd = line_width)

  plot(
    hc_obj,
    hang = -1,
    axes = FALSE,
    ann = FALSE,
    leaflab = "none"
  )

  usr <- par("usr")
  x_pos <- seq_along(hc_obj$order)
  lab_order <- hc_obj$labels[hc_obj$order]

  text(
    x = x_pos,
    y = usr[3] - 0.03 * (usr[4] - usr[3]),
    labels = lab_order,
    srt = 90,
    adj = 1,
    cex = cex_lab
  )

  mtext(main_letter, side = 3, adj = -0.05, line = -1.0, font = 2, cex = 1.8)
}

# =========================================================
# 11. 循环输出多种方法
# =========================================================
all_orders <- list()

for (mat_name in names(mat_list)) {

  mat_use <- mat_list[[mat_name]]

  for (cfg in cluster_configs) {

    cfg_name  <- cfg$name
    dist_name <- cfg$dist
    hc_name   <- cfg$hclust

    message("[RUN] matrix=", mat_name, " | dist=", dist_name, " | hclust=", hc_name)

    dist_mat <- get_dist(mat_use, dist_name)
    hc <- hclust(dist_mat, method = hc_name)

    ordered_diseases <- hc$labels[hc$order]

    # 保存顺序
    all_orders[[paste(mat_name, cfg_name, sep = "__")]] <- data.frame(
      matrix_type = mat_name,
      config = cfg_name,
      order = seq_along(ordered_diseases),
      disease = ordered_diseases,
      stringsAsFactors = FALSE
    )

    # 文件名
    base_name <- paste0("L1_disease_hclust_", mat_name, "_", cfg_name)

    # PDF
    pdf(file.path(outdir, paste0(base_name, ".pdf")),
        width = 6.6, height = 4.2)
    plot_dendrogram_only(
      hc_obj = hc,
      main_letter = "b",
      cex_lab = 1.05,
      line_width = 1.0
    )
    dev.off()

    # PNG
    png(file.path(outdir, paste0(base_name, ".png")),
        width = 6.6, height = 4.2, units = "in", res = 450)
    plot_dendrogram_only(
      hc_obj = hc,
      main_letter = "b",
      cex_lab = 1.05,
      line_width = 1.0
    )
    dev.off()

    # 距离矩阵
    write.table(
      as.matrix(dist_mat),
      file = file.path(outdir, paste0(base_name, "_distance_matrix.tsv")),
      sep = "\t", quote = FALSE, col.names = NA
    )
  }
}

# =========================================================
# 12. 保存全部聚类顺序
# =========================================================
order_df <- dplyr::bind_rows(all_orders)

write.table(
  order_df,
  file = file.path(outdir, "L1_disease_hclust_all_orders.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# 也保存原始矩阵
write.table(
  mat,
  file = file.path(outdir, "L1_disease_feature_matrix_raw.tsv"),
  sep = "\t", quote = FALSE, col.names = NA
)

if (run_zscore) {
  write.table(
    mat_list$zscore,
    file = file.path(outdir, "L1_disease_feature_matrix_zscore.tsv"),
    sep = "\t", quote = FALSE, col.names = NA
  )
}

cat("\n[DONE]\n")
cat("Output dir: ", outdir, "\n")
cat("Generated configs:\n")
for (cfg in cluster_configs) {
  cat(" - ", cfg$name, "\n", sep = "")
}
cat("Matrix types: ", paste(names(mat_list), collapse = ", "), "\n")