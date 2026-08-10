library(Seurat)
library(dplyr)
library(dendextend)
library(circlize)

## 0) 设定你要画 dendrogram 的“单位”
## group_col = 你要聚类的标签（例如 subCluster / majorCluster）
group_col <- "subCluster"     # 你也可以改成 "majorCluster"
color_col <- "majorCluster"   # 用哪个列来给标签上色（通常是 majorCluster）

## 1) 计算平均表达：group_col × gene
## 使用 SCT / RNA 都行：你用哪个 assay 就改 DefaultAssay
DefaultAssay(seurat_obj) <- "RNA"

# 返回 gene × group 的矩阵，转置成 group × gene
ave_exp <- AverageExpression(
  seurat_obj,
  group.by = group_col,
  assays = "RNA",
  slot = "data",
  verbose = FALSE
)$RNA
ave_exp <- t(ave_exp)  # group × gene

# 可选：只取变异基因让速度更快
# hvgs <- VariableFeatures(seurat_obj)
# ave_exp <- ave_exp[, hvgs, drop = FALSE]

## 2) log1p（模仿原文）
ave_exp <- log1p(ave_exp)

## 3) （可选）去掉某类，比如 Epi（按你原图逻辑）
# ave_exp <- ave_exp[substr(rownames(ave_exp), 1, 3) != "Epi", , drop = FALSE]

## 4) 准备 meta_sc：每个 group 的 majorCluster，用于上色
meta_sc <- seurat_obj@meta.data %>%
  dplyr::select(all_of(c(group_col, color_col))) %>%
  distinct() %>%
  arrange(.data[[group_col]])

# ave_exp 行顺序与 meta_sc 对齐
meta_sc <- meta_sc[match(rownames(ave_exp), meta_sc[[group_col]]), , drop = FALSE]


# ave_exp 必须是：行=cluster，列=gene，且全部是 numeric
ave_exp <- as.matrix(ave_exp)
storage.mode(ave_exp) <- "numeric"


## 5) 相关距离 + ward.D 聚类（完全模仿）
hc <- as.dist((1 - cor(t(ave_exp))) / 2) %>%
  hclust(method = "ward.D") %>%
  as.dendrogram() %>%
  dendextend::set("labels_cex", 0.5)

## 6) 给标签上色（用 majorCluster 颜色）
# 你可以用你现在那套：celltype_colors

celltype_colors <- c(
  "B"            = "#1F77B4",  # 蓝
  "Plasma"                = "#AEC7E8",  # 浅蓝
  "pDC"       = "#FF7F0E",  # 橙
  "DC"               = "#FFBB78",  # 浅橙
  "Mono"          = "#2CA02C",  # 绿
  "T CD4 naive"           = "#98DF8A",  # 浅绿
  "T CD4 non-naive"          = "#D62728",  # 红
  "T CD8 naive"   = "#FF9896",  # 粉红
  "T CD8 non-naive"            = "#9467BD" ,  # 紫

  "UTC"            = "#C5B0D5" ,  # 浅紫

  "ILC"            = "#8C564B"   # 棕


)


colors_to_use <- base::unname(celltype_colors[as.character(meta_sc[[color_col]])])

colors_to_use <- colors_to_use[order.dendrogram(hc)]

labels_colors(hc) <- colors_to_use




## 7) 画 circular dendrogram（模仿）
pdf("Figure1d_subCluster_circular.pdf", width = 6, height = 6)
circlize_dendrogram(
  hc,
  labels_track_height = 0.6,
  dend_track_height = 0.2
)
dev.off()
