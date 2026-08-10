setwd("data/analysis/atlas_objects")

library(dplyr)
library(ggplot2)
GSE120575<- readRDS("GSE120575/GSE120575.rds")
## 1. 读入 metadata
meta <- seurat_obj@meta.data

## 4. 构建 疾病 × 细胞类型 计数矩阵
celltype_mat <- table(meta$disease_short, meta$majorCluster)

## 5. 转为比例（每个疾病内部 = 1）
celltype_prop <- celltype_mat / rowSums(celltype_mat)





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


# clustering
hc <- hclust(dist(celltype_prop), method = "ward.D")
celltype_prop <- celltype_prop[hc$order, ]

celltype_df <- as.data.frame(celltype_prop)
colnames(celltype_df) <- c("Disease", "CellType", "Proportion")

# plot
ggplot(celltype_df,
       aes(x = Disease, y = Proportion, fill = CellType)) +
  geom_col(position = position_fill(reverse = TRUE)) +
  scale_fill_manual(values = celltype_colors) +
  coord_flip() +
  theme_classic()
ggsave(
  filename = "CellType_Composition_per_Disease.pdf",
  width = 8,
  height = 4
)




##########################################
### Fig. 1c (Cell type composition per disease)
##########################################

library(dplyr)
library(ggplot2)
library(dendextend)

# 1) load and filter
meta <- seurat_obj@meta.data



# 2) create mat
major_mat <- table(meta$disease_Level_3, meta$majorCluster)
tmp <- major_mat / rowSums(major_mat)
range(rowSums(tmp))  # 应该全是 1

# hclust（保持与你原代码一致：dist + ward.D，并转 dendrogram）
hc <- hclust(dist(tmp), method = "ward.D") %>%
  as.dendrogram() %>%
  dendextend::set("labels_cex", 0.5)

# dendrogram（单独输出，模仿原来的 majorCluster_dendrogram.pdf）
pdf("majorCluster_dendrogram_disease_Level_3.pdf", width = 2, height = 4.5)
plot_horiz.dendrogram(hc)
dev.off()

# 按聚类顺序重排（模仿原代码）
tmp <- tmp[order.dendrogram(hc), , drop = FALSE]

major_df <- tmp %>% as.data.frame()  # freq
colnames(major_df) <- c("Disease", "majorCluster", "count")


# 3) plot(horiz) —— 主题/字号/布局模仿你原图
ggplot(major_df, aes(x = Disease, y = count, fill = majorCluster)) +
  geom_col(position = position_fill(reverse = TRUE)) +
  # 颜色：优先用你文章统一的 color$majorCluster；如果你想用 celltype_colors 就替换这里
  scale_fill_manual(values = celltype_colors)+
  coord_flip() +
  labs(
    x = "", y = "Proportion of cells",
    title = "Cell type composition per disease (Level 3)",
    fill = "Cell type"
  ) +
  scale_x_discrete(limits = levels(major_df$Disease)) +
  theme_classic() +
  theme(
    line = element_line(linewidth = 0.3),
    text = element_text(size = 6,  color = "black"),
    axis.text = element_text(size = 8,  color = "black"),
    legend.key.size = unit(3, "mm"),
    plot.title = element_text(size = 8, hjust = 0.5),
    aspect.ratio = 2 / 1.5
  ) -> p

ggsave(
  plot = p, filename = "Figure1c_disease_Level_3.pdf",
  width = 6, height = 4
)
