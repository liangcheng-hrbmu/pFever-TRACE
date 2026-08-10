# ============================================================
# Final Fig3 b-e
#
# b Disease similarity based on inflammatory program profiles
# c Recurrent inflammatory axes
# d Immune compartment contribution to inflammatory programs
# e Inflammatory burden in Level 2 cell states
#
# Fig3a = existing global inflammatory heatmap
#
# Key:
#   IL6-JAK-STAT3 + JAK-STAT -> JAK STAT signaling
#   disease group colors and cell type colors match Fig3a
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(ggplot2)
  library(Seurat)
  library(scales)
  library(grid)
})

# ============================================================
# 0. Paths
# ============================================================

l1_path <- "data/analysis/inflammation_signatures/14_LMEM_Corr_L1_v2/LMEM_L1_all_results.tsv"
l2_path <- "data/analysis/inflammation_signatures/16_LMEM_Corr_L2_all/LMEM_L2_all_results.tsv"
seurat_rds <- "data/analysis/atlas_objects/integrated_pediatric.rds"

out_dir <- "data/analysis/inflammation_signatures/14_LMEM_Corr_L1_v2/Fig3_final_bcde"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 1. Settings
# ============================================================

p_cutoff <- 0.05

# FALSE = use p_value first, consistent with Fig3a
# TRUE  = use adj_p first
prefer_adj_p <- FALSE

# Fig3c
top_n_programs_c <- 10

# Fig3d: core programs only
core_programs_d <- c(
  "JAK STAT signaling",
  "Chemokines",
  "IFN type 1/2",
  "TNF via NFκB",
  "Antigen presentation",
  "MHC I II presentation",
  "Complement",
  "CD8 exhaustion",
  "TCR activation",
  "Th1"
)

# Fig3e
top_l2_n <- 12

# burden mode:
#   sig_up_beta_sum = sum(beta) for p < cutoff and beta > 0
#   sig_any_abs_sum = sum(abs(beta)) for p < cutoff
BURDEN_MODE_L2 <- "sig_up_beta_sum"

# ============================================================
# 2. Disease and color settings
# ============================================================

disease_short_map <- c(
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

disease_group_order <- c(
  "virus",
  "bacteria",
  "Sepsis",
  "Hyperinflam",
  "Autoimmune"
)

row_group_cols <- c(
  "virus"       = "#6699CC",
  "bacteria"    = "#C48EAD",
  "Sepsis"      = "#FF2A2A",
  "Hyperinflam" = "#E3A722",
  "Autoimmune"  = "#36A596"
)

group_cols <- row_group_cols

disease_group_labels <- c(
  "virus"       = "Virus",
  "bacteria"    = "Bacteria",
  "Sepsis"      = "Sepsis",
  "Hyperinflam" = "Hyperinflam.",
  "Autoimmune"  = "Autoimmune"
)

celltype_order <- c(
  "B", "Plasma",
  "pDC", "DC", "Mono",
  "T CD4 naive", "T CD4 non-naive",
  "T CD8 naive", "T CD8 non-naive",
  "UTC", "ILC",
  "Other"
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
  "ILC" = "#6f2141",
  "Other" = "#BDBDBD"
)

# ============================================================
# 3. Helper functions
# ============================================================

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

  x <- gsub(
    "cytokine_and receptors_ antiinflammatory",
    "cytokine_andreceptors_antiinflammatory",
    x,
    fixed = TRUE
  )

  x <- gsub("  +", " ", x)
  trimws(x)
}

parse_factor_info <- function(df) {
  sp <- strsplit(df$factor, "-X-", fixed = TRUE)

  df$factor_id <- vapply(
    sp,
    function(x) ifelse(length(x) >= 1, x[1], NA_character_),
    character(1)
  )

  df$factor_src <- vapply(
    sp,
    function(x) ifelse(length(x) >= 2, x[2], NA_character_),
    character(1)
  )

  df$factor_feature <- vapply(
    sp,
    function(x) ifelse(length(x) >= 3, x[3], NA_character_),
    character(1)
  )

  df
}

infer_program <- function(x) {
  x <- clean_factor_feature(x)

  dplyr::case_when(
    grepl("IL6-JAK-STAT3_signaling$", x) ~ "JAK STAT signaling",
    grepl("JAK-STAT_signaling$", x) ~ "JAK STAT signaling",

    grepl("IFN_Type_1_2_Lambda$", x) ~ "IFN type 1/2",
    grepl("IFN_response$", x) ~ "IFN response",
    grepl("IFNG_response$", x) ~ "IFNG response",
    grepl("IL4-IL13_response$", x) ~ "IL4 13 response",

    grepl("TNF-via-NFkB_signaling$", x) ~ "TNF via NFκB",
    grepl("TNF_receptors_ligands$", x) ~ "TNF receptor ligand",

    grepl("chemokine_receptors$", x) ~ "Chemokine receptors",
    grepl("chemokines$", x) ~ "Chemokines",

    grepl("cytokine_and_receptors_proinflammatory$", x) ~ "Pro inflammatory cytokines",
    grepl("cytokine_andreceptors_antiinflammatory$", x) ~ "Anti inflammatory cytokines",

    grepl("antigen_presentation_molecules$", x) ~ "Antigen presentation",
    grepl("MHC-II-I presentation$", x) ~ "MHC I II presentation",
    grepl("antigen-crosspresentation$", x) ~ "Antigen crosspresentation",

    grepl("adhesion_molecules$", x) ~ "Adhesion molecules",
    grepl("complement_production$", x) ~ "Complement",

    grepl("CD8T_tcr_activation$", x) ~ "TCR activation",
    grepl("CD8T_exhaustion$", x) ~ "CD8 exhaustion",
    grepl("effector$", x) ~ "Effector",

    grepl("CD4T_TH1_UP$", x) ~ "Th1",
    grepl("CD4T_TH2_UP$", x) ~ "Th2",
    grepl("CD4T_TH17_UP$", x) ~ "Th17",
    grepl("CD4T_TFH_UP$", x) ~ "Tfh",
    grepl("Tregs_FoxP3_stabilization$", x) ~ "FOXP3 stabilization",

    TRUE ~ NA_character_
  )
}

short_program <- function(x) {
  x <- as.character(x)

  dplyr::case_when(
    x == "JAK STAT signaling" ~ "JAK STAT",
    x == "IFN type 1/2" ~ "IFN type 1/2",
    x == "TNF via NFκB" ~ "TNF via NFκB",
    x == "TNF receptor ligand" ~ "TNF receptor/ligand",
    x == "Pro inflammatory cytokines" ~ "Pro inflam cytokines",
    x == "Anti inflammatory cytokines" ~ "Anti inflam cytokines",
    x == "Antigen presentation" ~ "Ag presentation",
    x == "MHC I II presentation" ~ "MHC I/II",
    x == "Antigen crosspresentation" ~ "Ag crosspresentation",
    x == "CD8 exhaustion" ~ "CD8 T exhaustion",
    x == "FOXP3 stabilization" ~ "FOXP3",
    TRUE ~ x
  )
}

get_p_col <- function(df, prefer_adj = FALSE) {
  if (prefer_adj) {
    candidates <- c(
      "adj_p", "FDR", "fdr", "padj", "p_adj", "p.adjust",
      "p_value_adj", "p_val_adj",
      "p_value", "pval", "p_val"
    )
  } else {
    candidates <- c(
      "p_value", "pval", "p_val",
      "adj_p", "FDR", "fdr", "padj", "p_adj", "p.adjust",
      "p_value_adj", "p_val_adj"
    )
  }

  hit <- candidates[candidates %in% colnames(df)]

  if (length(hit) == 0) {
    stop(
      "No p value column found. Available columns:\n",
      paste(colnames(df), collapse = ", ")
    )
  }

  hit[1]
}

pick_first_col <- function(df, candidates) {
  hit <- candidates[candidates %in% colnames(df)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

clean_parent_celltype <- function(x) {
  x <- as.character(x)
  x <- trimws(x)

  dplyr::case_when(
    x %in% c("B", "B cells", "B cell") ~ "B",
    x %in% c("Plasma", "Plasma cells", "Plasma cell") ~ "Plasma",
    x %in% c("pDC", "pDCs") ~ "pDC",
    x %in% c("DC", "Dendritic cells", "Dendritic cell") ~ "DC",
    x %in% c("Mono", "Monocytes", "Monocyte") ~ "Mono",

    x %in% c("T CD4 naive", "T_CD4_Naive", "CD4 naive", "T CD4 Naive") ~ "T CD4 naive",
    x %in% c("T CD4 non-naive", "T_CD4_non_naive", "CD4 non-naive", "T CD4 non naive") ~ "T CD4 non-naive",

    x %in% c("T CD8 naive", "T_CD8_Naive", "CD8 naive", "T CD8 Naive") ~ "T CD8 naive",
    x %in% c("T CD8 non-naive", "T_CD8_non_naive", "CD8 non-naive", "T CD8 non naive") ~ "T CD8 non-naive",

    x %in% c("UTC") ~ "UTC",
    x %in% c("ILC") ~ "ILC",

    TRUE ~ x
  )
}

meta_to_chr <- function(x) {
  # Robust conversion for Seurat metadata columns.
  # Handles factor / character / data.frame / list-column.

  if (is.data.frame(x)) {
    x <- x[[1]]
  }

  if (is.list(x) && !is.factor(x)) {
    x <- vapply(
      x,
      function(z) {
        if (length(z) == 0 || all(is.na(z))) {
          NA_character_
        } else {
          paste(as.character(z), collapse = ";")
        }
      },
      character(1)
    )
  }

  as.character(x)
}

calc_hc_from_matrix <- function(mat) {
  mat <- as.matrix(mat)

  mat_z <- scale(mat)
  mat_z[is.na(mat_z)] <- 0
  mat_z[!is.finite(mat_z)] <- 0

  cm <- suppressWarnings(
    cor(t(mat_z), method = "pearson", use = "pairwise.complete.obs")
  )
  cm[is.na(cm)] <- 0
  diag(cm) <- 1

  dist_mat <- as.dist(1 - cm)
  hclust(dist_mat, method = "average")
}

get_cluster_info <- function(hc) {
  n <- length(hc$labels)

  leaf_y <- numeric(n)
  leaf_y[hc$order] <- seq_len(n)

  node_pos <- vector("list", n - 1)
  segs <- list()

  get_child <- function(v) {
    if (v < 0) {
      idx <- -v
      list(
        x = 0,
        y = leaf_y[idx],
        leaves = hc$labels[idx]
      )
    } else {
      node_pos[[v]]
    }
  }

  for (i in seq_len(n - 1)) {
    left <- get_child(hc$merge[i, 1])
    right <- get_child(hc$merge[i, 2])

    x_parent <- hc$height[i]
    y_parent <- mean(c(left$y, right$y))
    leaves_i <- sort(c(left$leaves, right$leaves))

    segs[[length(segs) + 1]] <- data.frame(
      x = left$x,
      xend = x_parent,
      y = left$y,
      yend = left$y
    )

    segs[[length(segs) + 1]] <- data.frame(
      x = right$x,
      xend = x_parent,
      y = right$y,
      yend = right$y
    )

    segs[[length(segs) + 1]] <- data.frame(
      x = x_parent,
      xend = x_parent,
      y = left$y,
      yend = right$y
    )

    node_pos[[i]] <- list(
      x = x_parent,
      y = y_parent,
      leaves = leaves_i
    )
  }

  list(
    segments = dplyr::bind_rows(segs),
    labels = data.frame(
      disease = hc$labels,
      x = 0,
      y = leaf_y,
      stringsAsFactors = FALSE
    )
  )
}

# ============================================================
# 4. Load L1 LMEM result
# ============================================================

l1_raw <- fread(l1_path, data.table = FALSE)
colnames(l1_raw) <- trimws(colnames(l1_raw))

cat("\n[READ L1]\n")
print(l1_path)

cat("\n[Columns]\n")
print(colnames(l1_raw))

if (!"factor_feature" %in% colnames(l1_raw)) {
  l1_raw <- parse_factor_info(l1_raw)
}

p_col_l1 <- get_p_col(l1_raw, prefer_adj = prefer_adj_p)
message("[L1 using p column] ", p_col_l1)

required_cols_l1 <- c(
  "celltype", "group1", "group0",
  "beta", "factor", "factor_feature", p_col_l1
)

missing_cols_l1 <- setdiff(required_cols_l1, colnames(l1_raw))

if (length(missing_cols_l1) > 0) {
  stop(
    "Missing L1 required columns: ",
    paste(missing_cols_l1, collapse = ", ")
  )
}

l1 <- l1_raw %>%
  dplyr::filter(group0 == "Pediatric healthy") %>%
  dplyr::mutate(
    celltype = clean_parent_celltype(trimws(as.character(celltype))),
    disease = disease_short_map[group1],
    disease_group = disease_group_map[disease],
    program = infer_program(factor_feature),
    program_short = short_program(program),
    p_use = .data[[p_col_l1]],
    sig_up = !is.na(p_use) & p_use < p_cutoff & beta > 0,
    sig_any = !is.na(p_use) & p_use < p_cutoff
  ) %>%
  dplyr::filter(
    !is.na(disease),
    !is.na(disease_group),
    !is.na(program)
  ) %>%
  dplyr::filter(celltype %in% celltype_order)

l1$disease <- factor(l1$disease, levels = disease_order)
l1$disease_group <- factor(l1$disease_group, levels = disease_group_order)
l1$celltype <- factor(l1$celltype, levels = celltype_order)

cat("\n[L1 disease count]\n")
print(table(l1$disease))

cat("\n[L1 disease group count]\n")
print(table(l1$disease_group))

cat("\n[L1 program count after merging JAK terms]\n")
print(table(l1$program))

# ============================================================
# Fig3b. Simplified disease similarity tree
# ============================================================

fig3b_dir <- file.path(out_dir, "Fig3b")
dir.create(fig3b_dir, recursive = TRUE, showWarnings = FALSE)

tree_df <- l1 %>%
  dplyr::mutate(
    feature = paste(celltype, program, sep = " | ")
  ) %>%
  dplyr::group_by(disease, feature) %>%
  dplyr::summarise(
    beta = mean(beta, na.rm = TRUE),
    .groups = "drop"
  )

mat_df <- tree_df %>%
  tidyr::pivot_wider(
    names_from = feature,
    values_from = beta,
    values_fill = 0
  ) %>%
  dplyr::arrange(factor(disease, levels = disease_order))

mat <- mat_df %>%
  tibble::column_to_rownames("disease") %>%
  as.matrix()

feature_sd <- apply(mat, 2, sd, na.rm = TRUE)
mat <- mat[, feature_sd > 0, drop = FALSE]

if (nrow(mat) < 3) {
  stop("Too few diseases for clustering.")
}

if (ncol(mat) < 3) {
  stop("Too few informative inflammatory features for clustering.")
}

hc <- calc_hc_from_matrix(mat)
tree_info <- get_cluster_info(hc)

label_b <- tree_info$labels %>%
  dplyr::mutate(
    disease_group = disease_group_map[disease],
    disease_group = factor(disease_group, levels = disease_group_order)
  )

max_h <- max(hc$height, na.rm = TRUE)
if (!is.finite(max_h) || max_h <= 0) max_h <- 1

x_dot <- -0.035 * max_h
x_label <- -0.080 * max_h

p3b <- ggplot() +
  geom_segment(
    data = tree_info$segments,
    aes(x = x, y = y, xend = xend, yend = yend),
    linewidth = 0.55,
    color = "#2B2B2B"
  ) +
  geom_point(
    data = label_b,
    aes(x = x_dot, y = y, color = disease_group),
    size = 3.0
  ) +
  geom_text(
    data = label_b,
    aes(x = x_label, y = y, label = disease),
    hjust = 0,
    size = 3.2,
    color = "#111111"
  ) +
  scale_color_manual(
    values = group_cols,
    breaks = disease_group_order,
    labels = disease_group_labels,
    drop = FALSE,
    name = "Disease group"
  ) +
  scale_x_reverse(
    limits = c(max_h * 1.05, -max_h * 0.30),
    breaks = pretty(c(0, max_h), n = 4),
    name = "1 - Pearson correlation"
  ) +
  scale_y_continuous(
    limits = c(0.4, max(label_b$y) + 0.4),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  labs(
    title = "Disease similarity based on inflammatory program profiles",
    y = NULL
  ) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 9) +
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    axis.text.x = element_text(size = 8, color = "black"),
    axis.title.x = element_text(size = 8.5),
    legend.position = "right",
    legend.title = element_text(size = 8, face = "bold"),
    legend.text = element_text(size = 7.5),
    legend.key.size = unit(3.3, "mm"),
    plot.margin = margin(5, 15, 5, 5)
  )

ggsave(
  file.path(fig3b_dir, "Fig3b_disease_similarity_inflammatory_profiles.pdf"),
  p3b,
  width = 5.4,
  height = 3.8,
  device = cairo_pdf
)

ggsave(
  file.path(fig3b_dir, "Fig3b_disease_similarity_inflammatory_profiles.png"),
  p3b,
  width = 5.4,
  height = 3.8,
  dpi = 600
)

write.csv(
  data.frame(disease = rownames(mat), mat, check.names = FALSE),
  file.path(fig3b_dir, "Fig3b_inflammatory_profile_matrix.csv"),
  row.names = FALSE
)

# ============================================================
# Fig3c. Simplified recurrent inflammatory axes
# ============================================================

fig3c_dir <- file.path(out_dir, "Fig3c")
dir.create(fig3c_dir, recursive = TRUE, showWarnings = FALSE)

fig3c_df <- l1 %>%
  dplyr::filter(sig_up) %>%
  dplyr::group_by(program, program_short, disease_group) %>%
  dplyr::summarise(
    n_enriched = dplyr::n(),
    mean_beta = mean(beta, na.rm = TRUE),
    .groups = "drop"
  )

program_summary_c <- fig3c_df %>%
  dplyr::group_by(program, program_short) %>%
  dplyr::summarise(
    total_enriched = sum(n_enriched, na.rm = TRUE),
    shared_groups = dplyr::n_distinct(disease_group),
    .groups = "drop"
  ) %>%
  dplyr::arrange(
    dplyr::desc(shared_groups),
    dplyr::desc(total_enriched)
  ) %>%
  dplyr::slice_head(n = top_n_programs_c)

fig3c_plot <- fig3c_df %>%
  dplyr::filter(program %in% program_summary_c$program) %>%
  dplyr::left_join(
    program_summary_c,
    by = c("program", "program_short")
  ) %>%
  dplyr::mutate(
    disease_group = factor(as.character(disease_group), levels = disease_group_order),
    program_short = factor(
      as.character(program_short),
      levels = rev(program_summary_c$program_short)
    )
  )

lab_c <- program_summary_c %>%
  dplyr::mutate(
    program_short = factor(
      as.character(program_short),
      levels = rev(program_summary_c$program_short)
    ),
    label = paste0(shared_groups)
  )

x_max_c <- max(program_summary_c$total_enriched, na.rm = TRUE)
if (!is.finite(x_max_c) || x_max_c <= 0) x_max_c <- 1

p3c <- ggplot(
  fig3c_plot,
  aes(x = n_enriched, y = program_short, fill = disease_group)
) +
  geom_col(
    width = 0.72,
    color = "white",
    linewidth = 0.22
  ) +
  geom_text(
    data = lab_c,
    aes(
      x = total_enriched + x_max_c * 0.04,
      y = program_short,
      label = label
    ),
    inherit.aes = FALSE,
    hjust = 0.5,
    size = 3.0,
    color = "#333333"
  ) +
  annotate(
    "text",
    x = x_max_c * 1.04,
    y = length(levels(fig3c_plot$program_short)) + 0.75,
    label = "Groups",
    hjust = 0.5,
    size = 3.0,
    fontface = "bold"
  ) +
  scale_fill_manual(
    values = group_cols,
    breaks = disease_group_order,
    labels = disease_group_labels,
    drop = FALSE,
    name = "Disease group"
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.20))
  ) +
  labs(
    title = "Recurrent inflammatory axes",
    x = "Number of enriched inflammatory factors",
    y = NULL
  ) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 9) +
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0),
    axis.text.y = element_text(size = 8.4, color = "black"),
    axis.text.x = element_text(size = 8, color = "black"),
    axis.title.x = element_text(size = 9),
    axis.line = element_line(linewidth = 0.35),
    axis.ticks = element_line(linewidth = 0.3),
    legend.position = "right",
    legend.title = element_text(size = 8, face = "bold"),
    legend.text = element_text(size = 7.5),
    legend.key.size = unit(3.2, "mm"),
    plot.margin = margin(5, 18, 5, 5)
  )

ggsave(
  file.path(fig3c_dir, "Fig3c_recurrent_inflammatory_axes.pdf"),
  p3c,
  width = 6.2,
  height = 4.0,
  device = cairo_pdf
)

ggsave(
  file.path(fig3c_dir, "Fig3c_recurrent_inflammatory_axes.png"),
  p3c,
  width = 6.2,
  height = 4.0,
  dpi = 600
)

write.csv(
  fig3c_plot,
  file.path(fig3c_dir, "Fig3c_data.csv"),
  row.names = FALSE
)

# ============================================================
# Fig3d. Immune compartment contribution to core programs
# ============================================================

fig3d_dir <- file.path(out_dir, "Fig3d")
dir.create(fig3d_dir, recursive = TRUE, showWarnings = FALSE)

core_programs_present <- core_programs_d[
  core_programs_d %in% unique(l1$program)
]

if (length(core_programs_present) == 0) {
  stop("None of the selected core programs are present in L1 data.")
}

fig3d_df <- l1 %>%
  dplyr::filter(
    sig_up,
    program %in% core_programs_present
  ) %>%
  dplyr::group_by(program, program_short, celltype) %>%
  dplyr::summarise(
    n_enriched = dplyr::n(),
    burden = sum(beta, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::group_by(program, program_short) %>%
  dplyr::mutate(
    total_enriched = sum(n_enriched, na.rm = TRUE),
    prop = n_enriched / total_enriched * 100
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    program = factor(program, levels = core_programs_present),
    program_short = factor(
      program_short,
      levels = short_program(core_programs_present)
    ),
    celltype = factor(as.character(celltype), levels = celltype_order)
  )

p3d <- ggplot(
  fig3d_df,
  aes(x = program_short, y = prop, fill = celltype)
) +
  geom_col(
    width = 0.76,
    color = "white",
    linewidth = 0.18
  ) +
  scale_fill_manual(
    values = celltype_cols,
    breaks = celltype_order,
    drop = FALSE,
    name = "Cell type"
  ) +
  scale_y_continuous(
    limits = c(0, 100),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    title = "Immune compartment contribution to inflammatory programs",
    x = NULL,
    y = "Proportion of enriched factors (%)"
  ) +
  theme_classic(base_size = 9) +
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7.5, color = "black"),
    axis.text.y = element_text(size = 8, color = "black"),
    axis.title.y = element_text(size = 9),
    axis.line = element_line(linewidth = 0.32),
    axis.ticks = element_line(linewidth = 0.25),
    legend.position = "right",
    legend.title = element_text(size = 8, face = "bold"),
    legend.text = element_text(size = 7),
    legend.key.size = unit(3, "mm"),
    plot.margin = margin(5, 5, 5, 5)
  )

ggsave(
  file.path(fig3d_dir, "Fig3d_immune_compartment_contribution_to_programs.pdf"),
  p3d,
  width = 7.2,
  height = 3.8,
  device = cairo_pdf
)

ggsave(
  file.path(fig3d_dir, "Fig3d_immune_compartment_contribution_to_programs.png"),
  p3d,
  width = 7.2,
  height = 3.8,
  dpi = 600
)

write.csv(
  fig3d_df,
  file.path(fig3d_dir, "Fig3d_data.csv"),
  row.names = FALSE
)

# ============================================================
# Load L2 LMEM result for Fig3e
# ============================================================

l2_raw <- fread(l2_path, data.table = FALSE)
colnames(l2_raw) <- trimws(colnames(l2_raw))

cat("\n[READ L2]\n")
print(l2_path)

cat("\n[Columns]\n")
print(colnames(l2_raw))

if (!"factor_feature" %in% colnames(l2_raw)) {
  l2_raw <- parse_factor_info(l2_raw)
}

cell_col_l2 <- "celltype_L2"

if (!cell_col_l2 %in% colnames(l2_raw)) {
  stop(
    "Column celltype_L2 not found. Available columns:\n",
    paste(colnames(l2_raw), collapse = ", ")
  )
}

p_col_l2 <- get_p_col(l2_raw, prefer_adj = prefer_adj_p)
message("[L2 using p column] ", p_col_l2)

required_cols_l2 <- c(
  cell_col_l2, "group1", "group0",
  "beta", "factor", "factor_feature", p_col_l2
)

missing_cols_l2 <- setdiff(required_cols_l2, colnames(l2_raw))

if (length(missing_cols_l2) > 0) {
  stop(
    "Missing L2 required columns: ",
    paste(missing_cols_l2, collapse = ", ")
  )
}

l2 <- l2_raw %>%
  dplyr::filter(group0 == "Pediatric healthy") %>%
  dplyr::mutate(
    cell_state = trimws(as.character(.data[[cell_col_l2]])),
    disease = disease_short_map[group1],
    disease_group = disease_group_map[disease],
    program = infer_program(factor_feature),
    program_short = short_program(program),
    p_use = .data[[p_col_l2]],
    sig_up = !is.na(p_use) & p_use < p_cutoff & beta > 0,
    sig_any = !is.na(p_use) & p_use < p_cutoff
  ) %>%
  dplyr::filter(
    !is.na(disease),
    !is.na(disease_group),
    !is.na(program)
  ) %>%
  dplyr::filter(
    !is.na(cell_state),
    cell_state != ""
  )

l2$disease <- factor(l2$disease, levels = disease_order)
l2$disease_group <- factor(l2$disease_group, levels = disease_group_order)

cat("\n[L2 rows]\n")
print(nrow(l2))

cat("\n[L2 program count after merging JAK terms]\n")
print(table(l2$program))

# ============================================================
# Robust Seurat metadata mapping for L2 parent cell type
# ============================================================

obj <- readRDS(seurat_rds)
meta <- as.data.frame(obj@meta.data)

l1_candidates <- c(
  "majorCluster",
  "aifi_L1_paper",
  "aifi_L1_aligned",
  "aifi_L1",
  "celltype",
  "celltype_L1"
)

l2_candidates <- c(
  "subCluster",
  "aifi_L2_paper",
  "aifi_L2_aligned",
  "aifi_L2",
  "celltype_L2",
  "subtype"
)

l1_meta_col <- pick_first_col(meta, l1_candidates)
l2_meta_col <- pick_first_col(meta, l2_candidates)

if (is.na(l1_meta_col) || is.na(l2_meta_col)) {
  cat("\n[Metadata columns]\n")
  print(colnames(meta))
  stop("Cannot find Level 1 or Level 2 annotation columns in Seurat metadata.")
}

cat("\n[Using metadata columns]\n")
cat("Level 1:", l1_meta_col, "\n")
cat("Level 2:", l2_meta_col, "\n")

state_parent_input <- tibble::tibble(
  cell_state = trimws(meta_to_chr(meta[[l2_meta_col]])),
  parent_celltype = trimws(meta_to_chr(meta[[l1_meta_col]]))
)

cat("\n[Check state_parent_input]\n")
print(head(state_parent_input))

cat("\n[cell_state class]\n")
print(class(state_parent_input$cell_state))

cat("\n[parent_celltype class]\n")
print(class(state_parent_input$parent_celltype))

state_parent <- state_parent_input %>%
  dplyr::filter(
    !is.na(cell_state),
    cell_state != "",
    !is.na(parent_celltype),
    parent_celltype != ""
  ) %>%
  dplyr::mutate(
    parent_celltype = clean_parent_celltype(parent_celltype)
  ) %>%
  dplyr::count(cell_state, parent_celltype, name = "n") %>%
  dplyr::group_by(cell_state) %>%
  dplyr::slice_max(order_by = n, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::select(cell_state, parent_celltype)

cat("\n[Parent cell type mapping]\n")
print(head(state_parent, 20))

cat("\n[Parent cell type counts]\n")
print(table(state_parent$parent_celltype))

unmatched_states <- setdiff(unique(l2$cell_state), unique(state_parent$cell_state))

cat("\n[Unmatched L2 states between LMEM and Seurat metadata]\n")
print(head(unmatched_states, 50))
cat("Unmatched count:", length(unmatched_states), "\n")

# ============================================================
# Fig3e. Simplified Level 2 inflammatory cell states
# ============================================================

fig3e_dir <- file.path(out_dir, "Fig3e")
dir.create(fig3e_dir, recursive = TRUE, showWarnings = FALSE)

if (BURDEN_MODE_L2 == "sig_any_abs_sum") {
  state_burden <- l2 %>%
    dplyr::group_by(cell_state) %>%
    dplyr::summarise(
      inflammatory_burden = sum(abs(beta[sig_any]), na.rm = TRUE),
      n_sig_up = sum(sig_up, na.rm = TRUE),
      n_sig_any = sum(sig_any, na.rm = TRUE),
      .groups = "drop"
    )
} else {
  state_burden <- l2 %>%
    dplyr::group_by(cell_state) %>%
    dplyr::summarise(
      inflammatory_burden = sum(beta[sig_up], na.rm = TRUE),
      n_sig_up = sum(sig_up, na.rm = TRUE),
      n_sig_any = sum(sig_any, na.rm = TRUE),
      .groups = "drop"
    )
}

dominant_program_e <- l2 %>%
  dplyr::filter(sig_up) %>%
  dplyr::group_by(cell_state, program, program_short) %>%
  dplyr::summarise(
    program_burden = sum(beta, na.rm = TRUE),
    program_hits = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::group_by(cell_state) %>%
  dplyr::slice_max(order_by = program_burden, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::select(
    cell_state,
    dominant_program = program_short
  )

dominant_group_e <- l2 %>%
  dplyr::filter(sig_up) %>%
  dplyr::group_by(cell_state, disease_group) %>%
  dplyr::summarise(
    group_burden = sum(beta, na.rm = TRUE),
    group_hits = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::group_by(cell_state) %>%
  dplyr::slice_max(order_by = group_burden, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    dominant_group = as.character(disease_group),
    dominant_group_label = disease_group_labels[dominant_group]
  ) %>%
  dplyr::select(
    cell_state,
    dominant_group,
    dominant_group_label
  )

fig3e_df <- state_burden %>%
  dplyr::left_join(state_parent, by = "cell_state") %>%
  dplyr::left_join(dominant_program_e, by = "cell_state") %>%
  dplyr::left_join(dominant_group_e, by = "cell_state") %>%
  dplyr::mutate(
    parent_celltype = ifelse(is.na(parent_celltype), "Other", parent_celltype),
    parent_celltype = clean_parent_celltype(parent_celltype),
    parent_celltype = ifelse(parent_celltype %in% celltype_order, parent_celltype, "Other"),
    dominant_program = ifelse(is.na(dominant_program), "-", dominant_program),
    dominant_group = ifelse(is.na(dominant_group), NA_character_, dominant_group),
    dominant_group_label = ifelse(is.na(dominant_group_label), "-", dominant_group_label)
  ) %>%
  dplyr::filter(inflammatory_burden > 0) %>%
  dplyr::arrange(
    dplyr::desc(inflammatory_burden),
    dplyr::desc(n_sig_up)
  ) %>%
  dplyr::slice_head(n = top_l2_n) %>%
  dplyr::mutate(
    cell_state_label = stringr::str_wrap(cell_state, width = 30),
    cell_state_label = factor(cell_state_label, levels = rev(cell_state_label)),
    parent_celltype = factor(parent_celltype, levels = celltype_order),
    dominant_group = factor(dominant_group, levels = disease_group_order)
  )

xmax_e <- max(fig3e_df$inflammatory_burden, na.rm = TRUE)
if (!is.finite(xmax_e) || xmax_e <= 0) xmax_e <- 1

x_prog_e <- xmax_e * 1.05
x_group_e <- xmax_e * 1.40
x_factor_e <- xmax_e * 1.70
x_limit_e <- xmax_e * 1.88

p3e <- ggplot(
  fig3e_df,
  aes(x = inflammatory_burden, y = cell_state_label, fill = parent_celltype)
) +
  geom_col(
    width = 0.72,
    color = "white",
    linewidth = 0.25
  ) +
  geom_text(
    aes(label = dominant_program),
    x = x_prog_e,
    hjust = 0,
    size = 2.7,
    color = "#333333"
  ) +
  geom_point(
    aes(color = dominant_group),
    x = x_group_e,
    shape = 15,
    size = 3.2,
    na.rm = TRUE
  ) +
  geom_text(
    aes(label = dominant_group_label),
    x = x_group_e + xmax_e * 0.035,
    hjust = 0,
    size = 2.7,
    color = "#333333"
  ) +
  geom_text(
    aes(label = paste0("n=", n_sig_up)),
    x = x_factor_e,
    hjust = 0,
    size = 2.7,
    color = "#333333"
  ) +
  annotate(
    "text",
    x = x_prog_e,
    y = top_l2_n + 0.8,
    label = "Dominant program",
    hjust = 0,
    size = 2.9,
    fontface = "bold"
  ) +
  annotate(
    "text",
    x = x_group_e,
    y = top_l2_n + 0.8,
    label = "Dominant group",
    hjust = 0,
    size = 2.9,
    fontface = "bold"
  ) +
  annotate(
    "text",
    x = x_factor_e,
    y = top_l2_n + 0.8,
    label = "Factors",
    hjust = 0,
    size = 2.9,
    fontface = "bold"
  ) +
  scale_fill_manual(
    values = celltype_cols,
    breaks = celltype_order,
    drop = FALSE,
    name = "Parent cell type"
  ) +
  scale_color_manual(
    values = group_cols,
    breaks = disease_group_order,
    labels = disease_group_labels,
    drop = FALSE,
    name = "Dominant disease group"
  ) +
  scale_x_continuous(
    limits = c(0, x_limit_e),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    title = "Inflammatory burden in Level 2 cell states",
    x = "Inflammatory activation burden",
    y = NULL
  ) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 8.5) +
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0),
    axis.text.y = element_text(size = 7.5, color = "black", lineheight = 0.9),
    axis.text.x = element_text(size = 7.5, color = "black"),
    axis.title.x = element_text(size = 8.8),
    axis.line = element_line(linewidth = 0.32),
    axis.ticks = element_line(linewidth = 0.25),
    legend.position = "right",
    legend.title = element_text(size = 8, face = "bold"),
    legend.text = element_text(size = 7),
    legend.key.size = unit(3, "mm"),
    plot.margin = margin(5, 20, 5, 5)
  ) +
  guides(
    fill = guide_legend(order = 1, ncol = 1),
    color = guide_legend(order = 2, override.aes = list(shape = 15, size = 4))
  )

ggsave(
  file.path(fig3e_dir, "Fig3e_level2_inflammatory_cell_states.pdf"),
  p3e,
  width = 9.3,
  height = 4.7,
  device = cairo_pdf
)

ggsave(
  file.path(fig3e_dir, "Fig3e_level2_inflammatory_cell_states.png"),
  p3e,
  width = 9.3,
  height = 4.7,
  dpi = 600
)

write.csv(
  fig3e_df,
  file.path(fig3e_dir, "Fig3e_data.csv"),
  row.names = FALSE
)

write.csv(
  state_burden,
  file.path(fig3e_dir, "Fig3e_all_L2_burden.csv"),
  row.names = FALSE
)

# ============================================================
# Done
# ============================================================

cat("\n[DONE]\n")
cat("Output dir:\n", out_dir, "\n\n")

cat("Generated:\n")
cat("Fig3b/Fig3b_disease_similarity_inflammatory_profiles.png\n")
cat("Fig3c/Fig3c_recurrent_inflammatory_axes.png\n")
cat("Fig3d/Fig3d_immune_compartment_contribution_to_programs.png\n")
cat("Fig3e/Fig3e_level2_inflammatory_cell_states.png\n")

cat("\nJAK terms merged as: JAK STAT signaling\n")
cat("P column preference: ", ifelse(prefer_adj_p, "adj_p first", "p_value first"), "\n")
cat("L2 burden mode: ", BURDEN_MODE_L2, "\n")