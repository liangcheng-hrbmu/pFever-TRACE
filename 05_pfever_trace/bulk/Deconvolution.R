library(Matrix)
library(data.table)
library(InstaPrism)

ref_dir <- "data/analysis/09_bulk_deconvolution/03_instaprism_ref_input_L1type_L2state"
bulk_path <- "data/analysis/09_bulk_deconvolution/E-MTAB-11671_mixture_counts_filtered.tsv"
outdir <- "data/analysis/09_bulk_deconvolution/06_instaprism_run_L1type_L2state"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

read_expr_table <- function(path) {
  dt <- fread(path, data.table = FALSE)
  rn <- dt[[1]]
  mat <- as.matrix(dt[, -1, drop = FALSE])
  rownames(mat) <- rn
  storage.mode(mat) <- "numeric"
  mat
}

clean_bulk <- function(mat) {
  mat <- mat[!is.na(rownames(mat)) & rownames(mat) != "", , drop = FALSE]
  mat <- mat[!duplicated(rownames(mat)), , drop = FALSE]
  mat <- mat[rowSums(mat, na.rm = TRUE) > 0, , drop = FALSE]
  mat
}

get_frac_from_deconv <- function(deconv_res) {
  x <- try(t(deconv_res@Post.ini.ct@theta), silent = TRUE)
  if (!inherits(x, "try-error")) return(x)

  x <- try(t(slot(deconv_res, "theta")), silent = TRUE)
  if (!inherits(x, "try-error")) return(x)

  stop("Cannot find theta in deconv_res")
}

cat("[READ] single-cell reference inputs\n")
sc_counts <- readMM(file.path(ref_dir, "sc_counts_for_instaprism.mtx"))
sc_counts <- as(sc_counts, "dgCMatrix")

genes <- fread(file.path(ref_dir, "sc_genes_for_instaprism.tsv.gz"), data.table = FALSE)
cells <- fread(file.path(ref_dir, "sc_cells_for_instaprism.tsv.gz"), data.table = FALSE)
meta  <- fread(file.path(ref_dir, "sc_meta_for_instaprism.tsv.gz"), data.table = FALSE)

rownames(sc_counts) <- genes$gene
colnames(sc_counts) <- cells$cell_id

meta <- meta[match(colnames(sc_counts), meta$cell_id), , drop = FALSE]
stopifnot(all(colnames(sc_counts) == meta$cell_id))

cat("[INFO] sc_counts dim:", nrow(sc_counts), "genes x", ncol(sc_counts), "cells\n")
cat("[INFO] n cell types:", length(unique(meta$cell_type)), "\n")
cat("[INFO] n cell states:", length(unique(meta$cell_state)), "\n")

cat("[READ] bulk\n")
bulk <- read_expr_table(bulk_path)
bulk <- round(bulk)
bulk <- clean_bulk(bulk)

cat("[INFO] bulk dim:", nrow(bulk), "genes x", ncol(bulk), "samples\n")

common_genes <- intersect(rownames(sc_counts), rownames(bulk))
cat("[INFO] common genes:", length(common_genes), "\n")
if (length(common_genes) < 1000) {
  stop("Too few common genes: ", length(common_genes))
}

sc_counts_use <- sc_counts[common_genes, , drop = FALSE]
bulk_use <- bulk[common_genes, , drop = FALSE]

cell_type_labels  <- as.character(meta$cell_type)
cell_state_labels <- as.character(meta$cell_state)

cat("[RUN] refPrepare\n")
refPhi_obj <- tryCatch(
  {
    refPrepare(
      sc_Expr = sc_counts_use,
      cell.type.labels = cell_type_labels,
      cell.state.labels = cell_state_labels
    )
  },
  error = function(e) {
    message("[WARN] sparse matrix failed in refPrepare, retry with dense matrix")
    refPrepare(
      sc_Expr = as.matrix(sc_counts_use),
      cell.type.labels = cell_type_labels,
      cell.state.labels = cell_state_labels
    )
  }
)

saveRDS(refPhi_obj, file = file.path(outdir, "InstaPrism_refPhi_obj.rds"))

cat("[RUN] InstaPrism\n")
deconv_res <- InstaPrism(
  bulk_Expr = bulk_use,
  refPhi_cs = refPhi_obj
)

saveRDS(deconv_res, file = file.path(outdir, "InstaPrism_deconv_res_allgenes.rds"))

cat("[SAVE] fractions\n")
frac <- get_frac_from_deconv(deconv_res)
write.table(
  frac,
  file = file.path(outdir, "InstaPrism_cell_fraction_allgenes.tsv"),
  sep = "\t", quote = FALSE, col.names = NA
)

cat("[INFO] fraction row sums:\n")
print(summary(rowSums(frac)))

cat("[SAVE] Z array\n")
Z <- get_Z_array(deconv_res)
saveRDS(Z, file = file.path(outdir, "InstaPrism_Z_allgenes.rds"))

cat("[INFO] dim(Z):\n")
print(dim(Z))
cat("[INFO] cell types in Z:\n")
print(dimnames(Z)[[3]])

for (ct in dimnames(Z)[[3]]) {
  mat <- Z[, , ct, drop = FALSE]
  mat <- mat[, , 1]
  mat <- t(mat)
  safe_ct <- gsub("[^A-Za-z0-9_]+", "_", ct)

  write.table(
    mat,
    file = file.path(outdir, paste0("Z_", safe_ct, "_allgenes_genes_by_samples.tsv")),
    sep = "\t", quote = FALSE, col.names = NA
  )
}

summary_df <- data.frame(
  metric = c(
    "n_sc_genes",
    "n_sc_cells",
    "n_bulk_samples",
    "n_common_genes",
    "n_cell_types",
    "n_cell_states"
  ),
  value = c(
    nrow(sc_counts),
    ncol(sc_counts),
    ncol(bulk),
    length(common_genes),
    length(unique(cell_type_labels)),
    length(unique(cell_state_labels))
  )
)

write.table(
  summary_df,
  file = file.path(outdir, "InstaPrism_run_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

cat("[DONE]\n")
cat("[OUTDIR]", outdir, "\n")




frac <- read.delim(
  "data/analysis/09_bulk_deconvolution/06_instaprism_run_L1type_L2state/InstaPrism_cell_fraction_allgenes.tsv",
  check.names = FALSE,
  row.names = 1
)

Z <- readRDS(
  "data/analysis/09_bulk_deconvolution/06_instaprism_run_L1type_L2state/InstaPrism_Z_allgenes.rds"
)

summary(frac)
sort(apply(frac, 2, median, na.rm = TRUE), decreasing = TRUE)
dim(Z)
dimnames(Z)[[3]]




celltypes <- dimnames(Z)[[3]]

check_one_gene <- function(g) {
  x <- sapply(celltypes, function(ct) mean(Z[, g, ct], na.rm = TRUE))
  sort(x, decreasing = TRUE)
}

check_one_gene("FCER1A")
check_one_gene("CLEC10A")
check_one_gene("CD1C")

check_one_gene("NKG7")
check_one_gene("CCL5")
check_one_gene("GZMK")
check_one_gene("CTSW")
