#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(edgeR)
  library(optparse)
  library(data.table)
})

read_counts_gz <- function(path) {
  dt <- data.table::fread(cmd = paste("zcat", shQuote(path)), data.table = FALSE)
  gene <- dt[[1]]
  mat <- as.matrix(dt[, -1, drop=FALSE])
  rownames(mat) <- gene
  storage.mode(mat) <- "numeric"
  mat
}

# build a full-rank design; drop covariates that cause rank deficiency
build_design_fullrank <- function(disease_f, sex_f, age_f, dataset_f, include_dataset) {
  # droplevels to remove unused levels in THIS subset
  disease_f  <- droplevels(disease_f)
  sex_f      <- droplevels(sex_f)
  age_f      <- droplevels(age_f)
  dataset_f  <- droplevels(dataset_f)

  # start with candidate terms
  use_sex <- nlevels(sex_f) >= 2
  use_age <- nlevels(age_f) >= 2
  use_ds  <- include_dataset && nlevels(dataset_f) >= 2

  # helper to try build and check rank
  try_build <- function(use_sex, use_age, use_ds) {
    if (use_ds && use_sex && use_age) {
      design <- model.matrix(~ 0 + disease_f + sex_f + age_f + dataset_f)
    } else if (use_ds && use_sex && !use_age) {
      design <- model.matrix(~ 0 + disease_f + sex_f + dataset_f)
    } else if (use_ds && !use_sex && use_age) {
      design <- model.matrix(~ 0 + disease_f + age_f + dataset_f)
    } else if (use_ds && !use_sex && !use_age) {
      design <- model.matrix(~ 0 + disease_f + dataset_f)
    } else if (!use_ds && use_sex && use_age) {
      design <- model.matrix(~ 0 + disease_f + sex_f + age_f)
    } else if (!use_ds && use_sex && !use_age) {
      design <- model.matrix(~ 0 + disease_f + sex_f)
    } else if (!use_ds && !use_sex && use_age) {
      design <- model.matrix(~ 0 + disease_f + age_f)
    } else {
      design <- model.matrix(~ 0 + disease_f)
    }
    colnames(design) <- make.names(colnames(design))
    r <- qr(design)$rank
    if (r == ncol(design)) {
      return(list(ok=TRUE, design=design, use_sex=use_sex, use_age=use_age, use_ds=use_ds))
    } else {
      return(list(ok=FALSE, design=design, use_sex=use_sex, use_age=use_age, use_ds=use_ds))
    }
  }

  # attempt hierarchy: full -> drop dataset -> drop age -> drop sex -> disease only
  attempts <- list(
    c(use_sex, use_age, use_ds),
    c(use_sex, use_age, FALSE),
    c(use_sex, FALSE, FALSE),
    c(FALSE, use_age, FALSE),
    c(FALSE, FALSE, FALSE)
  )

  for (a in attempts) {
    res <- try_build(a[1], a[2], a[3])
    if (res$ok) return(res)
  }

  # last resort: force disease only
  res <- try_build(FALSE, FALSE, FALSE)
  return(res)
}

opt_list <- list(
  make_option("--counts", type="character", help="pseudobulk_counts.tsv.gz (gene x sample)"),
  make_option("--meta", type="character", help="pseudobulk_meta.tsv"),
  make_option("--out_dir", type="character", help="output directory"),
  make_option("--healthy_label", type="character", default="Pediatric healthy"),
  make_option("--disease_col", type="character", default="disease"),
  make_option("--celltype_col", type="character", default="Level1"),
  make_option("--sex_col", type="character", default="sex"),
  make_option("--binned_age_col", type="character", default="binned_age"),
  make_option("--dataset_col", type="character", default="dataset_id"),
  make_option("--min_samples_per_group", type="integer", default=3),
  make_option("--fdr", type="double", default=0.01),
  make_option("--logfc", type="double", default=1.5),
  make_option("--include_dataset_as_covariate", action="store_true", default=FALSE),
  make_option("--only_union", action="store_true", default=FALSE)
)
opt <- parse_args(OptionParser(option_list = opt_list))

stopifnot(file.exists(opt$counts), file.exists(opt$meta))
dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(opt$out_dir, "by_celltype"), recursive = TRUE, showWarnings = FALSE)

cat("[READ]", opt$counts, "\n")
counts <- read_counts_gz(opt$counts)

cat("[READ]", opt$meta, "\n")
meta <- data.table::fread(opt$meta, data.table = FALSE, sep = "\t")

pb_ids <- colnames(counts)
meta_idx <- match(pb_ids, meta$pb_id)
if (any(is.na(meta_idx))) {
  missing <- pb_ids[is.na(meta_idx)]
  stop("Some samples in counts not found in meta$pb_id, e.g.: ", paste(head(missing, 5), collapse=", "))
}
meta <- meta[meta_idx, , drop=FALSE]

# clean
for (cc in c(opt$disease_col, opt$celltype_col, opt$sex_col, opt$binned_age_col, opt$dataset_col)) {
  if (cc %in% colnames(meta)) meta[[cc]] <- trimws(as.character(meta[[cc]]))
}

if (!(opt$healthy_label %in% unique(meta[[opt$disease_col]]))) {
  cat("[WARN] healthy_label not found. Available:\n")
  print(sort(unique(meta[[opt$disease_col]])))
  stop("healthy_label not found: ", opt$healthy_label)
}

celltypes <- sort(unique(meta[[opt$celltype_col]]))
cat("[INFO] n_celltypes =", length(celltypes), "\n")

deg_union <- character(0)
union_records <- list()

for (ct in celltypes) {
  idx_ct <- which(meta[[opt$celltype_col]] == ct)
  meta_ct <- meta[idx_ct, , drop=FALSE]
  if (nrow(meta_ct) < (opt$min_samples_per_group * 2)) {
    cat("[SKIP] celltype=", ct, " too few samples:", nrow(meta_ct), "\n")
    next
  }
  counts_ct <- counts[, idx_ct, drop=FALSE]

  disease_all <- meta_ct[[opt$disease_col]]
  diseases <- sort(setdiff(unique(disease_all), opt$healthy_label))
  if (length(diseases) == 0) {
    cat("[SKIP] celltype=", ct, " has only healthy.\n")
    next
  }

  ct_dir <- file.path(opt$out_dir, "by_celltype", gsub("[^A-Za-z0-9._-]+", "_", ct))
  dir.create(ct_dir, recursive = TRUE, showWarnings = FALSE)

  for (dis in diseases) {
    keep <- which(disease_all %in% c(opt$healthy_label, dis))
    meta_pair <- meta_ct[keep, , drop=FALSE]
    counts_pair <- counts_ct[, keep, drop=FALSE]

    n_h <- sum(meta_pair[[opt$disease_col]] == opt$healthy_label, na.rm=TRUE)
    n_d <- sum(meta_pair[[opt$disease_col]] == dis, na.rm=TRUE)
    if (n_h < opt$min_samples_per_group || n_d < opt$min_samples_per_group) {
      cat("[SKIP] ", ct, " | ", dis, " group size too small (H=", n_h, ", D=", n_d, ")\n", sep="")
      next
    }

    # drop missing covariates (sex/age; dataset if requested)
    ok <- !is.na(meta_pair[[opt$sex_col]]) & meta_pair[[opt$sex_col]] != "" &
          !is.na(meta_pair[[opt$binned_age_col]]) & meta_pair[[opt$binned_age_col]] != ""
    if (opt$include_dataset_as_covariate) {
      ok <- ok & !is.na(meta_pair[[opt$dataset_col]]) & meta_pair[[opt$dataset_col]] != ""
    }
    if (sum(ok) < (opt$min_samples_per_group * 2)) {
      cat("[SKIP] ", ct, " | ", dis, " after dropping missing covariates too few samples.\n", sep="")
      next
    }
    meta_pair <- meta_pair[ok, , drop=FALSE]
    counts_pair <- counts_pair[, ok, drop=FALSE]

    # factors (droplevels happens inside build_design_fullrank)
    disease_f <- factor(meta_pair[[opt$disease_col]], levels=c(opt$healthy_label, dis))
    sex_f <- factor(meta_pair[[opt$sex_col]])
    age_f <- factor(meta_pair[[opt$binned_age_col]])
    dataset_f <- factor(meta_pair[[opt$dataset_col]])

    # DGEList
    y <- DGEList(counts = counts_pair, group = disease_f)

    keep_genes <- filterByExpr(y, group = disease_f)
    y <- y[keep_genes, , keep.lib.sizes = FALSE]
    y <- calcNormFactors(y, method="TMM", logratioTrim=0.3)

    des_res <- build_design_fullrank(disease_f, sex_f, age_f, dataset_f, opt$include_dataset_as_covariate)
    design <- des_res$design

    if (!des_res$use_ds && opt$include_dataset_as_covariate) {
      cat("[WARN] ", ct, " | ", dis, " dropped dataset covariate (not enough levels or rank issue)\n", sep="")
    }
    if (!des_res$use_age) {
      cat("[WARN] ", ct, " | ", dis, " dropped binned_age (not enough levels or rank issue)\n", sep="")
    }
    if (!des_res$use_sex) {
      cat("[WARN] ", ct, " | ", dis, " dropped sex (not enough levels or rank issue)\n", sep="")
    }

    y <- estimateDisp(y, design)
    fit <- glmQLFit(y, design, robust=TRUE)

    c_healthy <- make.names(paste0("disease_f", opt$healthy_label))
    c_disease <- make.names(paste0("disease_f", dis))
    if (!(c_healthy %in% colnames(design)) || !(c_disease %in% colnames(design))) {
      cat("[WARN] design missing disease cols for ", ct, " | ", dis, "\n", sep="")
      next
    }

    contrast_vec <- rep(0, ncol(design))
    names(contrast_vec) <- colnames(design)
    contrast_vec[c_disease] <- 1
    contrast_vec[c_healthy] <- -1

    qlf <- glmQLFTest(fit, contrast = contrast_vec)

    tab <- topTags(qlf, n = Inf)$table
    tab$gene <- rownames(tab)
    tab$celltype <- ct
    tab$disease_comp <- paste0(dis, "_vs_", opt$healthy_label)
    tab$n_healthy <- sum(disease_f == opt$healthy_label)
    tab$n_disease <- sum(disease_f == dis)

    deg <- tab[!is.na(tab$FDR) & tab$FDR < opt$fdr & abs(tab$logFC) > opt$logfc, , drop=FALSE]

    cat("[OK] ", ct, " | ", dis,
        " genes_tested=", nrow(tab),
        " DEGs=", nrow(deg),
        " (FDR<", opt$fdr, ", |logFC|>", opt$logfc, ")\n", sep="")

    if (!opt$only_union) {
      out_full <- file.path(ct_dir, paste0(gsub("[^A-Za-z0-9._-]+","_",dis), "_vs_healthy.edgeR.full.tsv"))
      out_deg  <- file.path(ct_dir, paste0(gsub("[^A-Za-z0-9._-]+","_",dis), "_vs_healthy.edgeR.DEG.tsv"))
      data.table::fwrite(tab, out_full, sep="\t")
      data.table::fwrite(deg, out_deg, sep="\t")
    }

    if (nrow(deg) > 0) {
      deg_union <- union(deg_union, deg$gene)
      union_records[[length(union_records) + 1]] <- deg
    }
  }
}

union_path <- file.path(opt$out_dir, sprintf("DEG_union_FDR%.3g_abslogFC%.3g.tsv", opt$fdr, opt$logfc))
union_list_path <- file.path(opt$out_dir, sprintf("DEG_union_FDR%.3g_abslogFC%.3g.txt", opt$fdr, opt$logfc))

if (length(union_records) > 0) {
  union_df <- do.call(rbind, union_records)
  union_df <- union_df[order(union_df$FDR, -abs(union_df$logFC)), ]
  union_df <- union_df[!duplicated(union_df$gene), ]
  data.table::fwrite(union_df, union_path, sep="\t")
} else {
  data.table::fwrite(data.frame(), union_path, sep="\t")
}

writeLines(sort(unique(deg_union)), con = union_list_path)

cat("\n[SUMMARY]\n")
cat("DEG union n_genes =", length(unique(deg_union)), "\n")
cat("[WRITE]", union_path, "\n")
cat("[WRITE]", union_list_path, "\n")
cat("[DONE]\n")
