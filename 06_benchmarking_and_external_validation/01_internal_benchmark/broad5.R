suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(glmnet)
})

# =========================================================
# Internal matched CV benchmark
# Med139 ridge multinomial baseline versus pFever internal CV
#
# Design:
#   Same held out samples as pFever internal split
#   Train Med139 ridge model on remaining E MTAB 11671 samples
#   Predict the same held out samples
#   Compare paired predictions with pFever majority vote
#
# Final internal benchmark:
#   Med139 imbalance only lambda.min
#   vs
#   pFever internal majority vote
# =========================================================

# -------------------------------
# 0. Paths
# -------------------------------

BASE <- "data/analysis/09_bulk_deconvolution"
EXT_ROOT <- file.path(BASE, "12_external_validation_multi_cohort_assets_v1")

train_expr_fp <- file.path(BASE, "E-MTAB-11671_mixture_counts_filtered.tsv")

pfever_l1_pred_fp <- file.path(
  BASE,
  "08_bulk_classifier/03_patient_classifier_broad5_bulkZ_L1x11_final_withLinearSVC/query_level1_predictions_all.tsv"
)

signature_genes_fp <- file.path(
  EXT_ROOT,
  "results/E-MTAB-12793/comparison_ours_vs_145_external/wholeblood_signature_genes_used.tsv"
)

out_dir <- file.path(
  BASE,
  "08_bulk_classifier/04_internal_matched_cv_med139_vs_pfever"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(train_expr_fp)) stop("Missing train expression file: ", train_expr_fp)
if (!file.exists(pfever_l1_pred_fp)) stop("Missing pFever L1 prediction file: ", pfever_l1_pred_fp)
if (!file.exists(signature_genes_fp)) stop("Missing signature gene file: ", signature_genes_fp)

# -------------------------------
# 1. Helper functions
# -------------------------------

norm_broad5 <- function(x) {
  x <- trimws(as.character(x))
  xl <- tolower(x)

  out <- rep(NA_character_, length(x))

  out[xl %in% c("viral infection", "viral", "virus", "viral infections", "viral inf")] <- "Viral infection"
  out[xl %in% c("bacterial infection", "bacterial", "bacteria", "bacterial infections", "bact.inf", "bact inf")] <- "Bacterial infection"
  out[xl %in% c("malaria", "plasmodium")] <- "Malaria"
  out[xl %in% c("hyperinflammatory disease", "hyperinflammatory", "hyperinflam", "hyperinflammation")] <- "Hyperinflammatory disease"
  out[xl %in% c("autoimmune disease", "autoimmune", "autoimmunity")] <- "Autoimmune disease"
  out[xl %in% c("healthy", "normal", "control", "pediatric healthy")] <- "Healthy"
  out[xl %in% c("tie", "not assigned", "unassigned", "unknown")] <- "Not assigned"

  out[is.na(out) & grepl(
    "influenza|iav|respiratory syncytial|rsv|adeno|adenovirus|rhinovirus|rhino|dengue|denv|sars|covid|hfmd|enterovirus",
    xl
  )] <- "Viral infection"

  out[is.na(out) & grepl(
    "mening|pneumo|staph|aureus|streptococcal|streptococcus|\\bgas\\b|escherichia|e\\.coli|e coli|ecoli|tuberculosis|\\btb\\b|sepsis",
    xl
  )] <- "Bacterial infection"

  out[is.na(out) & grepl(
    "kawasaki|\\bkd\\b|mis-c|misc|multisystem inflammatory",
    xl
  )] <- "Hyperinflammatory disease"

  out[is.na(out) & grepl(
    "juvenile idiopathic|juvenile idopathic|\\bjia\\b|arthritis|systemic lupus|\\bsle\\b",
    xl
  )] <- "Autoimmune disease"

  out[is.na(out) & grepl("malaria", xl)] <- "Malaria"

  out
}

read_expr <- function(fp) {
  expr <- data.table::fread(fp, data.table = FALSE, check.names = FALSE)

  gene_col <- colnames(expr)[1]
  genes <- trimws(as.character(expr[[gene_col]]))

  expr[[gene_col]] <- NULL

  mat <- as.matrix(expr)
  mode(mat) <- "numeric"
  rownames(mat) <- genes

  keep <- !is.na(rownames(mat)) & rownames(mat) != ""
  mat <- mat[keep, , drop = FALSE]

  if (anyDuplicated(rownames(mat)) > 0) {
    dt <- data.table::as.data.table(mat, keep.rownames = "gene")
    dt <- dt[, lapply(.SD, mean, na.rm = TRUE), by = gene]
    genes2 <- dt$gene
    dt$gene <- NULL

    mat <- as.matrix(dt)
    mode(mat) <- "numeric"
    rownames(mat) <- genes2
  }

  mat
}

auto_log1p <- function(mat) {
  q99 <- suppressWarnings(as.numeric(stats::quantile(mat, 0.99, na.rm = TRUE)))
  mx <- suppressWarnings(max(mat, na.rm = TRUE))

  if (is.finite(q99) && is.finite(mx) && (q99 > 100 || mx > 1000)) {
    message("[INFO] auto log1p applied to expression matrix.")
    return(log1p(mat))
  }

  message("[INFO] expression matrix kept as is.")
  mat
}

read_signature_genes <- function(fp) {
  sig <- data.table::fread(fp, data.table = FALSE)

  candidate_cols <- c(
    "gene", "genes", "gene_symbol", "Gene", "symbol", "Symbol",
    "feature", "features", "transcript", "probe", "gene_name"
  )

  gene_col <- intersect(candidate_cols, colnames(sig))[1]

  if (is.na(gene_col)) {
    gene_col <- colnames(sig)[1]
  }

  genes <- unique(trimws(as.character(sig[[gene_col]])))
  genes <- genes[!is.na(genes) & genes != ""]

  genes
}

majority_vote <- function(x) {
  x <- x[!is.na(x) & x != ""]

  if (length(x) == 0) {
    return("Not assigned")
  }

  tb <- sort(table(x), decreasing = TRUE)

  if (length(tb) >= 2 && tb[1] == tb[2]) {
    return("Not assigned")
  }

  names(tb)[1]
}

calc_metrics <- function(truth, pred, labels = NULL) {
  truth <- as.character(truth)
  pred <- as.character(pred)

  if (is.null(labels)) {
    labels <- sort(unique(truth))
  }

  labels <- labels[labels %in% unique(truth)]

  acc <- mean(truth == pred)

  per_class_list <- lapply(labels, function(cls) {
    tp <- sum(truth == cls & pred == cls)
    fp <- sum(truth != cls & pred == cls)
    fn <- sum(truth == cls & pred != cls)
    support <- sum(truth == cls)

    precision <- ifelse(tp + fp == 0, 0, tp / (tp + fp))
    recall <- ifelse(tp + fn == 0, 0, tp / (tp + fn))
    f1 <- ifelse(
      precision + recall == 0,
      0,
      2 * precision * recall / (precision + recall)
    )

    data.frame(
      class = cls,
      precision = precision,
      recall = recall,
      F1 = f1,
      support = support,
      stringsAsFactors = FALSE
    )
  })

  per_class <- dplyr::bind_rows(per_class_list)

  weighted_F1 <- sum(per_class$F1 * per_class$support) / sum(per_class$support)
  macro_F1 <- mean(per_class$F1)
  balanced_accuracy <- mean(per_class$recall)

  all_labs <- sort(unique(c(truth, pred)))

  cm <- table(
    factor(truth, levels = all_labs),
    factor(pred, levels = all_labs)
  )

  t_sum <- rowSums(cm)
  p_sum <- colSums(cm)
  s <- sum(cm)
  c_val <- sum(diag(cm))

  numerator <- c_val * s - sum(t_sum * p_sum)
  denominator <- sqrt((s^2 - sum(p_sum^2)) * (s^2 - sum(t_sum^2)))
  mcc <- ifelse(denominator == 0, NA_real_, numerator / denominator)

  list(
    metrics = data.frame(
      accuracy = acc,
      weighted_F1 = weighted_F1,
      macro_F1 = macro_F1,
      balanced_accuracy = balanced_accuracy,
      MCC = mcc,
      n = length(truth),
      n_classes = length(labels),
      stringsAsFactors = FALSE
    ),
    per_class = per_class
  )
}

mcnemar_exact <- function(a_correct, b_correct) {
  both_correct <- sum(a_correct & b_correct)
  both_wrong <- sum(!a_correct & !b_correct)
  a_correct_b_wrong <- sum(a_correct & !b_correct)
  a_wrong_b_correct <- sum(!a_correct & b_correct)

  discordant <- a_correct_b_wrong + a_wrong_b_correct

  if (discordant == 0) {
    p <- 1
  } else {
    p <- stats::binom.test(
      x = min(a_correct_b_wrong, a_wrong_b_correct),
      n = discordant,
      p = 0.5,
      alternative = "two.sided"
    )$p.value
  }

  data.frame(
    both_correct = both_correct,
    both_wrong = both_wrong,
    med_correct_pfever_wrong = a_correct_b_wrong,
    med_wrong_pfever_correct = a_wrong_b_correct,
    discordant_total = discordant,
    mcnemar_exact_p = p,
    stringsAsFactors = FALSE
  )
}

make_class_weights <- function(y) {
  tab <- table(y)
  raw <- 1 / as.numeric(tab)
  names(raw) <- names(tab)

  w <- raw[as.character(y)]
  w <- as.numeric(w)
  w <- w / mean(w)

  w
}

# -------------------------------
# 2. Read pFever internal predictions
# -------------------------------

message("[READ] pFever Level1 internal predictions: ", pfever_l1_pred_fp)

l1 <- data.table::fread(pfever_l1_pred_fp, data.table = FALSE)

required_cols <- c("sample_id", "Level1", "disease_broad", "pred_disease", "split")
missing_cols <- setdiff(required_cols, colnames(l1))

if (length(missing_cols) > 0) {
  stop("Missing required columns in pFever L1 prediction file: ", paste(missing_cols, collapse = ", "))
}

l1 <- l1 %>%
  dplyr::mutate(
    sample_id = trimws(as.character(sample_id)),
    Level1 = trimws(as.character(Level1)),
    truth = norm_broad5(disease_broad),
    pred_level1 = norm_broad5(pred_disease),
    split = as.integer(split)
  )

primary_labels <- c(
  "Viral infection",
  "Bacterial infection",
  "Malaria",
  "Hyperinflammatory disease",
  "Autoimmune disease"
)

l1 <- l1 %>%
  dplyr::filter(truth %in% primary_labels)

pfever_majority <- l1 %>%
  dplyr::group_by(split, sample_id, truth) %>%
  dplyr::summarise(
    pfever_prediction = majority_vote(pred_level1),
    n_level1_voted = dplyr::n_distinct(Level1),
    .groups = "drop"
  )

message("[INFO] pFever internal majority samples: ", nrow(pfever_majority))
message("[INFO] pFever splits: ", paste(sort(unique(pfever_majority$split)), collapse = ", "))

data.table::fwrite(
  pfever_majority,
  file.path(out_dir, "internal_pfever_majority_from_level1.tsv"),
  sep = "\t"
)

sample_label <- pfever_majority %>%
  dplyr::distinct(sample_id, truth)

if (anyDuplicated(sample_label$sample_id) > 0) {
  stop("Duplicated sample_id with inconsistent labels detected.")
}

# -------------------------------
# 3. Read expression and Med genes
# -------------------------------

message("[READ] training expression: ", train_expr_fp)

expr <- read_expr(train_expr_fp)
expr <- auto_log1p(expr)

message("[READ] Med signature genes: ", signature_genes_fp)

sig_genes <- read_signature_genes(signature_genes_fp)
use_genes <- intersect(sig_genes, rownames(expr))

message("[INFO] signature genes raw: ", length(sig_genes))
message("[INFO] signature genes available in training expression: ", length(use_genes))

if (length(use_genes) < 20) {
  stop("Too few Med signature genes available. Please check signature gene file.")
}

data.table::fwrite(
  data.frame(gene = use_genes, stringsAsFactors = FALSE),
  file.path(out_dir, "med_signature_genes_available_in_internal_training.tsv"),
  sep = "\t"
)

all_samples <- intersect(sample_label$sample_id, colnames(expr))

message("[INFO] labeled samples: ", nrow(sample_label))
message("[INFO] samples available in expression: ", length(all_samples))

if (length(all_samples) < 50) {
  stop("Too few samples matched between labels and expression.")
}

sample_label <- sample_label %>%
  dplyr::filter(sample_id %in% all_samples)

pfever_majority <- pfever_majority %>%
  dplyr::filter(sample_id %in% all_samples)

expr_use <- expr[use_genes, all_samples, drop = FALSE]

# -------------------------------
# 4. Internal matched CV for Med139 ridge
# -------------------------------

set.seed(1)

med_pred_rows <- list()
med_split_metrics_rows <- list()
med_per_class_rows <- list()
med_best_rows <- list()

splits <- sort(unique(pfever_majority$split))

for (sp in splits) {
  message("\n[RUN] split ", sp)

  test_samples <- pfever_majority %>%
    dplyr::filter(split == sp) %>%
    dplyr::pull(sample_id) %>%
    unique()

  train_samples <- setdiff(sample_label$sample_id, test_samples)

  train_lab <- sample_label %>%
    dplyr::filter(sample_id %in% train_samples) %>%
    dplyr::arrange(match(sample_id, train_samples))

  test_lab <- sample_label %>%
    dplyr::filter(sample_id %in% test_samples) %>%
    dplyr::arrange(match(sample_id, test_samples))

  X_train <- t(expr_use[, train_lab$sample_id, drop = FALSE])
  X_test <- t(expr_use[, test_lab$sample_id, drop = FALSE])

  y_train <- factor(train_lab$truth, levels = primary_labels)
  y_test <- factor(test_lab$truth, levels = primary_labels)

  keep_train <- !is.na(y_train)
  keep_test <- !is.na(y_test)

  X_train <- X_train[keep_train, , drop = FALSE]
  y_train <- droplevels(y_train[keep_train])
  train_lab <- train_lab[keep_train, , drop = FALSE]

  X_test <- X_test[keep_test, , drop = FALSE]
  y_test <- droplevels(y_test[keep_test])
  test_lab <- test_lab[keep_test, , drop = FALSE]

  w_train <- make_class_weights(y_train)

  class_counts <- table(y_train)
  nfolds_inner <- min(5, min(class_counts))
  nfolds_inner <- max(2, nfolds_inner)

  message("[INFO] split ", sp, " train n=", length(y_train), ", test n=", length(y_test))
  message("[INFO] split ", sp, " train class counts: ", paste(names(class_counts), class_counts, sep = "=", collapse = "; "))
  message("[INFO] split ", sp, " inner nfolds=", nfolds_inner)

  cvfit <- glmnet::cv.glmnet(
    x = X_train,
    y = y_train,
    family = "multinomial",
    alpha = 0,
    weights = w_train,
    type.measure = "class",
    nfolds = nfolds_inner,
    standardize = TRUE
  )

  pred <- predict(
    cvfit,
    newx = X_test,
    s = "lambda.min",
    type = "class"
  )

  pred <- as.character(pred[, 1])

  pred_df <- data.frame(
    split = sp,
    sample_id = test_lab$sample_id,
    truth = as.character(y_test),
    med139_prediction = pred,
    med139_lambda_min = as.numeric(cvfit$lambda.min),
    med139_lambda_1se = as.numeric(cvfit$lambda.1se),
    n_med_genes_used = length(use_genes),
    stringsAsFactors = FALSE
  )

  med_pred_rows[[as.character(sp)]] <- pred_df

  met <- calc_metrics(
    truth = pred_df$truth,
    pred = pred_df$med139_prediction,
    labels = primary_labels
  )

  split_met <- met$metrics %>%
    dplyr::mutate(
      split = sp,
      model = "med139_ridge_imbalance_only_lambda_min",
      n_med_genes_used = length(use_genes),
      .before = 1
    )

  med_split_metrics_rows[[as.character(sp)]] <- split_met

  per_cls <- met$per_class %>%
    dplyr::mutate(
      split = sp,
      model = "med139_ridge_imbalance_only_lambda_min",
      .before = 1
    )

  med_per_class_rows[[as.character(sp)]] <- per_cls

  med_best_rows[[as.character(sp)]] <- data.frame(
    split = sp,
    n_train = length(y_train),
    n_test = length(y_test),
    n_med_genes_used = length(use_genes),
    lambda_min = as.numeric(cvfit$lambda.min),
    lambda_1se = as.numeric(cvfit$lambda.1se),
    inner_nfolds = nfolds_inner,
    stringsAsFactors = FALSE
  )
}

med_pred_all <- dplyr::bind_rows(med_pred_rows)
med_split_metrics <- dplyr::bind_rows(med_split_metrics_rows)
med_per_class <- dplyr::bind_rows(med_per_class_rows)
med_best <- dplyr::bind_rows(med_best_rows)

data.table::fwrite(
  med_pred_all,
  file.path(out_dir, "internal_med139_ridge_predictions_all.tsv"),
  sep = "\t"
)

data.table::fwrite(
  med_split_metrics,
  file.path(out_dir, "internal_med139_ridge_split_metrics.tsv"),
  sep = "\t"
)

data.table::fwrite(
  med_per_class,
  file.path(out_dir, "internal_med139_ridge_per_class_metrics.tsv"),
  sep = "\t"
)

data.table::fwrite(
  med_best,
  file.path(out_dir, "internal_med139_ridge_best_params.tsv"),
  sep = "\t"
)

# -------------------------------
# 5. Paired comparison with pFever
# -------------------------------

paired <- med_pred_all %>%
  dplyr::inner_join(
    pfever_majority %>%
      dplyr::select(split, sample_id, truth, pfever_prediction, n_level1_voted),
    by = c("split", "sample_id", "truth")
  ) %>%
  dplyr::mutate(
    med139_correct = med139_prediction == truth,
    pfever_correct = pfever_prediction == truth
  )

if (nrow(paired) != nrow(med_pred_all)) {
  warning("Paired table row count differs from Med prediction row count.")
}

data.table::fwrite(
  paired,
  file.path(out_dir, "internal_med139_vs_pfever_paired_predictions.tsv"),
  sep = "\t"
)

med_overall <- calc_metrics(
  truth = paired$truth,
  pred = paired$med139_prediction,
  labels = primary_labels
)$metrics %>%
  dplyr::mutate(
    model = "med139_ridge_imbalance_only_lambda_min",
    n_med_genes_used = length(use_genes),
    .before = 1
  )

pfever_overall <- calc_metrics(
  truth = paired$truth,
  pred = paired$pfever_prediction,
  labels = primary_labels
)$metrics %>%
  dplyr::mutate(
    model = "pFever_internal_majority",
    n_med_genes_used = NA_integer_,
    .before = 1
  )

overall_metrics <- dplyr::bind_rows(med_overall, pfever_overall)

data.table::fwrite(
  overall_metrics,
  file.path(out_dir, "internal_med139_vs_pfever_overall_metrics.tsv"),
  sep = "\t"
)

pfever_split_metrics <- paired %>%
  dplyr::group_by(split) %>%
  dplyr::group_modify(~{
    calc_metrics(
      truth = .x$truth,
      pred = .x$pfever_prediction,
      labels = primary_labels
    )$metrics
  }) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    model = "pFever_internal_majority",
    .before = 1
  )

split_metrics_all <- dplyr::bind_rows(
  med_split_metrics %>%
    dplyr::select(model, split, dplyr::everything()),
  pfever_split_metrics %>%
    dplyr::select(model, split, dplyr::everything())
)

data.table::fwrite(
  split_metrics_all,
  file.path(out_dir, "internal_med139_vs_pfever_split_metrics.tsv"),
  sep = "\t"
)

med_per_class_overall <- calc_metrics(
  truth = paired$truth,
  pred = paired$med139_prediction,
  labels = primary_labels
)$per_class %>%
  dplyr::mutate(
    model = "med139_ridge_imbalance_only_lambda_min",
    .before = 1
  )

pfever_per_class_overall <- calc_metrics(
  truth = paired$truth,
  pred = paired$pfever_prediction,
  labels = primary_labels
)$per_class %>%
  dplyr::mutate(
    model = "pFever_internal_majority",
    .before = 1
  )

per_class_overall <- dplyr::bind_rows(
  med_per_class_overall,
  pfever_per_class_overall
)

data.table::fwrite(
  per_class_overall,
  file.path(out_dir, "internal_med139_vs_pfever_per_class_metrics.tsv"),
  sep = "\t"
)

mc <- mcnemar_exact(
  a_correct = paired$med139_correct,
  b_correct = paired$pfever_correct
) %>%
  dplyr::mutate(
    comparison = "med139_ridge_imbalance_only_lambda_min vs pFever_internal_majority",
    .before = 1
  )

data.table::fwrite(
  mc,
  file.path(out_dir, "internal_med139_vs_pfever_mcnemar.tsv"),
  sep = "\t"
)

# -------------------------------
# 6. Print summary
# -------------------------------

cat("\n[DONE] Internal matched CV Med139 vs pFever\n")
cat("[OUT DIR] ", out_dir, "\n\n")

cat("[Overall metrics]\n")
print(overall_metrics)

cat("\n[McNemar]\n")
print(mc)

cat("\n[Split metrics]\n")
print(split_metrics_all)

cat("\n[Files]\n")
cat("  internal_med139_ridge_predictions_all.tsv\n")
cat("  internal_pfever_majority_from_level1.tsv\n")
cat("  internal_med139_vs_pfever_paired_predictions.tsv\n")
cat("  internal_med139_vs_pfever_overall_metrics.tsv\n")
cat("  internal_med139_vs_pfever_split_metrics.tsv\n")
cat("  internal_med139_vs_pfever_per_class_metrics.tsv\n")
cat("  internal_med139_vs_pfever_mcnemar.tsv\n")