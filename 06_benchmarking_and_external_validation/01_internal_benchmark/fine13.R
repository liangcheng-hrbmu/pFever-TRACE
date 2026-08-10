suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(glmnet)
})

# =========================================================
# Internal matched CV benchmark, Fine 13
# Med139 ridge multinomial baseline versus pFever internal CV
#
# Design:
#   Same outer held out samples as pFever Fine 13 internal split
#   Train Med139 ridge model on remaining E MTAB 11671 samples
#   Select lambda by stratified inner CV within training samples only
#   Predict the same outer held out samples
#   Compare paired predictions with pFever Fine 13 majority vote
#
# Main comparison:
#   Med139 imbalance only lambda.min
#   vs
#   pFever Fine 13 internal majority vote
# =========================================================

# -------------------------------
# 0. Paths
# -------------------------------

BASE <- "data/analysis/09_bulk_deconvolution"
EXT_ROOT <- file.path(BASE, "12_external_validation_multi_cohort_assets_v1")

train_expr_fp <- file.path(BASE, "E-MTAB-11671_mixture_counts_filtered.tsv")

signature_genes_fp <- file.path(
  EXT_ROOT,
  "results/E-MTAB-12793/comparison_ours_vs_145_external/wholeblood_signature_genes_used.tsv"
)

out_dir <- file.path(
  BASE,
  "08_bulk_classifier/04_internal_matched_cv_med139_vs_pfever_fine13"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# 如果自动找错，就把这里改成真实 Fine 13 的 query_level1_predictions_all.tsv
pfever_l1_pred_fp_manual <- NA_character_

candidate_pfever_files <- c(
  file.path(BASE, "08_bulk_classifier/03_patient_classifier_fine13_bulkZ_L1x11_final_withLinearSVC/query_level1_predictions_all.tsv"),
  file.path(BASE, "08_bulk_classifier/03_patient_classifier_fine_bulkZ_L1x11_final_withLinearSVC/query_level1_predictions_all.tsv"),
  file.path(BASE, "08_bulk_classifier/03_patient_classifier_fine13_bulkZ_L1x11/query_level1_predictions_all.tsv"),
  file.path(BASE, "08_bulk_classifier/03_patient_classifier_fine_bulkZ_L1x11/query_level1_predictions_all.tsv")
)

if (!is.na(pfever_l1_pred_fp_manual) && file.exists(pfever_l1_pred_fp_manual)) {
  pfever_l1_pred_fp <- pfever_l1_pred_fp_manual
} else {
  hit <- candidate_pfever_files[file.exists(candidate_pfever_files)]

  if (length(hit) == 0) {
    all_query <- list.files(
      file.path(BASE, "08_bulk_classifier"),
      pattern = "query_level1_predictions_all.tsv$",
      recursive = TRUE,
      full.names = TRUE
    )

    hit <- all_query[
      grepl("fine|13", all_query, ignore.case = TRUE) &
        !grepl("broad5|broad", all_query, ignore.case = TRUE)
    ]

    cat("\n[DIAG] candidate Fine 13 query files found:\n")
    print(hit)
  }

  if (length(hit) == 0) {
    stop(
      "Cannot find Fine 13 pFever query_level1_predictions_all.tsv. Tried:\n",
      paste(candidate_pfever_files, collapse = "\n")
    )
  }

  pfever_l1_pred_fp <- hit[1]
}

if (!file.exists(train_expr_fp)) {
  stop("Missing train expression file: ", train_expr_fp)
}
if (!file.exists(signature_genes_fp)) {
  stop("Missing signature gene file: ", signature_genes_fp)
}
if (!file.exists(pfever_l1_pred_fp)) {
  stop("Missing pFever Fine 13 L1 file: ", pfever_l1_pred_fp)
}

# -------------------------------
# 1. Helper functions
# -------------------------------

norm_fine13 <- function(x) {
  x0 <- trimws(as.character(x))
  xl <- tolower(x0)
  out <- rep(NA_character_, length(x0))

  out[xl %in% c(
    "iav", "flu", "influenza", "influenza a",
    "influenza a virus", "influenza a virus (iav)"
  )] <- "IAV"

  out[xl %in% c(
    "rsv", "respiratory syncytial virus",
    "respiratory syncytial",
    "respiratory syncytial virus infection"
  )] <- "RSV"

  out[xl %in% c(
    "adenovirus", "adeno", "adv", "adenovirus infection"
  )] <- "Adenovirus"

  out[xl %in% c(
    "rhinovirus", "rhino", "rhinovirus infection"
  )] <- "Rhinovirus"

  out[xl %in% c(
    "meningococcus", "meningococcal", "meningococcal infection",
    "neisseria meningitidis"
  )] <- "Meningococcus"

  out[xl %in% c(
    "pneumococcus", "pneumococcal", "pneumococcal infection",
    "streptococcus pneumoniae"
  )] <- "Pneumococcus"

  out[xl %in% c(
    "s. aureus", "staph", "staphylococcus aureus",
    "staphylococcus aureus infection", "aureus"
  )] <- "S. aureus"

  out[xl %in% c(
    "gas", "group a streptococcus",
    "group a streptococcal",
    "group a streptococcal infection"
  )] <- "GAS"

  out[xl %in% c(
    "e. coli", "e coli", "ecoli", "escherichia coli",
    "escherichia coli infection"
  )] <- "E. coli"

  out[xl %in% c(
    "tuberculosis", "tb", "mycobacterium tuberculosis",
    "tuberculosis infection"
  )] <- "Tuberculosis"

  out[xl %in% c(
    "malaria", "plasmodium", "malaria infection"
  )] <- "Malaria"

  out[xl %in% c(
    "kd", "kawasaki", "kawasaki disease",
    "kawasaki disease (kd)"
  )] <- "KD"

  out[xl %in% c(
    "jia", "juvenile idiopathic arthritis",
    "juvenile idopathic arthritis",
    "juvenile idiopathic arthritis (jia)"
  )] <- "JIA"

  out[xl %in% c(
    "tie", "not assigned", "unassigned", "unknown"
  )] <- "Not assigned"

  out[is.na(out) & grepl("influenza|\\biav\\b|\\bflu\\b", xl)] <- "IAV"
  out[is.na(out) & grepl("respiratory syncytial|\\brsv\\b", xl)] <- "RSV"
  out[is.na(out) & grepl("adenovirus|\\badeno\\b", xl)] <- "Adenovirus"
  out[is.na(out) & grepl("rhinovirus|\\brhino\\b", xl)] <- "Rhinovirus"

  out[is.na(out) & grepl("meningococcus|meningococcal|meningitidis", xl)] <- "Meningococcus"
  out[is.na(out) & grepl("pneumococcus|pneumococcal|pneumoniae", xl)] <- "Pneumococcus"
  out[is.na(out) & grepl("staph|aureus", xl)] <- "S. aureus"
  out[is.na(out) & grepl("group a streptococcus|group a streptococcal|\\bgas\\b", xl)] <- "GAS"
  out[is.na(out) & grepl("escherichia|e\\.coli|e coli|ecoli", xl)] <- "E. coli"
  out[is.na(out) & grepl("tuberculosis|\\btb\\b", xl)] <- "Tuberculosis"

  out[is.na(out) & grepl("malaria|plasmodium", xl)] <- "Malaria"
  out[is.na(out) & grepl("kawasaki|\\bkd\\b", xl)] <- "KD"
  out[is.na(out) & grepl("juvenile idiopathic|juvenile idopathic|\\bjia\\b", xl)] <- "JIA"

  out
}

detect_col <- function(df, candidates, required = TRUE, label = "column") {
  cols <- colnames(df)
  lower_map <- stats::setNames(cols, tolower(cols))

  for (cc in candidates) {
    if (cc %in% cols) return(cc)
    if (tolower(cc) %in% names(lower_map)) return(lower_map[[tolower(cc)]])
  }

  if (required) {
    stop(
      "Cannot find ", label, ". Tried: ",
      paste(candidates, collapse = ", "),
      "\nCurrent columns: ",
      paste(cols, collapse = " | ")
    )
  }

  NA_character_
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
    message("[INFO] auto log1p applied.")
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
  if (is.na(gene_col)) gene_col <- colnames(sig)[1]

  genes <- unique(trimws(as.character(sig[[gene_col]])))
  genes <- genes[!is.na(genes) & genes != ""]

  genes
}

majority_vote <- function(x) {
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) return("Not assigned")

  tb <- sort(table(x), decreasing = TRUE)

  if (length(tb) >= 2 && tb[1] == tb[2]) {
    return("Not assigned")
  }

  names(tb)[1]
}

calc_metrics <- function(truth, pred, labels = NULL) {
  truth <- as.character(truth)
  pred <- as.character(pred)

  if (is.null(labels)) labels <- sort(unique(truth))
  labels <- labels[labels %in% unique(truth)]

  acc <- mean(truth == pred)

  per_class <- dplyr::bind_rows(lapply(labels, function(cls) {
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
  }))

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

make_stratified_foldid <- function(y, nfolds, seed = 1) {
  y <- as.character(y)
  foldid <- rep(NA_integer_, length(y))

  set.seed(seed)

  for (cls in unique(y)) {
    idx <- which(y == cls)
    idx <- sample(idx, length(idx))

    cls_fold <- rep(seq_len(nfolds), length.out = length(idx))
    cls_fold <- sample(cls_fold, length(cls_fold))

    foldid[idx] <- cls_fold
  }

  if (any(is.na(foldid))) {
    stop("Failed to create stratified foldid.")
  }

  foldid
}

# -------------------------------
# 2. Read pFever Fine 13 internal predictions
# -------------------------------

message("[READ] pFever Fine 13 Level1 internal predictions: ", pfever_l1_pred_fp)

l1 <- data.table::fread(pfever_l1_pred_fp, data.table = FALSE)

cat("\n[DIAG] pFever file used:\n")
cat(pfever_l1_pred_fp, "\n")

cat("\n[DIAG] columns in pFever file:\n")
print(colnames(l1))

cat("\n[DIAG] first rows:\n")
print(utils::head(l1, 5))

sample_col <- detect_col(
  l1,
  c("sample_id", "sample", "patient_id"),
  label = "sample column"
)

level1_col <- detect_col(
  l1,
  c("Level1", "level1", "L1", "celltype", "cell_type"),
  label = "Level1 column"
)

split_col <- detect_col(
  l1,
  c("split", "fold", "cv_split"),
  label = "split column"
)

truth_col <- detect_col(
  l1,
  c(
    "disease_clean",
    "disease_fine",
    "disease_fine13",
    "fine13",
    "fine_label",
    "disease_label",
    "disease_original",
    "truth_fine",
    "truth_fine13",
    "true_disease_fine",
    "true_disease",
    "truth",
    "label",
    "disease"
  ),
  required = TRUE,
  label = "Fine 13 truth column"
)

pred_col <- detect_col(
  l1,
  c("pred_disease", "pred_label", "prediction", "pred", "y_pred"),
  label = "prediction column"
)

cat("\n[DIAG] detected columns:\n")
cat("sample_col:", sample_col, "\n")
cat("level1_col:", level1_col, "\n")
cat("truth_col:", truth_col, "\n")
cat("pred_col:", pred_col, "\n")
cat("split_col:", split_col, "\n")

l1 <- l1 %>%
  dplyr::mutate(
    sample_id = trimws(as.character(.data[[sample_col]])),
    Level1 = trimws(as.character(.data[[level1_col]])),
    truth_raw = trimws(as.character(.data[[truth_col]])),
    pred_raw = trimws(as.character(.data[[pred_col]])),
    truth = norm_fine13(truth_raw),
    pred_level1 = norm_fine13(pred_raw),
    split = as.integer(.data[[split_col]])
  )

primary_labels <- c(
  "IAV", "RSV", "Adenovirus", "Rhinovirus",
  "Meningococcus", "Pneumococcus", "S. aureus", "GAS", "E. coli", "Tuberculosis",
  "Malaria", "KD", "JIA"
)

message("\n[DIAG] truth raw label counts:")
print(sort(table(l1$truth_raw, useNA = "ifany"), decreasing = TRUE))

message("\n[DIAG] truth label counts after Fine 13 mapping:")
print(sort(table(l1$truth, useNA = "ifany"), decreasing = TRUE))

message("\n[DIAG] pred raw label counts:")
print(sort(table(l1$pred_raw, useNA = "ifany"), decreasing = TRUE))

message("\n[DIAG] pred label counts after Fine 13 mapping:")
print(sort(table(l1$pred_level1, useNA = "ifany"), decreasing = TRUE))

l1 <- l1 %>%
  dplyr::filter(truth %in% primary_labels)

if (nrow(l1) == 0) {
  stop("No Fine 13 truth labels detected. Please check truth column in: ", pfever_l1_pred_fp)
}

pfever_majority <- l1 %>%
  dplyr::group_by(split, sample_id, truth) %>%
  dplyr::summarise(
    pfever_prediction = majority_vote(pred_level1),
    n_level1_voted = dplyr::n_distinct(Level1),
    .groups = "drop"
  )

message("[INFO] pFever Fine 13 internal majority samples: ", nrow(pfever_majority))
message("[INFO] pFever splits: ", paste(sort(unique(pfever_majority$split)), collapse = ", "))

data.table::fwrite(
  pfever_majority,
  file.path(out_dir, "internal_pfever_majority_from_level1.tsv"),
  sep = "\t"
)

sample_label <- pfever_majority %>%
  dplyr::distinct(sample_id, truth)

if (anyDuplicated(sample_label$sample_id) > 0) {
  stop("Duplicated sample_id with inconsistent Fine 13 labels detected.")
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
  stop("Too few Med signature genes available.")
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
# 4. Internal matched CV for Med139 ridge Fine 13
# -------------------------------

set.seed(1)

med_pred_rows <- list()
med_split_metrics_rows <- list()
med_per_class_rows <- list()
med_best_rows <- list()

splits <- sort(unique(pfever_majority$split))

for (sp in splits) {
  message("\n[RUN] Fine 13 split ", sp)

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

  if (any(class_counts < 3)) {
    stop(
      "At least one class has fewer than 3 training samples in split ",
      sp,
      ". This is too small for stable inner CV. Counts: ",
      paste(names(class_counts), class_counts, sep = "=", collapse = "; ")
    )
  }

  nfolds_inner <- min(5, min(class_counts))
  nfolds_inner <- max(2, nfolds_inner)

  foldid_inner <- make_stratified_foldid(
    y = y_train,
    nfolds = nfolds_inner,
    seed = 1000 + as.integer(sp)
  )

  message("[INFO] split ", sp, " train n=", length(y_train), ", test n=", length(y_test))
  message("[INFO] split ", sp, " train class counts: ", paste(names(class_counts), class_counts, sep = "=", collapse = "; "))
  message("[INFO] split ", sp, " inner nfolds=", nfolds_inner)
  message("[INFO] split ", sp, " inner fold class table:")
  print(table(foldid_inner, y_train))

  cvfit <- glmnet::cv.glmnet(
    x = X_train,
    y = y_train,
    family = "multinomial",
    alpha = 0,
    weights = w_train,
    type.measure = "class",
    foldid = foldid_inner,
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

cat("\n[DONE] Fine 13 internal matched CV Med139 vs pFever\n")
cat("[pFever L1 file] ", pfever_l1_pred_fp, "\n")
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