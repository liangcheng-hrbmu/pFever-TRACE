#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import glob
import json
import argparse
import warnings
import zipfile
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import joblib
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import f1_score, balanced_accuracy_score, accuracy_score
from xgboost import XGBClassifier

try:
    import optuna
except Exception:
    optuna = None

SUFFIX = "_adataMerged_SPECTRAgenes.log1p.h5ad"


def ensure_dir(p):
    os.makedirs(p, exist_ok=True)


def celltype_from_h5ad(path: str) -> str:
    base = os.path.basename(path)
    if base.endswith(SUFFIX):
        return base[:-len(SUFFIX)]
    return os.path.splitext(base)[0]


def clean_label(x):
    if pd.isna(x):
        return np.nan
    x = str(x).strip()
    if x in {"", "nan", "None", "NA"}:
        return np.nan
    return x


def normalize_log1p_per_cell(X, target_sum=1e4):
    lib = X.sum(axis=1, keepdims=True)
    lib[lib == 0] = 1.0
    Xn = X / lib * float(target_sum)
    return np.log1p(Xn)


def as_dense(x):
    if sp.issparse(x):
        return x.toarray()
    return np.asarray(x)


def evaluate_metrics_str(y_true_str, y_pred_str):
    y_true = pd.Series(y_true_str, dtype="object").astype(str)
    y_pred = pd.Series(y_pred_str, dtype="object").astype(str)
    labels = sorted(set(y_true.tolist()) | set(y_pred.tolist()))
    return {
        "WF1": float(f1_score(y_true, y_pred, average="weighted", labels=labels)),
        "BAS": float(balanced_accuracy_score(y_true, y_pred)),
        "ACC": float(accuracy_score(y_true, y_pred)),
    }


def suggest_params(trial, n_class, use_gpu=False):
    return {
        "n_estimators": trial.suggest_int("n_estimators", 150, 1500),
        "max_depth": trial.suggest_int("max_depth", 3, 10),
        "learning_rate": trial.suggest_float("learning_rate", 0.01, 0.3, log=True),
        "subsample": trial.suggest_float("subsample", 0.6, 1.0),
        "colsample_bytree": trial.suggest_float("colsample_bytree", 0.4, 1.0),
        "min_child_weight": trial.suggest_float("min_child_weight", 1.0, 10.0),
        "reg_alpha": trial.suggest_float("reg_alpha", 1e-8, 10.0, log=True),
        "reg_lambda": trial.suggest_float("reg_lambda", 1e-8, 10.0, log=True),
        "objective": "multi:softmax",
        "num_class": int(n_class),
        "eval_metric": "mlogloss",
        "random_state": 1234,
        "tree_method": "hist",
        "device": "cuda" if use_gpu else "cpu",
        "verbosity": 0,
        "n_jobs": 5,
    }


def fit_best_xgb(X_train, y_train, X_val, y_val_str, train_label_le, n_trials=50, use_gpu=False):
    n_class = len(np.unique(y_train))
    if n_class < 2:
        raise ValueError("Need at least 2 classes in training data.")

    def _score_model(model):
        pred_val_code = model.predict(X_val)
        pred_val_str = train_label_le.inverse_transform(pred_val_code.astype(int))
        return f1_score(pd.Series(y_val_str).astype(str), pd.Series(pred_val_str).astype(str), average="weighted")

    if optuna is None or n_trials <= 0:
        params = {
            "n_estimators": 500,
            "max_depth": 6,
            "learning_rate": 0.05,
            "subsample": 0.8,
            "colsample_bytree": 0.8,
            "min_child_weight": 1.0,
            "reg_alpha": 1e-6,
            "reg_lambda": 1.0,
            "objective": "multi:softmax",
            "num_class": int(n_class),
            "eval_metric": "mlogloss",
            "random_state": 1234,
            "tree_method": "hist",
            "device": "cuda" if use_gpu else "cpu",
            "verbosity": 0,
            "n_jobs": 5,
        }
        model = XGBClassifier(**params)
        model.fit(X_train, y_train, verbose=False)
        score = _score_model(model)
        return model, params, score

    def objective(trial):
        params = suggest_params(trial, n_class=n_class, use_gpu=use_gpu)
        model = XGBClassifier(**params)
        model.fit(X_train, y_train, verbose=False)
        return _score_model(model)

    study = optuna.create_study(direction="maximize")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        study.optimize(objective, n_trials=n_trials, show_progress_bar=False)

    best_params = suggest_params(study.best_trial, n_class=n_class, use_gpu=use_gpu)
    model = XGBClassifier(**best_params)
    model.fit(X_train, y_train, verbose=False)
    return model, best_params, float(study.best_value)


def remap_indices(indices, remap):
    return np.array([remap[i] for i in indices if i in remap], dtype=int)


def write_pred_zip(df, zip_fp, inner_name=None):
    if inner_name is None:
        inner_name = os.path.basename(zip_fp).replace(".zip", ".tsv")
    tmp_tsv = zip_fp.replace(".zip", ".tsv")
    df.to_csv(tmp_tsv, sep="\t", index=False)
    with zipfile.ZipFile(zip_fp, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        zf.write(tmp_tsv, arcname=inner_name)
    return tmp_tsv


def main():
    ap = argparse.ArgumentParser(description="Fit unintegrated XGBoost models using same splits as corrected models (coverage-safe v4).")
    ap.add_argument("--main_h5ad", required=True)
    ap.add_argument("--corrected_input_dir", required=True)
    ap.add_argument("--split_dir", required=True)
    ap.add_argument("--out_dir", required=True)
    ap.add_argument("--disease_key", default="disease_Level_3")
    ap.add_argument("--counts_layer", default="counts")
    ap.add_argument("--library_size", type=float, default=1e4)
    ap.add_argument("--n_trials", type=int, default=50)
    ap.add_argument("--use_gpu", action="store_true")
    ap.add_argument("--only", nargs="*", default=None)
    args = ap.parse_args()

    ensure_dir(args.out_dir)
    pred_dir = os.path.join(args.out_dir, "predictions")
    ensure_dir(pred_dir)

    print("[READ]", args.main_h5ad)
    main_adata = sc.read_h5ad(args.main_h5ad)
    if args.counts_layer not in main_adata.layers:
        raise ValueError(f"counts layer not found: {args.counts_layer}")
    if args.disease_key not in main_adata.obs.columns:
        raise ValueError(f"main_h5ad missing obs column: {args.disease_key}")

    main_adata.obs[args.disease_key] = main_adata.obs[args.disease_key].map(clean_label)

    h5ads = sorted(glob.glob(os.path.join(args.corrected_input_dir, f"*{SUFFIX}")))
    if not h5ads:
        raise FileNotFoundError(f"No corrected per-celltype h5ad found in {args.corrected_input_dir}")

    only_set = set(args.only) if args.only else None

    for corr_fp in h5ads:
        celltype = celltype_from_h5ad(corr_fp)
        if only_set is not None and celltype not in only_set:
            continue

        split_fp = os.path.join(args.split_dir, f"{celltype}.split_indices.npz")
        corr_model_fp = os.path.join(args.split_dir, f"{celltype}.xgb_model.joblib")
        if not os.path.exists(split_fp):
            print(f"[SKIP] {celltype}: missing split")
            continue
        if not os.path.exists(corr_model_fp):
            print(f"[SKIP] {celltype}: missing corrected model bundle")
            continue

        print(f"\n[RUN] {celltype}")
        corr = sc.read_h5ad(corr_fp)
        corr_bundle = joblib.load(corr_model_fp)
        corr_le = corr_bundle["label_encoder"]
        trained_classes = set(map(lambda x: str(x).strip(), corr_le.classes_.tolist()))

        cells = corr.obs_names.astype(str).tolist()
        genes = corr.var_names.astype(str).tolist()
        keep_cells = [x for x in cells if x in main_adata.obs_names]
        if len(keep_cells) == 0:
            print(f"[SKIP] {celltype}: no overlapping cells in main_h5ad")
            continue

        sub = main_adata[keep_cells, genes].copy()
        y_all = pd.Series(sub.obs[args.disease_key]).map(clean_label).astype("object").values

        class_mask = np.array([y in trained_classes for y in y_all], dtype=bool)
        if class_mask.sum() == 0:
            print(f"[SKIP] {celltype}: no cells remain after corrected-class filtering")
            continue

        orig_idx = np.where(class_mask)[0]
        remap = {old: new for new, old in enumerate(orig_idx)}

        X_counts = as_dense(sub.layers[args.counts_layer])[class_mask].astype(np.float64, copy=False)
        Xf = normalize_log1p_per_cell(X_counts, target_sum=args.library_size)
        yf = pd.Series(y_all[class_mask]).map(clean_label).astype(str).values
        kept_obs_names = np.asarray(sub.obs_names)[class_mask]

        split_obj = np.load(split_fp)
        train_idx = np.asarray(split_obj["train_idx"], dtype=int)
        val_idx = np.asarray(split_obj["val_idx"], dtype=int)
        test_idx = np.asarray(split_obj["test_idx"], dtype=int)

        train_idx_f = remap_indices(train_idx, remap)
        val_idx_f = remap_indices(val_idx, remap)
        test_idx_f = remap_indices(test_idx, remap)

        if len(train_idx_f) == 0 or len(val_idx_f) == 0 or len(test_idx_f) == 0:
            print(f"[SKIP] {celltype}: empty train/val/test after corrected-class filtering")
            continue

        train_label_le = LabelEncoder()
        y_train_str = yf[train_idx_f]
        train_label_le.fit(y_train_str)

        X_train = Xf[train_idx_f]
        y_train = train_label_le.transform(y_train_str)
        X_val = Xf[val_idx_f]
        y_val_str = yf[val_idx_f]
        X_test = Xf[test_idx_f]
        y_test_str = yf[test_idx_f]

        unseen_val = int(np.sum(~np.isin(y_val_str, train_label_le.classes_)))
        unseen_test = int(np.sum(~np.isin(y_test_str, train_label_le.classes_)))
        if unseen_val > 0 or unseen_test > 0:
            print(f"[WARN] {celltype}: keep unseen classes in eval/test (val={unseen_val}, test={unseen_test})")

        model, best_params, best_val_wf1 = fit_best_xgb(
            X_train, y_train, X_val, y_val_str,
            train_label_le=train_label_le,
            n_trials=args.n_trials,
            use_gpu=args.use_gpu,
        )

        pred_test_code = model.predict(X_test).astype(int)
        pred_test_str = train_label_le.inverse_transform(pred_test_code)
        metrics = evaluate_metrics_str(y_test_str, pred_test_str)
        metrics.update({
            "best_val_WF1": best_val_wf1,
            "celltype": celltype,
            "n_train": int(len(train_idx_f)),
            "n_val": int(len(val_idx_f)),
            "n_test": int(len(test_idx_f)),
            "n_classes_train": int(len(train_label_le.classes_)),
            "n_unseen_val_classes": unseen_val,
            "n_unseen_test_classes": unseen_test,
            "input_type": "unintegrated_lognorm_from_main_h5ad",
            "used_gpu": bool(args.use_gpu),
        })

        out_bundle = {
            "model": model,
            "label_encoder": train_label_le,
            "params": best_params,
            "corrected_classes": sorted(trained_classes),
        }
        model_fp = os.path.join(args.out_dir, f"{celltype}.xgb_model.joblib")
        metrics_fp = os.path.join(args.out_dir, f"{celltype}.metrics.json")
        split_out_fp = os.path.join(args.out_dir, f"{celltype}.split_indices.npz")
        pred_zip_fp = os.path.join(pred_dir, f"run1_{celltype}_pred_test.zip")
        pred_tsv_fp = os.path.join(pred_dir, f"run1_{celltype}_pred_test.tsv")
        test_npy_fp = os.path.join(pred_dir, f"run1_{celltype}_test_idxs.npy")

        joblib.dump(out_bundle, model_fp)
        with open(metrics_fp, "w", encoding="utf-8") as f:
            json.dump(metrics, f, ensure_ascii=False, indent=2)
        np.savez(split_out_fp, train_idx=train_idx_f, val_idx=val_idx_f, test_idx=test_idx_f)
        np.save(test_npy_fp, test_idx_f)

        pred_df = pd.DataFrame({
            "celltype": celltype,
            "cell_id": kept_obs_names[test_idx_f].astype(str),
            "true_label": pd.Series(y_test_str).astype(str).str.strip(),
            "pred_label": pd.Series(pred_test_str).astype(str).str.strip(),
            "pred_code": pred_test_code,
        })
        pred_df.to_csv(pred_tsv_fp, sep="\t", index=False)
        write_pred_zip(pred_df, pred_zip_fp)

        print(f"[WRITE] {model_fp}")
        print(f"[WRITE] {metrics_fp}")
        print(f"[WRITE] {pred_zip_fp}")

    print("\n[DONE]")
    print("out_dir:", args.out_dir)


if __name__ == "__main__":
    main()
