#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import json
import glob
import argparse
import warnings
from collections import Counter
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
import anndata as ad
import optuna

from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.svm import LinearSVC, SVC
from sklearn.neighbors import KNeighborsClassifier
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import (
    f1_score,
    balanced_accuracy_score,
    matthews_corrcoef,
)


def ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def clean_str(x):
    if pd.isna(x):
        return np.nan
    x = str(x).strip()
    if x in {"", "nan", "None", "NA"}:
        return np.nan
    return x


def parse_split_id(path: str) -> int:
    m = re.search(r"SPLIT_(\d+)", os.path.basename(path))
    if m is None:
        raise ValueError(f"Cannot parse split id from filename: {path}")
    return int(m.group(1))


def load_pb(path: str) -> ad.AnnData:
    adata = ad.read_h5ad(path)
    adata.obs_names_make_unique()
    for col in adata.obs.columns:
        if str(adata.obs[col].dtype) == "category" or adata.obs[col].dtype == object:
            adata.obs[col] = adata.obs[col].map(clean_str)
    return adata


def make_clf_from_trial(trial: optuna.Trial, classifier_space: List[str]):
    clf_type = trial.suggest_categorical("classifier_type", classifier_space)

    if clf_type == "LinearSVC":
        C = trial.suggest_float("C", 1e-3, 1e2, log=True)
        model = LinearSVC(C=C, dual="auto", max_iter=5000)
    elif clf_type == "SVC":
        C = trial.suggest_float("C", 1e-3, 1e2, log=True)
        kernel = trial.suggest_categorical("kernel", ["linear", "rbf"])
        if kernel == "linear":
            model = SVC(C=C, kernel="linear")
        else:
            gamma = trial.suggest_float("gamma", 1e-4, 1e1, log=True)
            model = SVC(C=C, kernel="rbf", gamma=gamma)
    elif clf_type == "KNN":
        n_neighbors = trial.suggest_int("n_neighbors", 3, 15)
        weights = trial.suggest_categorical("weights", ["uniform", "distance"])
        p = trial.suggest_categorical("p", [1, 2])
        model = KNeighborsClassifier(n_neighbors=n_neighbors, weights=weights, p=p)
    else:
        raise ValueError(f"Unsupported classifier type: {clf_type}")

    return model


def build_model_from_params(params: Dict):
    params = dict(params)
    clf_type = params.pop("classifier_type")

    if clf_type == "LinearSVC":
        return LinearSVC(C=float(params["C"]), dual="auto", max_iter=5000)
    elif clf_type == "SVC":
        kernel = params["kernel"]
        if kernel == "linear":
            return SVC(C=float(params["C"]), kernel="linear")
        return SVC(C=float(params["C"]), kernel="rbf", gamma=float(params["gamma"]))
    elif clf_type == "KNN":
        return KNeighborsClassifier(
            n_neighbors=int(params["n_neighbors"]),
            weights=params["weights"],
            p=int(params["p"]),
        )
    else:
        raise ValueError(f"Unsupported classifier_type: {clf_type}")


def optimize_classifier(
    X: np.ndarray,
    y_str: np.ndarray,
    n_trials: int,
    cv_splits: int,
    random_state: int,
    classifier_space: List[str],
):
    y_str = np.asarray(y_str).astype(str)
    X = np.asarray(X, dtype=np.float32)

    le = LabelEncoder()
    y = le.fit_transform(y_str)

    skf = StratifiedKFold(n_splits=cv_splits, shuffle=True, random_state=random_state)

    def objective(trial: optuna.Trial):
        model = make_clf_from_trial(trial, classifier_space)
        try:
            scores = cross_val_score(
                model,
                X,
                y,
                cv=skf,
                scoring="f1_weighted",
                n_jobs=1,
                error_score=np.nan,
            )
            score = float(np.nanmean(scores))
            if np.isnan(score):
                return 0.0
            return score
        except Exception:
            return 0.0

    study = optuna.create_study(direction="maximize")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        study.optimize(objective, n_trials=n_trials, show_progress_bar=False)

    return study.best_trial.params, float(study.best_value), le


def majority_vote(preds: List[str]) -> str:
    preds = [str(x) for x in preds if pd.notna(x)]
    if len(preds) == 0:
        return "NA"
    cnt = Counter(preds)
    top = cnt.most_common()
    if len(top) == 1:
        return top[0][0]
    if top[0][1] == top[1][1]:
        return "TIE"
    return top[0][0]


def compute_metrics(y_true: List[str], y_pred: List[str]) -> Dict[str, float]:
    y_true = pd.Series(y_true).astype(str)
    y_pred = pd.Series(y_pred).astype(str)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        wf1 = float(f1_score(y_true, y_pred, average="weighted"))
        bas = float(balanced_accuracy_score(y_true, y_pred))
        mcc = float(matthews_corrcoef(y_true, y_pred))
    return {"WF1": wf1, "BAS": bas, "MCC": mcc}


def main():
    ap = argparse.ArgumentParser(
        description="Train per-Level1 patient classifiers on pseudobulk embeddings and aggregate by majority vote (label-encoded, warning-safe version)."
    )
    ap.add_argument("--pseudobulk_dir", required=True)
    ap.add_argument("--out_dir", required=True)
    ap.add_argument("--level1_key", default="Level1")
    ap.add_argument("--sample_key", default="sample_id")
    ap.add_argument("--disease_key", default="disease_Level_3")
    ap.add_argument("--benchmark_key", default="benchmark_status")
    ap.add_argument("--primary_benchmark_label", default="primary_benchmark")
    ap.add_argument("--n_trials", type=int, default=30)
    ap.add_argument("--cv_splits", type=int, default=5)
    ap.add_argument("--random_state", type=int, default=1234)
    ap.add_argument("--classifier_space", nargs="*", default=["LinearSVC", "SVC", "KNN"])
    ap.add_argument("--only_splits", nargs="*", type=int, default=None)
    args = ap.parse_args()

    ensure_dir(args.out_dir)
    ensure_dir(os.path.join(args.out_dir, "best_params"))
    ensure_dir(os.path.join(args.out_dir, "predictions"))
    ensure_dir(os.path.join(args.out_dir, "metrics"))

    ref_files = sorted(glob.glob(os.path.join(args.pseudobulk_dir, "PSEUDOBULKs", "*reference_SPLIT_*.h5ad")))
    qry_files = sorted(glob.glob(os.path.join(args.pseudobulk_dir, "PSEUDOBULKs", "*query_SPLIT_*.h5ad")))

    if len(ref_files) == 0 or len(qry_files) == 0:
        raise FileNotFoundError("No pseudobulk reference/query files found under PSEUDOBULKs/")

    ref_map = {parse_split_id(fp): fp for fp in ref_files}
    qry_map = {parse_split_id(fp): fp for fp in qry_files}
    common_splits = sorted(set(ref_map).intersection(qry_map))

    if args.only_splits is not None:
        common_splits = [s for s in common_splits if s in set(args.only_splits)]
    if len(common_splits) == 0:
        raise RuntimeError("No common splits to process.")

    all_level1_metric_rows = []
    all_majority_rows = []
    split_summary_rows = []

    for split_id in common_splits:
        print(f"\n[RUN] split={split_id}")
        ref = load_pb(ref_map[split_id])
        qry = load_pb(qry_map[split_id])

        ref.obs[args.level1_key] = ref.obs[args.level1_key].map(clean_str)
        qry.obs[args.level1_key] = qry.obs[args.level1_key].map(clean_str)
        ref.obs[args.disease_key] = ref.obs[args.disease_key].map(clean_str)
        qry.obs[args.disease_key] = qry.obs[args.disease_key].map(clean_str)

        level1_list = sorted(set(ref.obs[args.level1_key].dropna().tolist()).intersection(
            set(qry.obs[args.level1_key].dropna().tolist())
        ))

        pred_rows = []
        best_params_rows = []
        level1_metric_rows = []

        for ct in level1_list:
            ref_ct = ref[ref.obs[args.level1_key] == ct].copy()
            qry_ct = qry[qry.obs[args.level1_key] == ct].copy()

            ref_tune = ref_ct[ref_ct.obs[args.benchmark_key] == args.primary_benchmark_label].copy()
            if ref_tune.n_obs < args.cv_splits:
                print(f"[SKIP] split={split_id} {ct}: too few reference rows for CV ({ref_tune.n_obs})")
                continue

            y_tune = ref_tune.obs[args.disease_key].astype(str).values
            class_counts = pd.Series(y_tune).value_counts()
            cv_splits = int(min(args.cv_splits, class_counts.min()))
            cv_splits = max(2, cv_splits)

            X_tune = np.asarray(ref_tune.X, dtype=np.float32)
            if np.isnan(X_tune).any():
                raise RuntimeError(f"NaN detected in tuning matrix for split={split_id}, Level1={ct}")

            best_params, best_cv_wf1, le = optimize_classifier(
                X=X_tune,
                y_str=y_tune,
                n_trials=args.n_trials,
                cv_splits=cv_splits,
                random_state=args.random_state,
                classifier_space=args.classifier_space,
            )

            model = build_model_from_params(best_params)

            X_ref_final = np.asarray(ref_ct.X, dtype=np.float32)
            y_ref_final_str = ref_ct.obs[args.disease_key].astype(str).values
            keep_mask = np.isin(y_ref_final_str, le.classes_)
            X_ref_final = X_ref_final[keep_mask]
            y_ref_final = le.transform(y_ref_final_str[keep_mask])

            if np.isnan(X_ref_final).any():
                raise RuntimeError(f"NaN detected in final reference matrix for split={split_id}, Level1={ct}")

            model.fit(X_ref_final, y_ref_final)

            X_q = np.asarray(qry_ct.X, dtype=np.float32)
            if np.isnan(X_q).any():
                raise RuntimeError(f"NaN detected in query matrix for split={split_id}, Level1={ct}")

            y_pred_code = model.predict(X_q)
            y_pred = le.inverse_transform(np.asarray(y_pred_code, dtype=int))

            pred_df = qry_ct.obs[[args.sample_key, args.level1_key, args.disease_key, args.benchmark_key]].copy()
            pred_df["pred_disease"] = y_pred
            pred_df["split"] = split_id
            pred_df["classifier_ct"] = ct
            pred_rows.append(pred_df)

            eval_mask = pred_df[args.benchmark_key].astype(str) == args.primary_benchmark_label
            if eval_mask.sum() > 0:
                met = compute_metrics(
                    y_true=pred_df.loc[eval_mask, args.disease_key].astype(str).tolist(),
                    y_pred=pred_df.loc[eval_mask, "pred_disease"].astype(str).tolist(),
                )
            else:
                met = {"WF1": np.nan, "BAS": np.nan, "MCC": np.nan}

            level1_metric_rows.append({
                "split": split_id,
                "Level1": ct,
                "n_ref_rows": ref_ct.n_obs,
                "n_ref_rows_tuned": ref_tune.n_obs,
                "n_qry_rows": qry_ct.n_obs,
                "best_cv_WF1": best_cv_wf1,
                "classifier_type": best_params.get("classifier_type"),
                "WF1": met["WF1"],
                "BAS": met["BAS"],
                "MCC": met["MCC"],
            })

            row = {"split": split_id, "Level1": ct, "best_cv_WF1": best_cv_wf1}
            row.update(best_params)
            best_params_rows.append(row)

        if len(pred_rows) == 0:
            raise RuntimeError(f"split={split_id}: no Level1 classifier produced predictions.")

        pred_all = pd.concat(pred_rows, axis=0, ignore_index=True)
        level1_metrics_df = pd.DataFrame(level1_metric_rows)
        best_params_df = pd.DataFrame(best_params_rows)

        maj_rows = []
        for sample_id, sub in pred_all.groupby(args.sample_key, dropna=False, observed=False):
            true_labels = sub[args.disease_key].astype(str).unique().tolist()
            bench_labels = sub[args.benchmark_key].astype(str).unique().tolist()
            pred_vote = majority_vote(sub["pred_disease"].astype(str).tolist())
            maj_rows.append({
                args.sample_key: sample_id,
                "split": split_id,
                "true_disease": true_labels[0] if len(true_labels) > 0 else np.nan,
                args.benchmark_key: bench_labels[0] if len(bench_labels) > 0 else np.nan,
                "pred_disease_majority": pred_vote,
                "n_level1_votes": sub[args.level1_key].nunique(),
            })
        maj_df = pd.DataFrame(maj_rows)

        eval_maj = maj_df[maj_df[args.benchmark_key].astype(str) == args.primary_benchmark_label].copy()
        split_metrics = compute_metrics(
            y_true=eval_maj["true_disease"].astype(str).tolist(),
            y_pred=eval_maj["pred_disease_majority"].astype(str).tolist(),
        )

        split_summary_rows.append({
            "split": split_id,
            "n_eval_samples_majority": eval_maj.shape[0],
            "WF1_majority": split_metrics["WF1"],
            "BAS_majority": split_metrics["BAS"],
            "MCC_majority": split_metrics["MCC"],
            "mean_level1_WF1": float(level1_metrics_df["WF1"].mean(skipna=True)),
            "mean_level1_BAS": float(level1_metrics_df["BAS"].mean(skipna=True)),
            "mean_level1_MCC": float(level1_metrics_df["MCC"].mean(skipna=True)),
        })

        pred_fp = os.path.join(args.out_dir, "predictions", f"query_level1_predictions_SPLIT_{split_id}.tsv")
        maj_fp = os.path.join(args.out_dir, "predictions", f"majority_vote_predictions_SPLIT_{split_id}.tsv")
        level1_met_fp = os.path.join(args.out_dir, "metrics", f"level1_metrics_SPLIT_{split_id}.tsv")
        best_fp = os.path.join(args.out_dir, "best_params", f"best_params_SPLIT_{split_id}.tsv")
        split_json_fp = os.path.join(args.out_dir, "metrics", f"majority_vote_metrics_SPLIT_{split_id}.json")

        pred_all.to_csv(pred_fp, sep="\t", index=False)
        maj_df.to_csv(maj_fp, sep="\t", index=False)
        level1_metrics_df.to_csv(level1_met_fp, sep="\t", index=False)
        best_params_df.to_csv(best_fp, sep="\t", index=False)
        with open(split_json_fp, "w", encoding="utf-8") as f:
            json.dump(split_summary_rows[-1], f, ensure_ascii=False, indent=2)

        all_level1_metric_rows.append(level1_metrics_df)
        all_majority_rows.append(maj_df)

        print(f"[OK] split={split_id} majority WF1={split_metrics['WF1']:.4f} BAS={split_metrics['BAS']:.4f}")

    all_level1_metrics = pd.concat(all_level1_metric_rows, axis=0, ignore_index=True)
    all_majority = pd.concat(all_majority_rows, axis=0, ignore_index=True)
    split_summary_df = pd.DataFrame(split_summary_rows)

    split_summary_fp = os.path.join(args.out_dir, "scenario1_split_summary.tsv")
    all_level1_fp = os.path.join(args.out_dir, "scenario1_level1_metrics_all.tsv")
    all_majority_fp = os.path.join(args.out_dir, "scenario1_majority_predictions_all.tsv")

    split_summary_df.to_csv(split_summary_fp, sep="\t", index=False)
    all_level1_metrics.to_csv(all_level1_fp, sep="\t", index=False)
    all_majority.to_csv(all_majority_fp, sep="\t", index=False)

    overall = {
        "n_splits": len(common_splits),
        "WF1_majority_mean": float(split_summary_df["WF1_majority"].mean()),
        "WF1_majority_std": float(split_summary_df["WF1_majority"].std(ddof=1)) if len(split_summary_df) > 1 else 0.0,
        "BAS_majority_mean": float(split_summary_df["BAS_majority"].mean()),
        "BAS_majority_std": float(split_summary_df["BAS_majority"].std(ddof=1)) if len(split_summary_df) > 1 else 0.0,
        "MCC_majority_mean": float(split_summary_df["MCC_majority"].mean()),
        "MCC_majority_std": float(split_summary_df["MCC_majority"].std(ddof=1)) if len(split_summary_df) > 1 else 0.0,
    }
    overall_fp = os.path.join(args.out_dir, "scenario1_overall_metrics.json")
    with open(overall_fp, "w", encoding="utf-8") as f:
        json.dump(overall, f, ensure_ascii=False, indent=2)

    print("\n[DONE]")
    print("split summary:", split_summary_fp, split_summary_df.shape)
    print("level1 metrics:", all_level1_fp, all_level1_metrics.shape)
    print("majority predictions:", all_majority_fp, all_majority.shape)
    print("overall:", overall_fp)
    print(json.dumps(overall, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
