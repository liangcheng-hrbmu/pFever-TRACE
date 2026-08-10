#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Fast panel-b style held-out test validation from existing outputs.

Compared with the original slower script, this version:
- does NOT use Optuna
- does NOT do nested tuning/retraining
- trains ONE fixed-parameter XGBoost model per gene set
- evaluates on held-out test split only
- plots a faceted panel-b style figure: one column per cell type, two rows (WF1/BAS)

Inputs:
- 16_SHAP_input_data/*.h5ad
- 17_xgboost_results/*.split_indices.npz
- 18_SHAP_results/*.allCells.dSHAP_sum.tsv

Outputs:
- fig3b_metrics_raw.tsv
- fig3b_metrics_tidy.tsv
- fig3b_gene_sets.tsv
- Fig3b_testsplit_validation.png / pdf
"""

from __future__ import annotations

import argparse
import glob
import json
import os
from typing import Dict, List, Sequence, Tuple

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from sklearn.metrics import balanced_accuracy_score, f1_score, accuracy_score
from sklearn.model_selection import GroupShuffleSplit
from sklearn.preprocessing import LabelEncoder
import matplotlib.pyplot as plt
import xgboost as xgb
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

H5AD_SUFFIX = "_adataMerged_SPECTRAgenes.log1p.h5ad"
plt.rcParams["figure.dpi"] = 160
plt.rcParams["savefig.bbox"] = "tight"


def ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def clean_str(x):
    if pd.isna(x):
        return np.nan
    x = str(x).strip()
    if x in {"", "nan", "None", "NA"}:
        return np.nan
    return x


def celltype_from_h5ad(path: str) -> str:
    base = os.path.basename(path)
    if base.endswith(H5AD_SUFFIX):
        return base[: -len(H5AD_SUFFIX)]
    return os.path.splitext(base)[0]


def read_h5ad_df(h5ad_path: str, disease_key: str, sample_key: str, study_key: str):
    adata = sc.read_h5ad(h5ad_path)
    for key in [disease_key, sample_key, study_key]:
        if key not in adata.obs.columns:
            raise ValueError(f"{h5ad_path}: obs missing required column '{key}'")

    obs = adata.obs.copy()
    obs[disease_key] = obs[disease_key].map(clean_str)
    obs[sample_key] = obs[sample_key].map(clean_str)
    obs[study_key] = obs[study_key].map(clean_str)

    keep = obs[disease_key].notna() & obs[sample_key].notna() & obs[study_key].notna()
    adata = adata[keep].copy()
    obs = adata.obs.copy()
    obs[disease_key] = obs[disease_key].map(clean_str)
    obs[sample_key] = obs[sample_key].map(clean_str)
    obs[study_key] = obs[study_key].map(clean_str)

    X = adata.X
    if sp.issparse(X):
        X = X.tocsr()
    else:
        X = np.asarray(X)

    df = pd.DataFrame(index=np.arange(adata.n_obs))
    df["disease"] = obs[disease_key].astype(str).values
    df["sampleID"] = obs[sample_key].astype(str).values
    df["studyID"] = obs[study_key].astype(str).values
    df["obs_name"] = adata.obs_names.astype(str)
    return {"X": X, "meta": df}, adata.var_names.astype(str).tolist()


def load_split_indices(npz_path: str, n_obs: int, groups: Sequence[str], y: Sequence[str], random_state: int = 1234):
    z = np.load(npz_path)
    keys = set(z.files)
    test_idx = z["test_idx"].astype(int)
    train_idx = z["train_idx"].astype(int) if "train_idx" in keys else None
    val_idx = z["val_idx"].astype(int) if "val_idx" in keys else None

    if train_idx is not None and val_idx is not None:
        return train_idx, val_idx, test_idx

    all_idx = np.arange(n_obs)
    trainval_mask = np.ones(n_obs, dtype=bool)
    trainval_mask[test_idx] = False
    trainval_idx = all_idx[trainval_mask]

    gss = GroupShuffleSplit(n_splits=1, test_size=0.2, random_state=random_state)
    tr_rel, va_rel = next(gss.split(trainval_idx, y=np.array(y)[trainval_idx], groups=np.array(groups)[trainval_idx]))
    train_idx = trainval_idx[tr_rel]
    val_idx = trainval_idx[va_rel]
    return train_idx.astype(int), val_idx.astype(int), test_idx.astype(int)


def positive_fraction(X) -> np.ndarray:
    if sp.issparse(X):
        return np.asarray((X > 0).mean(axis=0)).ravel()
    return (np.asarray(X) > 0).mean(axis=0)


def load_dshap_ranking(path: str) -> List[str]:
    df = pd.read_csv(path, sep="\t")
    cols = [str(c).strip() for c in df.columns]
    df.columns = cols
    if "gene" not in df.columns:
        raise ValueError(f"dSHAP file missing gene column: {path}")
    if "dSHAP" not in df.columns:
        numeric_cols = [c for c in df.columns if c != "gene"]
        if not numeric_cols:
            raise ValueError(f"dSHAP file missing numeric score column: {path}")
        score_col = numeric_cols[0]
    else:
        score_col = "dSHAP"
    df[score_col] = pd.to_numeric(df[score_col], errors="coerce")
    df = df.dropna(subset=[score_col]).sort_values(score_col, ascending=False)
    return df["gene"].astype(str).tolist()


def build_gene_sets(
    var_names: List[str],
    X,
    dshap_ranked: List[str],
    top_k_list: Sequence[int],
    min_expr_frac: float,
    n_random: int,
    seed: int,
    include_all_genes: bool,
) -> Tuple[Dict[str, List[str]], pd.DataFrame]:
    frac = positive_fraction(X)
    gs = [g for g, f in zip(var_names, frac) if f >= min_expr_frac]
    if len(gs) == 0:
        raise RuntimeError("No genes remain after min_expr_frac filter.")
    gs_set = set(gs)

    gene_sets: Dict[str, List[str]] = {}
    rows = []
    rng = np.random.default_rng(seed)

    if include_all_genes:
        gene_sets["all_genes"] = list(gs)
        rows.append({"method": "all_genes", "top_n": len(gs), "random_id": np.nan, "n_genes": len(gs)})

    for k in top_k_list:
        dgenes = [g for g in dshap_ranked if g in gs_set][:k]
        gene_sets[f"dshap_top{k}"] = dgenes
        rows.append({"method": "dshap", "top_n": k, "random_id": np.nan, "n_genes": len(dgenes)})

        pool = [g for g in gs if g not in set(dgenes)]
        if len(pool) < len(dgenes):
            pool = list(gs)

        for rid in range(n_random):
            if len(dgenes) == 0:
                rgenes = []
            else:
                take = min(len(pool), len(dgenes))
                rgenes = rng.choice(pool, size=take, replace=False).tolist()
            gene_sets[f"random_top{k}_rep{rid:02d}"] = rgenes
            rows.append({"method": "random", "top_n": k, "random_id": rid, "n_genes": len(rgenes)})

    return gene_sets, pd.DataFrame(rows)


def make_model(n_classes: int, use_gpu: bool, xgb_n_jobs: int):
    kwargs = dict(
        objective="multi:softprob",
        num_class=n_classes,
        eval_metric="mlogloss",
        n_estimators=120,
        max_depth=3,
        learning_rate=0.08,
        subsample=0.8,
        colsample_bytree=0.8,
        min_child_weight=1.0,
        reg_alpha=0.0,
        reg_lambda=1.0,
        random_state=1234,
        n_jobs=int(xgb_n_jobs),
        tree_method="hist",
    )
    if use_gpu:
        kwargs["device"] = "cuda"
    return xgb.XGBClassifier(**kwargs)


def evaluate_gene_set(
    X,
    y_str,
    groups,
    var_names,
    genes: List[str],
    train_idx,
    val_idx,
    test_idx,
    use_gpu: bool,
    xgb_n_jobs: int,
):
    gene_to_idx = {g: i for i, g in enumerate(var_names)}
    feat_idx = [gene_to_idx[g] for g in genes if g in gene_to_idx]
    if len(feat_idx) == 0:
        return None

    if sp.issparse(X):
        Xs = X[:, feat_idx].toarray()
    else:
        Xs = np.asarray(X[:, feat_idx])

    y_str = np.asarray(y_str).astype(str)

    train_classes = sorted(pd.Series(y_str[train_idx]).dropna().unique().tolist())
    le = LabelEncoder().fit(train_classes)

    def mask_known(idx):
        ys = y_str[idx]
        keep = np.isin(ys, le.classes_)
        return idx[keep]

    train_idx2 = mask_known(np.asarray(train_idx, dtype=int))
    val_idx2 = mask_known(np.asarray(val_idx, dtype=int))
    test_idx2 = mask_known(np.asarray(test_idx, dtype=int))

    if len(train_idx2) == 0 or len(test_idx2) == 0:
        return None
    if len(np.unique(y_str[train_idx2])) < 2:
        return None

    # Use train+val to fit once, then evaluate on held-out test.
    trainval_idx2 = np.concatenate([train_idx2, val_idx2]) if len(val_idx2) > 0 else train_idx2
    X_trainval = Xs[trainval_idx2]
    y_trainval = le.transform(y_str[trainval_idx2])
    X_test = Xs[test_idx2]
    y_test = le.transform(y_str[test_idx2])

    model = make_model(n_classes=len(le.classes_), use_gpu=use_gpu, xgb_n_jobs=xgb_n_jobs)
    model.fit(X_trainval, y_trainval, verbose=False)
    pred_test = model.predict(X_test)

    return {
        "WF1": float(f1_score(y_test, pred_test, average="weighted")),
        "BAS": float(balanced_accuracy_score(y_test, pred_test)),
        "ACC": float(accuracy_score(y_test, pred_test)),
        "best_val_WF1": np.nan,
        "n_test": int(len(test_idx2)),
        "n_classes_train": int(len(le.classes_)),
        "dropped_test_unseen_labels": int(len(test_idx) - len(test_idx2)),
        "used_genes": int(len(feat_idx)),
    }


def load_existing_disease_metrics(model_dir: str) -> pd.DataFrame:
    rows = []
    for fp in sorted(glob.glob(os.path.join(model_dir, "*.metrics.json"))):
        if fp.endswith(".studyID.metrics.json"):
            continue
        name = os.path.basename(fp).replace(".metrics.json", "")
        with open(fp, "r", encoding="utf-8") as f:
            obj = json.load(f)
        row = {
            "celltype": name,
            "WF1_main_model": obj.get("WF1", obj.get("wf1", np.nan)),
            "BAS_main_model": obj.get("BAS", obj.get("bas", np.nan)),
            "ACC_main_model": obj.get("ACC", obj.get("acc", np.nan)),
            "metrics_file": fp,
        }
        rows.append(row)
    return pd.DataFrame(rows)


def plot_panel_b(metrics_raw: pd.DataFrame, out_png: str, out_pdf: str, include_all_genes: bool = True):
    if metrics_raw.empty:
        return

    celltypes = metrics_raw["celltype"].drop_duplicates().tolist()
    if len(celltypes) == 0:
        return

    color_map = {
        "5":  "#7DA0F9",
        "10": "#F3A26B",
        "20": "#86D487",
    }

    metrics_order = ["WF1", "BAS"]
    ncols = len(celltypes)
    fig, axes = plt.subplots(
        2, ncols,
        figsize=(1.9 * ncols, 4.8),
        sharey="row"
    )

    if ncols == 1:
        axes = np.array(axes).reshape(2, 1)

    for col_i, ct in enumerate(celltypes):
        sub_ct = metrics_raw[metrics_raw["celltype"] == ct].copy()

        x_map = {"5": 0, "10": 1, "20": 2}
        if include_all_genes:
            x_map["GS"] = 3

        xticks = [x_map["5"], x_map["10"], x_map["20"]]
        xticklabels = ["5", "10", "20"]

        if include_all_genes:
            xticks.append(x_map["GS"])
            gs_vals = sub_ct.loc[sub_ct["method"] == "all_genes", "used_genes"].dropna().unique().tolist()
            gs_lab = str(int(gs_vals[0])) if len(gs_vals) else "GS"
            xticklabels.append(gs_lab)

        for row_i, metric in enumerate(metrics_order):
            ax = axes[row_i, col_i]

            for k in ["5", "10", "20"]:
                vals = sub_ct.loc[
                    (sub_ct["method"] == "random") &
                    (sub_ct["top_n_label"] == k),
                    metric
                ].dropna().values

                if len(vals) > 0:
                    ax.boxplot(
                        [vals],
                        positions=[x_map[k]],
                        widths=0.55,
                        patch_artist=True,
                        boxprops=dict(facecolor="white", edgecolor=color_map[k], linewidth=1.2),
                        medianprops=dict(color=color_map[k], linewidth=1.2),
                        whiskerprops=dict(color=color_map[k], linewidth=1.0),
                        capprops=dict(color=color_map[k], linewidth=1.0),
                        flierprops=dict(
                            marker='o', markersize=2.5,
                            markerfacecolor=color_map[k],
                            markeredgecolor=color_map[k],
                            alpha=0.6
                        ),
                    )

            for k in ["5", "10", "20"]:
                vals = sub_ct.loc[
                    (sub_ct["method"] == "dshap") &
                    (sub_ct["top_n_label"] == k),
                    metric
                ].dropna().values
                if len(vals) > 0:
                    ax.scatter(
                        [x_map[k]] * len(vals),
                        vals,
                        marker="+",
                        s=28,
                        c="black",
                        linewidths=1.1,
                        zorder=4
                    )

            if include_all_genes:
                vals = sub_ct.loc[
                    sub_ct["method"] == "all_genes",
                    metric
                ].dropna().values
                if len(vals) > 0:
                    ax.scatter(
                        [x_map["GS"]] * len(vals),
                        vals,
                        marker="x",
                        s=18,
                        c="forestgreen",
                        linewidths=1.1,
                        zorder=4
                    )

            if row_i == 0:
                ax.set_title(ct, fontsize=10, pad=10)

            if col_i == 0:
                ax.set_ylabel(metric, fontsize=11)
            else:
                ax.set_ylabel("")

            ax.set_xticks(xticks)
            if row_i == 1:
                ax.set_xticklabels(xticklabels, rotation=45, ha="right", fontsize=8)
            else:
                ax.set_xticklabels([])

            ax.set_ylim(0.15, 0.95)
            ax.grid(axis="y", alpha=0.18, linewidth=0.7)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.spines["left"].set_linewidth(0.8)
            ax.spines["bottom"].set_linewidth(0.8)
            ax.tick_params(axis="both", labelsize=8, width=0.8)

    legend_handles = [
        Patch(facecolor="white", edgecolor=color_map["5"],  label="Random selection (top 5 genes)"),
        Patch(facecolor="white", edgecolor=color_map["10"], label="Random selection (top 10 genes)"),
        Patch(facecolor="white", edgecolor=color_map["20"], label="Random selection (top 20 genes)"),
        Line2D([0], [0], marker="x", color="forestgreen", linestyle="None", markersize=5, label="No selection"),
        Line2D([0], [0], marker="+", color="black", linestyle="None", markersize=6, label="d-SHAP-based selection"),
    ]

    fig.legend(
        handles=legend_handles,
        loc="lower center",
        ncol=5,
        frameon=False,
        fontsize=9,
        bbox_to_anchor=(0.5, -0.02)
    )

    fig.subplots_adjust(wspace=0.10, hspace=0.03, bottom=0.23)
    fig.savefig(out_png, dpi=300)
    fig.savefig(out_pdf)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description="Fast panel-b style held-out test validation from existing outputs.")
    ap.add_argument("--input_dir", required=True, help="16_SHAP_input_data directory")
    ap.add_argument("--split_dir", required=True, help="17_xgboost_results directory containing *.split_indices.npz")
    ap.add_argument("--dshap_dir", required=True, help="18_SHAP_results directory containing *.allCells.dSHAP_sum.tsv")
    ap.add_argument("--out_dir", required=True)
    ap.add_argument("--disease_key", default="disease")
    ap.add_argument("--sample_key", default="sampleID")
    ap.add_argument("--study_key", default="studyID")
    ap.add_argument("--top_k", nargs="*", type=int, default=[5, 10, 20])
    ap.add_argument("--n_random", type=int, default=20)
    ap.add_argument("--min_expr_frac", type=float, default=0.05)
    ap.add_argument("--seed", type=int, default=1234)
    ap.add_argument("--use_gpu", action="store_true")
    ap.add_argument("--include_all_genes", action="store_true")
    ap.add_argument("--xgb_n_jobs", type=int, default=8)
    ap.add_argument("--only", nargs="*", default=None)
    args = ap.parse_args()

    ensure_dir(args.out_dir)
    np.random.seed(args.seed)

    h5ads = sorted(glob.glob(os.path.join(args.input_dir, f"*{H5AD_SUFFIX}")))
    if not h5ads:
        raise FileNotFoundError(f"No per-celltype h5ad found in {args.input_dir}")

    only_set = set(args.only) if args.only else None
    results = []
    geneset_rows = []

    for h5ad in h5ads:
        ct = celltype_from_h5ad(h5ad)
        if only_set is not None and ct not in only_set:
            continue

        split_fp = os.path.join(args.split_dir, f"{ct}.split_indices.npz")
        dshap_fp = os.path.join(args.dshap_dir, f"{ct}.allCells.dSHAP_sum.tsv")
        if not (os.path.exists(split_fp) and os.path.exists(dshap_fp)):
            print(f"[SKIP] {ct}: missing split or dSHAP file")
            continue

        print(f"\n[RUN] {ct}")
        data, var_names = read_h5ad_df(h5ad, args.disease_key, args.sample_key, args.study_key)
        X = data["X"]
        meta = data["meta"]

        train_idx, val_idx, test_idx = load_split_indices(
            split_fp,
            n_obs=len(meta),
            groups=meta["sampleID"].values,
            y=meta["disease"].values,
            random_state=args.seed
        )

        dshap_ranked = load_dshap_ranking(dshap_fp)
        gene_sets, gene_set_df = build_gene_sets(
            var_names=var_names,
            X=X,
            dshap_ranked=dshap_ranked,
            top_k_list=args.top_k,
            min_expr_frac=args.min_expr_frac,
            n_random=args.n_random,
            seed=args.seed,
            include_all_genes=args.include_all_genes,
        )
        gene_set_df.insert(0, "celltype", ct)
        geneset_rows.append(gene_set_df)

        for set_name, genes in gene_sets.items():
            rec = evaluate_gene_set(
                X=X,
                y_str=meta["disease"].values,
                groups=meta["sampleID"].values,
                var_names=var_names,
                genes=genes,
                train_idx=train_idx,
                val_idx=val_idx,
                test_idx=test_idx,
                use_gpu=args.use_gpu,
                xgb_n_jobs=args.xgb_n_jobs,
            )
            if rec is None:
                print(f"  [SKIP] {ct} {set_name}: no usable data after label filtering")
                continue

            if set_name == "all_genes":
                method = "all_genes"
                top_n = np.nan
                random_id = np.nan
                top_label = "GS"
            elif set_name.startswith("dshap_top"):
                method = "dshap"
                top_n = int(set_name.replace("dshap_top", ""))
                random_id = np.nan
                top_label = str(top_n)
            else:
                method = "random"
                tmp = set_name.replace("random_top", "")
                top_n = int(tmp.split("_rep")[0])
                random_id = int(tmp.split("_rep")[1])
                top_label = str(top_n)

            results.append({
                "celltype": ct,
                "gene_set_name": set_name,
                "method": method,
                "top_n": top_n,
                "top_n_label": top_label,
                "random_id": random_id,
                **rec,
            })
            print(f"  [OK] {ct} {set_name}: WF1={rec['WF1']:.3f} BAS={rec['BAS']:.3f} used_genes={rec['used_genes']}")

    metrics_raw = pd.DataFrame(results)
    genesets_df = pd.concat(geneset_rows, axis=0, ignore_index=True) if geneset_rows else pd.DataFrame()
    metrics_raw.to_csv(os.path.join(args.out_dir, "fig3b_metrics_raw.tsv"), sep="\t", index=False)
    genesets_df.to_csv(os.path.join(args.out_dir, "fig3b_gene_sets.tsv"), sep="\t", index=False)

    if not metrics_raw.empty:
        tidy = metrics_raw.melt(
            id_vars=[c for c in metrics_raw.columns if c not in ["WF1", "BAS", "ACC"]],
            value_vars=[c for c in ["WF1", "BAS", "ACC"] if c in metrics_raw.columns],
            var_name="metric",
            value_name="value",
        )
    else:
        tidy = pd.DataFrame()
    tidy.to_csv(os.path.join(args.out_dir, "fig3b_metrics_tidy.tsv"), sep="\t", index=False)

    main_df = load_existing_disease_metrics(args.split_dir)
    main_df.to_csv(os.path.join(args.out_dir, "fig3b_existing_main_model_metrics.tsv"), sep="\t", index=False)

    plot_panel_b(
        metrics_raw=metrics_raw,
        out_png=os.path.join(args.out_dir, "Fig3b_testsplit_validation.png"),
        out_pdf=os.path.join(args.out_dir, "Fig3b_testsplit_validation.pdf"),
        include_all_genes=args.include_all_genes,
    )

    print("\n[DONE]")
    print("metrics_raw:", os.path.join(args.out_dir, "fig3b_metrics_raw.tsv"))
    print("metrics_tidy:", os.path.join(args.out_dir, "fig3b_metrics_tidy.tsv"))
    print("figure_png:", os.path.join(args.out_dir, "Fig3b_testsplit_validation.png"))
    print("figure_pdf:", os.path.join(args.out_dir, "Fig3b_testsplit_validation.pdf"))


if __name__ == "__main__":
    main()
