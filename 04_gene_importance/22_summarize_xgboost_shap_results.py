#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
汇总多细胞类型 XGBoost + SHAP + dSHAP/sSHAP 比较结果。

输入目录默认对应你当前流程：
- 17_xgboost_results                 disease 分类模型
- 18_SHAP_results                    dSHAP 结果
- 19_studyID_xgboost_results        studyID 分类模型
- 20_SHAP_results_studyID            sSHAP 结果
- 21_compare_dshap_sshap_results     dSHAP vs sSHAP 比较结果

输出：
- 多个 summary tsv
- 多张 png 图
- 一个 summary_report.md
"""

from __future__ import annotations

import argparse
import glob
import json
import os
from typing import Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


plt.rcParams["figure.dpi"] = 150
plt.rcParams["savefig.bbox"] = "tight"


METRIC_BASE_COLS = [
    "celltype", "task", "WF1", "BAS", "ACC",
    "used_gpu", "best_iteration_on_val", "final_n_estimators", "metrics_file"
]


def ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)



def load_json(path: str) -> Dict:
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)



def strip_suffix(name: str, suffix: str) -> str:
    if name.endswith(suffix):
        return name[: -len(suffix)]
    return name



def safe_read_tsv(path: str) -> Optional[pd.DataFrame]:
    if not os.path.exists(path):
        return None
    try:
        return pd.read_csv(path, sep="\t")
    except Exception:
        return None



def safe_spearman(x: pd.Series, y: pd.Series) -> float:
    x = pd.to_numeric(x, errors="coerce")
    y = pd.to_numeric(y, errors="coerce")
    mask = x.notna() & y.notna()
    if mask.sum() < 2:
        return float("nan")
    xr = x[mask].rank(method="average")
    yr = y[mask].rank(method="average")
    return float(xr.corr(yr, method="pearson"))



def safe_pearson(x: pd.Series, y: pd.Series) -> float:
    x = pd.to_numeric(x, errors="coerce")
    y = pd.to_numeric(y, errors="coerce")
    mask = x.notna() & y.notna()
    if mask.sum() < 2:
        return float("nan")
    return float(x[mask].corr(y[mask], method="pearson"))



def save_df(df: pd.DataFrame, path: str) -> None:
    df.to_csv(path, sep="\t", index=False)
    print(f"[WRITE] {path} {df.shape}")



def metric_get(obj: Dict, key: str):
    aliases = {
        "WF1": ["WF1", "wf1", "weighted_f1", "weighted_F1"],
        "BAS": ["BAS", "bas", "balanced_accuracy", "balanced_acc"],
        "ACC": ["ACC", "acc", "accuracy"],
        "used_gpu": ["used_gpu"],
        "best_iteration_on_val": ["best_iteration_on_val", "best_iteration"],
        "final_n_estimators": ["final_n_estimators", "n_estimators"],
    }
    for k in aliases.get(key, [key]):
        if k in obj:
            return obj.get(k)
    if isinstance(obj.get("metrics"), dict):
        for k in aliases.get(key, [key]):
            if k in obj["metrics"]:
                return obj["metrics"].get(k)
    return None



def empty_metrics_df() -> pd.DataFrame:
    return pd.DataFrame(columns=METRIC_BASE_COLS)



def collect_metrics(metrics_dir: str, study_id: bool) -> pd.DataFrame:
    if study_id:
        files = sorted(glob.glob(os.path.join(metrics_dir, "*.studyID.metrics.json")))
        suffix = ".studyID.metrics.json"
        task = "studyID"
    else:
        files = sorted(glob.glob(os.path.join(metrics_dir, "*.metrics.json")))
        files = [f for f in files if not f.endswith(".studyID.metrics.json")]
        suffix = ".metrics.json"
        task = "disease"

    print(f"[INFO] collect_metrics(task={task}) matched files: {len(files)} in {metrics_dir}")
    if files:
        print(f"[INFO] first matched metrics file: {files[0]}")

    rows: List[Dict] = []
    for fp in files:
        name = os.path.basename(fp)
        celltype = strip_suffix(name, suffix)
        try:
            obj = load_json(fp)
        except Exception as e:
            print(f"[WARN] failed to read json: {fp} ({e})")
            continue

        row = {
            "celltype": celltype,
            "task": task,
            "WF1": metric_get(obj, "WF1"),
            "BAS": metric_get(obj, "BAS"),
            "ACC": metric_get(obj, "ACC"),
            "used_gpu": metric_get(obj, "used_gpu"),
            "best_iteration_on_val": metric_get(obj, "best_iteration_on_val"),
            "final_n_estimators": metric_get(obj, "final_n_estimators"),
            "metrics_file": fp,
        }
        rows.append(row)

    if not rows:
        print(f"[WARN] no usable metrics rows collected for task={task}")
        return empty_metrics_df()

    df = pd.DataFrame(rows)
    for col in METRIC_BASE_COLS:
        if col not in df.columns:
            df[col] = np.nan
    for col in ["WF1", "BAS", "ACC", "best_iteration_on_val", "final_n_estimators"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    df = df[METRIC_BASE_COLS].copy()
    sort_cols = [c for c in ["WF1", "BAS", "ACC"] if c in df.columns]
    if sort_cols:
        df = df.sort_values(sort_cols, ascending=[False] * len(sort_cols), na_position="last").reset_index(drop=True)
    return df



def select_metric_cols(df: pd.DataFrame, prefix: str) -> pd.DataFrame:
    base = pd.DataFrame(columns=["celltype", f"WF1_{prefix}", f"BAS_{prefix}", f"ACC_{prefix}"])
    if df is None or df.empty:
        return base
    out = df.copy()
    for col in ["celltype", "WF1", "BAS", "ACC"]:
        if col not in out.columns:
            out[col] = np.nan
    out = out[["celltype", "WF1", "BAS", "ACC"]].copy()
    out = out.rename(columns={"WF1": f"WF1_{prefix}", "BAS": f"BAS_{prefix}", "ACC": f"ACC_{prefix}"})
    return out



def collect_shap_top_genes(shap_dir: str, mode: str, kind: str, top_n: int) -> pd.DataFrame:
    if kind == "d":
        pattern = os.path.join(shap_dir, f"*.{mode}.dSHAP_sum.tsv")
        value_col = "dSHAP"
        suffix = f".{mode}.dSHAP_sum.tsv"
    elif kind == "s":
        pattern = os.path.join(shap_dir, f"*.{mode}.studyID.sSHAP_sum.tsv")
        value_col = "sSHAP"
        suffix = f".{mode}.studyID.sSHAP_sum.tsv"
    else:
        raise ValueError("kind must be 'd' or 's'")

    files = sorted(glob.glob(pattern))
    print(f"[INFO] collect_shap_top_genes(kind={kind}, mode={mode}) matched files: {len(files)} in {shap_dir}")
    rows: List[pd.DataFrame] = []
    for fp in files:
        df = safe_read_tsv(fp)
        if df is None or "gene" not in df.columns or value_col not in df.columns:
            continue
        celltype = strip_suffix(os.path.basename(fp), suffix)
        tmp = df[["gene", value_col]].copy()
        tmp[value_col] = pd.to_numeric(tmp[value_col], errors="coerce").fillna(0.0)
        tmp = tmp.sort_values(value_col, ascending=False).head(top_n).copy()
        tmp.insert(0, "celltype", celltype)
        tmp.insert(1, "rank", np.arange(1, len(tmp) + 1))
        tmp.insert(2, "mode", mode)
        rows.append(tmp)

    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()



def infer_compare_prefixes(compare_dir: str, mode: str) -> List[Tuple[str, str]]:
    files = sorted(glob.glob(os.path.join(compare_dir, f"*.{mode}.studyID.dSHAP_vs_sSHAP.tsv")))
    print(f"[INFO] infer_compare_prefixes(mode={mode}) matched files: {len(files)} in {compare_dir}")
    prefixes: List[Tuple[str, str]] = []
    for fp in files:
        base = os.path.basename(fp)
        prefix = strip_suffix(base, ".dSHAP_vs_sSHAP.tsv")
        celltype = strip_suffix(prefix, f".{mode}.studyID")
        prefixes.append((celltype, prefix))
    return prefixes



def collect_compare(compare_dir: str, mode: str, top_n: int) -> Tuple[pd.DataFrame, pd.DataFrame]:
    prefixes = infer_compare_prefixes(compare_dir, mode)
    summary_rows: List[Dict] = []
    top_rows: List[pd.DataFrame] = []

    for celltype, prefix in prefixes:
        all_fp = os.path.join(compare_dir, f"{prefix}.dSHAP_vs_sSHAP.tsv")
        top_fp = os.path.join(compare_dir, f"{prefix}.top20_priority_genes.tsv")
        high_fp = os.path.join(compare_dir, f"{prefix}.highD_lowS_genes.tsv")

        df = safe_read_tsv(all_fp)
        if df is None or not {"gene", "dSHAP", "sSHAP"}.issubset(df.columns):
            continue

        for col in ["dSHAP", "sSHAP", "d_over_s", "d_minus_s", "priority_zscore"]:
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors="coerce")

        if "priority_zscore" not in df.columns:
            eps = 1e-8
            d_z = (df["dSHAP"] - df["dSHAP"].mean()) / (df["dSHAP"].std(ddof=0) + eps)
            s_z = (df["sSHAP"] - df["sSHAP"].mean()) / (df["sSHAP"].std(ddof=0) + eps)
            df["priority_zscore"] = d_z - s_z

        pearson = safe_pearson(df["dSHAP"], df["sSHAP"])
        spearman = safe_spearman(df["dSHAP"], df["sSHAP"])
        n_high = 0
        if os.path.exists(high_fp):
            high_df = safe_read_tsv(high_fp)
            if high_df is not None:
                n_high = int(high_df.shape[0])

        best_row = df.sort_values(["priority_zscore", "d_over_s", "dSHAP"], ascending=[False, False, False]).iloc[0]
        summary_rows.append(
            {
                "celltype": celltype,
                "mode": mode,
                "n_genes": int(df.shape[0]),
                "pearson_dshap_sshap": pearson,
                "spearman_dshap_sshap": spearman,
                "mean_dSHAP": float(df["dSHAP"].mean()),
                "mean_sSHAP": float(df["sSHAP"].mean()),
                "median_d_over_s": float(pd.to_numeric(df.get("d_over_s", pd.Series(dtype=float)), errors="coerce").median()),
                "n_highD_lowS": n_high,
                "top_priority_gene": best_row.get("gene"),
                "top_priority_zscore": float(best_row.get("priority_zscore")),
                "top_gene_dSHAP": float(best_row.get("dSHAP")),
                "top_gene_sSHAP": float(best_row.get("sSHAP")),
                "compare_file": all_fp,
            }
        )

        if os.path.exists(top_fp):
            top_df = safe_read_tsv(top_fp)
        else:
            top_df = None
        if top_df is None or top_df.empty:
            top_df = df.sort_values(["priority_zscore", "d_over_s", "dSHAP"], ascending=[False, False, False]).head(top_n).copy()
        else:
            top_df = top_df.head(top_n).copy()

        keep_cols = [c for c in ["gene", "dSHAP", "sSHAP", "d_over_s", "d_minus_s", "priority_zscore"] if c in top_df.columns]
        top_df = top_df[keep_cols].copy()
        top_df.insert(0, "celltype", celltype)
        top_df.insert(1, "rank", np.arange(1, len(top_df) + 1))
        top_df.insert(2, "mode", mode)
        top_rows.append(top_df)

    summary_df = pd.DataFrame(summary_rows)
    if not summary_df.empty:
        summary_df = summary_df.sort_values(["spearman_dshap_sshap", "pearson_dshap_sshap"], ascending=[True, True]).reset_index(drop=True)

    top_df = pd.concat(top_rows, ignore_index=True) if top_rows else pd.DataFrame()
    return summary_df, top_df



def build_frequency_table(top_df: pd.DataFrame) -> pd.DataFrame:
    if top_df.empty:
        return pd.DataFrame()
    grp = top_df.groupby("gene", as_index=False).agg(
        celltype_count=("celltype", "nunique"),
        occurrences=("gene", "size"),
        mean_priority_zscore=("priority_zscore", "mean"),
        median_priority_zscore=("priority_zscore", "median"),
        mean_dSHAP=("dSHAP", "mean"),
        mean_sSHAP=("sSHAP", "mean"),
    )
    return grp.sort_values(["celltype_count", "mean_priority_zscore", "mean_dSHAP"], ascending=[False, False, False]).reset_index(drop=True)



def build_heatmap_matrix(top_df: pd.DataFrame, top_gene_list: List[str]) -> pd.DataFrame:
    if top_df.empty or not top_gene_list:
        return pd.DataFrame()
    sub = top_df[top_df["gene"].isin(top_gene_list)].copy()
    mat = sub.pivot_table(index="gene", columns="celltype", values="priority_zscore", aggfunc="max", fill_value=0.0)
    mat = mat.reindex(top_gene_list)
    return mat



def plot_hbar(df: pd.DataFrame, y_col: str, x_col: str, title: str, out_fp: str) -> None:
    if df.empty or x_col not in df.columns or y_col not in df.columns:
        return
    plot_df = df[[y_col, x_col]].dropna().sort_values(x_col, ascending=True)
    if plot_df.empty:
        return
    height = max(4, 0.35 * len(plot_df) + 1.5)
    plt.figure(figsize=(8, height))
    plt.barh(plot_df[y_col], plot_df[x_col])
    plt.xlabel(x_col)
    plt.ylabel(y_col)
    plt.title(title)
    plt.tight_layout()
    plt.savefig(out_fp)
    plt.close()
    print(f"[WRITE] {out_fp}")



def plot_scatter_merge(df: pd.DataFrame, x_col: str, y_col: str, label_col: str, title: str, out_fp: str) -> None:
    if df.empty or x_col not in df.columns or y_col not in df.columns:
        return
    sub = df[[label_col, x_col, y_col]].dropna().copy()
    if sub.empty:
        return
    plt.figure(figsize=(6.5, 6))
    plt.scatter(sub[x_col], sub[y_col])
    for _, row in sub.iterrows():
        plt.annotate(str(row[label_col]), (row[x_col], row[y_col]), fontsize=7, alpha=0.8)
    xmin = float(sub[x_col].min())
    xmax = float(sub[x_col].max())
    ymin = float(sub[y_col].min())
    ymax = float(sub[y_col].max())
    lo = min(xmin, ymin)
    hi = max(xmax, ymax)
    plt.plot([lo, hi], [lo, hi], linestyle="--")
    plt.xlabel(x_col)
    plt.ylabel(y_col)
    plt.title(title)
    plt.tight_layout()
    plt.savefig(out_fp)
    plt.close()
    print(f"[WRITE] {out_fp}")



def plot_gene_frequency(freq_df: pd.DataFrame, out_fp: str, top_k: int = 20) -> None:
    if freq_df.empty:
        return
    plot_df = freq_df.head(top_k).copy().sort_values(["celltype_count", "mean_priority_zscore"], ascending=[True, True])
    plt.figure(figsize=(8, max(4, 0.35 * len(plot_df) + 1.5)))
    plt.barh(plot_df["gene"], plot_df["celltype_count"])
    plt.xlabel("Number of cell types in top priority list")
    plt.ylabel("Gene")
    plt.title(f"Top recurrent priority genes (top {top_k})")
    plt.tight_layout()
    plt.savefig(out_fp)
    plt.close()
    print(f"[WRITE] {out_fp}")



def plot_heatmap(mat: pd.DataFrame, title: str, out_fp: str) -> None:
    if mat.empty:
        return
    arr = mat.values.astype(float)
    fig_w = max(7, 0.5 * max(1, arr.shape[1]) + 2)
    fig_h = max(5, 0.35 * max(1, arr.shape[0]) + 2)
    plt.figure(figsize=(fig_w, fig_h))
    im = plt.imshow(arr, aspect="auto")
    plt.colorbar(im, fraction=0.03, pad=0.02, label="priority_zscore")
    plt.xticks(np.arange(arr.shape[1]), mat.columns, rotation=90)
    plt.yticks(np.arange(arr.shape[0]), mat.index)
    plt.title(title)
    plt.tight_layout()
    plt.savefig(out_fp)
    plt.close()
    print(f"[WRITE] {out_fp}")



def write_report(
    out_fp: str,
    disease_df: pd.DataFrame,
    study_df: pd.DataFrame,
    merged_metrics: pd.DataFrame,
    compare_df: pd.DataFrame,
    freq_df: pd.DataFrame,
) -> None:
    lines: List[str] = []
    lines.append("# XGBoost + SHAP 结果汇总\n")

    lines.append("## 1. 模型层面\n")
    lines.append(f"- disease 模型细胞类型数: {int(disease_df.shape[0])}\n")
    lines.append(f"- studyID 模型细胞类型数: {int(study_df.shape[0])}\n")
    if not disease_df.empty and "WF1" in disease_df.columns:
        top_row = disease_df.dropna(subset=["WF1"]).head(1)
        if not top_row.empty:
            lines.append(f"- disease WF1 中位数: {disease_df['WF1'].median():.4f}; 最好: {top_row.iloc[0]['celltype']} ({top_row.iloc[0]['WF1']:.4f})\n")
    if not study_df.empty and "WF1" in study_df.columns:
        top_row = study_df.dropna(subset=["WF1"]).head(1)
        if not top_row.empty:
            lines.append(f"- studyID WF1 中位数: {study_df['WF1'].median():.4f}; 最好: {top_row.iloc[0]['celltype']} ({top_row.iloc[0]['WF1']:.4f})\n")

    if not merged_metrics.empty and {"WF1_disease", "WF1_studyID"}.issubset(merged_metrics.columns):
        corr = merged_metrics["WF1_disease"].corr(merged_metrics["WF1_studyID"])
        if corr is not None and not pd.isna(corr):
            lines.append(f"- disease vs studyID 的 WF1 Pearson 相关: {float(corr):.4f}\n")

    lines.append("\n## 2. dSHAP vs sSHAP 比较\n")
    lines.append(f"- compare 文件细胞类型数: {int(compare_df.shape[0])}\n")
    if not compare_df.empty and "spearman_dshap_sshap" in compare_df.columns:
        lines.append(f"- dSHAP vs sSHAP Spearman 中位数: {compare_df['spearman_dshap_sshap'].median():.4f}\n")
        low_corr = compare_df.sort_values("spearman_dshap_sshap", ascending=True).head(min(5, len(compare_df)))
        lines.append("- 相关性最低的细胞类型（越低说明 disease 信号和 studyID 信号越可分离）:\n")
        for _, row in low_corr.iterrows():
            lines.append(
                f"  - {row['celltype']}: spearman={row['spearman_dshap_sshap']:.4f}, top_priority_gene={row.get('top_priority_gene')}, highD_lowS={int(row.get('n_highD_lowS', 0))}\n"
            )

    lines.append("\n## 3. 高频 priority genes\n")
    lines.append(f"- priority gene 去重后总数: {int(freq_df.shape[0])}\n")
    if not freq_df.empty:
        lines.append("- 出现频率最高的前 10 个基因:\n")
        for _, row in freq_df.head(min(10, len(freq_df))).iterrows():
            lines.append(f"  - {row['gene']}: {int(row['celltype_count'])} 个 cell type, mean_priority_zscore={row['mean_priority_zscore']:.4f}\n")

    lines.append("\n## 4. 建议你正文里放的核心图\n")
    lines.append("1. disease/studyID 两套模型的 WF1 排名图\n")
    lines.append("2. disease WF1 vs studyID WF1 散点图\n")
    lines.append("3. dSHAP vs sSHAP 相关性条形图\n")
    lines.append("4. 高频 priority genes 条形图\n")
    lines.append("5. 高频 priority genes × cell type 热图\n")

    with open(out_fp, "w", encoding="utf-8") as f:
        f.writelines(lines)
    print(f"[WRITE] {out_fp}")



def main() -> None:
    ap = argparse.ArgumentParser(description="Summarize multi-celltype XGBoost + SHAP results.")
    ap.add_argument("--disease_metrics_dir", default="data/analysis/08_gene_importance/17_xgboost_results")
    ap.add_argument("--dshap_dir", default="data/analysis/08_gene_importance/18_SHAP_results")
    ap.add_argument("--study_metrics_dir", default="data/analysis/08_gene_importance/19_studyID_xgboost_results")
    ap.add_argument("--sshap_dir", default="data/analysis/08_gene_importance/20_SHAP_results_studyID")
    ap.add_argument("--compare_dir", default="data/analysis/08_gene_importance/21_compare_dshap_sshap_results")
    ap.add_argument("--out_dir", required=True)
    ap.add_argument("--mode", choices=["allCells", "testOnly"], default="allCells")
    ap.add_argument("--top_n", type=int, default=20, help="Top genes kept per cell type")
    ap.add_argument("--heatmap_top_genes", type=int, default=30, help="Number of recurrent genes in heatmap")
    args = ap.parse_args()

    ensure_dir(args.out_dir)

    disease_df = collect_metrics(args.disease_metrics_dir, study_id=False)
    study_df = collect_metrics(args.study_metrics_dir, study_id=True)

    dshap_top_df = collect_shap_top_genes(args.dshap_dir, mode=args.mode, kind="d", top_n=args.top_n)
    sshap_top_df = collect_shap_top_genes(args.sshap_dir, mode=args.mode, kind="s", top_n=args.top_n)
    compare_df, priority_top_df = collect_compare(args.compare_dir, mode=args.mode, top_n=args.top_n)
    freq_df = build_frequency_table(priority_top_df)

    merged_metrics = pd.merge(
        select_metric_cols(disease_df, prefix="disease"),
        select_metric_cols(study_df, prefix="studyID"),
        on="celltype",
        how="outer",
    )
    if not merged_metrics.empty:
        for col in ["WF1_disease", "WF1_studyID"]:
            if col not in merged_metrics.columns:
                merged_metrics[col] = np.nan
        merged_metrics = merged_metrics.sort_values(["WF1_disease", "WF1_studyID"], ascending=[False, False], na_position="last").reset_index(drop=True)

    save_df(disease_df, os.path.join(args.out_dir, "summary_disease_model_metrics.tsv"))
    save_df(study_df, os.path.join(args.out_dir, "summary_studyID_model_metrics.tsv"))
    save_df(merged_metrics, os.path.join(args.out_dir, "summary_model_metrics_merged.tsv"))
    save_df(dshap_top_df, os.path.join(args.out_dir, f"summary_dshap_top{args.top_n}_{args.mode}.tsv"))
    save_df(sshap_top_df, os.path.join(args.out_dir, f"summary_sshap_top{args.top_n}_{args.mode}.tsv"))
    save_df(compare_df, os.path.join(args.out_dir, f"summary_compare_by_celltype_{args.mode}.tsv"))
    save_df(priority_top_df, os.path.join(args.out_dir, f"summary_priority_top{args.top_n}_by_celltype_{args.mode}.tsv"))
    save_df(freq_df, os.path.join(args.out_dir, f"summary_priority_gene_frequency_{args.mode}.tsv"))

    plot_hbar(disease_df, y_col="celltype", x_col="WF1", title="Disease model WF1 by cell type", out_fp=os.path.join(args.out_dir, "fig1_disease_model_WF1.png"))
    plot_hbar(study_df, y_col="celltype", x_col="WF1", title="StudyID model WF1 by cell type", out_fp=os.path.join(args.out_dir, "fig2_studyID_model_WF1.png"))
    plot_scatter_merge(merged_metrics, x_col="WF1_disease", y_col="WF1_studyID", label_col="celltype", title="Disease WF1 vs StudyID WF1", out_fp=os.path.join(args.out_dir, "fig3_disease_vs_studyID_WF1.png"))
    if not compare_df.empty and "spearman_dshap_sshap" in compare_df.columns:
        plot_hbar(compare_df.sort_values("spearman_dshap_sshap", ascending=True), y_col="celltype", x_col="spearman_dshap_sshap", title=f"dSHAP vs sSHAP Spearman by cell type ({args.mode})", out_fp=os.path.join(args.out_dir, f"fig4_dshap_vs_sshap_spearman_{args.mode}.png"))
    plot_gene_frequency(freq_df, out_fp=os.path.join(args.out_dir, f"fig5_priority_gene_frequency_{args.mode}.png"), top_k=min(20, max(1, args.heatmap_top_genes)))

    top_genes = freq_df.head(args.heatmap_top_genes)["gene"].tolist() if not freq_df.empty else []
    mat = build_heatmap_matrix(priority_top_df, top_genes)
    if not mat.empty:
        mat.to_csv(os.path.join(args.out_dir, f"summary_priority_heatmap_matrix_{args.mode}.tsv"), sep="\t")
        print(f"[WRITE] {os.path.join(args.out_dir, f'summary_priority_heatmap_matrix_{args.mode}.tsv')}")
        plot_heatmap(mat, title=f"Priority genes across cell types ({args.mode})", out_fp=os.path.join(args.out_dir, f"fig6_priority_heatmap_{args.mode}.png"))

    write_report(
        os.path.join(args.out_dir, f"summary_report_{args.mode}.md"),
        disease_df=disease_df,
        study_df=study_df,
        merged_metrics=merged_metrics,
        compare_df=compare_df,
        freq_df=freq_df,
    )

    print("\n[DONE] summary completed.")


if __name__ == "__main__":
    main()
