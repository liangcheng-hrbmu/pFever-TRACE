#!/usr/bin/env python3
import os, glob, re, traceback
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import scrublet as scr
import celltypist




import re
from matplotlib import colors as mcolors

_HEX_RE = re.compile(r"^#(?:[0-9a-fA-F]{6}|[0-9a-fA-F]{8})$")

def _is_color(x):
    if x is None:
        return False
    s = str(x).strip()
    if _HEX_RE.match(s):
        return True
    # 支持 "r,g,b"（0-1 或 0-255）
    if re.match(r"^\s*\d*\.?\d+\s*,\s*\d*\.?\d+\s*,\s*\d*\.?\d+\s*$", s):
        return True
    try:
        mcolors.to_rgba(s)
        return True
    except Exception:
        return False

def _norm_color(x):
    """返回 matplotlib 可用的颜色：'#RRGGBB' 或 (r,g,b) 或 fallback 灰色"""
    s = str(x).strip()
    if _HEX_RE.match(s):
        return s
    if re.match(r"^\s*\d*\.?\d+\s*,\s*\d*\.?\d+\s*,\s*\d*\.?\d+\s*$", s):
        parts = [float(p) for p in re.split(r"\s*,\s*", s)]
        if max(parts) > 1.0:
            parts = [p/255.0 for p in parts]
        return tuple(parts[:3])
    try:
        mcolors.to_rgba(s)
        return s
    except Exception:
        return "#808080"












# =========================
# 配置区：你只需要改这里
# =========================
IN_ROOT   = "data/processed"              # 输入：每个GSE目录下的h5ad
OUT_ROOT  = "data/processed_celltypist"   # 输出：注释+图
MODEL_DIR = "data/models"                 # 你的pkl + colors csv

SAMPLE_COL = "sample_id"     # 你obs里区分每个上机库的列名
DATASET_COL = "dataset_id"   # obs里数据集列名（一般你已经有）
PATIENT_COL = "patient_id"   # 可选
GROUP_COLS_FOR_PLOTS = ["disease_Level_3", "Fever", "timepoint", "Age group", "sex"]  # 有就画，没有跳过

# 是否跑双细胞（建议True；sorted CD4 数据可以设 False 或只标记不滤）
RUN_DOUBLET = True
FILTER_DOUBLET = True        # True: 过滤 predicted_doublet；False: 只打标签不删除
MAX_CELLS_PER_DATASET_FOR_SCRUBLET = 250000  # 超大数据集时，Scrublet会慢；>这个数就只做标记不做（或你可改策略）

# CellTypist 预处理（推荐这样，简单稳）
MIN_CELLS_PER_GENE = 3
NORMALIZE_TARGET_SUM = 1e4

# UMAP 图（可选；若数据大，会抽样）
PLOT_UMAP = True
UMAP_MAX_CELLS = 50000

CONF_CUTOFF = 0.5  # 低置信度阈值（用于汇总表）

# =========================
# 工具函数
# =========================

def find_one(patterns, root):
    """在root下按多个pattern找文件，返回第一个匹配"""
    for pat in patterns:
        hits = sorted(glob.glob(os.path.join(root, pat)))
        if hits:
            return hits[0]
    return None

def load_colors_csv(csv_path):
    """
    自动识别：
      - label列：字符串列里 unique 多的那列
      - color列：能解析成颜色的比例最高的那列
      - order列：若存在 order/index/rank 则用于排序
    如果找不到合法颜色列：自动用 tab20 生成颜色
    """
    df = pd.read_csv(csv_path)

    # 1) order列（可选）
    lower = {c: c.lower() for c in df.columns}
    order_col = None
    for c in df.columns:
        if any(k in lower[c] for k in ["order", "index", "rank"]):
            order_col = c
            break
    if order_col is not None:
        df = df.sort_values(order_col)

    # 2) label列：优先找包含 type/label/name/cell 的列，否则选字符串列里 unique 多的
    label_col = None
    for c in df.columns:
        if any(k in lower[c] for k in ["type", "label", "name", "cell"]):
            label_col = c
            break
    if label_col is None:
        obj_cols = [c for c in df.columns if df[c].dtype == object]
        if not obj_cols:
            label_col = df.columns[0]
        else:
            label_col = max(obj_cols, key=lambda c: df[c].astype(str).nunique())

    # 3) color列：在所有列中找“合法颜色比例最高”的列
    best_color_col, best_score = None, -1
    for c in df.columns:
        vals = df[c].dropna().astype(str).head(300).tolist()
        if len(vals) == 0:
            continue
        score = sum(_is_color(v) for v in vals) / len(vals)
        if score > best_score:
            best_score, best_color_col = score, c

    labels = df[label_col].astype(str).tolist()

    # 若颜色列几乎都不是颜色（比如都是 1/2/3），则不用它，改用 tab20 生成
    if best_score < 0.2:
        cmap = plt.get_cmap("tab20")
        colors = [mcolors.to_hex(cmap(i % 20)) for i in range(len(labels))]
        color_map = dict(zip(labels, colors))
        return labels, color_map

    colors_raw = df[best_color_col].astype(str).tolist()
    colors = [_norm_color(x) for x in colors_raw]
    color_map = dict(zip(labels, colors))
    return labels, color_map

def ensure_counts_layer(adata):
    """保证 layers['counts'] 存在，并且是稀疏矩阵"""
    if "counts" in adata.layers:
        Xc = adata.layers["counts"]
    else:
        Xc = adata.X
        adata.layers["counts"] = Xc.copy()

    # 统一为CSR稀疏，避免dense触发奇怪行为/内存爆
    if not sp.issparse(adata.layers["counts"]):
        adata.layers["counts"] = sp.csr_matrix(adata.layers["counts"])
    else:
        adata.layers["counts"] = adata.layers["counts"].tocsr()

def preprocess_for_celltypist(adata):
    """filter genes + normalize_total + log1p；使用 counts layer 作为输入"""
    adata.X = adata.layers["counts"].copy()
    sc.pp.filter_genes(adata, min_cells=MIN_CELLS_PER_GENE)
    sc.pp.normalize_total(adata, target_sum=NORMALIZE_TARGET_SUM)
    sc.pp.log1p(adata)

def expected_doublet_rate(n_cells):
    # 10x经验：约0.8%/1000 cells；这里做个上下限夹逼
    r = 0.008 * (n_cells / 1000.0)
    return float(np.clip(r, 0.005, 0.20))

def run_scrublet_per_sample(adata):
    """按 sample_id 逐个跑 scrublet，写入 obs：doublet_score/predicted_doublet"""
    if SAMPLE_COL not in adata.obs.columns:
        adata.obs[SAMPLE_COL] = "sample1"

    adata.obs["doublet_score"] = np.nan
    adata.obs["predicted_doublet"] = False

    counts = adata.layers["counts"]

    # 超大数据集策略：只做标记/或跳过（这里选择跳过）
    if adata.n_obs > MAX_CELLS_PER_DATASET_FOR_SCRUBLET:
        print(f"  [WARN] cells={adata.n_obs} > {MAX_CELLS_PER_DATASET_FOR_SCRUBLET}, skip scrublet.")
        return adata

    for sid, idx in adata.obs.groupby(SAMPLE_COL).indices.items():
        idx = np.asarray(idx)
        X = counts[idx, :]

        rate = expected_doublet_rate(len(idx))
        s = scr.Scrublet(X, expected_doublet_rate=rate)
        scores, preds = s.scrub_doublets()

        adata.obs.iloc[idx, adata.obs.columns.get_loc("doublet_score")] = scores
        adata.obs.iloc[idx, adata.obs.columns.get_loc("predicted_doublet")] = preds

    return adata

def run_celltypist_levels(adata, model_L1, model_L2, model_L3):
    import pandas as pd
    import numpy as np

    def _pick_labels(pred):
        # pred 可能是 Series，也可能是 DataFrame（含 majority voting）
        if isinstance(pred, pd.DataFrame):
            for col in ["majority_voting", "majority vote", "majority_vote", "MV", "predicted_labels"]:
                if col in pred.columns:
                    return pred[col]
            return pred.iloc[:, -1]  # 兜底：最后一列
        return pred

    def _get_prob_matrix(res):
        # 不同版本可能叫 probability_matrix / prob_mat / decision_matrix / decision_mat
        for attr in ["probability_matrix", "prob_mat", "probability"]:
            if hasattr(res, attr):
                return getattr(res, attr), "prob"
        for attr in ["decision_matrix", "decision_mat"]:
            if hasattr(res, attr):
                return getattr(res, attr), "decision"
        return None, None

    for lvl, model_path in [("L1", model_L1), ("L2", model_L2), ("L3", model_L3)]:
        res = celltypist.annotate(adata, model=model_path, majority_voting=True)

        # labels
        labels = _pick_labels(res.predicted_labels)
        labels = pd.Series(np.asarray(labels).ravel(), index=adata.obs_names).astype(str)
        adata.obs[f"aifi_{lvl}"] = labels

        # confidence
        mat, kind = _get_prob_matrix(res)
        if mat is None:
            # 实在没有就给 NA
            adata.obs[f"aifi_{lvl}_conf"] = np.nan
        else:
            if isinstance(mat, pd.DataFrame):
                M = mat
            else:
                M = pd.DataFrame(mat, index=adata.obs_names)

            if kind == "decision":
                # decision -> probability: sigmoid
                M = 1.0 / (1.0 + np.exp(-M))

            conf = np.asarray(M.max(axis=1)).ravel()
            adata.obs[f"aifi_{lvl}_conf"] = conf

    return adata



def set_category_order_and_colors(adata, col, ordered_labels, color_map):
    """把obs[col]设置为有序category，并把颜色写入 uns"""
    if col not in adata.obs.columns:
        return None
    # 只保留本数据存在的类别，按官方顺序排列
    present = pd.unique(adata.obs[col].astype(str))
    cats = [x for x in ordered_labels if x in set(present)]
    if len(cats) == 0:
        # 如果没交集，就按出现顺序
        cats = list(present)

    adata.obs[col] = pd.Categorical(adata.obs[col].astype(str), categories=cats, ordered=True)
    # 写入scanpy用的颜色列表
    palette = [color_map.get(x, "#808080") for x in cats]
    adata.uns[f"{col}_colors"] = palette
    return cats, palette

def plot_qc(adata, out_png):
    """QC：n_counts/n_genes/percent_mt/doublet_score 分布"""
    # n_counts/n_genes
    if "n_counts" not in adata.obs.columns:
        X = adata.layers["counts"]
        adata.obs["n_counts"] = np.asarray(X.sum(axis=1)).ravel()
    if "n_genes" not in adata.obs.columns:
        X = adata.layers["counts"]
        adata.obs["n_genes"] = np.asarray((X > 0).sum(axis=1)).ravel()

    # percent.mt（人MT-；鼠mt-）
    if "percent_mt" not in adata.obs.columns:
        mt_mask = adata.var_names.str.startswith("MT-")
        if mt_mask.sum() == 0:
            mt_mask = adata.var_names.str.startswith("mt-")
        if mt_mask.sum() > 0:
            X = adata.layers["counts"]
            mt = np.asarray(X[:, mt_mask].sum(axis=1)).ravel()
            total = np.asarray(X.sum(axis=1)).ravel()
            adata.obs["percent_mt"] = (mt / np.maximum(total, 1)) * 100.0
        else:
            adata.obs["percent_mt"] = 0.0

    cols = ["n_counts", "n_genes", "percent_mt"]
    if "doublet_score" in adata.obs.columns:
        cols.append("doublet_score")

    fig, axes = plt.subplots(1, len(cols), figsize=(4*len(cols), 4))
    if len(cols) == 1:
        axes = [axes]
    for ax, c in zip(axes, cols):
        x = adata.obs[c].astype(float).values
        ax.hist(x, bins=60)
        ax.set_title(c)
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close(fig)

def plot_celltype_bar_by_sample(adata, label_col, out_png, max_samples=60):
    """按sample_id画堆叠比例条形图（太多样本就截取前max_samples）"""
    if SAMPLE_COL not in adata.obs.columns or label_col not in adata.obs.columns:
        return
    df = adata.obs[[SAMPLE_COL, label_col]].copy()
    df[SAMPLE_COL] = df[SAMPLE_COL].astype(str)
    df[label_col] = df[label_col].astype(str)

    tab = (df.groupby([SAMPLE_COL, label_col]).size()
             .reset_index(name="n"))
    # 转比例
    den = tab.groupby(SAMPLE_COL)["n"].sum().reset_index(name="tot")
    tab = tab.merge(den, on=SAMPLE_COL, how="left")
    tab["prop"] = tab["n"]/tab["tot"]

    # 样本太多就截取细胞数最多的前N个
    top_samples = (den.sort_values("tot", ascending=False)[SAMPLE_COL].head(max_samples).tolist())
    tab = tab[tab[SAMPLE_COL].isin(top_samples)].copy()

    pivot = tab.pivot(index=SAMPLE_COL, columns=label_col, values="prop").fillna(0.0)
    # 类别顺序用category顺序
    if isinstance(adata.obs[label_col].dtype, pd.CategoricalDtype):
        ordered = list(adata.obs[label_col].cat.categories)
        pivot = pivot.reindex(columns=[c for c in ordered if c in pivot.columns])

    colors = None
    if f"{label_col}_colors" in adata.uns and isinstance(adata.obs[label_col].dtype, pd.CategoricalDtype):
        cats = list(adata.obs[label_col].cat.categories)
        colmap = dict(zip(cats, adata.uns[f"{label_col}_colors"]))
        colors = [colmap.get(c, "#808080") for c in pivot.columns]

    fig, ax = plt.subplots(figsize=(max(10, 0.25*len(pivot.index)), 6))
    bottom = np.zeros(len(pivot.index))
    x = np.arange(len(pivot.index))
    for i, c in enumerate(pivot.columns):
        ax.bar(x, pivot[c].values, bottom=bottom, label=c, color=(colors[i] if colors else None), linewidth=0)
        bottom += pivot[c].values

    ax.set_xticks(x)
    ax.set_xticklabels(pivot.index.tolist(), rotation=90, fontsize=8)
    ax.set_ylabel("Proportion")
    ax.set_title(f"{label_col} composition by {SAMPLE_COL} (top {len(pivot.index)})")
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=7, frameon=False, ncol=1)
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close(fig)

def plot_umap(adata, color_col, out_png, title=None):
    """UMAP（抽样计算，避免太慢；不显式传palette，直接用 ad.uns['xxx_colors']）"""
    if not PLOT_UMAP or color_col not in adata.obs.columns:
        return

    # 抽样
    if adata.n_obs > UMAP_MAX_CELLS:
        idx = np.random.choice(adata.n_obs, UMAP_MAX_CELLS, replace=False)
        ad = adata[idx, :].copy()
    else:
        ad = adata.copy()

    # 保留并修正颜色：scanpy 会默认读 ad.uns[f"{color_col}_colors"]
    if f"{color_col}_colors" in ad.uns:
        ad.uns[f"{color_col}_colors"] = list(ad.uns[f"{color_col}_colors"])  # 关键：转成list避免numpy truth-value问题

    # 轻量高变+PCA+neighbors+umap（尽量从counts开始）
    ensure_counts_layer(ad)
    ad.X = ad.layers["counts"].copy()

    sc.pp.normalize_total(ad, target_sum=NORMALIZE_TARGET_SUM)
    sc.pp.log1p(ad)

    sc.pp.highly_variable_genes(ad, n_top_genes=3000, flavor="seurat_v3")
    ad = ad[:, ad.var["highly_variable"]].copy()

    sc.tl.pca(ad, n_comps=50)
    sc.pp.neighbors(ad, n_neighbors=15, n_pcs=30)
    sc.tl.umap(ad)

    fig = sc.pl.umap(ad, color=color_col, show=False, return_fig=True, title=title)
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    plt.close(fig)

# =========================
# 主流程
# =========================

def main():
    os.makedirs(OUT_ROOT, exist_ok=True)

    # 自动识别模型文件（按你截图的命名）
    model_L1 = find_one(["ref_pbmc*_AIFI_L1_*.pkl", "*AIFI*L1*.pkl"], MODEL_DIR)
    model_L2 = find_one(["ref_pbmc*_AIFI_L2_*.pkl", "*AIFI*L2*.pkl"], MODEL_DIR)
    model_L3 = find_one(["ref_pbmc*_AIFI_L3_*.pkl", "*AIFI*L3*.pkl"], MODEL_DIR)
    if not (model_L1 and model_L2 and model_L3):
        raise FileNotFoundError("Cannot find AIFI L1/L2/L3 model pkl in MODEL_DIR")

    # 自动识别颜色表
    colors_L1 = find_one(["AIFI_L1*_colors.csv", "*L1*colors*.csv"], MODEL_DIR)
    colors_L2 = find_one(["AIFI_L2*_colors.csv", "*L2*colors*.csv"], MODEL_DIR)
    colors_L3 = find_one(["AIFI_L3*_colors.csv", "*L3*colors*.csv"], MODEL_DIR)

    L1_order, L1_cmap = (None, None)
    L2_order, L2_cmap = (None, None)
    L3_order, L3_cmap = (None, None)
    if colors_L1: L1_order, L1_cmap = load_colors_csv(colors_L1)
    if colors_L2: L2_order, L2_cmap = load_colors_csv(colors_L2)
    if colors_L3: L3_order, L3_cmap = load_colors_csv(colors_L3)

    print("[MODEL] L1:", model_L1)
    print("[MODEL] L2:", model_L2)
    print("[MODEL] L3:", model_L3)
    print("[COLORS] L1:", colors_L1)
    print("[COLORS] L2:", colors_L2)
    print("[COLORS] L3:", colors_L3)

    files = sorted(glob.glob(os.path.join(IN_ROOT, "GSE*", "*.h5ad")))
    if len(files) == 0:
        raise RuntimeError(f"No h5ad found: {IN_ROOT}/GSE*/*.h5ad")

    summary = []
    errors = []

    for f in files:
        gse = os.path.basename(os.path.dirname(f))
        base = os.path.splitext(os.path.basename(f))[0]
        out_dir = os.path.join(OUT_ROOT, gse)
        os.makedirs(out_dir, exist_ok=True)

        out_h5ad = os.path.join(out_dir, f"{base}.AIFI_L1L2L3.h5ad")
        qc_png   = os.path.join(out_dir, f"{base}.QC.png")
        bar_png  = os.path.join(out_dir, f"{base}.L2_bar_by_sample.png")
        umap_png = os.path.join(out_dir, f"{base}.UMAP_L2.png")

        try:
            adata = sc.read_h5ad(f)
            n0 = adata.n_obs

            # 确保 dataset_id 存在
            if DATASET_COL not in adata.obs.columns:
                adata.obs[DATASET_COL] = gse

            ensure_counts_layer(adata)

            # 1) 双细胞
            if RUN_DOUBLET:
                adata = run_scrublet_per_sample(adata)
                dbl_rate = float(np.mean(adata.obs["predicted_doublet"])) if "predicted_doublet" in adata.obs else np.nan
                if FILTER_DOUBLET and "predicted_doublet" in adata.obs:
                    adata = adata[~adata.obs["predicted_doublet"]].copy()
            else:
                dbl_rate = np.nan

            n1 = adata.n_obs

            # 2) CellTypist 三层
            preprocess_for_celltypist(adata)
            adata = run_celltypist_levels(adata, model_L1, model_L2, model_L3)

            # 3) 设置类别顺序+颜色
            if L1_order and L1_cmap:
                set_category_order_and_colors(adata, "aifi_L1", L1_order, L1_cmap)
            if L2_order and L2_cmap:
                set_category_order_and_colors(adata, "aifi_L2", L2_order, L2_cmap)
            if L3_order and L3_cmap:
                set_category_order_and_colors(adata, "aifi_L3", L3_order, L3_cmap)

            # 4) 保存 h5ad
            adata.write_h5ad(out_h5ad)

            # 5) 画图：QC + L2组成 + UMAP(L2)
            plot_qc(adata, qc_png)
            plot_celltype_bar_by_sample(adata, "aifi_L2", bar_png, max_samples=60)
            plot_umap(adata, "aifi_L2", umap_png, title=f"{gse} | AIFI L2")

            # 6) 汇总
            row = {
                "dataset": gse,
                "file_in": f,
                "file_out": out_h5ad,
                "cells_before": int(n0),
                "cells_after_doublet_filter": int(n1),
                "doublet_rate_pred": dbl_rate
            }
            for lvl in ["L1", "L2", "L3"]:
                conf = adata.obs[f"aifi_{lvl}_conf"].to_numpy()
                row[f"{lvl}_lowconf_frac(<{CONF_CUTOFF})"] = float(np.mean(conf < CONF_CUTOFF))
                row[f"{lvl}_n_labels"] = int(pd.Series(adata.obs[f"aifi_{lvl}"]).nunique())
            summary.append(row)

            print(f"[OK] {gse} | {n0} -> {n1} cells | out={out_h5ad}")

        except Exception as e:
            errors.append({"dataset": gse, "file_in": f, "error": repr(e)})
            print(f"[FAIL] {gse} | {f}")
            traceback.print_exc()

    # 写 summary / errors
    sum_df = pd.DataFrame(summary).sort_values("dataset")
    sum_csv = os.path.join(OUT_ROOT, "AIFI_L1L2L3_batch_summary.csv")
    sum_df.to_csv(sum_csv, index=False)

    if errors:
        err_df = pd.DataFrame(errors)
        err_csv = os.path.join(OUT_ROOT, "AIFI_L1L2L3_errors.csv")
        err_df.to_csv(err_csv, index=False)
        print("[DONE] with errors:", err_csv)

    print("[DONE] summary:", sum_csv)

if __name__ == "__main__":
    main()
