#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os, argparse
import numpy as np
import pandas as pd

def ensure_dir(p): os.makedirs(p, exist_ok=True)

def hvgs_rank_by_study_frequency(
    adata,
    study_key="dataset_id",
    layer_counts="counts",
    min_mean=0.01,
    max_mean=4.0,
    min_disp=0.3,
    flavor="seurat",
):
    """
    参考 paper / 你 notebook 的思路：
    - 每个 study 独立 normalize_total + log1p
    - 每个 study 独立 highly_variable_genes(Seurat)
    - 统计每个基因在多少个 study 被标记为 HVG
    返回：hvg_count (pd.Series: gene -> n_studies_as_HVG)
    """
    import scanpy as sc
    import scipy.sparse as sp

    if study_key not in adata.obs.columns:
        raise ValueError(f"study_key '{study_key}' not found in adata.obs")
    if layer_counts not in adata.layers:
        raise ValueError(f"layer_counts '{layer_counts}' not found in adata.layers")

    tmp = adata.copy()
    if not sp.issparse(tmp.layers[layer_counts]):
        tmp.layers[layer_counts] = sp.csr_matrix(tmp.layers[layer_counts])

    # 统一做 log-normalized（仅用于 HVG 判定）
    sc.pp.normalize_total(tmp, target_sum=1e4, layer=layer_counts)
    sc.pp.log1p(tmp, layer=layer_counts)
    tmp.X = tmp.layers[layer_counts].copy()

    studies = tmp.obs[study_key].astype(str).values
    uniq = np.unique(studies)

    hvg_count = pd.Series(0, index=tmp.var_names, dtype=int)
    for s in uniq:
        idx = np.where(studies == s)[0]
        if idx.size == 0:
            continue
        sub = tmp[idx, :].copy()
        sc.pp.highly_variable_genes(
            sub,
            flavor=flavor,
            min_mean=min_mean,
            max_mean=max_mean,
            min_disp=min_disp,
            subset=False,
        )
        hv = sub.var["highly_variable"].reindex(tmp.var_names).fillna(False).values
        hvg_count += hv.astype(int)

    return hvg_count

def pick_hvgs_by_min_studies(hvg_count: pd.Series, min_hvgs=3000, start_min_studies=6):
    """
    paper 里的做法：找到“最低研究数量阈值”，让 HVG 列表大小 >= min_hvgs
    - 从 start_min_studies 往下放宽到 1
    - 若仍不足，则按频次排序补足到 min_hvgs
    返回：hvgs(list), used_min_studies(int)
    """
    used = start_min_studies
    hvgs = hvg_count.index[hvg_count.values >= used].tolist()

    if len(hvgs) >= min_hvgs:
        return hvgs, used

    for t in range(start_min_studies - 1, 0, -1):
        hvgs_t = hvg_count.index[hvg_count.values >= t].tolist()
        if len(hvgs_t) >= min_hvgs:
            return hvgs_t, t

    # 仍不足：按出现次数排序补足
    hvgs = hvg_count.sort_values(ascending=False).index[:min_hvgs].tolist()
    return hvgs, 1

def load_palette_csv(path):
    if (path is None) or (not os.path.exists(path)):
        return None
    df = pd.read_csv(path)
    if df.shape[1] < 2:
        return None
    label_candidates = [c for c in df.columns if any(k in c.lower() for k in ["type","cell","label","name"])]
    color_candidates = [c for c in df.columns if any(k in c.lower() for k in ["color","hex"])]
    label_col = label_candidates[0] if label_candidates else df.columns[0]
    color_col = color_candidates[0] if color_candidates else df.columns[1]
    pal = {}
    for _, r in df.iterrows():
        k = str(r[label_col])
        v = str(r[color_col]).strip()
        if v and (not v.startswith("#")) and len(v) == 6:
            v = "#" + v
        if v.startswith("#") and len(v) >= 7:
            pal[k] = v[:7]
    return pal if pal else None

def distinct_colors(n):
    import matplotlib.cm as cm
    import matplotlib.colors as mcolors
    cols = []
    base = cm.get_cmap("tab20")
    for i in range(min(n, base.N)):
        cols.append(mcolors.to_hex(base(i)))
    if len(cols) < n:
        k = n - len(cols)
        for i in range(k):
            h = i / max(1, k)
            cols.append(mcolors.to_hex(mcolors.hsv_to_rgb((h, 0.70, 0.92))))
    return cols[:n]

def palette_for_obs_categories(adata, col, pal_dict=None):
    if col not in adata.obs.columns:
        return None
    adata.obs[col] = adata.obs[col].astype("category")
    cats = list(adata.obs[col].cat.categories)
    auto = distinct_colors(len(cats))
    if pal_dict is None:
        return auto
    return [pal_dict.get(str(c), auto[i]) for i, c in enumerate(cats)]

def main():
    ap = argparse.ArgumentParser(
        description="Run scVI integration (HVG subset) and restore X_scVI/X_umap/leiden back to full-gene main object. Input already contains paper L1/L2."
    )

    ap.add_argument("--in_h5ad", required=True,
                    help="输入：merged_paperFiltered.paperL1L2.h5ad（已包含 aifi_L1_paper / aifi_L2_paper）")
    ap.add_argument("--out_dir", required=True)

    ap.add_argument("--counts_layer", default="counts")
    ap.add_argument("--batch_key", default="dataset_id")
    ap.add_argument("--study_key", default="dataset_id")

    # HVG 参数（按 paper / notebook 逻辑）
    ap.add_argument("--min_hvgs", type=int, default=3000, help="目标 HVG 数")
    ap.add_argument("--start_min_studies", type=int, default=6, help="从多少个 study 起步往下放宽（paper里常见>=6）")
    ap.add_argument("--min_mean", type=float, default=0.01)
    ap.add_argument("--max_mean", type=float, default=4.0)
    ap.add_argument("--min_disp", type=float, default=0.3)

    # scVI 参数
    ap.add_argument("--n_latent", type=int, default=30)
    ap.add_argument("--max_epochs", type=int, default=200)
    ap.add_argument("--use_gpu", action="store_true")
    ap.add_argument("--seed", type=int, default=0)

    # neighbors/UMAP/Leiden
    ap.add_argument("--n_neighbors", type=int, default=20)     # 你贴的 paper 默认 k=20
    ap.add_argument("--min_dist", type=float, default=0.5)
    ap.add_argument("--leiden_res", type=float, default=0.1)   # 你贴的 paper 示例 r=0.1

    # 画图（可关）
    ap.add_argument("--no_plot", action="store_true")
    ap.add_argument("--fig_size", type=float, default=6.0)
    ap.add_argument("--dpi", type=int, default=200)
    ap.add_argument("--pub_png_dpi", type=int, default=600)
    ap.add_argument("--pub_png_bg", default="white", choices=["white","transparent"])
    ap.add_argument("--pub_png_no_axes", action="store_true")
    ap.add_argument("--palette_L1", default=None)
    ap.add_argument("--palette_L2", default=None)

    # 你已有的列名（第二步已经生成好了）
    ap.add_argument("--paper_L1_col", default="aifi_L1_paper")
    ap.add_argument("--paper_L2_col", default="aifi_L2_paper")

    ap.add_argument("--out_main_h5ad", default="mainObject_fullGenes_scvi_umap.h5ad")

    args = ap.parse_args()

    ensure_dir(args.out_dir)
    list_dir = os.path.join(args.out_dir, "gene_lists"); ensure_dir(list_dir)
    qc_dir   = os.path.join(args.out_dir, "qc"); ensure_dir(qc_dir)
    fig_dir  = os.path.join(args.out_dir, "figs"); ensure_dir(fig_dir)
    pub_dir  = os.path.join(args.out_dir, "pub_png"); ensure_dir(pub_dir)

    import scanpy as sc
    import scipy.sparse as sp
    import scvi

    print("[READ]", args.in_h5ad)
    adata = sc.read_h5ad(args.in_h5ad)

    # 检查必需列
    for c in [args.paper_L1_col, args.paper_L2_col]:
        if c not in adata.obs.columns:
            raise SystemExit(f"[FATAL] missing obs['{c}'] in input. You said step2 already made it.")

    if args.counts_layer not in adata.layers:
        raise SystemExit(f"[FATAL] missing layers['{args.counts_layer}'] in input.")
    if not sp.issparse(adata.layers[args.counts_layer]):
        adata.layers[args.counts_layer] = sp.csr_matrix(adata.layers[args.counts_layer])

    # ===== 1) HVG（按 study 频次）=====
    print("[HVG] computing per-study HVGs and counting frequency ...")
    hvg_count = hvgs_rank_by_study_frequency(
        adata,
        study_key=args.study_key,
        layer_counts=args.counts_layer,
        min_mean=args.min_mean,
        max_mean=args.max_mean,
        min_disp=args.min_disp,
        flavor="seurat",
    )

    hvgs, used_min_studies = pick_hvgs_by_min_studies(
        hvg_count,
        min_hvgs=args.min_hvgs,
        start_min_studies=args.start_min_studies,
    )

    out_hvg = os.path.join(list_dir, "paperFiltered_HVGs.txt")
    with open(out_hvg, "w", encoding="utf-8") as f:
        f.write("\n".join(hvgs) + "\n")

    hvg_count.sort_values(ascending=False).to_csv(os.path.join(qc_dir, "hvg_n_studies.csv"))
    print(f"[HVG] n={len(hvgs)} used_min_studies={used_min_studies} -> {out_hvg}")

    # ===== 2) 训练对象 = HVG 子集（不落盘）=====
    adata_train = adata[:, hvgs].copy()
    if not sp.issparse(adata_train.layers[args.counts_layer]):
        adata_train.layers[args.counts_layer] = sp.csr_matrix(adata_train.layers[args.counts_layer])

    # ===== 3) scVI =====
    scvi.settings.seed = args.seed
    if args.batch_key not in adata_train.obs.columns:
        raise SystemExit(f"[FATAL] batch_key '{args.batch_key}' not found in obs.")

    scvi.model.SCVI.setup_anndata(
        adata_train,
        layer=args.counts_layer,
        batch_key=args.batch_key
    )
    model = scvi.model.SCVI(adata_train, n_latent=args.n_latent)

    train_kwargs = dict(max_epochs=args.max_epochs, early_stopping=True)
    if args.use_gpu:
        try:
            import torch
            if torch.cuda.is_available():
                train_kwargs.update(dict(accelerator="gpu", devices=1))
                print("[scVI] Using GPU.")
            else:
                print("[scVI] --use_gpu set but CUDA not available -> CPU.")
        except Exception:
            print("[scVI] torch check failed -> CPU.")

    print("[scVI] training ...")
    model.train(**train_kwargs)

    adata_train.obsm["X_scVI"] = model.get_latent_representation()

    # ===== 4) neighbors/UMAP/Leiden =====
    sc.pp.neighbors(adata_train, use_rep="X_scVI", n_neighbors=args.n_neighbors)
    sc.tl.umap(adata_train, min_dist=args.min_dist)
    sc.tl.leiden(adata_train, resolution=args.leiden_res, key_added="leiden_scvi")

    # ===== 5) restore 回 full-gene main object =====
    adata.obsm["X_scVI"] = adata_train.obsm["X_scVI"]
    adata.obsm["X_umap"] = adata_train.obsm["X_umap"]
    adata.obsp["connectivities"] = adata_train.obsp["connectivities"]
    adata.obsp["distances"]      = adata_train.obsp["distances"]
    adata.uns["neighbors"]       = adata_train.uns["neighbors"]
    adata.obs["leiden_scvi"]     = adata_train.obs["leiden_scvi"].astype(str).values

    # ===== 6) plot（可选）=====
    if not args.no_plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        pal_L1 = load_palette_csv(args.palette_L1)
        pal_L2 = load_palette_csv(args.palette_L2)

        def _format_pub_axes(ax):
            if args.pub_png_no_axes:
                ax.set_xlabel(""); ax.set_ylabel("")
                ax.set_xticks([]); ax.set_yticks([])
                for spine in ax.spines.values():
                    spine.set_visible(False)

        def save_umap(color, out_stem, pal_dict=None, title=None, legend=True):
            if color not in adata.obs.columns:
                print("[SKIP] missing obs:", color); return
            adata.obs[color] = adata.obs[color].astype("category")
            ncat = adata.obs[color].cat.categories.size
            legend_loc = "right margin" if legend else "none"
            if ncat > 25:
                legend_loc = "none"
            palette = palette_for_obs_categories(adata, color, pal_dict=pal_dict)

            fig = sc.pl.umap(
                adata,
                color=color,
                show=False,
                return_fig=True,
                title=title if title else color,
                palette=palette,
                legend_loc=legend_loc
            )
            fig.set_size_inches(args.fig_size, args.fig_size)

            out1 = os.path.join(fig_dir, f"{out_stem}.png")
            fig.savefig(out1, dpi=args.dpi, bbox_inches="tight", facecolor="white")
            print("[WRITE]", out1)

            ax = fig.axes[0] if fig.axes else None
            if ax is not None:
                _format_pub_axes(ax)

            face = "white" if args.pub_png_bg == "white" else "none"
            transparent = (args.pub_png_bg == "transparent")
            out2 = os.path.join(pub_dir, f"{out_stem}.pub.png")
            fig.savefig(out2, dpi=args.pub_png_dpi, bbox_inches="tight",
                        facecolor=face, transparent=transparent)
            plt.close(fig)
            print("[WRITE]", out2)

        save_umap(args.batch_key, f"umap_by_{args.batch_key}", title=f"UMAP (scVI) by {args.batch_key}")
        if "dataset_id" in adata.obs.columns:
            save_umap("dataset_id", "umap_by_dataset_id", title="UMAP (scVI) by dataset_id")
        save_umap(args.paper_L1_col, "umap_by_paper_Level1", pal_dict=pal_L1, title="Cell types (Level 1)")
        save_umap(args.paper_L2_col, "umap_by_paper_Level2", pal_dict=pal_L2, title="Cell types (Level 2)")
        save_umap("leiden_scvi", "umap_by_leiden_scvi", title="UMAP (scVI) by leiden_scvi")

    # ===== 7) save main object =====
    out_main = os.path.join(args.out_dir, args.out_main_h5ad)
    adata.write_h5ad(out_main, compression="gzip")
    print("[DONE] saved:", out_main)

    print("\n[SUMMARY]")
    print("Input:", args.in_h5ad)
    print("Output main:", out_main)
    print("HVG list:", out_hvg)
    print("HVG used_min_studies:", used_min_studies)
    print("Plots:", "OFF" if args.no_plot else (fig_dir, pub_dir))

if __name__ == "__main__":
    main()
