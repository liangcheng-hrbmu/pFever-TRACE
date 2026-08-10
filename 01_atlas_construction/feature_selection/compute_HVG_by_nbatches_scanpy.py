#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import argparse
import numpy as np
import pandas as pd

def ensure_dir(p: str):
    os.makedirs(p, exist_ok=True)

def pick_hvgs_by_nbatches(var: pd.DataFrame, min_studies: int, min_hvgs: int):
    """
    Scanpy highly_variable_genes(batch_key=...) 会产生：
      - highly_variable_nbatches: 每个基因在多少个 batch 中被判为 HVG
      - dispersions_norm / means 等统计
    策略：
      1) 先取 nbatches >= min_studies
      2) 若不足 min_hvgs：逐步降低阈值；仍不足则按 (nbatches desc, dispersions_norm desc) 补足到 min_hvgs
    """
    if "highly_variable_nbatches" not in var.columns:
        raise ValueError("Missing var['highly_variable_nbatches']. Did you run highly_variable_genes with batch_key?")

    nb = var["highly_variable_nbatches"].fillna(0).astype(int)

    hv = var.index[nb >= min_studies].tolist()

    used_thr = min_studies
    if len(hv) < min_hvgs:
        # 逐步放宽阈值
        for t in range(min_studies - 1, 0, -1):
            hv_t = var.index[nb >= t].tolist()
            if len(hv_t) >= min_hvgs:
                hv = hv_t
                used_thr = t
                break

        # 仍不足：按排序补足
        if len(hv) < min_hvgs:
            sort_cols = []
            if "dispersions_norm" in var.columns:
                sort_cols = ["highly_variable_nbatches", "dispersions_norm"]
                df = var.copy()
                df["highly_variable_nbatches"] = nb
                df["dispersions_norm"] = df["dispersions_norm"].fillna(-np.inf)
                df = df.sort_values(sort_cols, ascending=[False, False])
                hv = df.index[:min_hvgs].tolist()
            else:
                df = var.copy()
                df["highly_variable_nbatches"] = nb
                df = df.sort_values(["highly_variable_nbatches"], ascending=False)
                hv = df.index[:min_hvgs].tolist()
            used_thr = None  # 表示是“补足模式”

    return hv, used_thr

def main():
    ap = argparse.ArgumentParser(
        description="Select HVGs using scanpy highly_variable_genes with batch_key (nbatches strategy), like 04_HVGcomputation_inputObjects.ipynb"
    )
    ap.add_argument("--in_h5ad", required=True, help="Input h5ad (already merged + paperL1L2), e.g. merged_paperFiltered.paperL1L2.h5ad")
    ap.add_argument("--out_dir", required=True, help="Output directory")

    ap.add_argument("--counts_layer", default="counts", help="Raw counts layer name")
    ap.add_argument("--batch_key", default="libraryID", help="Batch key used in HVG: e.g. libraryID (paper) or dataset_id (your data)")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--max_cells", type=int, default=0, help="Optional subsample cells for speed (0 = use all)")

    # paper params
    ap.add_argument("--min_mean", type=float, default=0.01)
    ap.add_argument("--max_mean", type=float, default=4.0)
    ap.add_argument("--min_disp", type=float, default=0.3)

    ap.add_argument("--min_studies", type=int, default=5, help="Keep genes HVG in >=min_studies batches (nbatches)")
    ap.add_argument("--min_hvgs", type=int, default=3000, help="Ensure at least this many HVGs")

    # outputs
    ap.add_argument("--out_hvg_list", default=None, help="Output HVG list filename. Default: HVG_min{min_studies}_studies.txt")
    ap.add_argument("--out_stats_tsv", default="HVG_scanpy_nbatches_stats.tsv", help="Write var stats table (tsv)")
    ap.add_argument("--write_h5ad_with_hvgflag", action="store_true",
                    help="If set, write a copy of h5ad with var['highly_variable'] updated (no subsetting).")

    args = ap.parse_args()

    ensure_dir(args.out_dir)
    list_dir = os.path.join(args.out_dir, "gene_lists")
    ensure_dir(list_dir)

    import scanpy as sc
    import scipy.sparse as sp

    np.random.seed(args.seed)

    print("[READ]", args.in_h5ad)
    adata = sc.read_h5ad(args.in_h5ad)

    if args.batch_key not in adata.obs.columns:
        raise SystemExit(f"[FATAL] obs['{args.batch_key}'] not found. "
                         f"Your file has keys: {list(adata.obs.columns)[:30]} ... "
                         f"Please set --batch_key to the correct column (e.g. dataset_id / sample_id / libraryID).")

    if args.counts_layer not in adata.layers:
        raise SystemExit(f"[FATAL] layers['{args.counts_layer}'] not found in input h5ad.")

    # optional subsample cells
    if args.max_cells and args.max_cells > 0 and adata.n_obs > args.max_cells:
        keep = np.random.choice(adata.n_obs, size=args.max_cells, replace=False)
        adata = adata[keep, :].copy()
        print(f"[INFO] Subsampled to {adata.n_obs} cells for HVG selection.")

    # ensure sparse counts
    if not sp.issparse(adata.layers[args.counts_layer]):
        adata.layers[args.counts_layer] = sp.csr_matrix(adata.layers[args.counts_layer])

    # ---- Method (as in ipynb text) ----
    # 1) normalize_total(target_sum=1e4) on counts layer
    # 2) log1p on counts layer
    # 3) set X = log1p(counts) so scanpy HVG uses X
    # 4) scanpy.pp.highly_variable_genes(flavor="seurat", min_mean/min_disp/max_mean, batch_key=...)
    print("[STEP] normalize_total + log1p on layer:", args.counts_layer)
    sc.pp.normalize_total(adata, target_sum=1e4, layer=args.counts_layer)
    sc.pp.log1p(adata, layer=args.counts_layer)
    adata.X = adata.layers[args.counts_layer].copy()

    print("[STEP] scanpy.pp.highly_variable_genes with batch_key =", args.batch_key)
    sc.pp.highly_variable_genes(
        adata,
        flavor="seurat",
        min_mean=args.min_mean,
        max_mean=args.max_mean,
        min_disp=args.min_disp,
        batch_key=args.batch_key,
        subset=False
    )

    # pick HVGs by nbatches strategy
    hvgs, used_thr = pick_hvgs_by_nbatches(adata.var, args.min_studies, args.min_hvgs)

    # write stats table
    stats_cols = [c for c in [
        "highly_variable", "highly_variable_nbatches", "highly_variable_intersection",
        "means", "dispersions", "dispersions_norm"
    ] if c in adata.var.columns]
    df_stats = adata.var[stats_cols].copy()
    df_stats.index.name = "gene"
    out_stats = os.path.join(list_dir, args.out_stats_tsv)
    df_stats.to_csv(out_stats, sep="\t")

    # write HVG list
    out_list_name = args.out_hvg_list
    if out_list_name is None:
        out_list_name = f"HVG_min{args.min_studies}_studies.txt"
    out_list = os.path.join(list_dir, out_list_name)
    with open(out_list, "w", encoding="utf-8") as f:
        f.write("\n".join(hvgs) + "\n")

    # optionally write h5ad copy with var flags updated
    if args.write_h5ad_with_hvgflag:
        adata.var["highly_variable_final"] = adata.var_names.isin(set(hvgs))
        out_h5ad = os.path.join(args.out_dir, "with_HVGflag.h5ad")
        adata.write_h5ad(out_h5ad, compression="gzip")
        print("[WRITE]", out_h5ad)

    print("\n[DONE]")
    print("batch_key =", args.batch_key)
    print("min_mean/min_disp/max_mean =", args.min_mean, args.min_disp, args.max_mean)
    if used_thr is None:
        print(f"HVG selected by补足模式到 >= {args.min_hvgs}")
    else:
        print(f"HVG threshold used: nbatches >= {used_thr}")
    print("HVG n =", len(hvgs))
    print("HVG list:", out_list)
    print("Stats tsv:", out_stats)

if __name__ == "__main__":
    main()
