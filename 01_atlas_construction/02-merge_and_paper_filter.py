#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os, glob, argparse, re
import numpy as np
import pandas as pd

def ensure_dir(p):
    os.makedirs(p, exist_ok=True)

def infer_blacklist_mask(var_names):
    """
    论文式 blacklist：mt/ribo/TCR/BCR/MHC/HB/platelet
    """
    genes = np.array([str(x) for x in var_names])

    pat_mt   = re.compile(r"^(MT-|mt-)")
    pat_ribo = re.compile(r"^(RPL|RPS)")
    pat_tcr  = re.compile(r"^TR[ABDG][CVJ]")   # TRAC/TRBC1...
    pat_bcr  = re.compile(r"^(IGH|IGK|IGL)")   # IGHG1/IGKC...
    pat_mhc  = re.compile(r"^(HLA-|HLA)")      # HLA-A/HLA-DRA...
    pat_hb   = re.compile(r"^HB[ABDEGQZ]")     # HBB/HBA1...

    platelet_set = set(["PPBP", "PF4", "SDPR", "GP9", "ITGA2B", "GP1BA", "GP1BB", "NRGN"])

    is_mt   = np.array([bool(pat_mt.match(g)) for g in genes])
    is_ribo = np.array([bool(pat_ribo.match(g)) for g in genes])
    is_tcr  = np.array([bool(pat_tcr.match(g)) for g in genes])
    is_bcr  = np.array([bool(pat_bcr.match(g)) for g in genes])
    is_mhc  = np.array([bool(pat_mhc.match(g)) for g in genes])
    is_hb   = np.array([bool(pat_hb.match(g)) for g in genes])
    is_plat = np.array([g in platelet_set for g in genes])

    return (is_mt | is_ribo | is_tcr | is_bcr | is_mhc | is_hb | is_plat)

def main():
    ap = argparse.ArgumentParser(
        description="Merge *.withDisease.h5ad and apply paper-like gene filtering only (no HVG/scVI/UMAP)."
    )

    ap.add_argument("--in_root", required=True,
                    help="processed_celltypist 根目录：包含 GSE*/xxx.withDisease.h5ad")
    ap.add_argument("--out_dir", required=True)

    ap.add_argument("--counts_layer", default="counts",
                    help="必须是原始 counts（推荐 layers['counts']）")
    ap.add_argument("--study_key", default="dataset_id",
                    help="用来做 per-study 表达过滤的 obs 列（默认 dataset_id）")

    ap.add_argument("--do_paper_filter", action="store_true",
                    help="开启论文式过滤：去 mt/ribo/TCR/BCR/MHC/HB/platelet + counts>0表达过滤")
    ap.add_argument("--strict_each_study_expressed", action="store_true",
                    help="严格：每个study都至少表达一次；不加则任意study表达一次即可")

    ap.add_argument("--out_h5ad", default="merged_paperFiltered.h5ad",
                    help="输出 h5ad 文件名（默认 merged_paperFiltered.h5ad）")

    args = ap.parse_args()

    ensure_dir(args.out_dir)
    list_dir = os.path.join(args.out_dir, "gene_lists")
    ensure_dir(list_dir)

    import anndata as ad
    import scipy.sparse as sp

    # 1) 收集文件
    pattern = os.path.join(args.in_root, "GSE*", "*.withDisease.h5ad")
    files = sorted(glob.glob(pattern))
    if not files:
        raise SystemExit(f"[FATAL] No *.withDisease.h5ad found under --in_root (pattern={pattern})")

    # 2) 读入 + 最小修复
    adatas, keys = [], []
    for f in files:
        gse = os.path.basename(os.path.dirname(f))
        a = ad.read_h5ad(f)

        if "dataset_id" not in a.obs.columns:
            a.obs["dataset_id"] = gse

        adatas.append(a)
        keys.append(gse)

    # 3) 合并（join="inner"）
    adata = ad.concat(adatas, join="inner", label="dataset_merge_key", keys=keys, index_unique="-")
    del adatas
    print("[INFO] merged:", adata.shape)

    # 4) counts layer 检查 + 稀疏化
    if (not hasattr(adata, "layers")) or (args.counts_layer not in adata.layers.keys()):
        raise SystemExit(
            f"[FATAL] counts layer '{args.counts_layer}' not found. "
            f"Please ensure layers['{args.counts_layer}'] exists."
        )
    if not sp.issparse(adata.layers[args.counts_layer]):
        adata.layers[args.counts_layer] = sp.csr_matrix(adata.layers[args.counts_layer])

    # 5) 论文式过滤（可选）
    removed = []
    kept = list(adata.var_names)

    if args.do_paper_filter:
        if args.study_key not in adata.obs.columns:
            raise ValueError(f"[FATAL] study_key '{args.study_key}' not found in obs.")

        # (A) blacklist
        bl = infer_blacklist_mask(adata.var_names)
        keep0 = ~bl

        # (B) counts>0 表达过滤（按 study）
        X = adata.layers[args.counts_layer]
        studies = adata.obs[args.study_key].astype(str).values
        uniq = np.unique(studies)

        expr_masks = []
        for s in uniq:
            idx = np.where(studies == s)[0]
            if idx.size == 0:
                continue
            Xs = X[idx, :]
            m = np.array(Xs.getnnz(axis=0)).ravel() > 0
            expr_masks.append(m)

        if len(expr_masks) == 0:
            raise SystemExit("[FATAL] No studies found for expression filtering.")

        expr_masks = np.vstack(expr_masks)

        if args.strict_each_study_expressed:
            keep1 = expr_masks.all(axis=0)
        else:
            keep1 = expr_masks.any(axis=0)

        keep = keep0 & keep1

        removed = adata.var_names[~keep].tolist()
        kept = adata.var_names[keep].tolist()

        adata = adata[:, keep].copy()
        if not sp.issparse(adata.layers[args.counts_layer]):
            adata.layers[args.counts_layer] = sp.csr_matrix(adata.layers[args.counts_layer])

        # 保存 gene lists
        with open(os.path.join(list_dir, "removed_by_paper_filter.txt"), "w", encoding="utf-8") as f:
            f.write("\n".join(removed) + "\n")
        with open(os.path.join(list_dir, "kept_after_paper_filter.txt"), "w", encoding="utf-8") as f:
            f.write("\n".join(kept) + "\n")

        print(f"[INFO] paper_filter enabled: kept={len(kept)} removed={len(removed)}")
    else:
        print("[INFO] paper_filter disabled: keep all genes after concat (join=inner).")

    # 6) 写出主对象
    out_path = os.path.join(args.out_dir, args.out_h5ad)
    adata.write_h5ad(out_path, compression="gzip")
    print("[WRITE]", out_path)

    print("\n[SUMMARY]")
    print("n_obs =", adata.n_obs)
    print("n_vars =", adata.n_vars)
    print("counts_layer =", args.counts_layer)
    print("study_key =", args.study_key)
    print("paper_filter =", bool(args.do_paper_filter))
    if args.do_paper_filter:
        print("strict_each_study_expressed =", bool(args.strict_each_study_expressed))
        print("gene_lists dir =", list_dir)

if __name__ == "__main__":
    main()
