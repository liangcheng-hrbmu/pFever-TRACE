#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os, sys, math, copy
import argparse
import types

import numpy as np
import pandas as pd
import scanpy as sc

# ---- FIX: force tqdm.notebook -> tqdm (CLI-safe) ----
import tqdm as _tqdm
_nb = types.ModuleType("tqdm.notebook")
_nb.tqdm = _tqdm.tqdm
_nb.trange = _tqdm.trange
sys.modules["tqdm.notebook"] = _nb
# ----------------------------------------------------

import SEACells



def cumulative_explained_variance(adata, expVarThr=None, n_pc=None):
    cs = np.cumsum(adata.uns["pca"]["variance"] / np.sum(adata.uns["pca"]["variance"]))
    res = {"n_pcs": None, "expVar": None}

    if expVarThr is not None:
        nC = np.where(cs >= expVarThr)[0][0] + 1
        res["n_pcs"] = int(nC)
        res["expVar"] = float(cs[nC - 1])

    if n_pc is not None:
        res["n_pcs"] = int(n_pc)
        res["expVar"] = float(cs[n_pc - 1])

    return res


def get_soft_assignments_modified(model, top_k=5):
    """
    返回：
      soft_labels: (n_cells x top_k) 每个细胞 top-k 的 SEACell 标签（字符串）
      weights:     (n_cells x top_k) 对应权重
    """
    A = copy.deepcopy(model.A_.T)  # shape: (n_cells, n_seacells)
    labels = []
    weights = []

    for _ in range(top_k):
        l = A.argmax(1)  # best seacell index
        seacell_name = np.array([f"SEACell-{i}" for i in l], dtype=object)
        labels.append(seacell_name)
        weights.append(A[np.arange(A.shape[0]), l])
        A[np.arange(A.shape[0]), l] = -1

    weights = np.vstack(weights).T
    labels = np.vstack(labels).T

    soft_labels = pd.DataFrame(labels, index=model.ad.obs_names)
    return soft_labels, weights


def ensure_dir(p):
    os.makedirs(p, exist_ok=True)


def main():
    ap = argparse.ArgumentParser(description="Run SEACells per-sample h5ad and output seacell assignments (pkl).")
    ap.add_argument("--h5ad", required=True, help="Input per-sample .h5ad")
    ap.add_argument("--sample_key", default="sample_id", help="obs column for sample id (default: sample_id)")
    ap.add_argument("--library_key", default="libraryID", help="obs column for library/batch inside sample (default: libraryID)")
    ap.add_argument("--SCxSEACell", type=int, required=True, help="cells per SEACell (e.g. 50; smaller => more SEACells)")
    ap.add_argument("--outDir", required=True, help="Output directory for *.seacell.pkl")
    ap.add_argument("--minExpVar", type=float, default=0.90, help="PCA cumulative explained variance threshold (default 0.90)")
    ap.add_argument("--n_neighbors", type=int, default=10, help="SEACells n_neighbors (default 10)")
    ap.add_argument("--n_waypoint_eigs", type=int, default=7, help="SEACells n_waypoint_eigs (default 7)")
    ap.add_argument("--top_k", type=int, default=5, help="top-k soft assignments to compute (default 5)")
    ap.add_argument("--min_iter", type=int, default=10, help="SEACells min_iter (default 10)")
    ap.add_argument("--max_iter", type=int, default=10000, help="SEACells max_iter (default 10000)")
    ap.add_argument("--seed", type=int, default=42, help="Random seed")
    args = ap.parse_args()

    ensure_dir(args.outDir)

    np.random.seed(args.seed)

    # ---------------------------
    # Load
    # ---------------------------
    adata = sc.read_h5ad(args.h5ad)
    if args.sample_key not in adata.obs.columns:
        raise ValueError(f"[ERROR] obs missing sample_key='{args.sample_key}'. Available: {list(adata.obs.columns)[:30]} ...")

    sample_vals = adata.obs[args.sample_key].astype(str).str.strip().unique()
    if sample_vals.shape[0] != 1:
        raise ValueError(f"[ERROR] input h5ad must contain exactly 1 sample. Got {sample_vals.shape[0]} samples in {args.sample_key}.")
    sample_id = sample_vals[0]
    print(f"[INFO] Running SEACells for sample = {sample_id}")
    print(f"[INFO] n_cells = {adata.n_obs}")

    # ---------------------------
    # Preprocess (as your original)
    # ---------------------------
    adata.X = adata.X.astype(float)
    sc.pp.normalize_per_cell(adata)
    sc.pp.log1p(adata)

    # HVG (Seurat flavor, top 3000)
    batch_key = None
    if args.library_key in adata.obs.columns:
        if adata.obs[args.library_key].astype(str).nunique() > 1:
            batch_key = args.library_key

    sc.pp.highly_variable_genes(
        adata=adata,
        batch_key=batch_key,
        flavor="seurat",
        n_top_genes=3000
    )

    # PCA
    sc.tl.pca(adata, svd_solver="arpack", n_comps=50, use_highly_variable=True)
    n_pcs = cumulative_explained_variance(adata, expVarThr=args.minExpVar)["n_pcs"]
    if n_pcs is None or n_pcs < 2:
        n_pcs = 20  # fallback
        print(f"[WARN] PCA expVar threshold failed; fallback n_pcs={n_pcs}")

    adata.obsm["X_pca_npcs"] = adata.obsm["X_pca"][:, :n_pcs]
    build_kernel_on = "X_pca_npcs"

    # ---------------------------
    # Run SEACells
    # ---------------------------
    n_SEACells = int(math.ceil(adata.n_obs / args.SCxSEACell))
    print(f"[INFO] SCxSEACell={args.SCxSEACell} => n_SEACells={n_SEACells}")

    model = SEACells.core.SEACells(
        adata,
        n_neighbors=args.n_neighbors,
        build_kernel_on=build_kernel_on,
        n_SEACells=n_SEACells,
        n_waypoint_eigs=args.n_waypoint_eigs,
        convergence_epsilon=1e-5
    )

    model.construct_kernel_matrix()
    model.initialize_archetypes()
    model.fit(min_iter=args.min_iter, max_iter=args.max_iter)

    # ---------------------------
    # Extract assignments (FIXED)
    # ---------------------------
    soft_labels, weights = get_soft_assignments_modified(model, top_k=args.top_k)
    # 取 top-1 label & weight
    top1_label = soft_labels.iloc[:, 0].astype(str).values
    top1_weight = weights[:, 0].astype(float)

    # prefix with sample_id to make global unique
    seacell_full = np.array([f"{sample_id}_{x}" for x in top1_label], dtype=object)

    # non-trivial assignments count (A_ is (n_seacells x n_cells); transpose to (n_cells x n_seacells))
    non_trivial = (model.A_.T > 0.1).sum(axis=1).A1 if hasattr((model.A_.T > 0.1).sum(axis=1), "A1") else np.array((model.A_.T > 0.1).sum(axis=1)).ravel()

    adata.obs["SEACell"] = seacell_full
    adata.obs["SEACell_weight"] = top1_weight
    adata.obs["SEACell_non_trivial_assig"] = non_trivial.astype(int)

    # ---------------------------
    # QC metrics
    # ---------------------------
    compactness = SEACells.evaluate.compactness(adata, build_kernel_on)
    separation = SEACells.evaluate.separation(adata, build_kernel_on, nth_nbr=1)

    # compactness/separation are DataFrames keyed by 'SEACell'
    df = adata.obs[["SEACell", "SEACell_weight", "SEACell_non_trivial_assig"]].copy()
    df.index = adata.obs_names
    df = df.merge(compactness, on="SEACell", how="left")
    df = df.merge(separation, on="SEACell", how="left")
    df.index = adata.obs_names

    # ---------------------------
    # Save
    # ---------------------------
    out_pkl = os.path.join(args.outDir, f"{sample_id}.seacell.pkl")
    df.to_pickle(out_pkl)
    print(f"[WRITE] {out_pkl}")
    print("[DONE]")


if __name__ == "__main__":
    main()
