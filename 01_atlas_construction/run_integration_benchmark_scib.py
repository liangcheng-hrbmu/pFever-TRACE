#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Run integration benchmark for pediatric febrile PBMC atlas using scib-metrics.

This version is compatible with the scib-metrics API where BatchCorrection uses:
  bras, ilisi_knn, kbet_per_label, graph_connectivity, pcr_comparison

Main comparisons:
  Raw PCA
  scVI / scANVI / Harmony or any embedding stored in adata.obsm

This font42 version adds editable vector font settings for matplotlib outputs,
including the native scIB results table.
"""

import os
import argparse
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

import matplotlib as mpl
mpl.use("Agg")
mpl.rcParams.update({
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
    "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
    "font.family": "sans-serif",
    "text.usetex": False,
    "savefig.transparent": False,
    "figure.dpi": 300,
    "axes.unicode_minus": False
})

import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")


def ensure_dir(path):
    Path(path).mkdir(parents=True, exist_ok=True)


def sanitize_name(x):
    keep = []
    for ch in str(x):
        if ch.isalnum() or ch in ["_", "."]:
            keep.append(ch)
        else:
            keep.append("_")
    return "".join(keep)


def stratified_sample_obs(adata, batch_key, label_key, max_cells=80000, seed=5):
    if adata.n_obs <= max_cells:
        return adata.copy()

    rng = np.random.default_rng(seed)
    obs = adata.obs[[batch_key, label_key]].copy()
    obs[batch_key] = obs[batch_key].astype(str)
    obs[label_key] = obs[label_key].astype(str)

    frac = max_cells / adata.n_obs
    selected = []

    for _, sub in obs.groupby([batch_key, label_key], observed=True):
        n = int(np.ceil(sub.shape[0] * frac))
        n = max(20, n)
        n = min(n, sub.shape[0])
        idx = rng.choice(sub.index.to_numpy(), size=n, replace=False)
        selected.extend(idx.tolist())

    selected = pd.Index(selected).unique()

    if len(selected) > max_cells:
        selected = pd.Index(rng.choice(selected.to_numpy(), size=max_cells, replace=False))

    return adata[selected].copy()


def compute_raw_pca(adata, layer="counts", n_top_genes=3000, n_comps=50):
    print("[RAW PCA] computing raw log-normalized PCA")

    ad = adata.copy()

    if layer in ad.layers:
        ad.X = ad.layers[layer].copy()

    sc.pp.normalize_total(ad, target_sum=1e4)
    sc.pp.log1p(ad)

    try:
        sc.pp.highly_variable_genes(ad, n_top_genes=n_top_genes, flavor="seurat_v3")
    except Exception:
        sc.pp.highly_variable_genes(ad, n_top_genes=n_top_genes)

    if "highly_variable" in ad.var.columns and int(ad.var["highly_variable"].sum()) > 50:
        ad = ad[:, ad.var["highly_variable"].values].copy()

    sc.pp.scale(ad, max_value=10)
    sc.tl.pca(ad, n_comps=n_comps, svd_solver="arpack")

    adata.obsm["X_pca_unintegrated"] = ad.obsm["X_pca"]
    return adata


def attach_external_embedding(base, spec):
    """
    External embedding spec:
      Name::/path/to/file.h5ad::obsm_key

    Example:
      Harmony::/path/integrated_harmony_umap.h5ad::X_pca_harmony
    """
    parts = spec.split("::")
    if len(parts) != 3:
        raise ValueError(
            f"Bad --embedding format: {spec}\n"
            "Expected format: Name::/path/file.h5ad::obsm_key"
        )

    name, path, key = parts
    out_key = "X_bench_" + sanitize_name(name)

    print(f"[EMBEDDING] {name}: {path} :: {key}")
    ext = sc.read_h5ad(path)

    if key not in ext.obsm:
        raise KeyError(f"{key} not found in {path}. Available obsm keys: {list(ext.obsm.keys())}")

    common = base.obs_names.intersection(ext.obs_names)

    if len(common) == 0:
        raise ValueError(f"No overlapping cells between base h5ad and {path}")

    if len(common) < base.n_obs:
        print(f"[WARNING] {name} matched {len(common)} / {base.n_obs} cells. Subsetting base object.")
        base = base[common].copy()

    emb_df = pd.DataFrame(ext.obsm[key], index=ext.obs_names)
    base.obsm[out_key] = emb_df.loc[base.obs_names].values

    return base, out_key, name


def collect_base_embeddings(adata, requested_keys):
    """
    requested_keys:
      auto
      or comma-separated obsm keys, e.g. X_scANVI,X_scVI,X_pca_harmony
    """
    embedding_map = {}

    if requested_keys == "auto":
        candidates = [
            ("Raw PCA", "X_pca_unintegrated"),
            ("scANVI", "X_scANVI"),
            ("scANVI", "X_scanvi"),
            ("scVI", "X_scVI"),
            ("scVI", "X_scvi"),
            ("Harmony", "X_pca_harmony"),
            ("PCA", "X_pca"),
        ]

        used_keys = set()

        for name, key in candidates:
            if key in adata.obsm:
                out_key = "X_bench_" + sanitize_name(name)

                if out_key in used_keys:
                    continue

                adata.obsm[out_key] = adata.obsm[key]
                embedding_map[out_key] = name
                used_keys.add(out_key)

    else:
        for item in requested_keys.split(","):
            item = item.strip()
            if item == "":
                continue

            if item not in adata.obsm:
                raise KeyError(f"{item} not in base adata.obsm. Available: {list(adata.obsm.keys())}")

            out_key = "X_bench_" + sanitize_name(item)
            adata.obsm[out_key] = adata.obsm[item]
            embedding_map[out_key] = item

    return adata, embedding_map


def run_scib_benchmark(
    adata,
    batch_key,
    label_key,
    embedding_keys,
    pre_integrated_key,
    n_jobs,
    out_dir,
    disable_kbet=False,
):
    try:
        from scib_metrics.benchmark import Benchmarker, BioConservation, BatchCorrection
    except Exception as e:
        raise ImportError(
            "Cannot import scib_metrics. Install it with:\n"
            "  pip install scib-metrics\n"
            "or use a conda environment that already contains scib_metrics.\n"
            f"Original error: {e}"
        )

    bio_metrics = BioConservation(
        isolated_labels=False,
        nmi_ari_cluster_labels_leiden=True,
        nmi_ari_cluster_labels_kmeans=True,
        silhouette_label=True,
        clisi_knn=True,
    )

    batch_metrics = BatchCorrection(
        bras=True,
        ilisi_knn=True,
        kbet_per_label=(not disable_kbet),
        graph_connectivity=True,
        pcr_comparison=True,
    )

    print("[SCIB] metrics")
    print(f"  kBET per label: {not disable_kbet}")

    print("[SCIB] embeddings:")
    for k in embedding_keys:
        print("  ", k, adata.obsm[k].shape)

    kwargs = dict(
        adata=adata,
        batch_key=batch_key,
        label_key=label_key,
        embedding_obsm_keys=embedding_keys,
        bio_conservation_metrics=bio_metrics,
        batch_correction_metrics=batch_metrics,
        n_jobs=n_jobs,
        progress_bar=True,
    )

    if pre_integrated_key is not None and pre_integrated_key in adata.obsm:
        kwargs["pre_integrated_embedding_obsm_key"] = pre_integrated_key
        print("[SCIB] pre integrated embedding:", pre_integrated_key)

    print("[SCIB] benchmark started")
    bm = Benchmarker(**kwargs)
    bm.benchmark()

    print("[SCIB] collecting results")
    res_raw = bm.get_results(min_max_scale=False)
    res_scaled = bm.get_results(min_max_scale=True)

    raw_path = os.path.join(out_dir, "benchmark_results_raw.tsv")
    scaled_path = os.path.join(out_dir, "benchmark_results_scaled.tsv")

    res_raw.to_csv(raw_path, sep="\t")
    res_scaled.to_csv(scaled_path, sep="\t")

    print("[WRITE]", raw_path)
    print("[WRITE]", scaled_path)

    # Native scIB table plot.
    # This is the part that generates scib_results.svg in save_dir.
    # rc_context makes sure font settings are applied to this native plot.
    try:
        with mpl.rc_context({
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
           "font.family": "sans-serif",
            "text.usetex": False,
            "savefig.transparent": False,
            "figure.dpi": 300,
            "axes.unicode_minus": False
        }):
            bm.plot_results_table(min_max_scale=False, save_dir=out_dir)
    except Exception as e:
        print("[WARNING] Native scib plot failed:", e)

    return res_raw, res_scaled


def plot_simple_heatmap(df, out_prefix, title):
    numeric = df.copy()

    for col in numeric.columns:
        numeric[col] = pd.to_numeric(numeric[col], errors="coerce")

    numeric = numeric.dropna(axis=0, how="all").dropna(axis=1, how="all")

    if numeric.empty:
        print("[WARNING] No numeric values to plot for", title)
        return

    # Most scib result tables have metrics as rows and embeddings as columns.
    # Keep that orientation unless it is extremely tall.
    if numeric.shape[0] > numeric.shape[1] * 2:
        numeric = numeric.T

    fig_w = max(6, numeric.shape[1] * 0.7)
    fig_h = max(3.5, numeric.shape[0] * 0.42)

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    vmin = np.nanmin(numeric.values)
    vmax = np.nanmax(numeric.values)

    im = ax.imshow(numeric.values, aspect="auto", vmin=vmin, vmax=vmax)

    ax.set_xticks(np.arange(numeric.shape[1]))
    ax.set_xticklabels(numeric.columns.astype(str), rotation=60, ha="right", fontsize=8)

    ax.set_yticks(np.arange(numeric.shape[0]))
    ax.set_yticklabels(numeric.index.astype(str), fontsize=8)

    for i in range(numeric.shape[0]):
        for j in range(numeric.shape[1]):
            val = numeric.values[i, j]
            if np.isfinite(val):
                ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=6)

    ax.set_title(title, fontsize=11)

    cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cb.set_label("Score")

    fig.savefig(out_prefix + ".pdf", bbox_inches="tight")
    fig.savefig(out_prefix + ".svg", bbox_inches="tight")
    fig.savefig(out_prefix + ".png", bbox_inches="tight", dpi=300)
    plt.close(fig)

    print("[SAVE]", out_prefix + ".pdf")
    print("[SAVE]", out_prefix + ".svg")
    print("[SAVE]", out_prefix + ".png")


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--base_h5ad", required=True)
    parser.add_argument("--out_dir", required=True)

    parser.add_argument("--batch_key", default="dataset_id")
    parser.add_argument("--label_key", default="aifi_L1_paper")
    parser.add_argument("--counts_layer", default="counts")

    parser.add_argument("--max_cells", type=int, default=80000)
    parser.add_argument("--seed", type=int, default=5)

    parser.add_argument("--compute_raw_pca", action="store_true")
    parser.add_argument("--raw_pca_n_top_genes", type=int, default=3000)
    parser.add_argument("--raw_pca_n_comps", type=int, default=50)

    parser.add_argument(
        "--base_embedding_keys",
        default="auto",
        help="auto or comma-separated obsm keys in base_h5ad.",
    )

    parser.add_argument(
        "--embedding",
        action="append",
        default=[],
        help="External embedding spec: Name::/path/file.h5ad::obsm_key. Can be used multiple times.",
    )

    parser.add_argument("--n_jobs", type=int, default=16)

    parser.add_argument(
        "--disable_kbet",
        action="store_true",
        help="Disable kBET per label if benchmark is too slow or unstable.",
    )

    args = parser.parse_args()

    ensure_dir(args.out_dir)

    print("[READ]", args.base_h5ad)
    adata = sc.read_h5ad(args.base_h5ad)

    for col in [args.batch_key, args.label_key]:
        if col not in adata.obs.columns:
            raise KeyError(f"{col} not found in adata.obs. Available columns: {list(adata.obs.columns)}")

    adata.obs[args.batch_key] = adata.obs[args.batch_key].astype(str)
    adata.obs[args.label_key] = adata.obs[args.label_key].astype(str)

    print("[FILTER] remove missing batch or label")
    keep = (
        adata.obs[args.batch_key].notna()
        & adata.obs[args.label_key].notna()
        & (adata.obs[args.batch_key].astype(str) != "")
        & (adata.obs[args.label_key].astype(str) != "")
        & (adata.obs[args.batch_key].astype(str).str.lower() != "nan")
        & (adata.obs[args.label_key].astype(str).str.lower() != "nan")
    )

    adata = adata[keep].copy()

    print("[SAMPLE]")
    print("  before:", adata.n_obs)

    adata = stratified_sample_obs(
        adata,
        batch_key=args.batch_key,
        label_key=args.label_key,
        max_cells=args.max_cells,
        seed=args.seed,
    )

    print("  after:", adata.n_obs)

    if args.compute_raw_pca or "X_pca_unintegrated" not in adata.obsm:
        adata = compute_raw_pca(
            adata,
            layer=args.counts_layer,
            n_top_genes=args.raw_pca_n_top_genes,
            n_comps=args.raw_pca_n_comps,
        )

    print("[COLLECT BASE EMBEDDINGS]")
    adata, embedding_map = collect_base_embeddings(adata, args.base_embedding_keys)

    print("[ATTACH EXTERNAL EMBEDDINGS]")
    for spec in args.embedding:
        adata, key, name = attach_external_embedding(adata, spec)
        embedding_map[key] = name

    embedding_keys = []
    final_embedding_map = {}

    for key, name in embedding_map.items():
        if key not in adata.obsm:
            continue

        emb = np.asarray(adata.obsm[key])

        if emb.shape[0] != adata.n_obs:
            print(f"[WARNING] skip {key}: wrong n_obs {emb.shape[0]} vs {adata.n_obs}")
            continue

        if not np.isfinite(emb).all():
            print(f"[WARNING] skip {key}: contains NaN or Inf")
            continue

        if key not in embedding_keys:
            embedding_keys.append(key)
            final_embedding_map[key] = name

    if len(embedding_keys) < 2:
        raise ValueError(
            f"Need at least two embeddings. Found: {embedding_keys}\n"
            f"Available obsm keys: {list(adata.obsm.keys())}"
        )

    map_df = pd.DataFrame(
        [{"obsm_key": k, "method_name": v} for k, v in final_embedding_map.items()]
    )
    map_df.to_csv(os.path.join(args.out_dir, "benchmark_embedding_map.tsv"), sep="\t", index=False)

    meta_cols = [args.batch_key, args.label_key]
    extra_cols = [
        c for c in ["sample_id", "patient_id", "disease_Level_3", "sex", "binned_age"]
        if c in adata.obs.columns
    ]

    adata.obs[meta_cols + extra_cols].to_csv(
        os.path.join(args.out_dir, "benchmark_sampled_cell_metadata.tsv"),
        sep="\t",
    )

    sampled_h5ad = os.path.join(args.out_dir, "benchmark_input_sampled_embeddings.h5ad")
    adata.write(sampled_h5ad, compression="gzip")
    print("[WRITE]", sampled_h5ad)

    pre_key = None
    if "X_bench_Raw_PCA" in embedding_keys:
        pre_key = "X_bench_Raw_PCA"
    elif "X_bench_X_pca_unintegrated" in embedding_keys:
        pre_key = "X_bench_X_pca_unintegrated"

    res_raw, res_scaled = run_scib_benchmark(
        adata=adata,
        batch_key=args.batch_key,
        label_key=args.label_key,
        embedding_keys=embedding_keys,
        pre_integrated_key=pre_key,
        n_jobs=args.n_jobs,
        out_dir=args.out_dir,
        disable_kbet=args.disable_kbet,
    )

    plot_simple_heatmap(
        res_scaled,
        os.path.join(args.out_dir, "benchmark_results_scaled_heatmap"),
        "Integration benchmark scaled scores",
    )

    plot_simple_heatmap(
        res_raw,
        os.path.join(args.out_dir, "benchmark_results_raw_heatmap"),
        "Integration benchmark raw scores",
    )

    print("[DONE]")
    print("out_dir:", args.out_dir)


if __name__ == "__main__":
    main()
