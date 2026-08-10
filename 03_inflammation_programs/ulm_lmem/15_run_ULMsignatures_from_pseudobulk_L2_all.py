#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import argparse
import numpy as np
import pandas as pd
import anndata as ad
import decoupler as dc


def ensure_dir(p):
    os.makedirs(p, exist_ok=True)


def build_net_from_gene_weights(gene_weights_df: pd.DataFrame) -> pd.DataFrame:
    net = (
        gene_weights_df
        .stack()
        .reset_index()
        .rename(columns={"level_0": "source", "level_1": "target", 0: "weight"})
    )
    net["weight"] = pd.to_numeric(net["weight"], errors="coerce")
    net = net.dropna(subset=["weight"])
    net = net.loc[net["weight"] != 0].copy()
    return net


def zscore_by_gene(df: pd.DataFrame) -> pd.DataFrame:
    mu = df.mean(axis=0)
    sd = df.std(axis=0, ddof=0)
    sd = sd.replace(0, np.nan)
    out = (df - mu) / sd
    out = out.replace([np.inf, -np.inf], np.nan).fillna(0.0)
    return out


def main():
    ap = argparse.ArgumentParser(
        description="Run ULM signatures on corrected L2 pseudobulk tables across all L1 lineages."
    )
    ap.add_argument("--pseudobulk_expr", required=True)
    ap.add_argument("--pseudobulk_meta", required=True)
    ap.add_argument("--factor_gene_weights", required=True)
    ap.add_argument("--factor_best_match", required=True)
    ap.add_argument("--out_dir", required=True)

    ap.add_argument("--sample_key", default="sample_id")
    ap.add_argument("--celltype_key", default="celltype_L2")
    ap.add_argument("--lineage_key", default="celltype_L1_majority")
    ap.add_argument("--study_key", default="dataset_id")
    ap.add_argument("--disease_key", default="disease_Level_3")
    ap.add_argument("--sex_key", default="sex")
    ap.add_argument("--age_key", default="binned_age")
    ap.add_argument("--n_cells_key", default="n_cells")
    ap.add_argument("--min_n_cells", type=int, default=10)
    ap.add_argument("--min_n_genes", type=int, default=5)
    ap.add_argument("--min_pb_per_lineage", type=int, default=6,
                    help="Minimum number of retained L2 pseudobulks required to run one L1 lineage.")

    args = ap.parse_args()
    ensure_dir(args.out_dir)

    print("[READ]", args.pseudobulk_expr)
    pb_expr = pd.read_csv(args.pseudobulk_expr, sep="\t", index_col=0, compression="gzip")

    print("[READ]", args.pseudobulk_meta)
    pb_meta = pd.read_csv(args.pseudobulk_meta, sep="\t", index_col=0)

    print("[READ]", args.factor_gene_weights)
    gene_weights = pd.read_csv(args.factor_gene_weights, sep="\t", index_col=0, compression="gzip")

    print("[READ]", args.factor_best_match)
    factor_best = pd.read_csv(args.factor_best_match, sep="\t")

    needed_meta = [
        args.sample_key, args.celltype_key, args.lineage_key, args.study_key,
        args.disease_key, args.sex_key, args.age_key, args.n_cells_key
    ]
    miss = [x for x in needed_meta if x not in pb_meta.columns]
    if miss:
        raise ValueError(f"pseudobulk_meta missing columns: {miss}")

    if "factor" not in factor_best.columns or "Level1" not in factor_best.columns:
        raise ValueError("factor_best_match.tsv must contain columns: factor, Level1")

    common_rows = pb_expr.index.intersection(pb_meta.index)
    if len(common_rows) == 0:
        raise RuntimeError("No overlapping rows between pseudobulk_expr and pseudobulk_meta")
    pb_expr = pb_expr.loc[common_rows].copy()
    pb_meta = pb_meta.loc[common_rows].copy()

    pb_meta[args.n_cells_key] = pd.to_numeric(pb_meta[args.n_cells_key], errors="coerce")
    keep_pb = pb_meta[args.n_cells_key] >= args.min_n_cells
    print(f"[FILTER] keep pseudobulks with {args.n_cells_key} >= {args.min_n_cells}: {int(keep_pb.sum())}/{int(len(keep_pb))}")
    pb_meta = pb_meta.loc[keep_pb].copy()
    pb_expr = pb_expr.loc[pb_meta.index].copy()

    selected_factors = [f for f in factor_best["factor"].tolist() if f in gene_weights.index]
    if len(selected_factors) == 0:
        raise RuntimeError("No overlapping factors between factor_best_match and factor_gene_weights")
    gene_weights = gene_weights.loc[selected_factors].copy()

    common_genes = [g for g in gene_weights.columns if g in pb_expr.columns]
    if len(common_genes) == 0:
        raise RuntimeError("No overlapping genes between factor_gene_weights and pseudobulk_expr")
    gene_weights = gene_weights.loc[:, common_genes].copy()
    pb_expr = pb_expr.loc[:, common_genes].copy()

    factor2l1 = factor_best.set_index("factor")["Level1"].to_dict()

    adata_pb = ad.AnnData(
        X=pb_expr.values,
        obs=pb_meta.copy(),
        var=pd.DataFrame(index=pb_expr.columns)
    )
    adata_pb.obs_names = pb_expr.index
    adata_pb.var_names = pb_expr.columns

    lineages = sorted(adata_pb.obs[args.lineage_key].astype(str).unique().tolist())

    estimate_all = pd.DataFrame()
    pvals_all = pd.DataFrame()
    qc_rows = []

    for lineage in lineages:
        sub_idx = adata_pb.obs.index[adata_pb.obs[args.lineage_key].astype(str) == lineage]
        sub = adata_pb[sub_idx].copy()

        if sub.n_obs < args.min_pb_per_lineage:
            print(f"[SKIP] {lineage}: too few retained L2 pseudobulks ({sub.n_obs} < {args.min_pb_per_lineage})")
            qc_rows.append({
                "lineage": lineage,
                "n_pseudobulks": sub.n_obs,
                "n_l2_celltypes": int(sub.obs[args.celltype_key].astype(str).nunique()),
                "n_factors": 0,
                "status": "too_few_pseudobulks"
            })
            continue

        keep_f = [
            f for f in gene_weights.index
            if str(factor2l1.get(f, "")) == str(lineage) or str(factor2l1.get(f, "")) == "global"
        ]
        if len(keep_f) == 0:
            print(f"[SKIP] {lineage}: no matched factors")
            qc_rows.append({
                "lineage": lineage,
                "n_pseudobulks": sub.n_obs,
                "n_l2_celltypes": int(sub.obs[args.celltype_key].astype(str).nunique()),
                "n_factors": 0,
                "status": "no_matched_factors"
            })
            continue

        gw_sub = gene_weights.loc[keep_f].copy()
        net_df = build_net_from_gene_weights(gw_sub)
        if net_df.shape[0] == 0:
            print(f"[SKIP] {lineage}: empty network")
            qc_rows.append({
                "lineage": lineage,
                "n_pseudobulks": sub.n_obs,
                "n_l2_celltypes": int(sub.obs[args.celltype_key].astype(str).nunique()),
                "n_factors": len(keep_f),
                "status": "empty_network"
            })
            continue

        sub_expr = sub.to_df()
        sub_expr = np.log1p(sub_expr)
        sub_expr = zscore_by_gene(sub_expr)

        print(f"[RUN] {lineage}: pseudobulks={sub.n_obs}, l2_celltypes={sub.obs[args.celltype_key].astype(str).nunique()}, factors={len(keep_f)}, genes={sub.n_vars}")

        est, pv = dc.run_ulm(
            mat=sub_expr,
            net=net_df,
            source="source",
            target="target",
            weight="weight",
            verbose=False,
            use_raw=False,
            min_n=args.min_n_genes
        )

        est = pd.DataFrame(est)
        pv = pd.DataFrame(pv)
        est.index = sub.obs_names
        pv.index = sub.obs_names

        estimate_all = pd.concat([estimate_all, est], axis=0)
        pvals_all = pd.concat([pvals_all, pv], axis=0)
        qc_rows.append({
            "lineage": lineage,
            "n_pseudobulks": sub.n_obs,
            "n_l2_celltypes": int(sub.obs[args.celltype_key].astype(str).nunique()),
            "n_factors": len(keep_f),
            "status": "ok"
        })

    estimate_all = estimate_all.loc[~estimate_all.index.duplicated()].copy()
    pvals_all = pvals_all.loc[~pvals_all.index.duplicated()].copy()

    est_path = os.path.join(args.out_dir, "DecoupleR_ulmestimates_L2_Corr.csv")
    pval_path = os.path.join(args.out_dir, "DecoupleR_ulmpvals_L2_Corr.csv")
    meta_path = os.path.join(args.out_dir, "Pseudobulk_metadata_L2_Corr.tsv")
    h5ad_path = os.path.join(args.out_dir, "DecoupleR_ulmestimates_L2_Corr.h5ad")
    qc_path = os.path.join(args.out_dir, "ULM_L2_qc.tsv")

    estimate_all.to_csv(est_path)
    pvals_all.to_csv(pval_path)
    pb_meta.loc[estimate_all.index].to_csv(meta_path, sep="\t")
    pd.DataFrame(qc_rows).to_csv(qc_path, sep="\t", index=False)

    adata_out = ad.AnnData(
        X=estimate_all.fillna(0).values,
        obs=pb_meta.loc[estimate_all.index].copy(),
        var=pd.DataFrame(index=estimate_all.columns)
    )
    adata_out.obsm["ulm_estimate"] = estimate_all.copy()
    adata_out.obsm["ulm_pvals"] = pvals_all.reindex(index=estimate_all.index, columns=estimate_all.columns)
    adata_out.write_h5ad(h5ad_path)

    print("\n[DONE]")
    print("estimate:", est_path, estimate_all.shape)
    print("pvals:", pval_path, pvals_all.shape)
    print("meta:", meta_path, pb_meta.loc[estimate_all.index].shape)
    print("h5ad:", h5ad_path, adata_out.shape)
    print("qc:", qc_path, len(qc_rows))
    if estimate_all.size > 0:
        print("non_na_estimate_ratio:", float(estimate_all.notna().sum().sum()) / float(estimate_all.size))


if __name__ == "__main__":
    main()
