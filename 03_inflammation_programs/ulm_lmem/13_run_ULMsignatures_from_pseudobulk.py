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
    """
    Column-wise z-score (gene-wise across pseudobulks).
    Safe for zero-variance columns.
    """
    mu = df.mean(axis=0)
    sd = df.std(axis=0, ddof=0)
    sd = sd.replace(0, np.nan)
    out = (df - mu) / sd
    out = out.replace([np.inf, -np.inf], np.nan).fillna(0.0)
    return out


def main():
    ap = argparse.ArgumentParser(
        description="Run ULM signatures directly from corrected pseudobulk tables with filtering/log1p/z-score."
    )
    ap.add_argument("--pseudobulk_expr", required=True)
    ap.add_argument("--pseudobulk_meta", required=True)
    ap.add_argument("--factor_gene_weights", required=True)
    ap.add_argument("--factor_best_match", required=True)
    ap.add_argument("--out_dir", required=True)

    ap.add_argument("--sample_key", default="sample_id")
    ap.add_argument("--celltype_key", default="celltype")
    ap.add_argument("--study_key", default="dataset_id")
    ap.add_argument("--disease_key", default="disease_Level_3")
    ap.add_argument("--sex_key", default="sex")
    ap.add_argument("--age_key", default="binned_age")
    ap.add_argument("--min_n_genes", type=int, default=5)
    ap.add_argument("--n_cells_key", default="n_cells")
    ap.add_argument("--min_n_cells", type=int, default=10,
                    help="Minimum n_cells required for a pseudobulk row to be retained.")

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
        args.sample_key, args.celltype_key, args.study_key,
        args.disease_key, args.sex_key, args.age_key
    ]
    miss = [x for x in needed_meta if x not in pb_meta.columns]
    if miss:
        raise ValueError(f"pseudobulk_meta missing columns: {miss}")

    if args.n_cells_key not in pb_meta.columns:
        raise ValueError(f"pseudobulk_meta missing n_cells column: {args.n_cells_key}")

    if "factor" not in factor_best.columns or "Level1" not in factor_best.columns:
        raise ValueError("factor_best_match.tsv must contain columns: factor, Level1")

    common_rows = pb_expr.index.intersection(pb_meta.index)
    if len(common_rows) == 0:
        raise RuntimeError("No overlapping rows between pseudobulk_expr and pseudobulk_meta")
    pb_expr = pb_expr.loc[common_rows].copy()
    pb_meta = pb_meta.loc[common_rows].copy()

    # 1) filter low-cell-number pseudobulks
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

    celltypes = sorted(adata_pb.obs[args.celltype_key].astype(str).unique().tolist())

    estimate_all = pd.DataFrame()
    pvals_all = pd.DataFrame()
    qc_rows = []

    for ct in celltypes:
        sub_idx = adata_pb.obs.index[adata_pb.obs[args.celltype_key].astype(str) == ct]
        sub = adata_pb[sub_idx].copy()

        keep_f = [
            f for f in gene_weights.index
            if str(factor2l1.get(f, "")) == str(ct) or str(factor2l1.get(f, "")) == "global"
        ]
        if len(keep_f) == 0:
            print(f"[SKIP] {ct}: no matched factors")
            qc_rows.append({"celltype": ct, "n_pseudobulks": sub.n_obs, "n_factors": 0, "status": "no_matched_factors"})
            continue

        gw_sub = gene_weights.loc[keep_f].copy()
        net_df = build_net_from_gene_weights(gw_sub)
        if net_df.shape[0] == 0:
            print(f"[SKIP] {ct}: empty network")
            qc_rows.append({"celltype": ct, "n_pseudobulks": sub.n_obs, "n_factors": len(keep_f), "status": "empty_network"})
            continue

        # 2) log1p on corrected pseudobulk
        sub_expr = sub.to_df()
        sub_expr = np.log1p(sub_expr)

        # 3) z-score by gene within cell type
        sub_expr = zscore_by_gene(sub_expr)

        print(f"[RUN] {ct}: pseudobulks={sub.n_obs}, factors={len(keep_f)}, genes={sub.n_vars}")

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
        qc_rows.append({"celltype": ct, "n_pseudobulks": sub.n_obs, "n_factors": len(keep_f), "status": "ok"})

    estimate_all = estimate_all.loc[~estimate_all.index.duplicated()].copy()
    pvals_all = pvals_all.loc[~pvals_all.index.duplicated()].copy()

    est_path = os.path.join(args.out_dir, "DecoupleR_ulmestimates_L1_Corr.csv")
    pval_path = os.path.join(args.out_dir, "DecoupleR_ulmpvals_L1_Corr.csv")
    meta_path = os.path.join(args.out_dir, "Pseudobulk_metadata_L1_Corr.tsv")
    h5ad_path = os.path.join(args.out_dir, "DecoupleR_ulmestimates_L1_Corr.h5ad")
    qc_path = os.path.join(args.out_dir, "ULM_L1_qc.tsv")

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
    print("non_na_estimate_ratio:", float(estimate_all.notna().sum().sum()) / float(estimate_all.size))


if __name__ == "__main__":
    main()
