#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as ssp
import anndata as ad
import scvi


def ensure_dir(p):
    os.makedirs(p, exist_ok=True)


def clean_str(x):
    if pd.isna(x):
        return np.nan
    x = str(x).strip()
    if x in ("", "nan", "None", "NA"):
        return np.nan
    return x


def main():
    ap = argparse.ArgumentParser(
        description="Generate per-Level1 SHAP input h5ad files from scANVI corrected expression."
    )
    ap.add_argument("--in_h5ad", required=True)
    ap.add_argument("--scanvi_model_dir", required=True)
    ap.add_argument("--gene_selection_tsv", required=True)
    ap.add_argument("--out_dir", required=True)

    ap.add_argument("--celltype_key", default="aifi_L1_paper")
    ap.add_argument("--disease_key", default="disease_Level_3")
    ap.add_argument("--sample_key", default="sample_id")
    ap.add_argument("--patient_key", default="patient_id")
    ap.add_argument("--study_key", default="dataset_id")
    ap.add_argument("--sex_key", default="sex")
    ap.add_argument("--age_key", default="binned_age")
    ap.add_argument("--batch_key", default="dataset_id")

    ap.add_argument("--exclude_celltypes", nargs="*", default=["RBC", "Platelets", "Cycling cells", "Progenitors"])
    ap.add_argument("--unlabeled_category", default="Unknown")

    ap.add_argument("--library_size", type=float, default=1e4)
    ap.add_argument("--n_samples_norm", type=int, default=25)
    ap.add_argument("--batch_size", type=int, default=2048)

    args = ap.parse_args()

    ensure_dir(args.out_dir)

    print("[READ]", args.gene_selection_tsv)
    gs = pd.read_csv(args.gene_selection_tsv, sep="\t")
    if "gene" not in gs.columns:
        raise ValueError("geneSelection_forDownstream.tsv must contain a 'gene' column")
    selected_genes = gs["gene"].astype(str).dropna().drop_duplicates().tolist()
    print(f"[INFO] selected genes from file: {len(selected_genes)}")

    print("[READ]", args.in_h5ad)
    adata = sc.read_h5ad(args.in_h5ad)

    need_obs = [
        args.celltype_key, args.disease_key, args.sample_key,
        args.patient_key, args.study_key, args.sex_key,
        args.age_key, args.batch_key
    ]
    miss = [x for x in need_obs if x not in adata.obs.columns]
    if miss:
        raise ValueError(f"adata.obs missing columns: {miss}")

    adata.obs[args.celltype_key] = adata.obs[args.celltype_key].map(clean_str)
    adata.obs[args.disease_key] = adata.obs[args.disease_key].map(clean_str)
    adata.obs[args.sample_key] = adata.obs[args.sample_key].map(clean_str)
    adata.obs[args.patient_key] = adata.obs[args.patient_key].map(clean_str)
    adata.obs[args.study_key] = adata.obs[args.study_key].map(clean_str)
    adata.obs[args.sex_key] = adata.obs[args.sex_key].map(clean_str)
    adata.obs[args.age_key] = adata.obs[args.age_key].map(clean_str)
    adata.obs[args.batch_key] = adata.obs[args.batch_key].map(clean_str)

    keep_mask = ~adata.obs[args.celltype_key].isin(args.exclude_celltypes)
    if args.unlabeled_category is not None:
        keep_mask &= (adata.obs[args.celltype_key] != args.unlabeled_category)
    adata = adata[keep_mask].copy()

    genes_in_data = [g for g in selected_genes if g in adata.var_names]
    print(f"[INFO] selected genes present in adata: {len(genes_in_data)} / {len(selected_genes)}")
    if len(genes_in_data) == 0:
        raise RuntimeError("No selected genes overlap with adata.var_names")

    print("[LOAD MODEL]", args.scanvi_model_dir)
    model = scvi.model.SCANVI.load(args.scanvi_model_dir, adata=adata)

    transform_batches = sorted(pd.Series(adata.obs[args.batch_key]).dropna().unique().tolist())
    print(f"[INFO] transform_batch count: {len(transform_batches)}")
    print(transform_batches)

    celltypes = sorted(pd.Series(adata.obs[args.celltype_key]).dropna().unique().tolist())

    summary_rows = []

    for ct in celltypes:
        idx = np.where(adata.obs[args.celltype_key].values == ct)[0]
        if len(idx) == 0:
            continue

        # 关键修正：这里只筛细胞，不筛基因
        adata_ct = adata[idx, :].copy()
        print(f"\n[RUN] {ct}: cells={adata_ct.n_obs}, genes(model space)={adata_ct.n_vars}, genes(output)={len(genes_in_data)}")

        norm_expr = model.get_normalized_expression(
            adata=adata_ct,
            gene_list=genes_in_data,
            transform_batch=transform_batches,
            library_size=args.library_size,
            n_samples=args.n_samples_norm,
            return_mean=True,
            batch_size=args.batch_size
        )

        if isinstance(norm_expr, pd.DataFrame):
            X = norm_expr.loc[adata_ct.obs_names, genes_in_data].values
        else:
            X = np.asarray(norm_expr)

        X = np.log1p(X)

        obs_out = adata_ct.obs.copy()
        obs_out = obs_out.assign(
            disease=obs_out[args.disease_key].astype(str),
            sampleID=obs_out[args.sample_key].astype(str),
            patientID=obs_out[args.patient_key].astype(str),
            studyID=obs_out[args.study_key].astype(str)
        )

        out = ad.AnnData(
            X=ssp.csr_matrix(X),
            obs=obs_out,
            var=pd.DataFrame(index=genes_in_data)
        )
        out.uns["source_h5ad"] = args.in_h5ad
        out.uns["scanvi_model_dir"] = args.scanvi_model_dir
        out.uns["gene_selection_tsv"] = args.gene_selection_tsv
        out.uns["celltype_name"] = ct
        out.uns["is_scanvi_corrected"] = True
        out.uns["is_log1p"] = True

        out_fp = os.path.join(args.out_dir, f"{ct}_adataMerged_SPECTRAgenes.log1p.h5ad")
        out.write_h5ad(out_fp)
        print("[WRITE]", out_fp)

        summary_rows.append({
            "celltype": ct,
            "n_cells": out.n_obs,
            "n_genes": out.n_vars,
            "file": out_fp
        })

    summary_df = pd.DataFrame(summary_rows)
    summary_fp = os.path.join(args.out_dir, "SHAP_input_files_summary.tsv")
    summary_df.to_csv(summary_fp, sep="\t", index=False)

    print("\n[DONE]")
    print("summary:", summary_fp, summary_df.shape)
    print(summary_df.to_string(index=False))


if __name__ == "__main__":
    main()