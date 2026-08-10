#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import json
import pickle
import argparse
from typing import Dict, List, Tuple, Optional

import numpy as np
import pandas as pd
import anndata as ad


def ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def clean_str(x):
    if pd.isna(x):
        return np.nan
    x = str(x).strip()
    if x in {"", "nan", "None", "NA"}:
        return np.nan
    return x


def normalize_ct_name(x: str) -> str:
    """Normalize cell-type names parsed from file names or CLI arguments.

    InstaPrism Z files use safe names such as T_CD4_non_naive, whereas
    downstream metadata and command-line whitelists use names such as
    T CD4 non-naive.  This function keeps these names consistent so
    whitelist filtering does not accidentally drop non-naive T-cell types.
    """
    x = str(x).strip()
    x = x.replace("non naive", "non-naive")
    x = re.sub(r"\s+", " ", x)
    return x


def infer_level1_from_filename(fp: str) -> str:
    base = os.path.basename(fp)
    m = re.match(r"^Z_(.+?)_allgenes_genes_by_samples\.tsv$", base)
    if m is None:
        raise ValueError(f"Cannot parse Level1 from filename: {base}")

    safe_ct = m.group(1)

    special_map = {
        "T_CD4_non_naive": "T CD4 non-naive",
        "T_CD8_non_naive": "T CD8 non-naive",
        "T_CD4_naive": "T CD4 naive",
        "T_CD8_naive": "T CD8 naive",
    }
    if safe_ct in special_map:
        return special_map[safe_ct]

    return normalize_ct_name(safe_ct.replace("_", " "))


def read_expr_table(fp: str) -> pd.DataFrame:
    df = pd.read_csv(fp, sep="\t", index_col=0)
    df.index = df.index.astype(str).str.strip()
    df.columns = df.columns.astype(str).str.strip()
    df = df.loc[~df.index.duplicated(), :]
    df = df.loc[df.sum(axis=1) > 0, :]
    return df


def normalize_logcpm(expr_gxs: pd.DataFrame) -> pd.DataFrame:
    libsize = expr_gxs.sum(axis=0)
    libsize[libsize <= 0] = 1.0
    cpm = expr_gxs.div(libsize, axis=1) * 1e6
    return np.log1p(cpm)


def build_row_block(
    expr_gxs: pd.DataFrame,
    level1_name: str,
    sample_set: set,
    meta_df: pd.DataFrame,
    sample_key: str,
    patient_key: str,
    disease_key: str,
    benchmark_status_map: Optional[Dict[str, str]] = None,
) -> Tuple[np.ndarray, pd.DataFrame, List[str]]:
    samples_use = [s for s in expr_gxs.columns if s in sample_set]
    if len(samples_use) == 0:
        return np.empty((0, 0), dtype=np.float32), pd.DataFrame(), []

    expr_sub = expr_gxs.loc[:, samples_use]
    expr_sub = normalize_logcpm(expr_sub)

    meta_sub = meta_df.set_index(sample_key).loc[samples_use].reset_index()

    if benchmark_status_map is None:
        bench = ["unknown"] * meta_sub.shape[0]
    else:
        bench = meta_sub[disease_key].map(benchmark_status_map).fillna("unknown").tolist()

    obs_out = meta_sub[[sample_key, patient_key, disease_key]].copy()
    obs_out["Level1"] = level1_name
    obs_out["benchmark_status"] = bench
    obs_out["n_genes"] = expr_sub.shape[0]
    obs_out.index = (
        obs_out[sample_key].astype(str) + "|" +
        obs_out["Level1"].astype(str) + "|" +
        obs_out[disease_key].astype(str)
    )

    X = expr_sub.T.values.astype(np.float32)   # row = sample-level1
    var_names = expr_sub.index.tolist()

    return X, obs_out, var_names


def stack_to_anndata(blocks: List[Tuple[np.ndarray, pd.DataFrame]], var_names: List[str]) -> ad.AnnData:
    Xs = [x for x, obs in blocks if x.shape[0] > 0]
    Obs = [obs for x, obs in blocks if x.shape[0] > 0]

    if len(Xs) == 0:
        raise RuntimeError("No rows generated for this split.")

    X = np.vstack(Xs).astype(np.float32)
    obs = pd.concat(Obs, axis=0)
    var = pd.DataFrame(index=var_names)

    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.obs_names_make_unique()
    return adata


def main():
    ap = argparse.ArgumentParser(
        description="Generate bulk reference/query pseudobulk h5ad from InstaPrism Z matrices."
    )
    ap.add_argument("--z_dir", required=True, help="Directory with Z_<celltype>_allgenes_genes_by_samples.tsv")
    ap.add_argument("--meta_tsv", required=True, help="Bulk metadata TSV")
    ap.add_argument("--kfold_pkl", required=True, help="K_FOLD_sampleID.pkl from 52b")
    ap.add_argument("--out_dir", required=True)
    ap.add_argument("--sample_key", default="sample_id")
    ap.add_argument("--patient_key", default="patient_id")
    ap.add_argument("--disease_key", default="disease_broad5")
    ap.add_argument("--kfold_metadata_json", default=None)
    ap.add_argument("--allowed_level1", nargs="*", default=None,
                    help="Optional whitelist, e.g. Mono B Plasma 'T CD4 naive' ...")
    args = ap.parse_args()

    ensure_dir(args.out_dir)
    ensure_dir(os.path.join(args.out_dir, "reference"))
    ensure_dir(os.path.join(args.out_dir, "query"))
    ensure_dir(os.path.join(args.out_dir, "PSEUDOBULKs"))

    print("[READ META]", args.meta_tsv)
    meta = pd.read_csv(args.meta_tsv, sep="\t")
    need_cols = [args.sample_key, args.disease_key]
    if args.patient_key not in meta.columns:
        meta[args.patient_key] = meta[args.sample_key]
    need_cols.append(args.patient_key)

    miss = [c for c in need_cols if c not in meta.columns]
    if miss:
        raise ValueError(f"Metadata missing required columns: {miss}")

    for c in need_cols:
        meta[c] = meta[c].map(clean_str)

    meta = meta.dropna(subset=[args.sample_key, args.patient_key, args.disease_key]).copy()
    meta = meta.drop_duplicates(subset=[args.sample_key]).copy()

    benchmark_status_map = None
    if args.kfold_metadata_json and os.path.exists(args.kfold_metadata_json):
        with open(args.kfold_metadata_json, "r", encoding="utf-8") as f:
            meta_json = json.load(f)
        benchmark_status_map = {}
        for d in meta_json.get("primary_benchmark_diseases", []):
            benchmark_status_map[str(d)] = "primary_benchmark"
        for d in meta_json.get("reference_only_diseases", []):
            benchmark_status_map[str(d)] = "reference_only"

    print("[READ SPLITS]", args.kfold_pkl)
    with open(args.kfold_pkl, "rb") as handle:
        kfold = pickle.load(handle)

    z_files = sorted([
        os.path.join(args.z_dir, x) for x in os.listdir(args.z_dir)
        if re.match(r"^Z_.+_allgenes_genes_by_samples\.tsv$", x)
    ])
    if len(z_files) == 0:
        raise FileNotFoundError(f"No Z_*.tsv found in {args.z_dir}")

    ct_to_expr = {}
    gene_order = None

    allowed_set = None
    if args.allowed_level1 is not None:
        allowed_set = set(normalize_ct_name(x) for x in args.allowed_level1)
        print("[INFO] allowed_level1 normalized:", ";".join(sorted(allowed_set)))

    for fp in z_files:
        ct = normalize_ct_name(infer_level1_from_filename(fp))
        if allowed_set is not None and ct not in allowed_set:
            print(f"[SKIP] {os.path.basename(fp)} parsed as '{ct}' not in allowed_level1")
            continue

        expr = read_expr_table(fp)

        common_samples = [s for s in expr.columns if s in set(meta[args.sample_key].astype(str))]
        expr = expr.loc[:, common_samples]
        if expr.shape[1] == 0:
            continue

        if gene_order is None:
            gene_order = expr.index.tolist()
        else:
            expr = expr.loc[[g for g in gene_order if g in expr.index], :]
            if len(expr.index) != len(gene_order):
                gene_order = expr.index.tolist()
                for k in list(ct_to_expr.keys()):
                    ct_to_expr[k] = ct_to_expr[k].loc[gene_order, :]

        ct_to_expr[ct] = expr

    if len(ct_to_expr) == 0:
        raise RuntimeError("No usable Level1 Z matrices found after filtering.")

    print("[INFO] Level1 matrices loaded:", ";".join(sorted(ct_to_expr.keys())))
    print("[INFO] n_Level1:", len(ct_to_expr))

    split_summary_rows = []

    for split_id, (train_samples, test_samples) in enumerate(kfold):
        train_samples = set(pd.Index(np.asarray(train_samples, dtype=str)).tolist())
        test_samples = set(pd.Index(np.asarray(test_samples, dtype=str)).tolist())

        ref_blocks = []
        qry_blocks = []

        for ct, expr in ct_to_expr.items():
            expr = expr.loc[gene_order, :]

            X_ref, obs_ref, var_names = build_row_block(
                expr_gxs=expr,
                level1_name=ct,
                sample_set=train_samples,
                meta_df=meta,
                sample_key=args.sample_key,
                patient_key=args.patient_key,
                disease_key=args.disease_key,
                benchmark_status_map=benchmark_status_map,
            )
            if X_ref.shape[0] > 0:
                ref_blocks.append((X_ref, obs_ref))

            X_qry, obs_qry, var_names = build_row_block(
                expr_gxs=expr,
                level1_name=ct,
                sample_set=test_samples,
                meta_df=meta,
                sample_key=args.sample_key,
                patient_key=args.patient_key,
                disease_key=args.disease_key,
                benchmark_status_map=benchmark_status_map,
            )
            if X_qry.shape[0] > 0:
                qry_blocks.append((X_qry, obs_qry))

        ad_ref = stack_to_anndata(ref_blocks, gene_order)
        ad_qry = stack_to_anndata(qry_blocks, gene_order)

        ad_ref.obs["split"] = split_id
        ad_qry.obs["split"] = split_id
        ad_ref.obs["set"] = "reference"
        ad_qry.obs["set"] = "query"

        ref_pb_fp = os.path.join(args.out_dir, "PSEUDOBULKs", f"bulk_PSEUDOBULK_reference_SPLIT_{split_id}.h5ad")
        qry_pb_fp = os.path.join(args.out_dir, "PSEUDOBULKs", f"bulk_PSEUDOBULK_query_SPLIT_{split_id}.h5ad")
        ad_ref.write(ref_pb_fp, compression="gzip")
        ad_qry.write(qry_pb_fp, compression="gzip")

        ad_ref.obs.to_csv(
            os.path.join(args.out_dir, "reference", f"bulk_reference_samples_SPLIT_{split_id}.tsv"),
            sep="\t", index=True
        )
        ad_qry.obs.to_csv(
            os.path.join(args.out_dir, "query", f"bulk_query_samples_SPLIT_{split_id}.tsv"),
            sep="\t", index=True
        )

        split_summary_rows.append({
            "split": split_id,
            "reference_rows": ad_ref.n_obs,
            "query_rows": ad_qry.n_obs,
            "reference_samples": ad_ref.obs[args.sample_key].nunique(),
            "query_samples": ad_qry.obs[args.sample_key].nunique(),
            "n_genes": ad_ref.n_vars,
            "level1_used": ";".join(sorted(ad_ref.obs["Level1"].astype(str).unique().tolist())),
            "reference_fp": ref_pb_fp,
            "query_fp": qry_pb_fp,
        })

        print(f"[OK] split={split_id} | ref_pb={ad_ref.shape} | qry_pb={ad_qry.shape}")

    summary_df = pd.DataFrame(split_summary_rows)
    summary_df.to_csv(os.path.join(args.out_dir, "PSEUDOBULK_summary.tsv"), sep="\t", index=False)

    print("\n[DONE]")
    print(summary_df.to_string(index=False))


if __name__ == "__main__":
    main()