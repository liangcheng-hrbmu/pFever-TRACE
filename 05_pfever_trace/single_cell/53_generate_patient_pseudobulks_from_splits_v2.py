#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import pickle
import argparse
from typing import Dict, List, Tuple, Optional

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

try:
    import scvi
except Exception:
    scvi = None


def ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def clean_str(x):
    if pd.isna(x):
        return np.nan
    x = str(x).strip()
    if x in {"", "nan", "None", "NA"}:
        return np.nan
    return x


def aggregate_embedding_to_pseudobulk(
    emb: np.ndarray,
    obs: pd.DataFrame,
    sample_key: str,
    level1_key: str,
    disease_key: str,
    benchmark_status_map: Optional[Dict[str, str]] = None,
) -> ad.AnnData:
    # 只保留实际存在的观测组合，不产生 category 的笛卡尔积
    obs_tmp = obs[[sample_key, level1_key, disease_key]].copy()
    for c in [sample_key, level1_key, disease_key]:
        obs_tmp[c] = obs_tmp[c].map(clean_str).astype(str)

    df_emb = pd.DataFrame(
        np.asarray(emb, dtype=np.float32),
        index=obs.index,
        columns=[f"latent_{i}" for i in range(np.asarray(emb).shape[1])]
    )

    join_df = pd.concat([obs_tmp, df_emb], axis=1)

    pb = (
        join_df
        .groupby([sample_key, level1_key, disease_key], observed=True, dropna=False)
        .agg(
            n_cells=(sample_key, "size"),
            **{f"latent_{i}": (f"latent_{i}", "mean") for i in range(df_emb.shape[1])}
        )
        .reset_index()
    )

    # 理论上这里不该再有 NaN；留个保护
    latent_cols = [c for c in pb.columns if c.startswith("latent_")]
    pb = pb.dropna(subset=latent_cols, how="any").copy()

    if benchmark_status_map is None:
        pb["benchmark_status"] = "unknown"
    else:
        pb["benchmark_status"] = pb[disease_key].map(benchmark_status_map).fillna("unknown")

    obs_out = pb[[sample_key, level1_key, disease_key, "n_cells", "benchmark_status"]].copy()
    # 用 sample|Level1|disease 保证唯一
    obs_out.index = (
        pb[sample_key].astype(str) + "|" +
        pb[level1_key].astype(str) + "|" +
        pb[disease_key].astype(str)
    )

    X = pb[latent_cols].values.astype(np.float32)
    var_out = pd.DataFrame(index=latent_cols)

    adata_pb = ad.AnnData(X=X, obs=obs_out, var=var_out)
    adata_pb.obs_names_make_unique()
    return adata_pb


def get_embedding_matrix(
    adata,
    obs_sub: pd.DataFrame,
    embedding_obsm: Optional[str],
    scanvi_model_dir: Optional[str],
) -> np.ndarray:
    if embedding_obsm is not None and embedding_obsm != "NONE" and embedding_obsm in adata.obsm.keys():
        idx = adata.obs_names.get_indexer(obs_sub.index)
        X = adata.obsm[embedding_obsm][idx]
        return np.asarray(X, dtype=np.float32)

    if scanvi_model_dir is None:
        raise ValueError(
            f"embedding_obsm={embedding_obsm} not found in adata.obsm and no --scanvi_model_dir provided."
        )
    if scvi is None:
        raise ImportError("scvi-tools is required when embedding_obsm is unavailable.")

    adata_sub = adata[obs_sub.index].copy()

    model = None
    try:
        model = scvi.model.SCANVI.load(scanvi_model_dir, adata=adata)
    except Exception:
        try:
            model = scvi.model.SCVI.load(scanvi_model_dir, adata=adata)
        except Exception as e:
            raise RuntimeError(f"Failed to load scANVI/SCVI model from {scanvi_model_dir}: {e}")

    latent = model.get_latent_representation(adata=adata_sub)
    return np.asarray(latent, dtype=np.float32)


def main():
    ap = argparse.ArgumentParser(
        description="Generate reference/query patient pseudobulk embeddings from split cell IDs (NaN-safe version)."
    )
    ap.add_argument("--in_h5ad", required=True)
    ap.add_argument("--kfold_pkl", required=True)
    ap.add_argument("--out_dir", required=True)

    ap.add_argument("--embedding_obsm", default="X_scANVI")
    ap.add_argument("--scanvi_model_dir", default=None)
    ap.add_argument("--level1_key", default="aifi_L1_paper")
    ap.add_argument("--query_level1_key", default=None)
    ap.add_argument("--disease_key", default="disease_Level_3")
    ap.add_argument("--sample_key", default="sample_id")
    ap.add_argument("--patient_key", default="patient_id")
    ap.add_argument("--exclude_level1", nargs="*", default=["RBC", "Platelets", "Cycling cells", "Progenitors"])
    ap.add_argument("--kfold_metadata_json", default=None)
    args = ap.parse_args()

    ensure_dir(args.out_dir)
    ensure_dir(os.path.join(args.out_dir, "reference"))
    ensure_dir(os.path.join(args.out_dir, "query"))
    ensure_dir(os.path.join(args.out_dir, "PSEUDOBULKs"))

    q_level1_key = args.query_level1_key if args.query_level1_key else args.level1_key

    print("[READ]", args.in_h5ad)
    adata = sc.read_h5ad(args.in_h5ad)

    need_cols = [args.level1_key, args.disease_key, args.sample_key]
    if q_level1_key not in need_cols:
        need_cols.append(q_level1_key)
    if args.patient_key is not None:
        need_cols.append(args.patient_key)

    miss = [c for c in need_cols if c not in adata.obs.columns]
    if miss:
        raise ValueError(f"adata.obs missing required columns: {miss}")

    obs = adata.obs.copy()
    for c in set(need_cols):
        obs[c] = obs[c].map(clean_str)
    adata.obs = obs

    keep = ~adata.obs[args.level1_key].isin(args.exclude_level1)
    if q_level1_key != args.level1_key:
        keep &= ~adata.obs[q_level1_key].isin(args.exclude_level1)
    adata = adata[keep].copy()

    benchmark_status_map = None
    if args.kfold_metadata_json is not None and os.path.exists(args.kfold_metadata_json):
        import json
        with open(args.kfold_metadata_json, "r", encoding="utf-8") as f:
            meta = json.load(f)
        benchmark_status_map = {}
        for d in meta.get("primary_benchmark_diseases", []):
            benchmark_status_map[str(d)] = "primary_benchmark"
        for d in meta.get("reference_only_diseases", []):
            benchmark_status_map[str(d)] = "reference_only"

    print("[READ]", args.kfold_pkl)
    with open(args.kfold_pkl, "rb") as handle:
        kfold = pickle.load(handle)

    split_summary_rows = []

    for split_id, (train_cells, test_cells) in enumerate(kfold):
        train_cells = pd.Index(np.asarray(train_cells, dtype=str))
        test_cells = pd.Index(np.asarray(test_cells, dtype=str))

        ref_mask = adata.obs_names.astype(str).isin(train_cells)
        qry_mask = adata.obs_names.astype(str).isin(test_cells)

        adataR = adata[ref_mask].copy()
        adataQ = adata[qry_mask].copy()

        ref_fp = os.path.join(args.out_dir, "reference", f"scANVI_SPLIT_{split_id}_reference.h5ad")
        qry_fp = os.path.join(args.out_dir, "query", f"scANVI_SPLIT_{split_id}_query.h5ad")
        adataR.write(ref_fp, compression="gzip")
        adataQ.write(qry_fp, compression="gzip")

        ZR = get_embedding_matrix(adata, adataR.obs.copy(), args.embedding_obsm, args.scanvi_model_dir)
        ZQ = get_embedding_matrix(adata, adataQ.obs.copy(), args.embedding_obsm, args.scanvi_model_dir)

        adataPB_R = aggregate_embedding_to_pseudobulk(
            emb=ZR,
            obs=adataR.obs,
            sample_key=args.sample_key,
            level1_key=args.level1_key,
            disease_key=args.disease_key,
            benchmark_status_map=benchmark_status_map,
        )

        q_obs = adataQ.obs.copy().rename(columns={q_level1_key: "__query_level1__"})
        adataPB_Q = aggregate_embedding_to_pseudobulk(
            emb=ZQ,
            obs=q_obs,
            sample_key=args.sample_key,
            level1_key="__query_level1__",
            disease_key=args.disease_key,
            benchmark_status_map=benchmark_status_map,
        )
        adataPB_Q.obs.rename(columns={"__query_level1__": "Level1"}, inplace=True)
        adataPB_R.obs.rename(columns={args.level1_key: "Level1"}, inplace=True)

        adataPB_R.obs["split"] = split_id
        adataPB_Q.obs["split"] = split_id
        adataPB_R.obs["set"] = "reference"
        adataPB_Q.obs["set"] = "query"

        # 再做一次 NaN 保护
        if np.isnan(np.asarray(adataPB_R.X)).any():
            raise RuntimeError(f"NaN detected in reference pseudobulk matrix for split {split_id}")
        if np.isnan(np.asarray(adataPB_Q.X)).any():
            raise RuntimeError(f"NaN detected in query pseudobulk matrix for split {split_id}")

        ref_pb_fp = os.path.join(args.out_dir, "PSEUDOBULKs", f"scANVI_PSEUDOBULK_reference_SPLIT_{split_id}.h5ad")
        qry_pb_fp = os.path.join(args.out_dir, "PSEUDOBULKs", f"scANVI_PSEUDOBULK_query_SPLIT_{split_id}.h5ad")
        adataPB_R.write(ref_pb_fp, compression="gzip")
        adataPB_Q.write(qry_pb_fp, compression="gzip")

        split_summary_rows.append({
            "split": split_id,
            "reference_cells": adataR.n_obs,
            "query_cells": adataQ.n_obs,
            "reference_pseudobulk_rows": adataPB_R.n_obs,
            "query_pseudobulk_rows": adataPB_Q.n_obs,
            "reference_samples": adataPB_R.obs[args.sample_key].nunique(),
            "query_samples": adataPB_Q.obs[args.sample_key].nunique(),
            "reference_level1": ";".join(sorted(pd.unique(adataPB_R.obs["Level1"]).tolist())),
            "query_level1": ";".join(sorted(pd.unique(adataPB_Q.obs["Level1"]).tolist())),
            "reference_fp": ref_pb_fp,
            "query_fp": qry_pb_fp,
        })

        print(f"[OK] split={split_id} | ref_pb={adataPB_R.shape} | qry_pb={adataPB_Q.shape}")

    summary_df = pd.DataFrame(split_summary_rows)
    summary_fp = os.path.join(args.out_dir, "PSEUDOBULK_summary.tsv")
    summary_df.to_csv(summary_fp, sep="\t", index=False)

    print("\n[DONE]")
    print("summary:", summary_fp, summary_df.shape)
    print(summary_df.to_string(index=False))


if __name__ == "__main__":
    main()
