#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import pickle
import argparse
from typing import Dict, List, Optional, Sequence, Tuple

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


def require_complete_obs(adata: ad.AnnData, columns: Sequence[str], context: str) -> None:
    for col in columns:
        missing = int(adata.obs[col].isna().sum())
        if missing:
            raise ValueError(f"{context}: obs[{col!r}] contains {missing} missing values.")


def one_hot_encode_covariates(
    adata: ad.AnnData,
    covariate_keys: Sequence[str],
) -> Tuple[List[str], Dict[str, List[str]]]:
    encoded_keys: List[str] = []
    encoding_map: Dict[str, List[str]] = {}

    for key in covariate_keys:
        missing = int(adata.obs[key].isna().sum())
        if missing:
            raise ValueError(f"obs[{key!r}] contains {missing} missing values.")

        values = adata.obs[key].astype(str)
        categories = sorted(pd.unique(values).tolist())
        output_keys = []
        for category_id, category in enumerate(categories):
            output_key = f"__covariate_{key}_{category_id}"
            adata.obs[output_key] = (values == category).astype(np.float32)
            encoded_keys.append(output_key)
            output_keys.append(output_key)
        encoding_map[key] = output_keys

    return encoded_keys, encoding_map


def prepare_reference_and_query_obs(
    adata_ref: ad.AnnData,
    adata_query: ad.AnnData,
    labels_key: str,
    batch_key: str,
    categorical_covariate_keys: Sequence[str],
    continuous_covariate_keys: Sequence[str],
    unlabeled_category: str,
) -> None:
    categorical_keys = [batch_key, *categorical_covariate_keys]
    require_complete_obs(adata_ref, [labels_key, *categorical_keys], "reference")
    require_complete_obs(adata_query, categorical_keys, "query")
    require_complete_obs(adata_ref, continuous_covariate_keys, "reference")
    require_complete_obs(adata_query, continuous_covariate_keys, "query")

    reference_labels = pd.Categorical(adata_ref.obs[labels_key].astype(str))
    if unlabeled_category not in reference_labels.categories:
        reference_labels = reference_labels.add_categories([unlabeled_category])
    adata_ref.obs[labels_key] = reference_labels

    query_label_categories = list(reference_labels.categories)
    adata_query.obs[labels_key] = pd.Categorical(
        [unlabeled_category] * adata_query.n_obs,
        categories=query_label_categories,
    )

    for key in categorical_keys:
        adata_ref.obs[key] = adata_ref.obs[key].astype(str).astype("category")
        adata_query.obs[key] = adata_query.obs[key].astype(str).astype("category")

    for key in continuous_covariate_keys:
        adata_ref.obs[key] = pd.to_numeric(adata_ref.obs[key], errors="raise")
        adata_query.obs[key] = pd.to_numeric(adata_query.obs[key], errors="raise")


def train_and_map_fold_scanvi(
    adata_ref: ad.AnnData,
    adata_query: ad.AnnData,
    counts_layer: Optional[str],
    labels_key: str,
    batch_key: str,
    categorical_covariate_keys: Sequence[str],
    continuous_covariate_keys: Sequence[str],
    unlabeled_category: str,
    n_latent: int,
    max_epochs_scvi: int,
    max_epochs_scanvi: int,
    max_epochs_query: int,
    seed: int,
    model_dir: Optional[str],
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, Dict[str, str]]:
    if scvi is None:
        raise ImportError("scvi-tools is required for fold-specific scVI/scANVI fitting.")

    prepare_reference_and_query_obs(
        adata_ref=adata_ref,
        adata_query=adata_query,
        labels_key=labels_key,
        batch_key=batch_key,
        categorical_covariate_keys=categorical_covariate_keys,
        continuous_covariate_keys=continuous_covariate_keys,
        unlabeled_category=unlabeled_category,
    )

    scvi.settings.seed = seed
    scvi.model.SCVI.setup_anndata(
        adata_ref,
        layer=counts_layer,
        batch_key=batch_key,
        labels_key=labels_key,
        categorical_covariate_keys=list(categorical_covariate_keys) or None,
        continuous_covariate_keys=list(continuous_covariate_keys) or None,
    )

    reference_scvi = scvi.model.SCVI(adata_ref, n_latent=n_latent)
    reference_scvi.train(max_epochs=max_epochs_scvi)

    reference_scanvi = scvi.model.SCANVI.from_scvi_model(
        reference_scvi,
        unlabeled_category=unlabeled_category,
    )
    reference_scanvi.train(max_epochs=max_epochs_scanvi)

    reference_latent = np.asarray(
        reference_scanvi.get_latent_representation(),
        dtype=np.float32,
    )

    scvi.model.SCANVI.prepare_query_anndata(adata_query, reference_scanvi)
    query_scanvi = scvi.model.SCANVI.load_query_data(adata_query, reference_scanvi)
    query_scanvi.train(
        max_epochs=max_epochs_query,
        plan_kwargs={"weight_decay": 0.0},
    )

    query_latent = np.asarray(
        query_scanvi.get_latent_representation(),
        dtype=np.float32,
    )
    query_level1 = np.asarray(query_scanvi.predict(), dtype=str)

    model_paths: Dict[str, str] = {}
    if model_dir is not None:
        ensure_dir(model_dir)
        model_paths = {
            "reference_scvi_model": os.path.join(model_dir, "reference_scvi"),
            "reference_scanvi_model": os.path.join(model_dir, "reference_scanvi"),
            "query_scanvi_model": os.path.join(model_dir, "query_scanvi"),
        }
        reference_scvi.save(model_paths["reference_scvi_model"], overwrite=True)
        reference_scanvi.save(model_paths["reference_scanvi_model"], overwrite=True)
        query_scanvi.save(model_paths["query_scanvi_model"], overwrite=True)

    return reference_latent, query_latent, query_level1, model_paths


def main():
    ap = argparse.ArgumentParser(
        description=(
            "Fit fold-specific scVI/scANVI references, map held-out cells as unlabeled queries, "
            "and generate reference/query patient pseudobulk embeddings."
        )
    )
    ap.add_argument("--in_h5ad", required=True)
    ap.add_argument("--kfold_pkl", required=True)
    ap.add_argument("--out_dir", required=True)

    ap.add_argument("--level1_key", default="aifi_L1_paper")
    ap.add_argument("--disease_key", default="disease_Level_3")
    ap.add_argument("--sample_key", default="sample_id")
    ap.add_argument("--patient_key", default="patient_id")
    ap.add_argument("--counts_layer", default="counts",
                    help="Raw-count layer used to fit fold-specific scVI models; use NONE for adata.X")
    ap.add_argument("--batch_key", default="dataset_id")
    ap.add_argument("--categorical_covariate_keys", nargs="*", default=["sex", "binned_age"])
    ap.add_argument("--continuous_covariate_keys", nargs="*", default=[])
    ap.add_argument("--unlabeled_category", default="Unknown")
    ap.add_argument("--n_latent", type=int, default=30)
    ap.add_argument("--max_epochs_scvi", type=int, default=200)
    ap.add_argument("--max_epochs_scanvi", type=int, default=50)
    ap.add_argument("--max_epochs_query", type=int, default=100)
    ap.add_argument("--seed", type=int, default=5)
    ap.add_argument(
        "--save_models",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Save the fold-specific reference and query models (default: false)",
    )
    ap.add_argument("--exclude_level1", nargs="*", default=["RBC", "Platelets", "Cycling cells", "Progenitors"])
    ap.add_argument("--kfold_metadata_json", default=None)
    args = ap.parse_args()

    ensure_dir(args.out_dir)
    ensure_dir(os.path.join(args.out_dir, "reference"))
    ensure_dir(os.path.join(args.out_dir, "query"))
    ensure_dir(os.path.join(args.out_dir, "PSEUDOBULKs"))
    if args.save_models:
        ensure_dir(os.path.join(args.out_dir, "models"))

    print("[READ]", args.in_h5ad)
    adata = sc.read_h5ad(args.in_h5ad)

    need_cols = [
        args.level1_key,
        args.disease_key,
        args.sample_key,
        args.batch_key,
        *args.categorical_covariate_keys,
        *args.continuous_covariate_keys,
    ]
    if args.patient_key is not None:
        need_cols.append(args.patient_key)
    need_cols = list(dict.fromkeys(need_cols))

    miss = [c for c in need_cols if c not in adata.obs.columns]
    if miss:
        raise ValueError(f"adata.obs missing required columns: {miss}")

    obs = adata.obs.copy()
    for c in set(need_cols):
        obs[c] = obs[c].map(clean_str)
    adata.obs = obs

    counts_layer = None if str(args.counts_layer).upper() == "NONE" else args.counts_layer
    if counts_layer is not None and counts_layer not in adata.layers:
        raise ValueError(f"adata.layers missing counts layer: {counts_layer}")

    # scvi-tools 1.1.2 scArches cannot map queries with extra categorical
    # covariates, so encode them as fixed numeric inputs before splitting.
    encoded_covariate_keys, covariate_encoding_map = one_hot_encode_covariates(
        adata,
        args.categorical_covariate_keys,
    )
    model_continuous_covariate_keys = [
        *args.continuous_covariate_keys,
        *encoded_covariate_keys,
    ]
    if covariate_encoding_map:
        print("[INFO] one-hot encoded categorical covariates:", covariate_encoding_map)

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

        if adataR.n_obs == 0 or adataQ.n_obs == 0:
            raise RuntimeError(
                f"split={split_id}: empty reference/query after matching cell IDs "
                f"(reference={adataR.n_obs}, query={adataQ.n_obs})."
            )
        if adataR.obs_names.intersection(adataQ.obs_names).size:
            raise RuntimeError(f"split={split_id}: overlapping reference/query cells detected.")

        model_dir = (
            os.path.join(args.out_dir, "models", f"SPLIT_{split_id}")
            if args.save_models else None
        )
        ZR, ZQ, query_level1, model_paths = train_and_map_fold_scanvi(
            adata_ref=adataR,
            adata_query=adataQ,
            counts_layer=counts_layer,
            labels_key=args.level1_key,
            batch_key=args.batch_key,
            categorical_covariate_keys=[],
            continuous_covariate_keys=model_continuous_covariate_keys,
            unlabeled_category=args.unlabeled_category,
            n_latent=args.n_latent,
            max_epochs_scvi=args.max_epochs_scvi,
            max_epochs_scanvi=args.max_epochs_scanvi,
            max_epochs_query=args.max_epochs_query,
            seed=args.seed + split_id,
            model_dir=model_dir,
        )

        adataQ.obs["scanvi_predicted_level1"] = query_level1

        ref_core_mask = ~adataR.obs[args.level1_key].isin(args.exclude_level1).values
        qry_core_mask = ~adataQ.obs["scanvi_predicted_level1"].isin(args.exclude_level1).values
        adataR = adataR[ref_core_mask].copy()
        adataQ = adataQ[qry_core_mask].copy()
        ZR = ZR[ref_core_mask]
        ZQ = ZQ[qry_core_mask]

        adataR.obsm["X_fold_scANVI"] = ZR
        adataQ.obsm["X_fold_scANVI"] = ZQ

        ref_fp = os.path.join(args.out_dir, "reference", f"scANVI_SPLIT_{split_id}_reference.h5ad")
        qry_fp = os.path.join(args.out_dir, "query", f"scANVI_SPLIT_{split_id}_query.h5ad")
        adataR.write(ref_fp, compression="gzip")
        adataQ.write(qry_fp, compression="gzip")

        adataPB_R = aggregate_embedding_to_pseudobulk(
            emb=ZR,
            obs=adataR.obs,
            sample_key=args.sample_key,
            level1_key=args.level1_key,
            disease_key=args.disease_key,
            benchmark_status_map=benchmark_status_map,
        )

        q_obs = adataQ.obs.copy().rename(columns={"scanvi_predicted_level1": "__query_level1__"})
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
            "reference_scvi_model": model_paths.get("reference_scvi_model", ""),
            "reference_scanvi_model": model_paths.get("reference_scanvi_model", ""),
            "query_scanvi_model": model_paths.get("query_scanvi_model", ""),
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
