#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import json
import pickle
import argparse
from typing import List, Tuple, Dict

import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.model_selection import StratifiedGroupKFold


def ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def clean_str(x):
    if pd.isna(x):
        return np.nan
    x = str(x).strip()
    if x in {"", "nan", "None", "NA"}:
        return np.nan
    return x


def build_sample_table(
    obs: pd.DataFrame,
    disease_key: str,
    sample_key: str,
    patient_key: str | None = None,
) -> pd.DataFrame:
    keep_cols = [sample_key, disease_key]
    if patient_key is not None:
        keep_cols.append(patient_key)
    sample_df = obs[keep_cols].drop_duplicates().copy()

    # 防止同一个 sample 混了多个 disease
    n_per_sample = (
        sample_df.groupby(sample_key, dropna=False)[disease_key]
        .nunique()
        .reset_index(name="n_diseases_per_sample")
    )
    bad = n_per_sample.loc[n_per_sample["n_diseases_per_sample"] > 1, sample_key].tolist()
    if bad:
        raise ValueError(
            f"Found samples mapping to multiple diseases. First few problematic samples: {bad[:10]}"
        )

    sample_df = sample_df.drop_duplicates(subset=[sample_key]).copy()
    return sample_df


def main():
    ap = argparse.ArgumentParser(
        description=(
            "Generate 5-fold patient/sample grouped splits for patient-classifier workflow, "
            "while keeping underrepresented diseases in reference/train only."
        )
    )
    ap.add_argument("--in_h5ad", required=True, help="Input AnnData object")
    ap.add_argument("--out_dir", required=True, help="Output directory")
    ap.add_argument("--obs_name_key", default=None,
                    help="Optional obs column to use as cell IDs instead of adata.obs_names")
    ap.add_argument("--disease_key", default="disease_Level_3")
    ap.add_argument("--sample_key", default="sample_id")
    ap.add_argument("--patient_key", default="patient_id")
    ap.add_argument("--n_splits", type=int, default=5)
    ap.add_argument("--shuffle", action="store_true", help="Enable shuffle in StratifiedGroupKFold")
    ap.add_argument("--random_state", type=int, default=None)
    ap.add_argument("--exclude_diseases", nargs="*", default=None,
                    help="Diseases to remove entirely from reference and evaluation")
    ap.add_argument("--min_samples_for_eval", type=int, default=5,
                    help="Minimum number of samples required for a disease to enter primary benchmark")
    args = ap.parse_args()

    ensure_dir(args.out_dir)

    print("[READ]", args.in_h5ad)
    adata = sc.read_h5ad(args.in_h5ad, backed="r")

    need_cols = [args.disease_key, args.sample_key]
    if args.patient_key is not None:
        need_cols.append(args.patient_key)
    if args.obs_name_key is not None:
        need_cols.append(args.obs_name_key)

    miss = [c for c in need_cols if c not in adata.obs.columns]
    if miss:
        raise ValueError(f"adata.obs missing required columns: {miss}")

    obs = adata.obs.copy()
    obs[args.disease_key] = obs[args.disease_key].map(clean_str)
    obs[args.sample_key] = obs[args.sample_key].map(clean_str)
    if args.patient_key is not None:
        obs[args.patient_key] = obs[args.patient_key].map(clean_str)
    if args.obs_name_key is not None:
        obs[args.obs_name_key] = obs[args.obs_name_key].map(clean_str)

    keep = obs[args.disease_key].notna() & obs[args.sample_key].notna()
    if args.exclude_diseases:
        keep &= ~obs[args.disease_key].isin(args.exclude_diseases)
    obs = obs.loc[keep].copy()

    if obs.shape[0] == 0:
        raise RuntimeError("No cells left after filtering.")

    cell_ids = obs[args.obs_name_key].astype(str).values if args.obs_name_key else obs.index.astype(str).values
    obs["_cell_id_export"] = cell_ids

    sample_df = build_sample_table(
        obs=obs,
        disease_key=args.disease_key,
        sample_key=args.sample_key,
        patient_key=args.patient_key,
    )

    disease_support = (
        sample_df.groupby(args.disease_key, dropna=False)[args.sample_key]
        .nunique()
        .reset_index(name="n_samples")
        .sort_values(["n_samples", args.disease_key], ascending=[True, True])
        .reset_index(drop=True)
    )
    disease_support["benchmark_status"] = np.where(
        disease_support["n_samples"] >= args.min_samples_for_eval,
        "primary_benchmark",
        "reference_only",
    )

    primary_diseases = disease_support.loc[
        disease_support["benchmark_status"] == "primary_benchmark", args.disease_key
    ].astype(str).tolist()
    reference_only_diseases = disease_support.loc[
        disease_support["benchmark_status"] == "reference_only", args.disease_key
    ].astype(str).tolist()

    print("[INFO] diseases kept in primary benchmark:", primary_diseases)
    print("[INFO] diseases kept as reference-only:", reference_only_diseases)

    if len(primary_diseases) == 0:
        raise RuntimeError("No diseases satisfy min_samples_for_eval. Lower the threshold.")

    # 只对 primary benchmark diseases 做 5-fold
    primary_sample_df = sample_df[sample_df[args.disease_key].isin(primary_diseases)].copy()
    n_primary_samples = primary_sample_df.shape[0]
    if n_primary_samples < args.n_splits:
        raise RuntimeError(
            f"Only {n_primary_samples} primary-benchmark samples available, smaller than n_splits={args.n_splits}."
        )

    splitter = StratifiedGroupKFold(
        n_splits=args.n_splits,
        shuffle=args.shuffle,
        random_state=args.random_state,
    )

    primary_sample_ids = primary_sample_df[args.sample_key].astype(str).values
    primary_labels = primary_sample_df[args.disease_key].astype(str).values
    primary_groups = primary_sample_df[args.sample_key].astype(str).values  # same as sample-level groups

    split_indices = list(
        splitter.split(
            X=np.arange(primary_sample_df.shape[0]),
            y=primary_labels,
            groups=primary_groups,
        )
    )

    # rare/reference-only samples are ALWAYS assigned to train/reference
    reference_only_samples = (
        sample_df.loc[sample_df[args.disease_key].isin(reference_only_diseases), args.sample_key]
        .astype(str).tolist()
    )
    reference_only_sample_set = set(reference_only_samples)

    kfold_cell_ids: List[Tuple[np.ndarray, np.ndarray]] = []
    summary_rows: List[Dict] = []
    sample_rows: List[pd.DataFrame] = []

    for split_id, (train_idx, test_idx) in enumerate(split_indices):
        train_primary_samples = set(primary_sample_ids[train_idx].tolist())
        test_primary_samples = set(primary_sample_ids[test_idx].tolist())

        # train/reference = primary train + all reference_only
        train_sample_set = set(train_primary_samples) | reference_only_sample_set
        test_sample_set = set(test_primary_samples)

        train_mask = obs[args.sample_key].astype(str).isin(train_sample_set).values
        test_mask = obs[args.sample_key].astype(str).isin(test_sample_set).values

        # sanity: train/test should not overlap
        if np.any(train_mask & test_mask):
            raise RuntimeError(f"Split {split_id}: overlapping train/test cells detected.")

        train_cells = obs.loc[train_mask, "_cell_id_export"].astype(str).values
        test_cells = obs.loc[test_mask, "_cell_id_export"].astype(str).values
        kfold_cell_ids.append((np.array(train_cells, dtype=object), np.array(test_cells, dtype=object)))

        train_samples_sorted = sorted(train_sample_set)
        test_samples_sorted = sorted(test_sample_set)

        train_diseases = sorted(
            sample_df.loc[sample_df[args.sample_key].isin(train_samples_sorted), args.disease_key]
            .astype(str).unique().tolist()
        )
        test_diseases = sorted(
            sample_df.loc[sample_df[args.sample_key].isin(test_samples_sorted), args.disease_key]
            .astype(str).unique().tolist()
        )

        summary_rows.append({
            "split": split_id,
            "n_train_cells": int(train_mask.sum()),
            "n_test_cells": int(test_mask.sum()),
            "n_train_samples": len(train_samples_sorted),
            "n_test_samples": len(test_samples_sorted),
            "n_reference_only_samples_in_train": len(reference_only_sample_set),
            "primary_benchmark_diseases": ";".join(primary_diseases),
            "reference_only_diseases": ";".join(reference_only_diseases),
            "train_diseases": ";".join(train_diseases),
            "test_diseases": ";".join(test_diseases),
        })

        split_sample_df = sample_df.copy()
        split_sample_df["split"] = split_id
        split_sample_df["set"] = np.where(
            split_sample_df[args.sample_key].astype(str).isin(test_sample_set),
            "test",
            np.where(
                split_sample_df[args.disease_key].astype(str).isin(reference_only_diseases),
                "reference_only_train",
                "train",
            ),
        )
        split_sample_df["benchmark_status"] = np.where(
            split_sample_df[args.disease_key].astype(str).isin(reference_only_diseases),
            "reference_only",
            "primary_benchmark",
        )
        sample_rows.append(split_sample_df)

    # write pickle compatible with original notebook (train/test cell IDs)
    pkl_fp = os.path.join(args.out_dir, "K_FOLD_cellID.pkl")
    with open(pkl_fp, "wb") as handle:
        pickle.dump(kfold_cell_ids, handle, protocol=4)

    summary_df = pd.DataFrame(summary_rows)
    summary_fp = os.path.join(args.out_dir, "K_FOLD_summary.tsv")
    summary_df.to_csv(summary_fp, sep="\t", index=False)

    sample_split_df = pd.concat(sample_rows, axis=0, ignore_index=True).drop_duplicates()
    sample_split_fp = os.path.join(args.out_dir, "K_FOLD_sample_split.tsv")
    sample_split_df.to_csv(sample_split_fp, sep="\t", index=False)

    disease_count_rows = []
    for split_id in range(args.n_splits):
        sub = sample_split_df[sample_split_df["split"] == split_id].copy()
        tab = (
            sub.groupby(["set", args.disease_key, "benchmark_status"], dropna=False, observed=False)[args.sample_key]
            .nunique()
            .reset_index(name="n_samples")
        )
        tab["split"] = split_id
        disease_count_rows.append(tab)
    disease_count_df = pd.concat(disease_count_rows, axis=0, ignore_index=True)
    disease_count_fp = os.path.join(args.out_dir, "K_FOLD_sample_disease_counts.tsv")
    disease_count_df.to_csv(disease_count_fp, sep="\t", index=False)

    disease_support_fp = os.path.join(args.out_dir, "K_FOLD_disease_support.tsv")
    disease_support.to_csv(disease_support_fp, sep="\t", index=False)

    metadata = {
        "n_splits": args.n_splits,
        "shuffle": bool(args.shuffle),
        "random_state": args.random_state,
        "disease_key": args.disease_key,
        "sample_key": args.sample_key,
        "patient_key": args.patient_key,
        "min_samples_for_eval": args.min_samples_for_eval,
        "primary_benchmark_diseases": primary_diseases,
        "reference_only_diseases": reference_only_diseases,
        "excluded_diseases": args.exclude_diseases or [],
    }
    meta_fp = os.path.join(args.out_dir, "K_FOLD_metadata.json")
    with open(meta_fp, "w", encoding="utf-8") as f:
        json.dump(metadata, f, ensure_ascii=False, indent=2)

    print("\n[DONE]")
    print("pickle:", pkl_fp)
    print("summary:", summary_fp, summary_df.shape)
    print("sample split:", sample_split_fp, sample_split_df.shape)
    print("sample disease counts:", disease_count_fp, disease_count_df.shape)
    print("disease support:", disease_support_fp, disease_support.shape)
    print("metadata:", meta_fp)
    print(summary_df.to_string(index=False))


if __name__ == "__main__":
    main()
