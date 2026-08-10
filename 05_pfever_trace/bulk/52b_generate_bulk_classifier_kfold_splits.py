#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import json
import pickle
import argparse
from typing import List, Dict, Tuple

import numpy as np
import pandas as pd
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


def build_sample_table(meta: pd.DataFrame,
                       disease_key: str,
                       sample_key: str,
                       patient_key: str | None = None) -> pd.DataFrame:
    keep_cols = [sample_key, disease_key]
    if patient_key is not None and patient_key in meta.columns:
        keep_cols.append(patient_key)

    sample_df = meta[keep_cols].drop_duplicates().copy()

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
        description="Generate 5-fold grouped splits for bulk patient-classifier workflow."
    )
    ap.add_argument("--meta_tsv", required=True, help="Bulk sample metadata TSV")
    ap.add_argument("--out_dir", required=True)
    ap.add_argument("--sample_key", default="sample_id")
    ap.add_argument("--patient_key", default="patient_id")
    ap.add_argument("--disease_key", default="disease_broad5")
    ap.add_argument("--n_splits", type=int, default=5)
    ap.add_argument("--shuffle", action="store_true")
    ap.add_argument("--random_state", type=int, default=1234)
    ap.add_argument("--exclude_diseases", nargs="*", default=None)
    ap.add_argument("--min_samples_for_eval", type=int, default=5)
    args = ap.parse_args()

    ensure_dir(args.out_dir)

    print("[READ]", args.meta_tsv)
    meta = pd.read_csv(args.meta_tsv, sep="\t")

    need_cols = [args.sample_key, args.disease_key]
    if args.patient_key is not None:
        need_cols.append(args.patient_key)

    miss = [c for c in need_cols if c not in meta.columns]
    if miss:
        raise ValueError(f"Metadata missing required columns: {miss}")

    for c in need_cols:
        meta[c] = meta[c].map(clean_str)

    keep = meta[args.sample_key].notna() & meta[args.disease_key].notna()
    if args.patient_key is not None:
        keep &= meta[args.patient_key].notna()

    if args.exclude_diseases:
        keep &= ~meta[args.disease_key].isin(args.exclude_diseases)

    meta = meta.loc[keep].copy()
    if meta.shape[0] == 0:
        raise RuntimeError("No samples left after filtering.")

    sample_df = build_sample_table(
        meta=meta,
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

    print("[INFO] primary benchmark diseases:", primary_diseases)
    print("[INFO] reference-only diseases:", reference_only_diseases)

    if len(primary_diseases) == 0:
        raise RuntimeError("No diseases satisfy min_samples_for_eval.")

    primary_sample_df = sample_df[sample_df[args.disease_key].isin(primary_diseases)].copy()
    if primary_sample_df.shape[0] < args.n_splits:
        raise RuntimeError(
            f"Only {primary_sample_df.shape[0]} primary-benchmark samples available, "
            f"smaller than n_splits={args.n_splits}."
        )

    if args.patient_key is not None and args.patient_key in primary_sample_df.columns:
        primary_groups = primary_sample_df[args.patient_key].astype(str).values
    else:
        primary_groups = primary_sample_df[args.sample_key].astype(str).values

    splitter = StratifiedGroupKFold(
        n_splits=args.n_splits,
        shuffle=args.shuffle,
        random_state=args.random_state,
    )

    primary_sample_ids = primary_sample_df[args.sample_key].astype(str).values
    primary_labels = primary_sample_df[args.disease_key].astype(str).values

    split_indices = list(
        splitter.split(
            X=np.arange(primary_sample_df.shape[0]),
            y=primary_labels,
            groups=primary_groups,
        )
    )

    reference_only_samples = (
        sample_df.loc[sample_df[args.disease_key].isin(reference_only_diseases), args.sample_key]
        .astype(str)
        .tolist()
    )
    reference_only_sample_set = set(reference_only_samples)

    kfold_sample_ids: List[Tuple[np.ndarray, np.ndarray]] = []
    summary_rows: List[Dict] = []
    sample_rows: List[pd.DataFrame] = []

    for split_id, (train_idx, test_idx) in enumerate(split_indices):
        train_primary_samples = set(primary_sample_ids[train_idx].tolist())
        test_primary_samples = set(primary_sample_ids[test_idx].tolist())

        train_sample_set = set(train_primary_samples) | reference_only_sample_set
        test_sample_set = set(test_primary_samples)

        if len(train_sample_set.intersection(test_sample_set)) > 0:
            raise RuntimeError(f"Split {split_id}: overlapping train/test samples detected.")

        kfold_sample_ids.append((
            np.array(sorted(train_sample_set), dtype=object),
            np.array(sorted(test_sample_set), dtype=object),
        ))

        summary_rows.append({
            "split": split_id,
            "n_train_samples": len(train_sample_set),
            "n_test_samples": len(test_sample_set),
            "n_reference_only_samples_in_train": len(reference_only_sample_set),
            "primary_benchmark_diseases": ";".join(primary_diseases),
            "reference_only_diseases": ";".join(reference_only_diseases),
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

    pkl_fp = os.path.join(args.out_dir, "K_FOLD_sampleID.pkl")
    with open(pkl_fp, "wb") as handle:
        pickle.dump(kfold_sample_ids, handle, protocol=4)

    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(os.path.join(args.out_dir, "K_FOLD_summary.tsv"), sep="\t", index=False)

    sample_split_df = pd.concat(sample_rows, axis=0, ignore_index=True).drop_duplicates()
    sample_split_df.to_csv(os.path.join(args.out_dir, "K_FOLD_sample_split.tsv"), sep="\t", index=False)

    disease_support.to_csv(os.path.join(args.out_dir, "K_FOLD_disease_support.tsv"), sep="\t", index=False)

    metadata = {
        "n_splits": args.n_splits,
        "shuffle": bool(args.shuffle),
        "random_state": args.random_state,
        "sample_key": args.sample_key,
        "patient_key": args.patient_key,
        "disease_key": args.disease_key,
        "min_samples_for_eval": args.min_samples_for_eval,
        "primary_benchmark_diseases": primary_diseases,
        "reference_only_diseases": reference_only_diseases,
        "excluded_diseases": args.exclude_diseases or [],
    }
    with open(os.path.join(args.out_dir, "K_FOLD_metadata.json"), "w", encoding="utf-8") as f:
        json.dump(metadata, f, ensure_ascii=False, indent=2)

    print("\n[DONE]")
    print("pickle:", pkl_fp)
    print(summary_df.to_string(index=False))


if __name__ == "__main__":
    main()