#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
from pathlib import Path


def clean_str(x):
    if pd.isna(x):
        return None
    x = str(x).strip()
    if x == "" or x == "-" or x.lower() == "nan":
        return None
    return x


def first_nonnull(series):
    vals = [x for x in series if pd.notna(x) and str(x).strip() != ""]
    return vals[0] if len(vals) > 0 else None


def normalize_disease_name(x: str):
    if x is None:
        return None

    s = str(x).strip()
    s_lower = s.lower()

    if "healthy" in s_lower or "control" in s_lower:
        return "Healthy"

    if "kawasaki" in s_lower or s_lower == "kd":
        return "Kawasaki disease (KD)"
    if "mis-c" in s_lower or "multisystem inflammatory syndrome" in s_lower:
        return "Multisystem inflammatory syndrome in children (MIS-C)"

    if "juvenile idiopathic arthritis" in s_lower or s_lower == "jia":
        return "juvenile idiopathic arthritis (JIA)"
    if "lupus" in s_lower or s_lower == "sle":
        return "human systemic lupus erythematosus (SLE)"

    if "influenza" in s_lower:
        return "Influenza A virus (IAV)"
    if "respiratory syncytial virus" in s_lower or "rsv" in s_lower:
        return "Respiratory Syncytial Virus Infection"
    if "adenovirus" in s_lower:
        return "adenovirus infection"
    if "rhinovirus" in s_lower:
        return "rhinovirus infection"

    if "meningococcal" in s_lower:
        return "meningococcal infection"
    if "pneumococcal" in s_lower:
        return "pneumococcal infection"
    if "staphylococcus aureus" in s_lower:
        return "Staphylococcus aureus infection"
    if "streptococcal" in s_lower:
        return "group A streptococcal infection"
    if "escherichia coli" in s_lower or "e. coli" in s_lower:
        return "escherichia coli infection"
    if "tuberculosis" in s_lower or s_lower == "tb":
        return "Tuberculosis"
    if "bacterial" in s_lower:
        return "bacterial infection"
    if "sepsis" in s_lower:
        return "Pediatric Sepsis"

    if "malaria" in s_lower:
        return "malaria"

    return s


def map_to_broad(x: str):
    if x is None:
        return None

    s = str(x).strip().lower()

    if s == "healthy":
        return "Healthy"

    if (
        "respiratory syncytial virus" in s or
        "adenovirus" in s or
        "rhinovirus" in s or
        "influenza" in s or
        s == "viral infection"
    ):
        return "Viral infection"

    if (
        "meningococcal" in s or
        "pneumococcal" in s or
        "staphylococcus aureus" in s or
        "streptococcal" in s or
        "escherichia coli" in s or
        "tuberculosis" in s or
        "sepsis" in s or
        s == "bacterial infection"
    ):
        return "Bacterial infection"

    if "kawasaki" in s or "mis-c" in s or "multisystem inflammatory syndrome" in s:
        return "Hyperinflammatory disease"

    if (
        "juvenile idiopathic arthritis" in s or
        s == "jia" or
        "lupus" in s or
        s == "sle"
    ):
        return "Autoimmune disease"

    if "malaria" in s:
        return "Malaria"

    return None


def main():
    in_csv = "data/analysis/09_bulk_deconvolution/E-MTAB-11671_samples.csv"
    out_dir = Path("data/analysis/09_bulk_deconvolution")
    out_dir.mkdir(parents=True, exist_ok=True)

    out_basic = out_dir / "bulk_sample_metadata.tsv"
    out_full = out_dir / "bulk_sample_metadata_full.tsv"
    out_mapped = out_dir / "bulk_sample_metadata_mapped.tsv"

    print(f"[READ] {in_csv}")
    df = pd.read_csv(in_csv)

    print("[INFO] columns:")
    print(df.columns.tolist())

    possible_cols = [
        "Title",
        "factorvalue[disease]",
        "characteristics[disease]",
        "characteristics[age]",
        "characteristics[sex]",
        "characteristics[developmental stage]"
    ]
    for c in possible_cols:
        if c in df.columns:
            df[c] = df[c].map(clean_str)

    disease = df["factorvalue[disease]"] if "factorvalue[disease]" in df.columns else pd.Series([None] * len(df))
    if "characteristics[disease]" in df.columns:
        disease = disease.fillna(df["characteristics[disease]"])
    df["disease"] = disease

    if "Title" not in df.columns:
        raise ValueError("Column 'Title' not found in input CSV.")

    df["sample_id"] = df["Title"].map(clean_str)

    keep_cols = ["sample_id", "disease"]
    extra_cols = [
        c for c in [
            "characteristics[age]",
            "characteristics[sex]",
            "characteristics[developmental stage]"
        ] if c in df.columns
    ]

    out_df = df[keep_cols + extra_cols].copy()
    out_df = out_df[out_df["sample_id"].notna()].copy()
    out_df = out_df.drop_duplicates()

    out_df = (
        out_df.groupby("sample_id", as_index=False)
        .agg({c: first_nonnull for c in out_df.columns if c != "sample_id"})
    )

    out_df[["sample_id", "disease"]].to_csv(out_basic, sep="\t", index=False)
    out_df.to_csv(out_full, sep="\t", index=False)

    mapped_df = out_df.copy()
    mapped_df["disease_clean"] = mapped_df["disease"].map(normalize_disease_name)
    mapped_df["disease_broad"] = mapped_df["disease_clean"].map(map_to_broad)

    mapped_df.to_csv(out_mapped, sep="\t", index=False)

    print("\n[SAVED]")
    print(out_basic)
    print(out_full)
    print(out_mapped)

    print("\n[INFO] unique sample count:")
    print(mapped_df.shape[0])

    print("\n[INFO] raw disease counts:")
    print(mapped_df["disease"].value_counts(dropna=False))

    print("\n[INFO] mapped disease_clean counts:")
    print(mapped_df["disease_clean"].value_counts(dropna=False))

    print("\n[INFO] mapped disease_broad counts:")
    print(mapped_df["disease_broad"].value_counts(dropna=False))

    print("\n[DEBUG] unique disease_clean -> disease_broad:")
    tmp = (
        mapped_df[["disease_clean", "disease_broad"]]
        .drop_duplicates()
        .sort_values(["disease_broad", "disease_clean"], na_position="last")
    )
    print(tmp.to_string(index=False))

    print("\n[PREVIEW]")
    print(mapped_df.head(10))


if __name__ == "__main__":
    main()