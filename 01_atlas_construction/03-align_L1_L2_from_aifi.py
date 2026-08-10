#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Align AIFI hierarchical labels for a merged AnnData object:

Input:
  data/analysis/03_integrated_umap_scvi_v4/merged_paperFiltered.h5ad

Goal:
  1) Fix small inconsistencies between aifi_L2 and aifi_L3:
     - Force each L3 to map to a single L2 by majority vote (L3 -> L2 mode).
  2) Define "Level 1" (major lineages) exactly as in your figure legend:
     B, Plasma, pDC, DC, Mono,
     T CD4 naive, T CD4 non-naive,
     T CD8 naive, T CD8 non-naive,
     UTC, ILC, Cycling cells,
     Progenitors, Platelets, RBC
     - Map from (aligned) aifi_L2 to L1 via a precise dictionary.
  3) Final outputs:
     - aifi_L1_paper : final Level1 (forced to be unique per L3)
     - aifi_L2_paper : final Level2 (use aifi_L3)
     - also keep intermediate: aifi_L2_aligned, aifi_L1_aligned
  4) Write a new .h5ad and export QC tables to qc_align_paper/
"""

import os
import argparse
import pandas as pd
import scanpy as sc

# =========================
# L2 -> Level1 mapping (STRICT, matches your legend names)
# =========================
L2_TO_L1 = {
    # Mono
    "CD14 monocyte": "Mono",
    "CD16 monocyte": "Mono",

    # B lineage
    "Naive B cell": "B",
    "Memory B cell": "B",
    "Transitional B cell": "B",
    "Effector B cell": "B",

    # Plasma
    "Plasma cell": "Plasma",

    # DC / pDC
    "cDC2": "DC",
    "pDC": "pDC",

    # ILC / NK
    "ILC": "ILC",
    "CD56bright NK cell": "ILC",
    "CD56dim NK cell": "ILC",

    # Cycling cells (separate L1 class, per your legend)
    "Proliferating NK cell": "Cycling cells",
    "Proliferating T cell": "Cycling cells",

    # T cells (split naive vs non-naive)
    "Naive CD4 T cell": "T CD4 naive",
    "Memory CD4 T cell": "T CD4 non-naive",
    "Treg": "T CD4 non-naive",

    "Naive CD8 T cell": "T CD8 naive",
    "Memory CD8 T cell": "T CD8 non-naive",
    "CD8aa": "T CD8 non-naive",

    # UTC (unconventional T)
    "MAIT": "UTC",
    "gdT": "UTC",
    "DN T cell": "UTC",

    # Others
    "Progenitor cell": "Progenitors",
    "Erythrocyte": "RBC",
    "Platelet": "Platelets",
}

FINAL_L1_ORDER = [
    "B",
    "Plasma",
    "pDC",
    "DC",
    "Mono",
    "T CD4 naive",
    "T CD4 non-naive",
    "T CD8 naive",
    "T CD8 non-naive",
    "UTC",
    "ILC",
    "Cycling cells",
    "Progenitors",
    "Platelets",
    "RBC",
]


def majority_map(source: pd.Series, target: pd.Series) -> pd.Series:
    """
    Majority vote mapping:
      for each category in source -> the most frequent category in target.
    Returns:
      Series(index=source_category, value=target_mode)
    """
    tab = pd.crosstab(source.astype(str), target.astype(str))
    return tab.idxmax(axis=1)


def main():
    ap = argparse.ArgumentParser()

    ap.add_argument(
        "--in_h5ad",
        default="data/analysis/03_integrated_umap_scvi_v4/merged_paperFiltered.h5ad",
        help="Input h5ad (paperFiltered merged object).",
    )
    ap.add_argument(
        "--out_h5ad",
        default="data/analysis/03_integrated_umap_scvi_v4/merged_paperFiltered.paperL1L2.h5ad",
        help="Output h5ad with aligned Level1/Level2 labels.",
    )

    ap.add_argument("--col_L2", default="aifi_L2", help="Existing L2 column in obs.")
    ap.add_argument("--col_L3", default="aifi_L3", help="Existing L3 column in obs.")

    # Final columns you will use downstream
    ap.add_argument("--out_L1", default="aifi_L1_paper", help="Final Level1 column name.")
    ap.add_argument("--out_L2", default="aifi_L2_paper", help="Final Level2 column name (will be aifi_L3).")

    ap.add_argument(
        "--qc_dir",
        default="data/analysis/03_integrated_umap_scvi_v4/qc_align_paper",
        help="Directory to write QC tables.",
    )

    args = ap.parse_args()
    os.makedirs(args.qc_dir, exist_ok=True)

    print("[READ]", args.in_h5ad)
    adata = sc.read_h5ad(args.in_h5ad)

    if args.col_L2 not in adata.obs.columns:
        raise SystemExit(f"[FATAL] missing obs['{args.col_L2}']")
    if args.col_L3 not in adata.obs.columns:
        raise SystemExit(f"[FATAL] missing obs['{args.col_L3}']")

    L2_raw = adata.obs[args.col_L2].astype(str).fillna("NA")
    L3 = adata.obs[args.col_L3].astype(str).fillna("NA")

    # =========================
    # 1) Align L2 by L3 (L3 -> L2 majority vote)
    # =========================
    map_L3_to_L2 = majority_map(L3, L2_raw)
    L2_aligned = L3.map(map_L3_to_L2)

    adata.obs["aifi_L2_aligned"] = pd.Categorical(L2_aligned)

    changed = int((L2_aligned != L2_raw).sum())
    print(f"[INFO] L2 aligned by L3: changed {changed} / {adata.n_obs}")

    # =========================
    # 2) Map aligned L2 -> paper Level1 (strict mapping)
    # =========================
    L2a = adata.obs["aifi_L2_aligned"].astype(str)

    unknown_L2 = sorted(set(L2a.unique()) - set(L2_TO_L1.keys()))
    if unknown_L2:
        print("[WARN] Unknown L2 labels not in mapping (will set to UTC):")
        print(unknown_L2)

    L1_aligned = L2a.map(lambda x: L2_TO_L1.get(x, "UTC"))
    adata.obs["aifi_L1_aligned"] = pd.Categorical(L1_aligned, categories=FINAL_L1_ORDER, ordered=True)

    # =========================
    # 3) Force L3 -> L1 uniqueness (L3 -> L1 majority vote), final Level1
    # =========================
    map_L3_to_L1 = majority_map(L3, adata.obs["aifi_L1_aligned"].astype(str))
    L1_final = L3.map(map_L3_to_L1)

    adata.obs[args.out_L1] = pd.Categorical(L1_final, categories=FINAL_L1_ORDER, ordered=True)
    adata.obs[args.out_L2] = pd.Categorical(L3)

    # =========================
    # 4) QC outputs
    # =========================
    # Mapping tables
    pd.Series(map_L3_to_L2, name="L2_mode").to_csv(os.path.join(args.qc_dir, "map_L3_to_L2_mode.csv"))
    pd.Series(map_L3_to_L1, name="L1_mode").to_csv(os.path.join(args.qc_dir, "map_L3_to_L1_mode.csv"))

    # Crosstabs
    pd.crosstab(L3, L2_raw).to_csv(os.path.join(args.qc_dir, "crosstab_L3_to_L2_raw.csv"))
    pd.crosstab(L3, adata.obs["aifi_L2_aligned"].astype(str)).to_csv(
        os.path.join(args.qc_dir, "crosstab_L3_to_L2_aligned.csv")
    )
    pd.crosstab(L3, adata.obs["aifi_L1_aligned"].astype(str)).to_csv(
        os.path.join(args.qc_dir, "crosstab_L3_to_L1_aligned.csv")
    )
    pd.crosstab(L3, adata.obs[args.out_L1].astype(str)).to_csv(
        os.path.join(args.qc_dir, "crosstab_L3_to_L1_final.csv")
    )

    # Counts
    adata.obs[args.out_L1].value_counts(dropna=False).to_csv(os.path.join(args.qc_dir, "L1_final_counts.csv"))
    adata.obs[args.out_L2].value_counts(dropna=False).to_csv(os.path.join(args.qc_dir, "L2_final_counts.csv"))

    # Identify L3 that originally map to multiple L2 (inconsistency sources)
    tab_before = pd.crosstab(L3, L2_raw)
    n_L2_per_L3 = (tab_before > 0).sum(axis=1).sort_values(ascending=False)
    n_L2_per_L3.to_csv(os.path.join(args.qc_dir, "n_L2_per_L3_before.csv"))

    # Unknown L2 list (if any)
    pd.Series(unknown_L2, name="unknown_L2").to_csv(os.path.join(args.qc_dir, "unknown_L2_labels.csv"), index=False)

    # =========================
    # 5) Write output h5ad
    # =========================
    print("[WRITE]", args.out_h5ad)
    adata.write_h5ad(args.out_h5ad, compression="gzip")

    print("\n[SUMMARY]")
    print("Input:", args.in_h5ad)
    print("Output:", args.out_h5ad)
    print("QC dir:", args.qc_dir)
    print("Final cols:", args.out_L1, args.out_L2)
    print("Added cols:", "aifi_L2_aligned", "aifi_L1_aligned", args.out_L1, args.out_L2)


if __name__ == "__main__":
    main()
