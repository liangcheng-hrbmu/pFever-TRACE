#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os, argparse, re
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import scipy.sparse as sp


def ensure_dir(p):
    os.makedirs(p, exist_ok=True)


def clean_string_series(x: pd.Series) -> pd.Series:
    """strip + collapse whitespaces; return object dtype for old anndata compatibility"""
    s = x.astype(object)

    def _clean(v):
        if v is None or (isinstance(v, float) and np.isnan(v)):
            return np.nan
        v = str(v)
        v = re.sub(r"\s+", " ", v).strip()
        return v

    return s.map(_clean).astype(object)


def sanitize_obs_for_h5ad(adata: ad.AnnData):
    """
    Make obs safe for old anndata/h5py:
    - numeric/bool/category keep
    - object -> scalar str/NA
    """
    for c in list(adata.obs.columns):
        s = adata.obs[c]

        if pd.api.types.is_bool_dtype(s) or pd.api.types.is_numeric_dtype(s) or isinstance(s.dtype, pd.CategoricalDtype):
            continue

        if pd.api.types.is_object_dtype(s):
            def _to_scalar(v):
                if v is None or (isinstance(v, float) and np.isnan(v)):
                    return np.nan
                if isinstance(v, (list, dict, tuple, set)):
                    return str(v)
                if hasattr(v, "item") and callable(v.item):
                    try:
                        return v.item()
                    except Exception:
                        pass
                return str(v)

            adata.obs[c] = s.map(_to_scalar).astype(object)
        else:
            adata.obs[c] = s.astype(str).astype(object)

    return adata


def read_gene_list(path):
    genes = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            g = line.strip()
            if not g:
                continue
            genes.append(g)
    # keep order unique
    seen = set()
    out = []
    for g in genes:
        if g not in seen:
            out.append(g)
            seen.add(g)
    return out


def main():
    ap = argparse.ArgumentParser(description="Prepare downstream h5ad: drop celltypes + subset union genes + keep clinical keys")
    ap.add_argument("--in_h5ad", required=True)
    ap.add_argument("--union_genes", required=True, help="gene_universe_union.txt")
    ap.add_argument("--out_h5ad", required=True)

    ap.add_argument("--counts_layer", default="counts")  # if exists, will be subset & sparsified
    ap.add_argument("--celltype_key", default="aifi_L1_paper")
    ap.add_argument("--drop_celltypes", nargs="*", default=["RBC", "Platelets", "Cycling_cells", "Progenitors"])

    ap.add_argument("--required_obs", nargs="*", default=[
        "dataset_id", "sex", "binned_age", "aifi_L1_paper", "sample_id", "patient_id"
    ])

    ap.add_argument("--write_reports", action="store_true",
                    help="write dropped celltypes counts + gene missing list next to out_h5ad")

    args = ap.parse_args()
    out_dir = os.path.dirname(os.path.abspath(args.out_h5ad)) or "."
    ensure_dir(out_dir)

    print("[READ]", args.in_h5ad)
    adata = sc.read_h5ad(args.in_h5ad)

    # ---- basic checks ----
    missing_obs = [k for k in args.required_obs if k not in adata.obs.columns]
    if missing_obs:
        raise SystemExit(f"[FATAL] missing required obs columns: {missing_obs}")

    # clean strings for key obs
    for k in args.required_obs:
        if pd.api.types.is_object_dtype(adata.obs[k]) or str(adata.obs[k].dtype) == "category":
            adata.obs[k] = clean_string_series(adata.obs[k])

    # ---- drop celltypes ----
    ct = clean_string_series(adata.obs[args.celltype_key])
    drop_set = set(args.drop_celltypes)
    mask_keep = ~ct.isin(drop_set)
    n_before = adata.n_obs
    adata = adata[mask_keep].copy()
    n_after = adata.n_obs
    print(f"[INFO] drop_celltypes by {args.celltype_key}: removed {n_before - n_after} cells, kept {n_after}")

    if args.write_reports:
        # counts dropped/kept summary
        df_ct = pd.DataFrame({args.celltype_key: ct})
        df_ct["kept"] = mask_keep
        summ = (
            df_ct.groupby([args.celltype_key, "kept"])
            .size()
            .reset_index(name="n_cells")
            .sort_values(["kept", "n_cells"], ascending=[False, False])
        )
        summ_path = os.path.join(out_dir, "report_celltype_drop_summary.tsv")
        summ.to_csv(summ_path, sep="\t", index=False)
        print("[WRITE]", summ_path)

    # ---- subset genes to union ----
    union = read_gene_list(args.union_genes)
    union_set = set(union)
    genes_in = list(adata.var_names)

    present = [g for g in union if g in adata.var_names]
    missing = [g for g in union if g not in adata.var_names]

    if len(present) == 0:
        raise SystemExit("[FATAL] none of union genes are present in adata.var_names. Check gene symbols / IDs.")

    print(f"[INFO] union genes: total={len(union)} present={len(present)} missing={len(missing)}")
    adata = adata[:, present].copy()

    # ---- keep layers consistent ----
    if args.counts_layer in adata.layers:
        X = adata.layers[args.counts_layer]
        if not sp.issparse(X):
            adata.layers[args.counts_layer] = sp.csr_matrix(X)

    # ---- final sanitize & write ----
    adata = sanitize_obs_for_h5ad(adata)
    print("[WRITE]", args.out_h5ad, "shape=", adata.shape)
    adata.write_h5ad(args.out_h5ad, compression="gzip")

    if args.write_reports:
        miss_path = os.path.join(out_dir, "report_union_genes_missing.txt")
        with open(miss_path, "w", encoding="utf-8") as f:
            f.write("\n".join(missing) + "\n")
        print("[WRITE]", miss_path)

        keep_path = os.path.join(out_dir, "report_union_genes_kept.txt")
        with open(keep_path, "w", encoding="utf-8") as f:
            f.write("\n".join(present) + "\n")
        print("[WRITE]", keep_path)

    print("[DONE]")


if __name__ == "__main__":
    main()
