#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import glob
import argparse
import zipfile
import numpy as np
import pandas as pd
import scanpy as sc
import joblib

SUFFIX = "_adataMerged_SPECTRAgenes.log1p.h5ad"


def ensure_dir(p):
    os.makedirs(p, exist_ok=True)


def celltype_from_h5ad(path: str) -> str:
    base = os.path.basename(path)
    if base.endswith(SUFFIX):
        return base[:-len(SUFFIX)]
    return os.path.splitext(base)[0]


def normalize_split_indices(split_obj):
    if isinstance(split_obj, np.lib.npyio.NpzFile):
        if "test_idx" in split_obj.files:
            return np.asarray(split_obj["test_idx"], dtype=int)
        raise KeyError("split npz missing test_idx")
    return np.asarray(split_obj, dtype=int)


def main():
    ap = argparse.ArgumentParser(description="Export corrected test predictions from fitted disease models.")
    ap.add_argument("--input_dir", required=True)
    ap.add_argument("--model_dir", required=True)
    ap.add_argument("--out_dir", required=True)
    ap.add_argument("--disease_key", default="disease")
    args = ap.parse_args()

    ensure_dir(args.out_dir)
    h5ads = sorted(glob.glob(os.path.join(args.input_dir, f"*{SUFFIX}")))
    if not h5ads:
        raise FileNotFoundError(f"No files found in {args.input_dir}")

    for h5ad_fp in h5ads:
        celltype = celltype_from_h5ad(h5ad_fp)
        model_fp = os.path.join(args.model_dir, f"{celltype}.xgb_model.joblib")
        split_fp = os.path.join(args.model_dir, f"{celltype}.split_indices.npz")
        if not (os.path.exists(model_fp) and os.path.exists(split_fp)):
            print(f"[SKIP] {celltype}: missing model or split")
            continue

        print(f"[RUN] {celltype}")
        adata = sc.read_h5ad(h5ad_fp)
        bundle = joblib.load(model_fp)
        model = bundle["model"]
        le = bundle["label_encoder"]

        if args.disease_key not in adata.obs.columns:
            raise ValueError(f"{celltype}: missing obs column {args.disease_key}")

        y_all = adata.obs[args.disease_key].astype(str).values
        keep_classes = set(map(str, le.classes_.tolist()))
        mask = np.array([y in keep_classes for y in y_all], dtype=bool)
        if mask.sum() == 0:
            print(f"[SKIP] {celltype}: no cells match trained classes")
            continue

        adata_sub = adata[mask].copy()
        orig_idx = np.where(mask)[0]
        remap = {old: new for new, old in enumerate(orig_idx)}

        split_obj = np.load(split_fp)
        test_idx_raw = normalize_split_indices(split_obj)
        test_idx = np.array([remap[i] for i in test_idx_raw if i in remap], dtype=int)
        if test_idx.size == 0:
            print(f"[SKIP] {celltype}: no test indices remain after class filtering")
            continue

        X = adata_sub.X
        if hasattr(X, "toarray"):
            X = X.toarray()
        X = np.asarray(X)
        y_true = adata_sub.obs[args.disease_key].astype(str).values

        X_test = X[test_idx]
        y_true_test = y_true[test_idx]
        y_true_code = le.transform(y_true_test)
        y_pred_code = model.predict(X_test)
        y_pred = le.inverse_transform(y_pred_code)

        out_df = pd.DataFrame({
            "celltype": celltype,
            "test_index_in_filtered": test_idx,
            "true_label": y_true_test,
            "pred_label": y_pred,
            "true_code": y_true_code,
            "pred_code": y_pred_code,
        })

        out_zip = os.path.join(args.out_dir, f"run1_{celltype}_pred_test.zip")
        out_npy = os.path.join(args.out_dir, f"run1_{celltype}_test_idxs.npy")
        out_tsv = os.path.join(args.out_dir, f"run1_{celltype}_pred_test.tsv")

        out_df.to_csv(out_tsv, sep="\t", index=False)
        with zipfile.ZipFile(out_zip, "w", compression=zipfile.ZIP_DEFLATED) as zf:
            zf.write(out_tsv, arcname=os.path.basename(out_tsv))
        np.save(out_npy, test_idx)
        os.remove(out_tsv)

        print(f"[WRITE] {out_zip}")
        print(f"[WRITE] {out_npy}")

    print("\n[DONE]")


if __name__ == "__main__":
    main()
