#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import argparse
import scanpy as sc
import scvi
import numpy as np
import torch

def ensure_dir(p):
    os.makedirs(p, exist_ok=True)

def main():
    ap = argparse.ArgumentParser()

    ap.add_argument("--in_h5ad", required=True)
    ap.add_argument("--out_dir", required=True)

    ap.add_argument("--counts_layer", default="counts")
    ap.add_argument("--batch_key", default="sample_id")      # 你坚持用 sample_id
    ap.add_argument("--labels_key", default="aifi_L1_paper") # scANVI label
    ap.add_argument("--unlabeled_category", default="Unknown")

    ap.add_argument("--cat_cov", nargs="*", default=["dataset_id", "sex", "binned_age"])
    ap.add_argument("--cont_cov", nargs="*", default=[])

    ap.add_argument("--n_latent", type=int, default=30)
    ap.add_argument("--max_epochs_scvi", type=int, default=200)
    ap.add_argument("--max_epochs_scanvi", type=int, default=50)

    ap.add_argument("--seed", type=int, default=5)
    ap.add_argument("--num_workers", type=int, default=8)

    args = ap.parse_args()
    ensure_dir(args.out_dir)

    scvi.settings.seed = args.seed
    scvi.settings.dl_num_workers = args.num_workers
    torch.set_float32_matmul_precision("high")

    print("[READ]", args.in_h5ad)
    adata = sc.read_h5ad(args.in_h5ad)

    # basic checks
    for k in [args.batch_key, args.labels_key]:
        if k not in adata.obs.columns:
            raise ValueError(f"[ERROR] obs missing: {k}")

    for k in args.cat_cov:
        if k not in adata.obs.columns:
            raise ValueError(f"[ERROR] obs missing categorical covariate: {k}")

    # ensure categories
    for k in [args.batch_key, args.labels_key] + args.cat_cov:
        adata.obs[k] = adata.obs[k].astype("category")

    # setup anndata
    scvi.model.SCVI.setup_anndata(
        adata,
        layer=args.counts_layer,
        batch_key=args.batch_key,
        labels_key=args.labels_key,
        categorical_covariate_keys=args.cat_cov if len(args.cat_cov) else None,
        continuous_covariate_keys=args.cont_cov if len(args.cont_cov) else None,
    )

    # 1) Train SCVI
    print("[SCVI] training ...")
    scvi_model = scvi.model.SCVI(adata, n_latent=args.n_latent)
    scvi_model.train(max_epochs=args.max_epochs_scvi)

    scvi_dir = os.path.join(args.out_dir, "scvi_model")
    ensure_dir(scvi_dir)
    scvi_model.save(scvi_dir, overwrite=True)
    print("[SCVI] saved:", scvi_dir)

    # store latent in adata for quick QC
    adata.obsm["X_scVI"] = scvi_model.get_latent_representation()
    adata.write_h5ad(os.path.join(args.out_dir, "adata_with_scvi_latent.h5ad"), compression="gzip")

    # 2) Train SCANVI from SCVI
    print("[SCANVI] training ...")
    scanvi_model = scvi.model.SCANVI.from_scvi_model(
        scvi_model,
        unlabeled_category=args.unlabeled_category,
    )
    scanvi_model.train(max_epochs=args.max_epochs_scanvi)

    scanvi_dir = os.path.join(args.out_dir, "scanvi_model")
    ensure_dir(scanvi_dir)
    scanvi_model.save(scanvi_dir, overwrite=True)
    print("[SCANVI] saved:", scanvi_dir)

    # store scanvi latent as well
    adata.obsm["X_scANVI"] = scanvi_model.get_latent_representation()
    adata.write_h5ad(os.path.join(args.out_dir, "adata_with_scanvi_latent.h5ad"), compression="gzip")

    print("[DONE]")

if __name__ == "__main__":
    main()
