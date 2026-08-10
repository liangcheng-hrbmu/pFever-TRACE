#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os, argparse
import scanpy as sc

def ensure_dir(p): os.makedirs(p, exist_ok=True)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in_h5ad", required=True)
    ap.add_argument("--out_dir", required=True)
    ap.add_argument("--sample_key", default="sample_id")
    args = ap.parse_args()

    ensure_dir(args.out_dir)
    adata = sc.read_h5ad(args.in_h5ad)

    if args.sample_key not in adata.obs.columns:
        raise ValueError(f"obs 缺少 sample_key: {args.sample_key}")

    s = adata.obs[args.sample_key].astype(str)
    for sid in sorted(s.unique()):
        sub = adata[s == sid].copy()
        out = os.path.join(args.out_dir, f"{sid}.h5ad")
        sub.write_h5ad(out, compression="gzip")
        print("[WRITE]", out, "n_cells=", sub.n_obs, "n_genes=", sub.n_vars)

    print("[DONE] n_samples=", s.nunique(), "out_dir=", args.out_dir)

if __name__ == "__main__":
    main()
