#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os, argparse, gzip
import numpy as np
import pandas as pd

import scanpy as sc
import scipy.sparse as sp
import decoupler


def ensure_dir(p):
    os.makedirs(p, exist_ok=True)


def clean_str(s: pd.Series) -> pd.Series:
    return s.astype(str).str.strip()


def norm_ct_label(x: str) -> str:
    """统一细胞类型写法：小写 + 空格->下划线，便于匹配 drop list。"""
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return ""
    x = str(x).strip().lower()
    x = x.replace(" ", "_")
    return x


def mode1(x: pd.Series):
    """取众数（若全 NA 返回 NA）"""
    x = x.dropna()
    if x.empty:
        return np.nan
    vc = x.value_counts()
    return vc.index[0]


def write_counts_tsv_gz_gene_by_sample(X_gene_by_sample, genes, sample_ids, out_path):
    """
    写 edgeR 友好格式：行=gene，列=sample。
    X_gene_by_sample: CSR (n_genes x n_samples) 或 ndarray
    """
    if sp.issparse(X_gene_by_sample):
        X = X_gene_by_sample.tocsr()
        get_row = lambda i: X.getrow(i).toarray().ravel()
    else:
        X = np.asarray(X_gene_by_sample)
        get_row = lambda i: X[i, :]

    with gzip.open(out_path, "wt", encoding="utf-8") as f:
        # header
        f.write("gene\t" + "\t".join(sample_ids) + "\n")
        for i, g in enumerate(genes):
            row = get_row(i)
            # counts 用整数更稳
            row = np.asarray(row).ravel()
            # 可能是 float（sum 后），转 int64
            row = np.round(row).astype(np.int64)
            f.write(str(g) + "\t" + "\t".join(map(str, row.tolist())) + "\n")


def main():
    ap = argparse.ArgumentParser(
        description="Make pseudobulk for edgeR: patient_id x Level1, min_cells=20, sum counts. Output TSV (gene x sample) + meta."
    )

    ap.add_argument("--in_h5ad", required=True, help="Input h5ad (建议用 *.clinical.h5ad)")
    ap.add_argument("--out_dir", required=True)

    ap.add_argument("--counts_layer", default="counts")
    ap.add_argument("--patient_key", default="patient_id", help="sample_col in decoupler.get_pseudobulk")
    ap.add_argument("--celltype_key", default="aifi_L1_paper", help="groups_col in decoupler.get_pseudobulk")
    ap.add_argument("--disease_key", default="disease_Level_3", help="for meta; optional but recommended")
    ap.add_argument("--study_key", default="dataset_id", help="for meta; optional")

    ap.add_argument("--sex_key", default="sex")
    ap.add_argument("--binned_age_key", default="binned_age")
    ap.add_argument("--age_key", default="age_years")

    ap.add_argument("--min_cells", type=int, default=20, help="min cells per (patient, celltype) group")

    # 你要求排除的四类（兼容空格/下划线）
    ap.add_argument(
        "--drop_celltypes",
        nargs="*",
        default=["RBC", "Platelets", "Cycling_cells", "Progenitors"],
        help="Cell types to exclude before pseudobulk"
    )

    ap.add_argument("--out_counts", default="pseudobulk_counts.tsv.gz")
    ap.add_argument("--out_meta", default="pseudobulk_meta.tsv")

    args = ap.parse_args()
    ensure_dir(args.out_dir)

    print("[READ]", args.in_h5ad)
    adata = sc.read_h5ad(args.in_h5ad)

    # 基础检查
    for k in [args.patient_key, args.celltype_key]:
        if k not in adata.obs.columns:
            raise SystemExit(f"[FATAL] missing obs['{k}'] in h5ad")

    if args.counts_layer not in adata.layers:
        raise SystemExit(f"[FATAL] missing layers['{args.counts_layer}'] in h5ad")

    # counts 稀疏化
    if not sp.issparse(adata.layers[args.counts_layer]):
        adata.layers[args.counts_layer] = sp.csr_matrix(adata.layers[args.counts_layer])

    # 清洗关键列
    adata.obs[args.patient_key] = clean_str(adata.obs[args.patient_key])
    adata.obs[args.celltype_key] = clean_str(adata.obs[args.celltype_key])

    # 1) 排除指定细胞类型（按你给的四类）
    drop_norm = set(norm_ct_label(x) for x in args.drop_celltypes)
    ct_norm = adata.obs[args.celltype_key].map(norm_ct_label)

    keep_mask = ~ct_norm.isin(drop_norm)
    adata2 = adata[keep_mask].copy()

    print(f"[INFO] drop_celltypes={args.drop_celltypes}")
    print(f"[INFO] cells before={adata.n_obs} after_drop={adata2.n_obs} (dropped={adata.n_obs-adata2.n_obs})")

    # 2) 统计每个 (patient, celltype) 的 n_cells（用于 meta）
    obs = adata2.obs.copy()
    obs["_pid"] = clean_str(obs[args.patient_key])
    obs["_ct"] = clean_str(obs[args.celltype_key])

    # meta 需要的列（没有就填 NA）
    for k in [args.disease_key, args.study_key, args.sex_key, args.binned_age_key, args.age_key]:
        if k not in obs.columns:
            obs[k] = np.nan

    grp_cols = ["_pid", "_ct"]
    g = obs.groupby(grp_cols, dropna=False)

    meta_agg = g.agg(
        n_cells=(args.patient_key, "size"),
        disease=(args.disease_key, mode1),
        dataset_id=(args.study_key, mode1),
        sex=(args.sex_key, mode1),
        binned_age=(args.binned_age_key, mode1),
        age_years=(args.age_key, "median"),
    ).reset_index()

    # 3) decoupler pseudobulk：patient x celltype
    #    关键：min_cells=args.min_cells, mode='sum', layer=counts
    print("[STEP] decoupler.get_pseudobulk ...")
    pb = decoupler.get_pseudobulk(
        adata=adata2,
        min_cells=args.min_cells,
        sample_col=args.patient_key,
        groups_col=args.celltype_key,
        layer=args.counts_layer,
        mode="sum",
    )

    # pb.obs 通常含 sample_col / groups_col（列名可能是 sample / group 或你传入的 key）
    # 为稳妥：统一造一个 pb_id = patient|celltype
    pb_obs = pb.obs.copy()

    # 找到 patient / celltype 列（优先用你传入的 key）
    if args.patient_key in pb_obs.columns:
        pb_pid = clean_str(pb_obs[args.patient_key])
    elif "sample" in pb_obs.columns:
        pb_pid = clean_str(pb_obs["sample"])
    else:
        # 最差：用 index 解析
        pb_pid = pd.Series(pb_obs.index.astype(str)).str.split("|").str[0]

    if args.celltype_key in pb_obs.columns:
        pb_ct = clean_str(pb_obs[args.celltype_key])
    elif "group" in pb_obs.columns:
        pb_ct = clean_str(pb_obs["group"])
    else:
        pb_ct = pd.Series(pb_obs.index.astype(str)).str.split("|").str[-1]

    pb_id = (pb_pid.astype(str) + "|" + pb_ct.astype(str)).tolist()
    pb.obs_names = pb_id  # 让列名直接可用

    # 4) 合并 meta（n_cells / disease / sex / binned_age / dataset_id 等）
    meta_agg["_pb_id"] = meta_agg["_pid"].astype(str) + "|" + meta_agg["_ct"].astype(str)
    meta_agg = meta_agg.set_index("_pb_id")

    # 对齐 pb.obs_names
    meta_final = meta_agg.reindex(pb.obs_names).copy()

    # lib_size（pseudobulk library size）加上，edgeR 常用
    Xpb = pb.X
    if sp.issparse(Xpb):
        lib_size = np.asarray(Xpb.sum(axis=1)).ravel()
    else:
        lib_size = np.asarray(Xpb).sum(axis=1)
    meta_final["lib_size"] = lib_size.astype(np.float64)

    # 再把 patient/celltype 列写出来（便于 R 过滤）
    meta_final["patient_id"] = [x.split("|")[0] for x in pb.obs_names]
    meta_final["Level1"] = [x.split("|", 1)[1] for x in pb.obs_names]

    # 输出 meta
    out_meta = os.path.join(args.out_dir, args.out_meta)
    meta_final.reset_index(names="pb_id").to_csv(out_meta, sep="\t", index=False)
    print("[WRITE]", out_meta, "n=", meta_final.shape[0])

    # 5) 输出 counts：gene x sample（gzip）
    # pb.X: (n_pb x n_genes) -> 转成 (n_genes x n_pb)
    if sp.issparse(pb.X):
        X_gene_by_sample = pb.X.T.tocsr()
    else:
        X_gene_by_sample = np.asarray(pb.X).T

    out_counts = os.path.join(args.out_dir, args.out_counts)
    print("[STEP] writing counts (gene x sample) ...")
    write_counts_tsv_gz_gene_by_sample(
        X_gene_by_sample=X_gene_by_sample,
        genes=pb.var_names.tolist(),
        sample_ids=pb.obs_names.tolist(),
        out_path=out_counts
    )
    print("[WRITE]", out_counts)

    # 汇总
    print("\n[SUMMARY]")
    print("pseudobulk n_samples =", pb.n_obs)
    print("pseudobulk n_genes   =", pb.n_vars)
    print("min_cells per (patient, Level1) =", args.min_cells)
    print("counts file:", out_counts)
    print("meta   file:", out_meta)
    print("[DONE]")


if __name__ == "__main__":
    main()
