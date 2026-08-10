#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import argparse
import pickle

import numpy as np
import pandas as pd
import anndata as ad

import Spectra


def ensure_dir(p):
    os.makedirs(p, exist_ok=True)


def clean_str(x):
    if x is None:
        return None
    x = str(x)
    x = re.sub(r"\s+", " ", x).strip()
    if x in ("", "nan", "None", "NA"):
        return None
    return x


def normalize_group_name(x):
    """
    统一原始 gmt 第2列里的 group 名字
    例如：
      group=GLOBAL
      group=Lymphocyte T-CD8 Non Naive
    """
    x = clean_str(x)
    if x is None:
        return None

    # 去掉前缀 group=
    x = re.sub(r"^group\s*=\s*", "", x, flags=re.IGNORECASE)

    # 统一空格
    x = re.sub(r"\s+", " ", x).strip()

    # 统一一些常见写法
    x = x.replace("Non Naive", "NonNaive")
    x = x.replace("T CD4", "T-CD4")
    x = x.replace("T CD8", "T-CD8")

    return x


def source_group_to_target_celltype(group_name):
    """
    原始 gmt 里的 group -> 你当前 metacell_meta 的 celltype
    """
    g = normalize_group_name(group_name)
    if g is None:
        return None

    mapping = {
        "GLOBAL": "global",
        "Lymphocyte B": "B",
        "Monocytes": "Mono",
        "Lymphocyte T-CD4 Naive": "T CD4 naive",
        "Lymphocyte T-CD4 NonNaive": "T CD4 non-naive",
        "Lymphocyte T-CD8 Naive": "T CD8 naive",
        "Lymphocyte T-CD8 NonNaive": "T CD8 non-naive",
        "Plasma": "Plasma",
        "UTC": "UTC",
        "pDC": "pDC",
        "ILC": "ILC",
        "DC": "DC",
    }

    if g in mapping:
        return mapping[g]

    gl = g.lower()

    # 容错
    if gl == "global":
        return "global"
    if "lymphocyte b" in gl:
        return "B"
    if "monocyte" in gl:
        return "Mono"
    if "t-cd4 naive" in gl:
        return "T CD4 naive"
    if "t-cd4 nonnaive" in gl:
        return "T CD4 non-naive"
    if "t-cd8 naive" in gl:
        return "T CD8 naive"
    if "t-cd8 nonnaive" in gl:
        return "T CD8 non-naive"
    if gl == "plasma":
        return "Plasma"
    if gl == "utc":
        return "UTC"
    if gl == "pdc":
        return "pDC"
    if gl == "ilc":
        return "ILC"
    if gl == "dc":
        return "DC"

    # 未识别：并入 global，避免丢失
    return "global"


def read_original_gmt_as_nested_dict(gmt_path, expr_genes):
    """
    读取原始 gmt：
      geneset_name \t group=XXX \t gene1 \t gene2 ...

    返回：
      gene_set_dict: Spectra 用 nested dict
      net_df: ULM 用长表
    """
    expr_genes = set(expr_genes)

    gene_set_dict = {"global": {}}
    net_rows = []

    total_sets = 0
    kept_sets = 0

    with open(gmt_path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue

            parts = line.split("\t")
            if len(parts) < 3:
                continue

            total_sets += 1

            gs_name = clean_str(parts[0])
            raw_group = clean_str(parts[1])
            genes = [clean_str(x) for x in parts[2:]]
            genes = [g for g in genes if g is not None]

            if gs_name is None or raw_group is None or len(genes) == 0:
                continue

            target_group = source_group_to_target_celltype(raw_group)

            # 交到表达矩阵基因空间
            genes = [g for g in genes if g in expr_genes]
            genes = sorted(set(genes))

            # 每个 geneset 至少 2 个基因
            if len(genes) < 2:
                continue

            gene_set_dict.setdefault(target_group, {})
            gene_set_dict[target_group][gs_name] = genes

            for g in genes:
                net_rows.append((gs_name, g, target_group, normalize_group_name(raw_group)))

            kept_sets += 1

    gene_set_dict = {k: v for k, v in gene_set_dict.items() if len(v) > 0 or k == "global"}
    if "global" not in gene_set_dict:
        gene_set_dict["global"] = {}

    net_df = pd.DataFrame(
        net_rows,
        columns=["gene_set", "gene", "Level1", "source_group"]
    )

    return gene_set_dict, net_df, total_sets, kept_sets


def clean_celltype_series(s):
    s = s.astype(str).str.strip()
    s = s.replace({
        "nan": np.nan,
        "None": np.nan,
        "NA": np.nan,
        "": np.nan
    })
    return s


def main():
    ap = argparse.ArgumentParser(
        description="Run Spectra from metacells using ORIGINAL GMT (2nd column = group=...)."
    )
    ap.add_argument("--metacell_expr", required=True)
    ap.add_argument("--metacell_meta", required=True)
    ap.add_argument("--gmt", required=True)
    ap.add_argument("--out_dir", required=True)

    ap.add_argument("--celltype_key", default="celltype")

    ap.add_argument("--lam", type=float, default=0.001)
    ap.add_argument("--overlap_threshold", type=float, default=0.3)
    ap.add_argument("--delta", type=float, default=0.001)
    ap.add_argument("--rho", type=float, default=0.001)
    ap.add_argument("--n_top_vals", type=int, default=200)
    ap.add_argument("--min_gs_num", type=int, default=2)
    ap.add_argument("--num_epochs", type=int, default=15000)

    args = ap.parse_args()
    ensure_dir(args.out_dir)

    print("[READ]", args.metacell_expr)
    expr = pd.read_csv(args.metacell_expr, sep="\t", index_col=0, compression="gzip")

    print("[READ]", args.metacell_meta)
    meta = pd.read_csv(args.metacell_meta, sep="\t", index_col=0)

    common = expr.index.intersection(meta.index)
    if len(common) == 0:
        raise RuntimeError("No overlapping rows between metacell_expr and metacell_meta.")

    expr = expr.loc[common].copy()
    meta = meta.loc[common].copy()

    if args.celltype_key not in meta.columns:
        raise ValueError(f"metacell_meta missing column: {args.celltype_key}")

    # 清洗 celltype
    meta[args.celltype_key] = clean_celltype_series(meta[args.celltype_key])
    bad_mask = meta[args.celltype_key].isna()
    if bad_mask.any():
        print(f"[WARN] dropping {bad_mask.sum()} metacells with invalid '{args.celltype_key}'")
        expr = expr.loc[~bad_mask].copy()
        meta = meta.loc[~bad_mask].copy()

    meta["spectra_categories"] = meta[args.celltype_key].astype(str)
    meta_celltypes = sorted(meta["spectra_categories"].unique().tolist())
    print("[INFO] celltypes in metacell_meta =", meta_celltypes)

    # 读取原始 gmt
    print("[READ]", args.gmt)
    gene_set_dict_raw, net_df, total_sets, kept_sets = read_original_gmt_as_nested_dict(
        args.gmt,
        expr.columns
    )

    print("[INFO] genesets in original GMT =", total_sets)
    print("[INFO] genesets kept after expr intersection =", kept_sets)

    # 只保留在 meta 中真实存在的 celltype + global
    valid_ct = set(meta["spectra_categories"].unique().tolist())

    gene_set_dict = {
        k: v for k, v in gene_set_dict_raw.items()
        if (k == "global") or (k in valid_ct)
    }
    if "global" not in gene_set_dict:
        gene_set_dict["global"] = {}

    if net_df.shape[0] > 0:
        net_df = net_df[net_df["Level1"].isin(valid_ct | {"global"})].copy()

    # 只保留“有非空先验”的细胞类型
    nonempty_ct = sorted([k for k, v in gene_set_dict.items() if k != "global" and len(v) > 0])

    if len(nonempty_ct) == 0:
        raise RuntimeError("No non-empty cell-type-specific genesets after parsing original GMT.")

    keep_cells = meta["spectra_categories"].isin(nonempty_ct)
    n_drop = (~keep_cells).sum()
    if n_drop > 0:
        print(f"[WARN] dropping {n_drop} metacells whose celltypes have no non-empty genesets")

    expr = expr.loc[keep_cells].copy()
    meta = meta.loc[keep_cells].copy()

    valid_ct_after = sorted(meta["spectra_categories"].unique().tolist())

    gene_set_dict = {
        k: v for k, v in gene_set_dict.items()
        if (k == "global") or (k in valid_ct_after)
    }

    if net_df.shape[0] > 0:
        net_df = net_df[net_df["Level1"].isin(set(valid_ct_after) | {"global"})].copy()

    dict_keys_no_global = sorted([k for k in gene_set_dict.keys() if k != "global"])

    print("[INFO] celltypes kept for Spectra =", valid_ct_after)
    print("[CHECK] gene_set_dict keys(no global) =", dict_keys_no_global)

    if set(valid_ct_after) != set(dict_keys_no_global):
        raise ValueError(
            f"Mismatch before AnnData build: meta celltypes={set(valid_ct_after)}, "
            f"gene_set_dict keys(no global)={set(dict_keys_no_global)}"
        )

    n_gs = sum(len(v) for v in gene_set_dict.values())
    print("[INFO] metacells kept =", expr.shape[0])
    print("[INFO] genes =", expr.shape[1])
    print("[INFO] gene sets kept =", n_gs)

    # 构建 AnnData
    obs = meta.copy()
    obs.index = expr.index.astype(str)

    adata = ad.AnnData(
        X=expr.values.astype(np.float32),
        obs=obs,
        var=pd.DataFrame(index=expr.columns.astype(str))
    )
    adata.obs_names = expr.index.astype(str)
    adata.var_names = expr.columns.astype(str)
    adata.obs["spectra_categories"] = adata.obs["spectra_categories"].astype("category")

    obs_key = "spectra_categories"
    print("[CHECK] adata labels =", list(adata.obs[obs_key].cat.categories))

    adata_labels = set(adata.obs[obs_key].cat.categories.tolist())
    dict_labels = set([k for k in gene_set_dict.keys() if k != "global"])
    if adata_labels != dict_labels:
        raise ValueError(
            f"Final mismatch: adata labels={adata_labels}, gene_set_dict keys(no global)={dict_labels}"
        )

    # 保存解析后的 geneset 清单，便于核查
    parsed_out = os.path.join(args.out_dir, "parsed_genesets_used.tsv.gz")
    net_df.to_csv(parsed_out, sep="\t", index=False, compression="gzip")

    # 跑 Spectra
    print("[STEP] Running Spectra.est_spectra(...)")
    model = Spectra.est_spectra(
        adata=adata,
        gene_set_dictionary=gene_set_dict,
        use_highly_variable=False,
        cell_type_key=obs_key,
        use_weights=True,
        lam=args.lam,
        delta=args.delta,
        kappa=None,
        rho=args.rho,
        use_cell_types=True,
        n_top_vals=args.n_top_vals,
        label_factors=True,
        overlap_threshold=args.overlap_threshold,
        clean_gs=True,
        verbose=True,
        min_gs_num=args.min_gs_num,
        num_epochs=args.num_epochs,
    )

    # 保存模型
    model_pkl = os.path.join(args.out_dir, "spectra_model.pkl")
    with open(model_pkl, "wb") as f:
        pickle.dump(model, f)

    adata_h5ad = os.path.join(args.out_dir, "adata_spectra.h5ad")
    adata.write_h5ad(adata_h5ad, compression="gzip")

    net_out = os.path.join(args.out_dir, "gene_set_network_long.tsv")
    net_df.to_csv(net_out, sep="\t", index=False)

    # 提取结果
    if "SPECTRA_overlap" not in adata.uns or "SPECTRA_factors" not in adata.uns:
        raise RuntimeError("SPECTRA outputs missing in adata.uns. Check Spectra version/result.")

    factor_index = adata.uns["SPECTRA_overlap"].index

    if "spectra_vocab" in adata.var.columns:
        spectra_genes = adata.var_names[adata.var["spectra_vocab"].values]
    else:
        spectra_genes = adata.var_names

    gene_weights = pd.DataFrame(
        adata.uns["SPECTRA_factors"],
        index=factor_index,
        columns=spectra_genes
    )

    if "SPECTRA_cell_scores" in adata.obsm:
        factor_scores = pd.DataFrame(
            adata.obsm["SPECTRA_cell_scores"],
            index=adata.obs_names,
            columns=factor_index
        )
    else:
        factor_scores = pd.DataFrame(index=adata.obs_names)

    gene_weights_out = os.path.join(args.out_dir, "factor_gene_weights.tsv.gz")
    gene_weights.to_csv(gene_weights_out, sep="\t", compression="gzip")

    factor_scores_out = os.path.join(args.out_dir, "factor_scores.tsv.gz")
    factor_scores.to_csv(factor_scores_out, sep="\t", compression="gzip")

    overlap_out = os.path.join(args.out_dir, "SPECTRA_overlap.tsv")
    adata.uns["SPECTRA_overlap"].to_csv(overlap_out, sep="\t")

    print("\n[DONE]")
    print("adata_spectra:", adata_h5ad)
    print("spectra_model:", model_pkl)
    print("parsed_genesets_used:", parsed_out)
    print("gene_set_network_long:", net_out)
    print("factor_gene_weights:", gene_weights_out, gene_weights.shape)
    print("factor_scores:", factor_scores_out, factor_scores.shape)
    print("SPECTRA_overlap:", overlap_out)


if __name__ == "__main__":
    main()