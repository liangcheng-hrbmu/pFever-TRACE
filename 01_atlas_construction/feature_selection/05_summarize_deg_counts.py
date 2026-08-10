#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os, argparse
import pandas as pd
import numpy as np

def ensure_dir(p): os.makedirs(p, exist_ok=True)

def pick_col(cols, candidates):
    lower = {c.lower(): c for c in cols}
    for cand in candidates:
        if cand.lower() in lower:
            return lower[cand.lower()]
    # fuzzy contains
    for cand in candidates:
        for c in cols:
            if cand.lower() in c.lower():
                return c
    return None

def norm_str(x):
    if pd.isna(x):
        return np.nan
    return str(x).strip()

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--deg_union_tsv", required=True, help="DEG_union_*.tsv (edgeR union output)")
    ap.add_argument("--out_dir", required=True)
    ap.add_argument("--healthy_label", default="Pediatric healthy",
                    help="disease column 里健康的名字（默认 Pediatric healthy）")
    ap.add_argument("--gene_col", default=None, help="可选：显式指定基因列名")
    ap.add_argument("--disease_col", default=None, help="可选：显式指定疾病列名")
    ap.add_argument("--celltype_col", default=None, help="可选：显式指定细胞类型列名（用于诊断，不强依赖）")
    ap.add_argument("--count_mode", default="unique_genes", choices=["unique_genes", "rows"],
                    help="unique_genes=每个比较统计unique gene数（推荐）；rows=按行数计数")
    ap.add_argument("--min_abslogfc", type=float, default=None,
                    help="可选：再次过滤 |logFC| >= 阈值（如果 union 表没过滤干净）")
    ap.add_argument("--max_fdr", type=float, default=None,
                    help="可选：再次过滤 FDR <= 阈值（如果 union 表没过滤干净）")
    args = ap.parse_args()

    ensure_dir(args.out_dir)

    print("[READ]", args.deg_union_tsv)
    df = pd.read_csv(args.deg_union_tsv, sep="\t")
    df.columns = [str(c) for c in df.columns]

    gene_col = args.gene_col or pick_col(df.columns, ["gene", "genes", "symbol"])
    disease_col = args.disease_col or pick_col(df.columns, ["disease", "disease_level_3", "condition", "group"])
    celltype_col = args.celltype_col or pick_col(df.columns, ["celltype", "level1", "aifi_l1", "ct"])

    if gene_col is None or disease_col is None:
        raise SystemExit(f"[FATAL] cannot infer gene_col or disease_col. "
                         f"gene_col={gene_col}, disease_col={disease_col}. "
                         f"Please pass --gene_col/--disease_col explicitly.")

    # normalize
    df[gene_col] = df[gene_col].map(norm_str)
    df[disease_col] = df[disease_col].map(norm_str)
    if celltype_col and celltype_col in df.columns:
        df[celltype_col] = df[celltype_col].map(norm_str)

    # optional refilter
    if args.min_abslogfc is not None:
        logfc_col = pick_col(df.columns, ["logfc", "log_fc", "log2fc", "log2_fc"])
        if logfc_col is None:
            print("[WARN] --min_abslogfc set but cannot find logFC col; skip.")
        else:
            df = df[pd.to_numeric(df[logfc_col], errors="coerce").abs() >= args.min_abslogfc].copy()

    if args.max_fdr is not None:
        fdr_col = pick_col(df.columns, ["fdr", "padj", "qvalue", "adj_pval"])
        if fdr_col is None:
            print("[WARN] --max_fdr set but cannot find FDR col; skip.")
        else:
            df = df[pd.to_numeric(df[fdr_col], errors="coerce") <= args.max_fdr].copy()

    # build comparison label: healthy_vs_disease
    # union 表通常是 healthy vs disease 的结果汇总；这里按 disease 列做比较名
    # 排除 healthy 自己（如果 union 表里包含）
    df = df[df[disease_col].notna() & df[gene_col].notna()].copy()
    df = df[df[disease_col] != args.healthy_label].copy()

    df["comparison"] = "healthy_vs_" + df[disease_col].astype(str)

    # count
    if args.count_mode == "unique_genes":
        counts = df.groupby("comparison")[gene_col].nunique().sort_values(ascending=False)
    else:
        counts = df.groupby("comparison").size().sort_values(ascending=False)

    out1 = os.path.join(args.out_dir, "deg_counts_by_comparison.tsv")
    counts.reset_index().rename(columns={0: "n_degs", gene_col: "n_degs"}).to_csv(out1, sep="\t", index=False)
    print("[WRITE]", out1)

    # two-column layout like your example
    comp = counts.index.tolist()
    n = counts.values.tolist()
    left = list(zip(comp[0::2], n[0::2]))
    right = list(zip(comp[1::2], n[1::2]))

    rows = max(len(left), len(right))
    left += [("", "")] * (rows - len(left))
    right += [("", "")] * (rows - len(right))

    pretty = pd.DataFrame({
        "comparison_1": [x[0] for x in left],
        "DEG_n_1":      [x[1] for x in left],
        "comparison_2": [x[0] for x in right],
        "DEG_n_2":      [x[1] for x in right],
    })

    out2 = os.path.join(args.out_dir, "deg_counts_two_col_table.tsv")
    pretty.to_csv(out2, sep="\t", index=False)
    print("[WRITE]", out2)

    # small print
    print("\n[TOP 20 comparisons]")
    print(counts.head(20).to_string())

    if celltype_col and celltype_col in df.columns:
        print("\n[INFO] celltypes in union:", df[celltype_col].nunique())

if __name__ == "__main__":
    main()
