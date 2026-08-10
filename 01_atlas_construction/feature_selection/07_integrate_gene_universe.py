#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os, argparse, re
import pandas as pd

def ensure_dir(p): os.makedirs(p, exist_ok=True)

def read_gene_list(path: str):
    """
    支持：
    - .txt: 一行一个 gene（允许空行/注释）
    - .tsv/.csv: 自动找第一列含 gene 的列；否则用第一列
    """
    if not os.path.exists(path):
        raise FileNotFoundError(path)

    ext = os.path.splitext(path)[1].lower()

    if ext in [".txt", ".list", ".genes"]:
        genes = []
        with open(path, "r", encoding="utf-8") as f:
            for line in f:
                s = line.strip()
                if not s or s.startswith("#"):
                    continue
                genes.append(s)
        return genes

    # tsv/csv/anything tabular
    sep = "\t" if ext in [".tsv"] else None
    df = pd.read_csv(path, sep=sep, engine="python")
    if df.shape[0] == 0:
        return []
    # 找 gene 列
    cand = [c for c in df.columns if re.search(r"(gene|symbol)", str(c), flags=re.I)]
    if cand:
        col = cand[0]
    else:
        col = df.columns[0]
    genes = (
        df[col]
        .dropna()
        .astype(str)
        .map(lambda x: re.sub(r"\s+", " ", x).strip())
        .tolist()
    )
    genes = [g for g in genes if g and g.lower() not in ("nan", "none")]
    return genes

def normalize_gene_names(genes):
    # 不强行 upper，避免你数据里有大小写混乱；只做 strip/去空格
    out = []
    for g in genes:
        g = re.sub(r"\s+", " ", str(g)).strip()
        if g and g.lower() not in ("nan", "none"):
            out.append(g)
    return out

def overlap_counts(A, B, C):
    A, B, C = set(A), set(B), set(C)
    AB = A & B
    AC = A & C
    BC = B & C
    ABC = A & B & C
    onlyA = A - (B | C)
    onlyB = B - (A | C)
    onlyC = C - (A | B)
    return {
        "only_HVG": len(onlyA),
        "only_DEG": len(onlyB),
        "only_CURATED": len(onlyC),
        "HVG_and_DEG": len(AB - C),
        "HVG_and_CURATED": len(AC - B),
        "DEG_and_CURATED": len(BC - A),
        "HVG_and_DEG_and_CURATED": len(ABC),
        "union_total": len(A | B | C),
        "HVG_total": len(A),
        "DEG_total": len(B),
        "CURATED_total": len(C),
    }

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hvg", required=True, help="HVG list txt")
    ap.add_argument("--deg", required=True, help="DEG union list txt/tsv")
    ap.add_argument("--curated", required=True, help="Curated union list txt/tsv (from Supplementary Table 4)")
    ap.add_argument("--out_dir", required=True)

    ap.add_argument("--out_union_txt", default="gene_universe_union.txt")
    ap.add_argument("--out_union_tsv", default="gene_universe_union.tsv")
    ap.add_argument("--out_overlap_stats", default="overlap_stats.tsv")
    ap.add_argument("--out_source_sizes", default="source_sizes.tsv")

    args = ap.parse_args()
    ensure_dir(args.out_dir)

    hvg = normalize_gene_names(read_gene_list(args.hvg))
    deg = normalize_gene_names(read_gene_list(args.deg))
    cur = normalize_gene_names(read_gene_list(args.curated))

    # 去重（保留原样）
    hvg_set = set(hvg)
    deg_set = set(deg)
    cur_set = set(cur)
    uni = sorted(hvg_set | deg_set | cur_set)

    # 1) union gene list
    out_union_txt = os.path.join(args.out_dir, args.out_union_txt)
    with open(out_union_txt, "w", encoding="utf-8") as f:
        f.write("\n".join(uni) + "\n")

    # 2) union with provenance
    rows = []
    for g in uni:
        rows.append({
            "gene": g,
            "in_HVG": int(g in hvg_set),
            "in_DEG": int(g in deg_set),
            "in_CURATED": int(g in cur_set),
        })
    df = pd.DataFrame(rows)
    out_union_tsv = os.path.join(args.out_dir, args.out_union_tsv)
    df.to_csv(out_union_tsv, sep="\t", index=False)

    # 3) overlap stats
    stats = overlap_counts(hvg_set, deg_set, cur_set)
    out_overlap = os.path.join(args.out_dir, args.out_overlap_stats)
    pd.DataFrame([stats]).to_csv(out_overlap, sep="\t", index=False)

    # 4) sizes
    out_sizes = os.path.join(args.out_dir, args.out_source_sizes)
    pd.DataFrame([
        {"set": "HVG", "n_genes": len(hvg_set)},
        {"set": "DEG_union", "n_genes": len(deg_set)},
        {"set": "CURATED_union", "n_genes": len(cur_set)},
        {"set": "UNION_all", "n_genes": len(uni)},
    ]).to_csv(out_sizes, sep="\t", index=False)

    print("[DONE]")
    print("HVG n =", len(hvg_set))
    print("DEG n =", len(deg_set))
    print("CURATED n =", len(cur_set))
    print("UNION n =", len(uni))
    print("[WRITE]", out_union_txt)
    print("[WRITE]", out_union_tsv)
    print("[WRITE]", out_overlap)
    print("[WRITE]", out_sizes)

if __name__ == "__main__":
    main()
