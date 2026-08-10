#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os, argparse, re
import pandas as pd

def ensure_dir(p): os.makedirs(p, exist_ok=True)

def _clean(x):
    if pd.isna(x):
        return None
    s = re.sub(r"\s+", " ", str(x)).strip()
    if s == "" or s.lower() in ("nan", "none"):
        return None
    return s

def parse_input_genesets_xlsx(
    xlsx: str,
    sheet: str,
    group_row: int,
    geneset_row: int,
    genes_start_row: int,
):
    """
    解析 Supplementary_Table_4.xlsx 的 Input GeneSets（真实结构）：
      - group_row: 细胞类型大类所在行（如 GLOBAL, Lymphocyte B ...）
      - geneset_row: geneset 名称所在行（如 effector, B_chemokines ...）
      - genes_start_row: 基因起始行（往下直到空）
    返回：
      genesets: dict[geneset_name] = genes(list)
      geneset_to_group: dict[geneset_name] = group_name
    """
    df = pd.read_excel(xlsx, sheet_name=sheet, header=None)
    nrows, ncols = df.shape

    if group_row >= nrows or geneset_row >= nrows or genes_start_row >= nrows:
        raise SystemExit("[FATAL] group_row/geneset_row/genes_start_row out of range.")

    # 1) 每列所属 group（用“最近的非空group”填充，模拟合并单元格）
    current_group = None
    col_group = {}
    for j in range(ncols):
        g = _clean(df.iat[group_row, j])
        if g is not None:
            current_group = g
        col_group[j] = current_group  # 允许为 None（无group列会被跳过）

    # 2) geneset 名称
    genesets = {}
    geneset_to_group = {}

    for j in range(ncols):
        grp = col_group.get(j, None)
        gs = _clean(df.iat[geneset_row, j])
        if grp is None or gs is None:
            continue

        # 3) genes：从 genes_start_row 往下取非空
        genes = []
        for i in range(genes_start_row, nrows):
            v = _clean(df.iat[i, j])
            if v is None:
                continue
            # 一格里可能是 "GENE1, GENE2" 这种，拆一下更稳
            parts = re.split(r"[,\s]+", v)
            parts = [p for p in (re.sub(r"\s+", "", x) for x in parts) if p]
            genes.extend(parts)

        genes = sorted(set(genes))
        if len(genes) == 0:
            continue

        genesets[gs] = genes
        geneset_to_group[gs] = grp

    return genesets, geneset_to_group

def write_outputs(genesets, geneset_to_group, out_dir):
    ensure_dir(out_dir)
    by_group_dir = os.path.join(out_dir, "curated_mapped_by_group")
    ensure_dir(by_group_dir)

    # 1) GMT（description写group）
    out_gmt = os.path.join(out_dir, "curated_mapped.gmt")
    with open(out_gmt, "w", encoding="utf-8") as f:
        for gs_name, genes in genesets.items():
            grp = geneset_to_group.get(gs_name, "NA")
            f.write(gs_name + "\t" + f"group={grp}" + "\t" + "\t".join(genes) + "\n")

    # 2) 按 group 输出：每行 geneset_name \t gene1,gene2,...
    group_map = {}
    for gs_name, genes in genesets.items():
        grp = geneset_to_group.get(gs_name, "NA")
        group_map.setdefault(grp, {})
        group_map[grp][gs_name] = genes

    for grp, d in group_map.items():
        safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", grp)
        fn = os.path.join(by_group_dir, f"{safe}.txt")
        with open(fn, "w", encoding="utf-8") as f:
            for gs_name, genes in d.items():
                f.write(gs_name + "\t" + ",".join(genes) + "\n")

    # 3) union genes
    union = sorted(set(g for genes in genesets.values() for g in genes))
    out_union = os.path.join(out_dir, "curated_mapped_union_genes.txt")
    with open(out_union, "w", encoding="utf-8") as f:
        f.write("\n".join(union) + "\n")

    # 4) 简单统计表
    stat_rows = []
    for gs_name, genes in genesets.items():
        stat_rows.append({
            "group": geneset_to_group.get(gs_name, "NA"),
            "geneset": gs_name,
            "n_genes": len(genes),
        })
    stat = pd.DataFrame(stat_rows).sort_values(["group", "geneset"])
    out_stat = os.path.join(out_dir, "curated_mapped_geneset_stats.tsv")
    stat.to_csv(out_stat, sep="\t", index=False)

    print("[DONE] total genesets:", len(genesets))
    print("[DONE] mapped gmt:", out_gmt)
    print("[DONE] by-group dir:", by_group_dir)
    print("[DONE] union genes:", out_union, "n=", len(union))
    print("[DONE] stats:", out_stat)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--xlsx", required=True)
    ap.add_argument("--sheet", default="Input GeneSets")

    # ✅ 对应你这个文件的真实结构（建议就用默认）
    ap.add_argument("--group_row", type=int, default=0, help="大类标题行（如 GLOBAL / Lymphocyte B）")
    ap.add_argument("--geneset_row", type=int, default=2, help="geneset 名称行（如 effector / B_chemokines）")
    ap.add_argument("--genes_start_row", type=int, default=3, help="基因起始行")

    ap.add_argument("--out_dir", required=True)
    args = ap.parse_args()

    genesets, gs2grp = parse_input_genesets_xlsx(
        args.xlsx, args.sheet,
        group_row=args.group_row,
        geneset_row=args.geneset_row,
        genes_start_row=args.genes_start_row,
    )
    write_outputs(genesets, gs2grp, args.out_dir)

if __name__ == "__main__":
    main()
