#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import argparse
import warnings
import numpy as np
import pandas as pd

from sklearn.feature_selection import mutual_info_classif
from statsmodels.stats.multitest import multipletests
import decoupler as dc

try:
    from kneed import KneeLocator
    HAS_KNEE = True
except Exception:
    HAS_KNEE = False


def ensure_dir(p):
    os.makedirs(p, exist_ok=True)


def clean_str(x):
    if pd.isna(x):
        return np.nan
    x = str(x)
    x = re.sub(r"\s+", " ", x).strip()
    if x in ("", "nan", "None", "NA"):
        return np.nan
    return x


def normalize_ct(x):
    x = clean_str(x)
    if pd.isna(x):
        return np.nan
    mapping = {
        "B cell": "B",
        "B cells": "B",
        "B": "B",
        "Monocyte": "Mono",
        "Monocytes": "Mono",
        "Mono": "Mono",
        "T CD4 naive": "T CD4 naive",
        "T CD4 non-naive": "T CD4 non-naive",
        "T CD8 naive": "T CD8 naive",
        "T CD8 non-naive": "T CD8 non-naive",
        "Plasma": "Plasma",
        "UTC": "UTC",
        "pDC": "pDC",
        "ILC": "ILC",
        "DC": "DC",
        "global": "global",
        "GLOBAL": "global",
    }
    return mapping.get(x, x)


def infer_factor_level1(factor_name, valid_cts):
    """
    兼容：
      0-X-global-X-...
      12-X-B-X-...
      35-X-T CD4 naive-X-...
    """
    f = str(factor_name)

    if "-X-" in f:
        parts = f.split("-X-")
        for cand in parts[1:3]:
            cand2 = normalize_ct(cand)
            if cand2 in valid_cts:
                return cand2

    ordered = sorted(valid_cts, key=lambda x: len(str(x)), reverse=True)
    fl = f.lower()
    for ct in ordered:
        if str(ct).lower() in fl:
            return ct

    if "global" in fl:
        return "global"

    return np.nan


def bh_fdr(pvals):
    pvals = pd.to_numeric(pd.Series(pvals), errors="coerce")
    out = np.full(len(pvals), np.nan, dtype=float)
    mask = pvals.notna().values
    if mask.sum() > 0:
        out[mask] = multipletests(pvals.loc[mask].values, method="fdr_bh")[1]
    return out


def safe_mutual_info(x, y):
    try:
        val = mutual_info_classif(x.reshape(-1, 1), y, discrete_features=False)[0]
        return float(val)
    except Exception:
        return np.nan


def get_knee_n(sorted_abs_weights, min_genes=5):
    """
    对一个 factor 的绝对权重降序向量，用 knee 找拐点。
    若失败则返回 fallback。
    """
    n = len(sorted_abs_weights)
    if n == 0:
        return 0

    if HAS_KNEE and n >= 10:
        try:
            x = np.arange(1, n + 1)
            y = np.asarray(sorted_abs_weights, dtype=float)
            kl = KneeLocator(
                x, y,
                curve="convex",
                direction="decreasing"
            )
            if kl.knee is not None:
                return int(max(min_genes, kl.knee))
        except Exception:
            pass

    mx = float(np.nanmax(sorted_abs_weights)) if n > 0 else 0.0
    if mx > 0:
        keep = np.sum(np.asarray(sorted_abs_weights) >= 0.1 * mx)
        return int(max(min_genes, keep))

    return int(min_genes)


def factor_top_genes_by_knee(gene_weights_df, factor, min_genes=5, max_genes=8):
    """
    gene_weights_df: 行=factor, 列=gene
    返回该 factor 被 knee 保留下来的基因列表，并加上上限 cap
    """
    w = gene_weights_df.loc[factor].astype(float)
    w = w.replace([np.inf, -np.inf], np.nan).dropna()
    if w.shape[0] == 0:
        return []

    absw = w.abs().sort_values(ascending=False)
    k = get_knee_n(absw.values, min_genes=min_genes)

    if max_genes is not None:
        k = min(int(k), int(max_genes))

    genes = absw.index[:k].tolist()
    return genes


def main():
    ap = argparse.ArgumentParser(
        description="Spectra factor annotation with old fallback logic + capped gene selection."
    )

    ap.add_argument("--factor_gene_weights", required=True,
                    help="11_run_spectra 输出的 factor_gene_weights.tsv.gz")
    ap.add_argument("--factor_scores", required=True,
                    help="11_run_spectra 输出的 factor_scores.tsv.gz")
    ap.add_argument("--metacell_meta", required=True,
                    help="metacell_meta.tsv")
    ap.add_argument("--gene_set_network", required=True,
                    help="gene_set_network_long.tsv 或 parsed_genesets_used.tsv.gz")
    ap.add_argument("--out_dir", required=True)

    ap.add_argument("--celltype_key", default="celltype")
    ap.add_argument("--disease_key", default="disease_Level_3")
    ap.add_argument("--healthy_label", default="Pediatric healthy")

    ap.add_argument("--fdr_cutoff", type=float, default=0.05)
    ap.add_argument("--estimate_cutoff", type=float, default=0.0)

    ap.add_argument("--S_sel", type=float, default=0.6,
                    help="最终采用的 factor 重要性阈值比例")
    ap.add_argument("--max_n_genes", type=int, default=1100,
                    help="下游分析允许的最大基因数")
    ap.add_argument("--min_genes_per_factor", type=int, default=5)
    ap.add_argument("--max_genes_per_factor", type=int, default=8,
                    help="每个factor最多保留的代表基因数")
    ap.add_argument("--s_grid", default="0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0",
                    help="扫描 S 的网格")
    args = ap.parse_args()

    ensure_dir(args.out_dir)

    print("[READ]", args.factor_gene_weights)
    gene_weights = pd.read_csv(args.factor_gene_weights, sep="\t", index_col=0, compression="gzip")

    print("[READ]", args.factor_scores)
    factor_scores = pd.read_csv(args.factor_scores, sep="\t", index_col=0, compression="gzip")

    print("[READ]", args.metacell_meta)
    meta = pd.read_csv(args.metacell_meta, sep="\t", index_col=0)

    comp = "gzip" if args.gene_set_network.endswith(".gz") else None
    print("[READ]", args.gene_set_network)
    net_df = pd.read_csv(args.gene_set_network, sep="\t", compression=comp)

    if args.celltype_key not in meta.columns:
        raise ValueError(f"Missing metacell_meta column: {args.celltype_key}")
    if args.disease_key not in meta.columns:
        raise ValueError(f"Missing metacell_meta column: {args.disease_key}")

    meta[args.celltype_key] = meta[args.celltype_key].map(normalize_ct)
    meta[args.disease_key] = meta[args.disease_key].map(clean_str)

    valid_cts = sorted(pd.Series(meta[args.celltype_key]).dropna().unique().tolist())
    if "global" not in valid_cts:
        valid_cts = ["global"] + valid_cts

    common_rows = factor_scores.index.intersection(meta.index)
    if len(common_rows) == 0:
        raise RuntimeError("No overlapping rows between factor_scores and metacell_meta.")
    factor_scores = factor_scores.loc[common_rows].copy()
    meta = meta.loc[common_rows].copy()

    gene_weights["Level1"] = [infer_factor_level1(f, valid_cts) for f in gene_weights.index]
    bad = gene_weights["Level1"].isna()
    if bad.any():
        warnings.warn(f"{bad.sum()} factors could not infer Level1 and will be dropped.")
        gene_weights = gene_weights.loc[~bad].copy()

    factor_scores = factor_scores.loc[:, gene_weights.index].copy()

    need_cols = {"Level1", "gene_set", "gene"}
    miss = need_cols - set(net_df.columns)
    if len(miss) > 0:
        raise ValueError(f"gene_set_network missing columns: {sorted(miss)}")
    net_df["Level1"] = net_df["Level1"].map(normalize_ct)

    # =========================================================
    # 1) ULM annotation
    # =========================================================
    ulm_res_list = []

    level1_list = [x for x in pd.Series(gene_weights["Level1"]).dropna().unique().tolist()]
    level1_list = sorted(level1_list, key=lambda x: (x != "global", str(x)))

    print("[STEP 1] ULM annotation by Level1 ...")
    for l1 in level1_list:
        gw_i = gene_weights.loc[gene_weights["Level1"] == l1].drop(columns=["Level1"])
        if gw_i.shape[0] == 0:
            continue

        net_i = net_df.query("Level1 == @l1 or Level1 == 'global'").copy()
        if net_i.shape[0] == 0:
            print(f"[SKIP] {l1} has no network rows")
            continue

        common_genes = [g for g in gw_i.columns if g in set(net_i["gene"])]
        if len(common_genes) < 2:
            print(f"[SKIP] {l1} common genes < 2")
            continue

        gw_i = gw_i.loc[:, common_genes].copy()
        net_i = net_i[net_i["gene"].isin(common_genes)].copy()

        estimate, pvalues = dc.run_ulm(
            net=net_i,
            mat=gw_i,
            source="gene_set",
            target="gene",
            weight=None,
            use_raw=False
        )

        ulm_i = (
            estimate.melt(ignore_index=False, value_name="estimate")
            .reset_index()
            .rename(columns={"index": "factor", "variable": "gene_set"})
            .merge(
                pvalues.melt(ignore_index=False, value_name="p_value")
                .reset_index()
                .rename(columns={"index": "factor", "variable": "gene_set"}),
                on=["factor", "gene_set"]
            )
        )
        ulm_i["Level1"] = l1
        ulm_res_list.append(ulm_i)
        print(f"[OK] {l1}: factors={gw_i.shape[0]} genes={len(common_genes)} tests={ulm_i.shape[0]}")

    if len(ulm_res_list) == 0:
        raise RuntimeError("No ULM results generated.")

    ulm_res = pd.concat(ulm_res_list, axis=0, ignore_index=True)
    ulm_res["estimate"] = pd.to_numeric(ulm_res["estimate"], errors="coerce")
    ulm_res["p_value"] = pd.to_numeric(ulm_res["p_value"], errors="coerce")
    ulm_res["adj_p"] = bh_fdr(ulm_res["p_value"])
    ulm_res["estimate_is_na"] = ulm_res["estimate"].isna()

    ulm_res["annotated"] = (
        ulm_res["adj_p"].notna() &
        ulm_res["estimate"].notna() &
        (ulm_res["adj_p"] < args.fdr_cutoff) &
        (ulm_res["estimate"] > args.estimate_cutoff)
    )

    ulm_out = os.path.join(args.out_dir, "ULMres.tsv")
    ulm_res.to_csv(ulm_out, sep="\t", index=False)

    # =========================================================
    # 1.1) old logic: global strict first, fallback only if no strict at all
    # =========================================================
    ann_strict = ulm_res.loc[
        ulm_res["adj_p"].notna() &
        ulm_res["estimate"].notna() &
        (ulm_res["adj_p"] < args.fdr_cutoff) &
        (ulm_res["estimate"] > args.estimate_cutoff)
    ].copy()

    ann_fallback = ulm_res.loc[
        ulm_res["estimate"].notna() &
        (ulm_res["estimate"] > args.estimate_cutoff)
    ].copy()

    if ann_fallback.shape[0] == 0:
        ann_fallback = ulm_res.loc[ulm_res["estimate"].notna()].copy()

    if ann_strict.shape[0] > 0:
        use_df = ann_strict.copy()
        default_mode = "strict"
    else:
        use_df = ann_fallback.copy()
        default_mode = "fallback"

    if use_df.shape[0] == 0:
        raise RuntimeError("No usable ULM matches found after strict/fallback filtering.")

    idx = use_df.groupby("factor")["estimate"].idxmax()
    factor_best = use_df.loc[idx].copy()
    factor_best["selection_mode"] = default_mode
    factor_best = factor_best.sort_values(["Level1", "estimate"], ascending=[True, False])

    best_out = os.path.join(args.out_dir, "factor_best_match.tsv")
    factor_best.to_csv(best_out, sep="\t", index=False)

    selection_mode_df = pd.DataFrame({
        "factor": factor_best["factor"].values,
        "selection_mode": factor_best["selection_mode"].values
    })
    selection_mode_out = os.path.join(args.out_dir, "factor_selection_modes.tsv")
    selection_mode_df.to_csv(selection_mode_out, sep="\t", index=False)

    selected_factors = factor_best["factor"].drop_duplicates().tolist()
    strict_df = ann_strict.copy()

    # =========================================================
    # 2) factor-level metrics: MI_HvsD / MI_DvsD
    # =========================================================
    print("[STEP 2] Calculating MI metrics ...")
    metrics = pd.DataFrame(index=factor_scores.columns)

    disease = meta[args.disease_key]
    keep = ~disease.isna()
    disease = disease.loc[keep]
    fs = factor_scores.loc[keep].copy()

    if disease.nunique() >= 2:
        y_hvd = (disease != args.healthy_label).astype(int).values
        mi_hvd = {}
        mi_dvd = {}
        for f in fs.columns:
            x = fs[f].astype(float).values
            mi_hvd[f] = safe_mutual_info(x, y_hvd)
            mi_dvd[f] = safe_mutual_info(x, disease.values)
        metrics["MI_HvsD"] = pd.Series(mi_hvd)
        metrics["MI_DvsD"] = pd.Series(mi_dvd)
    else:
        metrics["MI_HvsD"] = np.nan
        metrics["MI_DvsD"] = np.nan

    metrics["selected"] = metrics.index.isin(set(selected_factors))
    metrics["annotated_strict"] = metrics.index.isin(set(strict_df["factor"].unique()))

    metrics = metrics.join(
        factor_best.set_index("factor")[["gene_set", "Level1", "estimate", "adj_p", "selection_mode"]],
        how="left"
    )

    metrics["MI_sum"] = metrics[["MI_HvsD", "MI_DvsD"]].sum(axis=1, skipna=True)

    max_mi = np.nanmax(metrics["MI_sum"].values) if metrics["MI_sum"].notna().any() else np.nan
    if np.isfinite(max_mi) and max_mi > 0:
        metrics["MI_sum_scaled"] = metrics["MI_sum"] / max_mi
    else:
        metrics["MI_sum_scaled"] = np.nan

    metrics_out = os.path.join(args.out_dir, "factor_metrics.tsv")
    metrics.to_csv(metrics_out, sep="\t")

    # =========================================================
    # 3) knee-based gene selection per factor
    # =========================================================
    print("[STEP 3] Knee-based gene selection per factor ...")
    gw_only = gene_weights.drop(columns=["Level1"]).copy()

    factor_gene_rows = []
    factor_gene_sets = {}
    factor_best_map = factor_best.set_index("factor")

    for f in selected_factors:
        top_genes = factor_top_genes_by_knee(
            gw_only,
            f,
            min_genes=args.min_genes_per_factor,
            max_genes=args.max_genes_per_factor
        )
        factor_gene_sets[f] = set(top_genes)

        if len(top_genes) == 0:
            continue

        tmp = pd.DataFrame({
            "factor": f,
            "gene": top_genes,
            "abs_weight_rank": np.arange(1, len(top_genes) + 1),
            "Level1": factor_best_map.loc[f, "Level1"] if f in factor_best_map.index else np.nan,
            "gene_set": factor_best_map.loc[f, "gene_set"] if f in factor_best_map.index else np.nan,
            "selection_mode": factor_best_map.loc[f, "selection_mode"] if f in factor_best_map.index else np.nan,
            "MI_HvsD": metrics.loc[f, "MI_HvsD"] if f in metrics.index else np.nan,
            "MI_DvsD": metrics.loc[f, "MI_DvsD"] if f in metrics.index else np.nan,
            "MI_sum": metrics.loc[f, "MI_sum"] if f in metrics.index else np.nan,
        })
        factor_gene_rows.append(tmp)

    if len(factor_gene_rows) == 0:
        raise RuntimeError("No genes selected by knee procedure.")

    factor_gene_selection = pd.concat(factor_gene_rows, axis=0, ignore_index=True)
    factor_gene_selection_out = os.path.join(args.out_dir, "factor_gene_selection_by_knee.tsv")
    factor_gene_selection.to_csv(factor_gene_selection_out, sep="\t", index=False)

    # =========================================================
    # 4) Scan S and choose downstream genes
    # =========================================================
    print("[STEP 4] Scanning S values ...")
    try:
        s_grid = [float(x) for x in str(args.s_grid).split(",")]
    except Exception:
        raise ValueError("--s_grid 格式应为逗号分隔数字，如 0.0,0.1,0.2,0.3")

    scan_rows = []
    for s in s_grid:
        keep_f = metrics.index[
            metrics["selected"] &
            metrics["MI_sum_scaled"].notna() &
            (metrics["MI_sum_scaled"] >= s)
        ].tolist()

        genes_union = set()
        for f in keep_f:
            genes_union.update(factor_gene_sets.get(f, set()))

        scan_rows.append({
            "S": s,
            "n_factors": len(keep_f),
            "n_genes": len(genes_union),
            "n_strict_factors": int(metrics.loc[keep_f, "annotated_strict"].sum()) if len(keep_f) > 0 else 0
        })

    s_scan = pd.DataFrame(scan_rows).sort_values("S")
    s_scan_out = os.path.join(args.out_dir, "S_scan.tsv")
    s_scan.to_csv(s_scan_out, sep="\t", index=False)

    keep_f_final = metrics.index[
        metrics["selected"] &
        metrics["MI_sum_scaled"].notna() &
        (metrics["MI_sum_scaled"] >= args.S_sel)
    ].tolist()

    genes_union_final = set()
    for f in keep_f_final:
        genes_union_final.update(factor_gene_sets.get(f, set()))

    if len(genes_union_final) > args.max_n_genes:
        ranked_factors = (
            metrics.loc[keep_f_final]
            .sort_values("MI_sum_scaled", ascending=False)
            .index
            .tolist()
        )
        genes_union_trim = set()
        keep_f_trim = []
        for f in ranked_factors:
            gset = factor_gene_sets.get(f, set())
            new_union = genes_union_trim.union(gset)
            if len(new_union) <= args.max_n_genes or len(keep_f_trim) == 0:
                genes_union_trim = new_union
                keep_f_trim.append(f)
            else:
                break
        keep_f_final = keep_f_trim
        genes_union_final = genes_union_trim

    selected_factor_final = pd.DataFrame({
        "factor": keep_f_final
    }).merge(
        metrics.reset_index().rename(columns={"index": "factor"}),
        on="factor",
        how="left"
    )
    selected_factor_final_out = os.path.join(args.out_dir, "selected_factors_final.tsv")
    selected_factor_final.to_csv(selected_factor_final_out, sep="\t", index=False)

    geneSelection_forDownstream = (
        factor_gene_selection.loc[
            factor_gene_selection["factor"].isin(keep_f_final),
            ["gene", "factor", "Level1", "gene_set", "selection_mode", "MI_sum"]
        ]
        .drop_duplicates()
        .sort_values(["MI_sum", "factor"], ascending=[False, True])
        .reset_index(drop=True)
    )
    geneSelection_forDownstream_out = os.path.join(args.out_dir, "geneSelection_forDownstream.tsv")
    geneSelection_forDownstream.to_csv(geneSelection_forDownstream_out, sep="\t", index=False)

    selected_scores = factor_scores.loc[:, keep_f_final].copy()
    selected_scores_out = os.path.join(args.out_dir, "selected_factor_scores.tsv.gz")
    selected_scores.to_csv(selected_scores_out, sep="\t", compression="gzip")

    selected_txt = os.path.join(args.out_dir, "selected_factors.txt")
    with open(selected_txt, "w", encoding="utf-8") as f:
        f.write("\n".join(keep_f_final) + "\n")

    print("\n[DONE]")
    print("ULMres:", ulm_out, ulm_res.shape)
    print("factor_best_match:", best_out, factor_best.shape)
    print("factor_selection_modes:", selection_mode_out, selection_mode_df.shape)
    print("factor_metrics:", metrics_out, metrics.shape)
    print("factor_gene_selection_by_knee:", factor_gene_selection_out, factor_gene_selection.shape)
    print("S_scan:", s_scan_out, s_scan.shape)
    print("selected_factors_final:", selected_factor_final_out, selected_factor_final.shape)
    print("geneSelection_forDownstream:", geneSelection_forDownstream_out, geneSelection_forDownstream.shape)
    print("selected_factors:", selected_txt, len(keep_f_final))
    print("selected_factor_scores:", selected_scores_out, selected_scores.shape)


if __name__ == "__main__":
    main()