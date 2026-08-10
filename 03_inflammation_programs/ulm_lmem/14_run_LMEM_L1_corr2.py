#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import argparse
import warnings
import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings("ignore")


def ensure_dir(p):
    os.makedirs(p, exist_ok=True)


def clean_str(x):
    if pd.isna(x):
        return np.nan
    x = str(x).strip()
    if x in ("", "nan", "None", "NA"):
        return np.nan
    return x


def bh_fdr(series):
    s = pd.to_numeric(series, errors="coerce")
    out = np.full(len(s), np.nan, dtype=float)
    mask = s.notna().values
    if mask.sum() > 0:
        out[mask] = multipletests(s.loc[mask].values, method="fdr_bh")[1]
    return out


def fit_pair_model(df_sub, factor, study_key, g1, g0, disease_key, min_per_group=3):
    """
    两组比较：
      grp = 1 for g1 (disease)
      grp = 0 for g0 (healthy)

    优先 MixedLM:
      y ~ grp, random intercept by study_key
    若失败则退化 OLS:
      y ~ grp
    """
    tmp = df_sub[[factor, disease_key, study_key]].copy().dropna()
    tmp = tmp.loc[tmp[disease_key].isin([g1, g0])].copy()
    if tmp.shape[0] == 0:
        return None

    tmp["grp"] = (tmp[disease_key] == g1).astype(int)

    vc = tmp["grp"].value_counts()
    n1 = int(vc.get(1, 0))
    n0 = int(vc.get(0, 0))
    if min(n1, n0) < min_per_group:
        return None

    y = pd.to_numeric(tmp[factor], errors="coerce")
    tmp = tmp.loc[y.notna()].copy()
    y = pd.to_numeric(tmp[factor], errors="coerce").astype(float)

    if y.shape[0] == 0:
        return None

    vc2 = tmp["grp"].value_counts()
    n1_eff = int(vc2.get(1, 0))
    n0_eff = int(vc2.get(0, 0))
    if min(n1_eff, n0_eff) < min_per_group:
        return None

    X = pd.DataFrame(
        {
            "const": 1.0,
            "grp": tmp["grp"].astype(float).values
        },
        index=tmp.index
    )

    fit = None
    model_type = None

    if tmp[study_key].nunique() >= 2:
        try:
            md = sm.MixedLM(endog=y.values, exog=X.values, groups=tmp[study_key].values)
            fit = md.fit(reml=False, method="lbfgs", maxiter=500, disp=False)
            model_type = "LMEM"
        except Exception:
            try:
                md = sm.MixedLM(endog=y.values, exog=X.values, groups=tmp[study_key].values)
                fit = md.fit(reml=False, method="powell", maxiter=500, disp=False)
                model_type = "LMEM"
            except Exception:
                fit = None

    if fit is None:
        try:
            od = sm.OLS(y.values, X.values)
            fit = od.fit()
            model_type = "OLS"
        except Exception:
            return None

    try:
        beta = float(fit.params[1])
        pval = float(fit.pvalues[1])
    except Exception:
        return None

    return {
        "beta": beta,
        "p_value": pval,
        "n_group1": n1_eff,
        "n_group0": n0_eff,
        "n_total": int(tmp.shape[0]),
        "n_studies": int(tmp[study_key].nunique()),
        "model_type": model_type
    }


def main():
    ap = argparse.ArgumentParser(description="Run L1 disease-vs-healthy comparison on ULM factor activities.")
    ap.add_argument("--ulm_estimate_csv", required=True,
                    help="DecoupleR_ulmestimates_L1_Corr.csv")
    ap.add_argument("--meta_tsv", required=True,
                    help="Pseudobulk_metadata_L1_Corr.tsv")
    ap.add_argument("--factor_best_match", required=True,
                    help="factor_best_match.tsv")
    ap.add_argument("--out_dir", required=True)

    ap.add_argument("--celltype_key", default="celltype")
    ap.add_argument("--study_key", default="dataset_id")
    ap.add_argument("--disease_key", default="disease_Level_3")
    ap.add_argument("--healthy_label", default="Pediatric healthy")
    ap.add_argument("--min_per_group", type=int, default=3,
                    help="Minimum sample size required in each group for formal modeling.")
    args = ap.parse_args()

    ensure_dir(args.out_dir)

    print("[READ]", args.ulm_estimate_csv)
    est = pd.read_csv(args.ulm_estimate_csv, index_col=0)

    print("[READ]", args.meta_tsv)
    meta = pd.read_csv(args.meta_tsv, sep="\t", index_col=0)

    print("[READ]", args.factor_best_match)
    factor_best = pd.read_csv(args.factor_best_match, sep="\t")

    common_rows = est.index.intersection(meta.index)
    if len(common_rows) == 0:
        raise RuntimeError("No overlapping rows between estimate and metadata.")
    est = est.loc[common_rows].copy()
    meta = meta.loc[common_rows].copy()

    needed = [args.celltype_key, args.study_key, args.disease_key]
    miss = [x for x in needed if x not in meta.columns]
    if miss:
        raise ValueError(f"metadata missing columns: {miss}")

    meta[args.celltype_key] = meta[args.celltype_key].map(clean_str)
    meta[args.study_key] = meta[args.study_key].map(clean_str)
    meta[args.disease_key] = meta[args.disease_key].map(clean_str)

    if "factor" not in factor_best.columns or "Level1" not in factor_best.columns:
        raise ValueError("factor_best_match.tsv must contain columns: factor, Level1")

    factor2ct = factor_best.set_index("factor")["Level1"].to_dict()
    dat = meta.join(est, how="inner")

    all_results = []
    diag_rows = []

    celltypes = sorted(dat[args.celltype_key].dropna().unique().tolist())

    for ct in celltypes:
        sub = dat.loc[dat[args.celltype_key] == ct].copy()
        if sub.shape[0] == 0:
            continue

        diseases = sorted(sub[args.disease_key].dropna().unique().tolist())
        if args.healthy_label not in diseases:
            diag_rows.append({
                "celltype": ct,
                "group1": np.nan,
                "group0": args.healthy_label,
                "comparison": np.nan,
                "n_group1": np.nan,
                "n_group0": int((sub[args.disease_key] == args.healthy_label).sum()),
                "n_total": int(sub.shape[0]),
                "n_studies": int(sub[args.study_key].nunique()),
                "pair_pass_min_n": False,
                "pair_skip_reason": "healthy label absent in this celltype"
            })
            continue

        disease_list = [d for d in diseases if d != args.healthy_label]
        if len(disease_list) == 0:
            continue

        factors_ct = [
            f for f in est.columns
            if factor2ct.get(f, None) == ct or factor2ct.get(f, None) == "global"
        ]
        if len(factors_ct) == 0:
            continue

        print(f"[RUN] {ct}: pseudobulks={sub.shape[0]}, factors={len(factors_ct)}, disease_vs_healthy={len(disease_list)}")

        for g1 in disease_list:
            g0 = args.healthy_label
            sub_pair = sub.loc[sub[args.disease_key].isin([g1, g0])].copy()
            vc = sub_pair[args.disease_key].value_counts()

            pair_n1 = int(vc.get(g1, 0))
            pair_n0 = int(vc.get(g0, 0))
            pair_pass = min(pair_n1, pair_n0) >= args.min_per_group

            diag_rows.append({
                "celltype": ct,
                "group1": g1,
                "group0": g0,
                "comparison": f"{g1}_vs_{g0}",
                "n_group1": pair_n1,
                "n_group0": pair_n0,
                "n_total": int(sub_pair.shape[0]),
                "n_studies": int(sub_pair[args.study_key].nunique()),
                "pair_pass_min_n": pair_pass,
                "pair_skip_reason": "" if pair_pass else f"min(n_group1,n_group0) < {args.min_per_group}"
            })

            if not pair_pass:
                continue

            for factor in factors_ct:
                if factor not in sub_pair.columns:
                    continue

                res = fit_pair_model(
                    df_sub=sub_pair,
                    factor=factor,
                    study_key=args.study_key,
                    g1=g1,
                    g0=g0,
                    disease_key=args.disease_key,
                    min_per_group=args.min_per_group
                )
                if res is None:
                    continue

                all_results.append({
                    "celltype": ct,
                    "factor": factor,
                    "comparison": f"{g1}_vs_{g0}",
                    "group1": g1,
                    "group0": g0,
                    **res
                })

    diag_df = pd.DataFrame(diag_rows)
    diag_out = os.path.join(args.out_dir, "LMEM_L1_pair_diagnostics.tsv")
    diag_df.to_csv(diag_out, sep="\t", index=False)

    if len(all_results) == 0:
        raise RuntimeError(f"No model results generated. Check diagnostics: {diag_out}")

    res_df = pd.DataFrame(all_results)
    res_df["adj_p"] = bh_fdr(res_df["p_value"])
    res_df["direction"] = np.where(res_df["beta"] > 0, "up_in_group1", "down_in_group1")

    out_all = os.path.join(args.out_dir, "LMEM_L1_all_results.tsv")
    res_df.to_csv(out_all, sep="\t", index=False)

    sig = res_df.loc[res_df["adj_p"] < 0.05].copy()
    out_sig = os.path.join(args.out_dir, "LMEM_L1_significant.tsv")
    sig.to_csv(out_sig, sep="\t", index=False)

    print("\n[DONE]")
    print("diagnostics:", diag_out, diag_df.shape)
    print("all_results:", out_all, res_df.shape)
    print("significant:", out_sig, sig.shape)
    print("model_type_counts:", res_df["model_type"].value_counts().to_dict())

    fail_df = diag_df.loc[~diag_df["pair_pass_min_n"]].copy()
    if fail_df.shape[0] > 0:
        print("\n[PAIR FAILED min_per_group]")
        print(fail_df.groupby(["group1", "group0"]).size().reset_index(name="n_celltypes").to_string(index=False))


if __name__ == "__main__":
    main()
