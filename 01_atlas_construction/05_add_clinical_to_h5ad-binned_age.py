#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os, argparse
import numpy as np
import pandas as pd


def ensure_dir(p):
    if p:
        os.makedirs(p, exist_ok=True)


def normalize_patient_id(s: pd.Series) -> pd.Series:
    return s.astype(str).str.strip()


def normalize_sex_value(v):
    """统一成 'M'/'F' 或 NaN；兼容 0/1、男/女、male/female 等。"""
    if v is None:
        return np.nan
    try:
        if pd.isna(v):
            return np.nan
    except Exception:
        pass

    s = str(v).strip()
    if s == "" or s.lower() in {"nan", "none", "na", "null"}:
        return np.nan

    s_low = s.lower()
    if s_low in {"1", "m", "male", "man", "boy", "男"}:
        return "M"
    if s_low in {"0", "2", "f", "female", "woman", "girl", "女"}:
        return "F"
    if s.upper() in {"M", "F"}:
        return s.upper()
    return s  # 兜底保留字符串


def bin_age_years(age: float) -> str:
    """
    你说“没具体分箱规则，你看着来”：
    给一个常用儿童分箱（尽量平衡、可解释）：
      <1, 1-2, 3-5, 6-11, 12-17, >=18
    """
    if age is None or pd.isna(age):
        return np.nan
    try:
        a = float(age)
    except Exception:
        return np.nan

    if a < 1:
        return "<1"
    elif a < 3:
        return "1-2"
    elif a < 6:
        return "3-5"
    elif a < 12:
        return "6-11"
    elif a < 18:
        return "12-17"
    else:
        return ">=18"


def make_obs_h5ad_safe(df: pd.DataFrame) -> pd.DataFrame:
    """
    写 h5ad 前，把 obs 清洗成稳定可写：
    - numeric/bool/category 保留
    - 其余列：逐元素转为 python str 或 np.nan
    """
    out = df.copy()
    for c in out.columns:
        s = out[c]

        if isinstance(s.dtype, pd.CategoricalDtype) or pd.api.types.is_bool_dtype(s) or pd.api.types.is_numeric_dtype(s):
            continue

        def _to_str_or_nan(v):
            if v is None:
                return np.nan
            try:
                if pd.isna(v):
                    return np.nan
            except Exception:
                pass
            # 防止 list/dict/tuple
            if isinstance(v, (list, dict, tuple, set)):
                return str(v)
            return str(v)

        out[c] = s.map(_to_str_or_nan).astype(object)

    return out


def main():
    ap = argparse.ArgumentParser(
        description="Add clinical columns (age_years/sex/binned_age) into h5ad.obs by patient_id, h5ad-safe."
    )
    ap.add_argument("--in_h5ad", required=True)
    ap.add_argument("--clinical_xlsx", required=True)
    ap.add_argument("--sheet", default=None, help="Excel sheet name (default: first sheet)")

    ap.add_argument("--patient_key_h5ad", default="patient_id")
    ap.add_argument("--patient_key_table", default="patient_id")

    ap.add_argument("--age_col", default="age_years")
    ap.add_argument("--sex_col", default="sex")

    ap.add_argument("--write_binned_age", action="store_true",
                    help="If set, create obs['binned_age'] from age_years with built-in bins.")

    ap.add_argument("--out_h5ad", required=True)

    args = ap.parse_args()
    ensure_dir(os.path.dirname(args.out_h5ad))

    import scanpy as sc

    print("[READ]", args.in_h5ad)
    adata = sc.read_h5ad(args.in_h5ad)

    if args.patient_key_h5ad not in adata.obs.columns:
        raise SystemExit(f"[FATAL] h5ad missing obs['{args.patient_key_h5ad}']")

    print("[READ]", args.clinical_xlsx)
    clin = pd.read_excel(args.clinical_xlsx, sheet_name=args.sheet)

    # 必需列检查
    need_cols = [args.patient_key_table, args.age_col, args.sex_col]
    missing = [c for c in need_cols if c not in clin.columns]
    if missing:
        raise SystemExit(f"[FATAL] clinical table missing columns: {missing}")

    # 规范化 patient_id
    clin = clin.copy()
    clin["_pid"] = normalize_patient_id(clin[args.patient_key_table])
    clin = clin.drop_duplicates("_pid", keep="first")

    # age -> numeric
    clin["_age"] = pd.to_numeric(clin[args.age_col], errors="coerce")

    # sex -> normalized
    clin["_sex"] = clin[args.sex_col].map(normalize_sex_value)

    # h5ad pid
    pid = normalize_patient_id(adata.obs[args.patient_key_h5ad])
    adata.obs["_pid"] = pid.values  # 临时列，后面删

    # 建映射
    age_map = pd.Series(clin["_age"].values, index=clin["_pid"].values).to_dict()
    sex_map = pd.Series(clin["_sex"].values, index=clin["_pid"].values).to_dict()

    # 写入
    adata.obs["age_years"] = adata.obs["_pid"].map(age_map).astype(float)
    adata.obs["sex"] = adata.obs["_pid"].map(sex_map).astype(object)
    adata.obs["sex"] = pd.Categorical(adata.obs["sex"], categories=["M", "F"])

    if args.write_binned_age:
        adata.obs["binned_age"] = adata.obs["age_years"].map(bin_age_years).astype(object)
        adata.obs["binned_age"] = adata.obs["binned_age"].astype("category")

    # 统计
    n_cells = adata.n_obs
    n_age = int(pd.notna(adata.obs["age_years"]).sum())
    n_sex = int(pd.notna(adata.obs["sex"]).sum())
    print(f"[INFO] cells total={n_cells}")
    print(f"[INFO] age_years filled cells={n_age} ({n_age/n_cells*100:.3f}%)")
    print(f"[INFO] sex filled cells={n_sex} ({n_sex/n_cells*100:.3f}%)")

    # 按 patient 粒度统计（用 obs 去重）
    pat_df = adata.obs[["_pid", "age_years", "sex"]].drop_duplicates("_pid")
    print(f"[INFO] patients total={pat_df.shape[0]}")
    print(f"[INFO] patients with age={int(pd.notna(pat_df['age_years']).sum())}")
    print(f"[INFO] patients with sex={int(pd.notna(pat_df['sex']).sum())}")

    # 写 h5ad 前清洗 obs，避免你遇到的 h5py 报错
    # 1) 去掉临时列
    adata.obs.drop(columns=["_pid"], inplace=True, errors="ignore")
    # 2) 整体安全化
    adata.obs = make_obs_h5ad_safe(adata.obs)

    print("[WRITE]", args.out_h5ad)
    adata.write_h5ad(args.out_h5ad, compression="gzip")
    print("[DONE]")


if __name__ == "__main__":
    main()
