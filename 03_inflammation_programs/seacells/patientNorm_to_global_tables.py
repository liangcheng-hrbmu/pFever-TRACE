#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import glob
import argparse

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

import torch
torch.set_float32_matmul_precision("high")
torch.set_printoptions(precision=3, sci_mode=False, edgeitems=7)
torch.multiprocessing.set_sharing_strategy("file_system")

import scvi


def ensure_dir(p):
    os.makedirs(p, exist_ok=True)


def str2bool(x):
    if isinstance(x, bool):
        return x
    x = str(x).strip().lower()
    if x in {"1", "true", "t", "yes", "y"}:
        return True
    if x in {"0", "false", "f", "no", "n"}:
        return False
    raise argparse.ArgumentTypeError(f"Cannot parse boolean from: {x}")


def load_all_seacell_pkls(seacell_pkl_dir):
    pkls = sorted(glob.glob(os.path.join(seacell_pkl_dir, "*.seacell.pkl")))
    if len(pkls) == 0:
        raise FileNotFoundError(f"No *.seacell.pkl found in: {seacell_pkl_dir}")

    dfs = []
    for p in pkls:
        df = pd.read_pickle(p)
        if "SEACell" not in df.columns:
            raise ValueError(f"{p} missing column 'SEACell'")
        tmp = df.copy()
        tmp.index = tmp.index.astype(str)
        dfs.append(tmp)

    out = pd.concat(dfs, axis=0)
    out = out[~out.index.duplicated(keep="first")]
    return out


def mode_or_first(x):
    vc = x.value_counts(dropna=True)
    if vc.shape[0] == 0:
        return np.nan
    return vc.index[0]


def aggregate_dense_matrix(expr_df, group_series, meta_df=None, how="mean"):
    """
    expr_df: cells x genes
    group_series: aligned with expr_df.index
    meta_df: aligned with expr_df.index
    returns:
      agg_expr_df: groups x genes
      agg_meta_df: groups x meta
    """
    group_series = group_series.astype(str)
    uniq = pd.Index(pd.unique(group_series))
    idx_map = {k: i for i, k in enumerate(uniq)}

    X = expr_df.values
    out = np.zeros((len(uniq), X.shape[1]), dtype=np.float64)
    n = np.zeros(len(uniq), dtype=np.int64)

    gids = np.array([idx_map[g] for g in group_series.values], dtype=np.int64)

    for i in range(X.shape[0]):
        g = gids[i]
        out[g] += X[i]
        n[g] += 1

    if how == "mean":
        out = out / n[:, None]
    elif how == "sum":
        pass
    else:
        raise ValueError("how must be mean or sum")

    agg_expr_df = pd.DataFrame(out, index=uniq, columns=expr_df.columns)

    if meta_df is not None:
        tmp = meta_df.copy()
        tmp["_group"] = group_series.values
        agg_meta_df = tmp.groupby("_group").agg(mode_or_first).reindex(uniq)
        agg_meta_df["n_cells"] = n
    else:
        agg_meta_df = pd.DataFrame(index=uniq)
        agg_meta_df["n_cells"] = n

    return agg_expr_df, agg_meta_df


def compute_group_purity(obs_df, group_col, label_col):
    purity = obs_df.groupby(group_col)[label_col].apply(
        lambda x: x.astype(str).value_counts(normalize=True).max()
    )
    return purity


def subset_layer_to_dense(adata, indices, layer_name):
    X = adata.layers[layer_name][indices]
    if sp.issparse(X):
        X = X.toarray()
    else:
        X = np.asarray(X)
    return X.astype(np.float64, copy=False)


def normalize_log1p_per_cell(X, target_sum=1e4):
    lib = X.sum(axis=1, keepdims=True)
    lib[lib == 0] = 1.0
    Xn = X / lib * float(target_sum)
    return np.log1p(Xn)


def get_corrected_expr(scanvi_model, obs_sub, transform_batches, genes, n_samples_norm, library_size):
    indices = obs_sub.index.to_list()
    expr = scanvi_model.get_normalized_expression(
        n_samples=n_samples_norm,
        transform_batch=transform_batches,
        indices=indices,
        library_size=library_size,
        return_numpy=True,
    )
    return pd.DataFrame(
        expr,
        index=obs_sub["cell_id"].astype(str).values,
        columns=genes
    )


def get_uncorrected_expr(adata, obs_sub, genes, counts_layer, library_size, mode="lognorm"):
    indices = obs_sub.index.to_list()
    X = subset_layer_to_dense(adata, indices, counts_layer)

    if mode == "lognorm":
        X = normalize_log1p_per_cell(X, target_sum=library_size)
    elif mode == "counts":
        pass
    else:
        raise ValueError("--uncorrected_mode must be one of: lognorm, counts")

    return pd.DataFrame(
        X,
        index=obs_sub["cell_id"].astype(str).values,
        columns=genes
    )


def variant_output_dirs(root_out, variant_name):
    base = os.path.join(root_out, variant_name)
    dirs = {
        "base": base,
        "pb_l1": os.path.join(base, "pseudobulk", "L1"),
        "pb_l2": os.path.join(base, "pseudobulk", "L2"),
        "mc": os.path.join(base, "metacell"),
    }
    for p in dirs.values():
        ensure_dir(p)
    return dirs


def export_variant_outputs(variant_name, out_root,
                           pb_expr_l1_list, pb_meta_l1_list,
                           pb_expr_l2_list, pb_meta_l2_list,
                           mc_expr_list, mc_meta_list):
    dirs = variant_output_dirs(out_root, variant_name)

    # L1 pseudobulk
    if len(pb_expr_l1_list) == 0:
        raise RuntimeError(f"[{variant_name}] No L1 pseudobulk rows generated.")

    pseudobulk_expr_l1 = pd.concat(pb_expr_l1_list, axis=0)
    pseudobulk_meta_l1 = pd.concat(pb_meta_l1_list, axis=0)

    pseudobulk_expr_l1 = pseudobulk_expr_l1[~pseudobulk_expr_l1.index.duplicated(keep="first")]
    pseudobulk_meta_l1 = pseudobulk_meta_l1.loc[pseudobulk_expr_l1.index]

    pseudobulk_expr_l1.to_csv(os.path.join(dirs["pb_l1"], "pseudobulk_expr.tsv.gz"), sep="\t", compression="gzip")
    pseudobulk_meta_l1.to_csv(os.path.join(dirs["pb_l1"], "pseudobulk_meta.tsv"), sep="\t")

    print("[WRITE]", os.path.join(dirs["pb_l1"], "pseudobulk_expr.tsv.gz"), pseudobulk_expr_l1.shape)
    print("[WRITE]", os.path.join(dirs["pb_l1"], "pseudobulk_meta.tsv"), pseudobulk_meta_l1.shape)

    # L2 pseudobulk
    if len(pb_expr_l2_list) > 0:
        pseudobulk_expr_l2 = pd.concat(pb_expr_l2_list, axis=0)
        pseudobulk_meta_l2 = pd.concat(pb_meta_l2_list, axis=0)

        pseudobulk_expr_l2 = pseudobulk_expr_l2[~pseudobulk_expr_l2.index.duplicated(keep="first")]
        pseudobulk_meta_l2 = pseudobulk_meta_l2.loc[pseudobulk_expr_l2.index]

        pseudobulk_expr_l2.to_csv(os.path.join(dirs["pb_l2"], "pseudobulk_expr.tsv.gz"), sep="\t", compression="gzip")
        pseudobulk_meta_l2.to_csv(os.path.join(dirs["pb_l2"], "pseudobulk_meta.tsv"), sep="\t")

        print("[WRITE]", os.path.join(dirs["pb_l2"], "pseudobulk_expr.tsv.gz"), pseudobulk_expr_l2.shape)
        print("[WRITE]", os.path.join(dirs["pb_l2"], "pseudobulk_meta.tsv"), pseudobulk_meta_l2.shape)
    else:
        print(f"[INFO] {variant_name}: No L2 pseudobulk generated.")

    # metacell
    if len(mc_expr_list) == 0:
        raise RuntimeError(f"[{variant_name}] No metacell rows generated.")

    metacell_expr = pd.concat(mc_expr_list, axis=0)
    metacell_meta = pd.concat(mc_meta_list, axis=0)

    metacell_expr = metacell_expr[~metacell_expr.index.duplicated(keep="first")]
    metacell_meta = metacell_meta.loc[metacell_expr.index]

    metacell_expr.to_csv(os.path.join(dirs["mc"], "metacell_expr.tsv.gz"), sep="\t", compression="gzip")
    metacell_meta.to_csv(os.path.join(dirs["mc"], "metacell_meta.tsv"), sep="\t")

    print("[WRITE]", os.path.join(dirs["mc"], "metacell_expr.tsv.gz"), metacell_expr.shape)
    print("[WRITE]", os.path.join(dirs["mc"], "metacell_meta.tsv"), metacell_meta.shape)

    return {
        "pb_expr_l1": pseudobulk_expr_l1,
        "pb_meta_l1": pseudobulk_meta_l1,
        "pb_expr_l2": None if len(pb_expr_l2_list) == 0 else pseudobulk_expr_l2,
        "pb_meta_l2": None if len(pb_meta_l2_list) == 0 else pseudobulk_meta_l2,
        "mc_expr": metacell_expr,
        "mc_meta": metacell_meta,
    }


def main():
    ap = argparse.ArgumentParser(
        description="Export corrected and/or uncorrected pseudobulk + metacell tables with L1/L2 annotations."
    )

    ap.add_argument("--in_h5ad", required=True)
    ap.add_argument("--scanvi_model_dir", default=None,
                    help="Required only when --export_corrected true")
    ap.add_argument("--seacell_pkl_dir", required=True)
    ap.add_argument("--out_dir", required=True)

    # keys
    ap.add_argument("--counts_layer", default="counts")
    ap.add_argument("--batch_key", default="dataset_id")
    ap.add_argument("--labels_key", default="aifi_L1_paper")
    ap.add_argument("--labels_key_l2", default="aifi_L2_paper")
    ap.add_argument("--unlabeled_category", default="Unknown")
    ap.add_argument("--sample_key", default="sample_id")
    ap.add_argument("--patient_key", default="patient_id")
    ap.add_argument("--disease_key", default="disease_Level_3")
    ap.add_argument("--sex_key", default="sex")
    ap.add_argument("--age_key", default="binned_age")

    ap.add_argument("--cat_cov", nargs="*", default=["sex", "binned_age"])
    ap.add_argument("--cont_cov", nargs="*", default=[])

    # filtering
    ap.add_argument("--exclude_celltypes", nargs="*", default=["RBC", "Platelets", "Cycling_cells", "Progenitors"])
    ap.add_argument("--min_seacell_purity", type=float, default=0.75)
    ap.add_argument("--purity_level", choices=["L1", "L2"], default="L1")

    # corrected expression
    ap.add_argument("--n_samples_norm", type=int, default=25)
    ap.add_argument("--library_size", type=float, default=10000)

    # export switches
    ap.add_argument("--export_corrected", type=str2bool, default=True)
    ap.add_argument("--export_uncorrected", type=str2bool, default=True)
    ap.add_argument("--uncorrected_mode", choices=["lognorm", "counts"], default="lognorm")
    ap.add_argument("--write_backward_compatible", type=str2bool, default=False,
                    help="If true, also write corrected L1 outputs to out_dir/backward_compatible/")

    ap.add_argument("--seed", type=int, default=5)
    ap.add_argument("--num_workers", type=int, default=0)

    args = ap.parse_args()

    ensure_dir(args.out_dir)

    if (not args.export_corrected) and (not args.export_uncorrected):
        raise ValueError("At least one of --export_corrected or --export_uncorrected must be true.")

    if args.export_corrected and (args.scanvi_model_dir is None):
        raise ValueError("--scanvi_model_dir is required when --export_corrected true.")

    scvi.settings.seed = args.seed
    scvi.settings.dl_num_workers = args.num_workers

    print("[READ]", args.in_h5ad)
    adata = sc.read_h5ad(args.in_h5ad)

    # basic checks
    required = [
        args.batch_key, args.labels_key, args.sample_key,
        args.patient_key, args.disease_key, args.sex_key, args.age_key
    ]
    if args.labels_key_l2 is not None:
        required.append(args.labels_key_l2)

    if args.export_uncorrected and args.counts_layer not in adata.layers.keys():
        raise ValueError(f"[ERROR] counts layer not found: {args.counts_layer}")

    for k in required:
        if k not in adata.obs.columns:
            raise ValueError(f"[ERROR] obs missing: {k}")

    # categories
    cat_keys_all = [
        args.batch_key, args.labels_key, args.sample_key,
        args.sex_key, args.age_key
    ] + list(args.cat_cov)

    if args.labels_key_l2 is not None:
        cat_keys_all.append(args.labels_key_l2)

    for k in cat_keys_all:
        if k in adata.obs.columns:
            adata.obs[k] = adata.obs[k].astype("category")

    adata.obs[args.patient_key] = adata.obs[args.patient_key].astype(str)
    adata.obs[args.disease_key] = adata.obs[args.disease_key].astype(str)

    scanvi_model = None
    transform_batches = None

    if args.export_corrected:
        # 确保 unlabeled category 存在
        if args.unlabeled_category not in adata.obs[args.labels_key].cat.categories:
            adata.obs[args.labels_key] = adata.obs[args.labels_key].cat.add_categories([args.unlabeled_category])

        scvi.model.SCANVI.setup_anndata(
            adata,
            layer=args.counts_layer,
            batch_key=args.batch_key,
            labels_key=args.labels_key,
            unlabeled_category=args.unlabeled_category,
            categorical_covariate_keys=args.cat_cov if len(args.cat_cov) else None,
            continuous_covariate_keys=args.cont_cov if len(args.cont_cov) else None,
        )

        print("[LOAD MODEL]", args.scanvi_model_dir)
        scanvi_model = scvi.model.SCANVI.load(args.scanvi_model_dir, adata=adata)

        transform_batches = list(adata.obs[args.batch_key].cat.categories)
        print("[INFO] transform_batch =", transform_batches)

    print("[READ SEACELLS]", args.seacell_pkl_dir)
    seacell_df = load_all_seacell_pkls(args.seacell_pkl_dir)

    sample_ids = list(adata.obs[args.sample_key].astype(str).unique())
    print("[INFO] n_samples =", len(sample_ids))

    genes = adata.var_names.astype(str).tolist()

    # 每个 variant 都单独存 list
    variant_buffers = {}
    if args.export_corrected:
        variant_buffers["corrected"] = {
            "pb_expr_l1": [], "pb_meta_l1": [],
            "pb_expr_l2": [], "pb_meta_l2": [],
            "mc_expr": [], "mc_meta": []
        }
    if args.export_uncorrected:
        variant_buffers[f"uncorrected_{args.uncorrected_mode}"] = {
            "pb_expr_l1": [], "pb_meta_l1": [],
            "pb_expr_l2": [], "pb_meta_l2": [],
            "mc_expr": [], "mc_meta": []
        }

    for sample_id in sample_ids:
        print(f"\n[PROCESS] {sample_id}")

        obs = adata.obs.copy()
        obs_reset = obs.reset_index().rename(columns={"index": "cell_id"})

        # 仍以 L1 排除大类
        mask = (
            (~obs_reset[args.labels_key].astype(str).isin(args.exclude_celltypes))
            & (obs_reset[args.sample_key].astype(str) == str(sample_id))
        )
        obs_sub = obs_reset.loc[mask].copy()

        if obs_sub.shape[0] == 0:
            print("[SKIP] no cells after exclude_celltypes")
            continue

        # base meta
        base_meta_cols = [
            args.sample_key, args.patient_key, args.disease_key,
            args.sex_key, args.age_key, args.labels_key, args.batch_key
        ]
        if args.labels_key_l2 is not None:
            base_meta_cols.append(args.labels_key_l2)

        base_meta = obs_sub[base_meta_cols].copy()

        # --------------------------
        # 准备各 variant 的表达矩阵
        # --------------------------
        expr_variant_map = {}

        if args.export_corrected:
            expr_variant_map["corrected"] = get_corrected_expr(
                scanvi_model=scanvi_model,
                obs_sub=obs_sub,
                transform_batches=transform_batches,
                genes=genes,
                n_samples_norm=args.n_samples_norm,
                library_size=args.library_size
            )

        if args.export_uncorrected:
            expr_variant_map[f"uncorrected_{args.uncorrected_mode}"] = get_uncorrected_expr(
                adata=adata,
                obs_sub=obs_sub,
                genes=genes,
                counts_layer=args.counts_layer,
                library_size=args.library_size,
                mode=args.uncorrected_mode
            )

        # --------------------------
        # 先计算 SEACell + purity（与 variant 无关）
        # --------------------------
        common_cells = pd.Index(obs_sub["cell_id"].astype(str).values).intersection(seacell_df.index)
        obs_mc = obs_sub.set_index("cell_id").loc[common_cells].copy()

        if len(common_cells) > 0:
            obs_mc = obs_mc.merge(
                seacell_df[["SEACell"]],
                left_index=True,
                right_index=True,
                how="left"
            )
            obs_mc = obs_mc[~obs_mc["SEACell"].isna()].copy()

        if obs_mc.shape[0] > 0:
            purity_l1 = compute_group_purity(obs_mc, "SEACell", args.labels_key)
            if args.labels_key_l2 is not None:
                purity_l2 = compute_group_purity(obs_mc, "SEACell", args.labels_key_l2)
            else:
                purity_l2 = pd.Series(index=purity_l1.index, data=np.nan)

            purity_filter_series = purity_l1 if args.purity_level == "L1" else purity_l2
            keep_seacells = purity_filter_series.index[purity_filter_series >= args.min_seacell_purity]
            obs_mc = obs_mc[obs_mc["SEACell"].isin(keep_seacells)].copy()
        else:
            purity_l1 = pd.Series(dtype=float)
            purity_l2 = pd.Series(dtype=float)

        # --------------------------
        # 对每个 variant 分别导出
        # --------------------------
        for variant_name, expr_df in expr_variant_map.items():
            print(f"[VARIANT] {variant_name}")

            # pseudobulk L1
            pb_group_l1 = obs_sub[args.sample_key].astype(str) + "|" + obs_sub[args.labels_key].astype(str)

            pb_meta_l1 = base_meta.copy()
            pb_how = "mean"
            if variant_name == "uncorrected_counts":
                pb_how = "sum"

            pb_expr_l1, pb_meta_l1_out = aggregate_dense_matrix(
                expr_df=expr_df,
                group_series=pb_group_l1,
                meta_df=pb_meta_l1,
                how=pb_how
            )

            pb_meta_l1_out["celltype"] = pb_meta_l1_out[args.labels_key].astype(str)
            pb_meta_l1_out["celltype_L1"] = pb_meta_l1_out[args.labels_key].astype(str)

            if args.labels_key_l2 is not None and args.labels_key_l2 in pb_meta_l1_out.columns:
                pb_meta_l1_out["celltype_L2_majority"] = pb_meta_l1_out[args.labels_key_l2].astype(str)

            variant_buffers[variant_name]["pb_expr_l1"].append(pb_expr_l1)
            variant_buffers[variant_name]["pb_meta_l1"].append(pb_meta_l1_out)

            print(f"[OK] {variant_name} pseudobulk L1 rows += {pb_expr_l1.shape[0]}")

            # pseudobulk L2
            if args.labels_key_l2 is not None:
                pb_group_l2 = obs_sub[args.sample_key].astype(str) + "|" + obs_sub[args.labels_key_l2].astype(str)

                pb_meta_l2 = base_meta.copy()
                pb_expr_l2, pb_meta_l2_out = aggregate_dense_matrix(
                    expr_df=expr_df,
                    group_series=pb_group_l2,
                    meta_df=pb_meta_l2,
                    how=pb_how
                )

                pb_meta_l2_out["celltype"] = pb_meta_l2_out[args.labels_key_l2].astype(str)
                pb_meta_l2_out["celltype_L2"] = pb_meta_l2_out[args.labels_key_l2].astype(str)
                pb_meta_l2_out["celltype_L1_majority"] = pb_meta_l2_out[args.labels_key].astype(str)

                variant_buffers[variant_name]["pb_expr_l2"].append(pb_expr_l2)
                variant_buffers[variant_name]["pb_meta_l2"].append(pb_meta_l2_out)

                print(f"[OK] {variant_name} pseudobulk L2 rows += {pb_expr_l2.shape[0]}")

            # metacell
            if obs_mc.shape[0] == 0:
                print(f"[SKIP] {variant_name} no SEACells pass purity / overlap")
                continue

            expr_mc = expr_df.loc[obs_mc.index].copy()
            mc_group = obs_mc["SEACell"].astype(str)

            mc_meta_cols = [
                args.sample_key, args.patient_key, args.disease_key,
                args.sex_key, args.age_key, args.labels_key, args.batch_key
            ]
            if args.labels_key_l2 is not None:
                mc_meta_cols.append(args.labels_key_l2)

            mc_meta = obs_mc[mc_meta_cols].copy()

            mc_expr, mc_meta_out = aggregate_dense_matrix(
                expr_df=expr_mc,
                group_series=mc_group,
                meta_df=mc_meta,
                how="mean"
            )

            # 维持 corrected 旧逻辑：mean 后 log1p
            if variant_name == "corrected":
                mc_expr = np.log1p(mc_expr)
            elif variant_name == "uncorrected_counts":
                mc_expr = np.log1p(mc_expr)
            # uncorrected_lognorm: 已经是 log1p 过的每细胞表达，均值后不再 log1p

            mc_meta_out["celltype"] = mc_meta_out[args.labels_key].astype(str)
            mc_meta_out["celltype_L1"] = mc_meta_out[args.labels_key].astype(str)

            if args.labels_key_l2 is not None and args.labels_key_l2 in mc_meta_out.columns:
                mc_meta_out["celltype_L2"] = mc_meta_out[args.labels_key_l2].astype(str)

            mc_meta_out["purity_L1"] = purity_l1.reindex(mc_meta_out.index).values
            mc_meta_out["purity_L2"] = purity_l2.reindex(mc_meta_out.index).values

            variant_buffers[variant_name]["mc_expr"].append(mc_expr)
            variant_buffers[variant_name]["mc_meta"].append(mc_meta_out)

            print(f"[OK] {variant_name} metacell rows += {mc_expr.shape[0]}")

    # --------------------------
    # 导出全部 variant
    # --------------------------
    exported = {}
    for variant_name, buf in variant_buffers.items():
        exported[variant_name] = export_variant_outputs(
            variant_name=variant_name,
            out_root=args.out_dir,
            pb_expr_l1_list=buf["pb_expr_l1"],
            pb_meta_l1_list=buf["pb_meta_l1"],
            pb_expr_l2_list=buf["pb_expr_l2"],
            pb_meta_l2_list=buf["pb_meta_l2"],
            mc_expr_list=buf["mc_expr"],
            mc_meta_list=buf["mc_meta"],
        )

    # backward-compatible copy: corrected L1
    if args.write_backward_compatible and ("corrected" in exported):
        bdir = os.path.join(args.out_dir, "backward_compatible")
        ensure_dir(bdir)

        exported["corrected"]["pb_expr_l1"].to_csv(
            os.path.join(bdir, "pseudobulk_expr.tsv.gz"),
            sep="\t", compression="gzip"
        )
        exported["corrected"]["pb_meta_l1"].to_csv(
            os.path.join(bdir, "pseudobulk_meta.tsv"),
            sep="\t"
        )

        print("[WRITE]", os.path.join(bdir, "pseudobulk_expr.tsv.gz"), exported["corrected"]["pb_expr_l1"].shape)
        print("[WRITE]", os.path.join(bdir, "pseudobulk_meta.tsv"), exported["corrected"]["pb_meta_l1"].shape)

    print("\n[DONE]")
    print("Output root:", args.out_dir)
    for variant_name in exported:
        print(" -", variant_name, "=>", os.path.join(args.out_dir, variant_name))


if __name__ == "__main__":
    main()