#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import glob
import argparse
import subprocess
import sys


SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
H5AD_SUFFIX = "_adataMerged_SPECTRAgenes.log1p.h5ad"


def celltype_from_h5ad(path):
    base = os.path.basename(path)
    if base.endswith(H5AD_SUFFIX):
        return base[:-len(H5AD_SUFFIX)]
    return os.path.splitext(base)[0]


def required_model_files(model_dir, celltype):
    model_fp = os.path.join(model_dir, f"{celltype}.xgb_model.joblib")
    split_fp = os.path.join(model_dir, f"{celltype}.split_indices.npz")
    return model_fp, split_fp


def expected_shap_outputs(out_dir, celltype, use_test_only=False):
    tag = "testOnly" if use_test_only else "allCells"
    return [
        os.path.join(out_dir, f"{celltype}.{tag}.dSHAP_by_disease.tsv"),
        os.path.join(out_dir, f"{celltype}.{tag}.dSHAP_sum.tsv"),
        os.path.join(out_dir, f"{celltype}.{tag}.cells_metadata.tsv"),
        os.path.join(out_dir, f"{celltype}.{tag}.shap_values.npy"),
        os.path.join(out_dir, f"{celltype}.{tag}.classes.json"),
    ]


def already_done(out_dir, celltype, use_test_only=False):
    return all(os.path.exists(fp) for fp in expected_shap_outputs(out_dir, celltype, use_test_only))


def main():
    ap = argparse.ArgumentParser(
        description="Run 18_compute_shap_one_celltype.py for all cell types."
    )
    ap.add_argument(
        "--input_dir",
        default="data/analysis/08_gene_importance/16_SHAP_input_data",
        help="Directory containing *_adataMerged_SPECTRAgenes.log1p.h5ad files"
    )
    ap.add_argument(
        "--compute_script",
        default=os.path.join(SCRIPT_DIR, "18_compute_shap_one_celltype.py"),
        help="Path to the one-celltype SHAP script"
    )
    ap.add_argument(
        "--model_dir",
        default="data/analysis/08_gene_importance/17_xgboost_results",
        help="Directory containing *.xgb_model.joblib and *.split_indices.npz"
    )
    ap.add_argument(
        "--out_dir",
        default="data/analysis/08_gene_importance/18_SHAP_results",
        help="Output directory for SHAP results"
    )
    ap.add_argument("--disease_key", default="disease")
    ap.add_argument("--sample_key", default="sampleID")
    ap.add_argument("--python", dest="python_exec", default=sys.executable)
    ap.add_argument(
        "--exclude",
        nargs="*",
        default=[],
        help="Cell types to skip, e.g. Mono B Plasma"
    )
    ap.add_argument(
        "--only",
        nargs="*",
        default=None,
        help="Only run these cell types"
    )
    ap.add_argument(
        "--use_test_only",
        action="store_true",
        help="Pass --use_test_only to the one-celltype SHAP script"
    )
    ap.add_argument(
        "--force",
        action="store_true",
        help="Rerun even if all expected SHAP outputs already exist"
    )
    ap.add_argument(
        "--dry_run",
        action="store_true",
        help="Only print commands, do not execute"
    )
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    pattern = os.path.join(args.input_dir, f"*{H5AD_SUFFIX}")
    h5ad_files = sorted(glob.glob(pattern))
    if not h5ad_files:
        raise FileNotFoundError(f"No h5ad files found: {pattern}")

    exclude_set = set(args.exclude)
    only_set = set(args.only) if args.only else None

    todo = []
    skipped = []

    for h5ad in h5ad_files:
        celltype = celltype_from_h5ad(h5ad)

        if only_set is not None and celltype not in only_set:
            skipped.append((celltype, "not_in_only"))
            continue

        if celltype in exclude_set:
            skipped.append((celltype, "excluded"))
            continue

        model_fp, split_fp = required_model_files(args.model_dir, celltype)
        if not (os.path.exists(model_fp) and os.path.exists(split_fp)):
            skipped.append((celltype, "missing_model_or_split"))
            continue

        if (not args.force) and already_done(args.out_dir, celltype, args.use_test_only):
            skipped.append((celltype, "already_done"))
            continue

        todo.append((celltype, h5ad, model_fp, split_fp))

    print(f"[INFO] found {len(h5ad_files)} total cell types")
    print(f"[INFO] to run: {len(todo)}")
    print(f"[INFO] skipped: {len(skipped)}")

    if skipped:
        for celltype, reason in skipped:
            print(f"[SKIP] {celltype}: {reason}")

    if not todo:
        print("[DONE] nothing to run")
        return

    failed = []

    for i, (celltype, h5ad, model_fp, split_fp) in enumerate(todo, start=1):
        cmd = [
            args.python_exec,
            args.compute_script,
            "--h5ad", h5ad,
            "--model_joblib", model_fp,
            "--split_npz", split_fp,
            "--out_dir", args.out_dir,
            "--disease_key", args.disease_key,
            "--sample_key", args.sample_key,
        ]
        if args.use_test_only:
            cmd.append("--use_test_only")

        print("\n" + "=" * 80)
        print(f"[RUN {i}/{len(todo)}] {celltype}")
        print("[CMD]", " ".join(cmd))

        if args.dry_run:
            continue

        ret = subprocess.run(cmd)
        if ret.returncode != 0:
            failed.append(celltype)
            print(f"[ERROR] failed: {celltype} (return code {ret.returncode})")

    print("\n" + "=" * 80)
    print(f"[SUMMARY] total={len(h5ad_files)} run={len(todo)} skipped={len(skipped)} failed={len(failed)}")
    if failed:
        print("[FAILED]", ", ".join(failed))
        raise SystemExit(1)
    print("[DONE] all requested cell types finished")


if __name__ == "__main__":
    main()
