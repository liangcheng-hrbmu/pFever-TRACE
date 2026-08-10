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


def expected_outputs(out_dir, celltype):
    prefix = f"{celltype}.studyID"
    return [
        os.path.join(out_dir, f"{prefix}.metrics.json"),
        os.path.join(out_dir, f"{prefix}.xgb_model.joblib"),
        os.path.join(out_dir, f"{prefix}.split_indices.npz"),
    ]


def already_done(out_dir, celltype):
    return all(os.path.exists(fp) for fp in expected_outputs(out_dir, celltype))


def main():
    ap = argparse.ArgumentParser(
        description="Run 19_xgboost_fit_one_celltype_studyID2.py for all cell types."
    )
    ap.add_argument(
        "--input_dir",
        default="data/analysis/08_gene_importance/16_SHAP_input_data",
        help="Directory containing *_adataMerged_SPECTRAgenes.log1p.h5ad files"
    )
    ap.add_argument(
        "--fit_script",
        default=os.path.join(SCRIPT_DIR, "19_xgboost_fit_one_celltype_studyID2.py"),
        help="Path to the one-celltype studyID fitting script"
    )
    ap.add_argument(
        "--out_dir",
        default="data/analysis/08_gene_importance/19_xgboost_results_studyID",
        help="Output directory"
    )
    ap.add_argument("--label_key", default="studyID")
    ap.add_argument("--group_key", default="sampleID")
    ap.add_argument("--random_state", type=int, default=1234)
    ap.add_argument("--n_trials", type=int, default=50)
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
        "--force",
        action="store_true",
        help="Rerun even if expected outputs already exist"
    )
    ap.add_argument(
        "--use_gpu",
        action="store_true",
        help="Pass --use_gpu to the one-celltype script"
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

        if (not args.force) and already_done(args.out_dir, celltype):
            skipped.append((celltype, "already_done"))
            continue

        todo.append((celltype, h5ad))

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

    for i, (celltype, h5ad) in enumerate(todo, start=1):
        cmd = [
            args.python_exec,
            args.fit_script,
            "--h5ad", h5ad,
            "--out_dir", args.out_dir,
            "--label_key", args.label_key,
            "--group_key", args.group_key,
            "--random_state", str(args.random_state),
            "--n_trials", str(args.n_trials),
        ]
        if args.use_gpu:
            cmd.append("--use_gpu")

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
