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


def input_sum_files(dshap_dir, sshap_dir, celltype, use_test_only=False):
    tag = "testOnly" if use_test_only else "allCells"
    dshap_fp = os.path.join(dshap_dir, f"{celltype}.{tag}.dSHAP_sum.tsv")
    sshap_fp = os.path.join(sshap_dir, f"{celltype}.{tag}.studyID.sSHAP_sum.tsv")
    return dshap_fp, sshap_fp


def compare_prefix(celltype, use_test_only=False):
    tag = "testOnly" if use_test_only else "allCells"
    return f"{celltype}.{tag}.studyID"


def expected_outputs(out_dir, celltype, use_test_only=False):
    prefix = compare_prefix(celltype, use_test_only)
    return [
        os.path.join(out_dir, f"{prefix}.dSHAP_vs_sSHAP.tsv"),
        os.path.join(out_dir, f"{prefix}.top20_priority_genes.tsv"),
        os.path.join(out_dir, f"{prefix}.highD_lowS_genes.tsv"),
    ]


def already_done(out_dir, celltype, use_test_only=False):
    return all(os.path.exists(fp) for fp in expected_outputs(out_dir, celltype, use_test_only))


def main():
    ap = argparse.ArgumentParser(
        description="Run 21_compare_dshap_sshap_one_celltype.py for all cell types."
    )
    ap.add_argument(
        "--input_dir",
        default="data/analysis/08_gene_importance/16_SHAP_input_data",
        help="Directory containing *_adataMerged_SPECTRAgenes.log1p.h5ad files"
    )
    ap.add_argument(
        "--compare_script",
        default=os.path.join(SCRIPT_DIR, "21_compare_dshap_sshap_one_celltype.py"),
        help="Path to the one-celltype compare script"
    )
    ap.add_argument(
        "--dshap_dir",
        default="data/analysis/08_gene_importance/18_SHAP_results",
        help="Directory containing disease dSHAP summary files"
    )
    ap.add_argument(
        "--sshap_dir",
        default="data/analysis/08_gene_importance/20_SHAP_results_studyID",
        help="Directory containing studyID sSHAP summary files"
    )
    ap.add_argument(
        "--out_dir",
        default="data/analysis/08_gene_importance/21_compare_dshap_sshap_results",
        help="Output directory for comparison results"
    )
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
        help="Compare testOnly outputs instead of allCells outputs"
    )
    ap.add_argument(
        "--force",
        action="store_true",
        help="Rerun even if expected outputs already exist"
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

        dshap_fp, sshap_fp = input_sum_files(
            args.dshap_dir,
            args.sshap_dir,
            celltype,
            use_test_only=args.use_test_only,
        )
        if not (os.path.exists(dshap_fp) and os.path.exists(sshap_fp)):
            skipped.append((celltype, "missing_dshap_or_sshap_sum"))
            continue

        if (not args.force) and already_done(args.out_dir, celltype, args.use_test_only):
            skipped.append((celltype, "already_done"))
            continue

        todo.append((celltype, dshap_fp, sshap_fp))

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

    for i, (celltype, dshap_fp, sshap_fp) in enumerate(todo, start=1):
        celltype_name = compare_prefix(celltype, use_test_only=args.use_test_only)
        cmd = [
            args.python_exec,
            args.compare_script,
            "--dshap_sum", dshap_fp,
            "--sshap_sum", sshap_fp,
            "--out_dir", args.out_dir,
            "--celltype_name", celltype_name,
        ]

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
