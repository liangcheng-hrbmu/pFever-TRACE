#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import glob
import argparse
import zipfile
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import confusion_matrix

SUFFIX = '_adataMerged_SPECTRAgenes.log1p.h5ad'


def ensure_dir(p):
    os.makedirs(p, exist_ok=True)


def celltype_from_h5ad(path: str) -> str:
    base = os.path.basename(path)
    if base.endswith(SUFFIX):
        return base[:-len(SUFFIX)]
    return os.path.splitext(base)[0]


def read_pred_zip(zip_fp: str) -> pd.DataFrame:
    with zipfile.ZipFile(zip_fp, 'r') as zf:
        names = [n for n in zf.namelist() if n.endswith('.tsv') or n.endswith('.csv')]
        if not names:
            raise FileNotFoundError(f'No tsv/csv found inside {zip_fp}')
        name = names[0]
        with zf.open(name) as f:
            if name.endswith('.tsv'):
                df = pd.read_csv(f, sep='\t')
            else:
                df = pd.read_csv(f)
    return df


def normalize_pred_columns(df: pd.DataFrame, celltype: str) -> pd.DataFrame:
    cols = {c.strip(): c for c in df.columns}
    true_col = None
    pred_col = None
    for c in ['true_label', 'y_true', 'true', 'label_true']:
        if c in cols:
            true_col = cols[c]
            break
    for c in ['pred_label', 'disease_pred', 'y_pred_label', 'y_pred', 'pred', 'label_pred']:
        if c in cols:
            pred_col = cols[c]
            break
    if true_col is None or pred_col is None:
        raise ValueError(f'{celltype}: cannot infer true/pred columns from {list(df.columns)}')
    out = pd.DataFrame({
        'celltype': celltype,
        'true_label': df[true_col].astype(str),
        'pred_label': df[pred_col].astype(str),
    })
    return out


def aggregate_predictions(pred_dir: str, corrected_input_dir: str | None = None) -> pd.DataFrame:
    rows = []
    if corrected_input_dir is not None:
        h5s = sorted(glob.glob(os.path.join(corrected_input_dir, f'*{SUFFIX}')))
        celltypes = [celltype_from_h5ad(x) for x in h5s]
    else:
        zips = sorted(glob.glob(os.path.join(pred_dir, 'run1_*_pred_test.zip')))
        celltypes = [os.path.basename(z).replace('run1_', '').replace('_pred_test.zip', '') for z in zips]
    for ct in celltypes:
        zip_fp = os.path.join(pred_dir, f'run1_{ct}_pred_test.zip')
        if not os.path.exists(zip_fp):
            continue
        df = read_pred_zip(zip_fp)
        rows.append(normalize_pred_columns(df, ct))
    if not rows:
        raise RuntimeError(f'No prediction zip files found/parsed in {pred_dir}')
    out = pd.concat(rows, ignore_index=True)
    out['true_label'] = out['true_label'].astype(str)
    out['pred_label'] = out['pred_label'].astype(str)
    return out


def cm_col_norm(df: pd.DataFrame, label_order: list[str]) -> pd.DataFrame:
    present = [x for x in label_order if x in set(df['true_label']) or x in set(df['pred_label'])]
    cm = confusion_matrix(df['true_label'], df['pred_label'], labels=present)
    cm = pd.DataFrame(cm, index=present, columns=present).T
    cm = cm.div(cm.sum(axis=0).replace(0, np.nan), axis=1)
    return cm


def plot_pair(cm1: pd.DataFrame, cm2: pd.DataFrame, label_order: list[str], out_fp: str):
    present = [x for x in label_order if x in cm1.index or x in cm2.index or x in cm1.columns or x in cm2.columns]
    cm1p = cm1.reindex(index=present, columns=present)
    cm2p = cm2.reindex(index=present, columns=present)

    fig, axes = plt.subplots(1, 2, figsize=(14.5, 6), constrained_layout=True)
    titles = ['XGBoost trained on\nscANVI-corrected expression', 'XGBoost trained on\nnon-corrected expression']
    for i, (ax, cm, title) in enumerate(zip(axes, [cm1p, cm2p], titles)):
        sns.heatmap(
            cm,
            ax=ax,
            cmap='coolwarm',
            vmin=0,
            vmax=1,
            square=True,
            cbar=(i == 1),
            cbar_kws={'label': 'Proportion\n(col. norm.)'},
            linewidths=0.2,
            linecolor='white'
        )
        ax.set_title(title, fontsize=12)
        ax.set_xlabel('True conditions')
        ax.set_ylabel('Predicted conditions')
        ax.tick_params(axis='x', rotation=90, labelsize=8)
        ax.tick_params(axis='y', rotation=0, labelsize=8)
    fig.savefig(out_fp, dpi=300, bbox_inches='tight')
    print('[WRITE]', out_fp)


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description='Plot Fig. 3a from corrected and unintegrated pred_test.zip files.')
    ap.add_argument('--corrected_pred_dir', required=True)
    ap.add_argument('--raw_pred_dir', required=True)
    ap.add_argument('--out_dir', required=True)
    ap.add_argument('--corrected_input_dir', default=None)
    ap.add_argument('--label_order', default='Pediatric healthy,Dengue human infection(DENV-1),HFMD(EV-A71)enterovirus,Influenza A virus (IAV),severe acute respiratory syndrome coronavirus 2 (SARS-CoV-2),bacterial infection,Pediatric Sepsis,Kawasaki disease (KD),Multisystem inflammatory syndrome in children (MIS-C),human systemic lupus erythematosus (SLE),juvenile idopathic arthritis (JIA)')
    args = ap.parse_args()

    ensure_dir(args.out_dir)
    label_order = [x.strip() for x in args.label_order.split(',') if x.strip()]

    corr = aggregate_predictions(args.corrected_pred_dir, args.corrected_input_dir)
    raw = aggregate_predictions(args.raw_pred_dir, args.corrected_input_dir)

    corr_cm = cm_col_norm(corr, label_order)
    raw_cm = cm_col_norm(raw, label_order)

    corr.to_csv(os.path.join(args.out_dir, 'aggregated_predictions_corrected.tsv'), sep='\t', index=False)
    raw.to_csv(os.path.join(args.out_dir, 'aggregated_predictions_unintegrated.tsv'), sep='\t', index=False)
    corr_cm.to_csv(os.path.join(args.out_dir, 'confusion_matrix_corrected_colnorm.tsv'), sep='\t')
    raw_cm.to_csv(os.path.join(args.out_dir, 'confusion_matrix_unintegrated_colnorm.tsv'), sep='\t')

    plot_pair(corr_cm, raw_cm, label_order, os.path.join(args.out_dir, 'Fig3a_corrected_vs_unintegrated_confusion.png'))
    plot_pair(corr_cm, raw_cm, label_order, os.path.join(args.out_dir, 'Fig3a_corrected_vs_unintegrated_confusion.pdf'))

    print('\n[DONE]')
