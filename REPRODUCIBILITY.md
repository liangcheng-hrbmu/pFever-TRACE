# Reproducibility notes

## Scope

This repository contains the analysis scripts supplied for the manuscript. It does not contain the public source datasets, controlled participant metadata, large intermediate objects, fitted models, or generated figures.

## Required local assets

- Public single-cell and bulk RNA-sequencing datasets listed in the manuscript supplementary tables.
- AIFI/CellTypist reference models and color maps used during annotation.
- Curated inflammation gene sets used by Spectra.
- Intermediate AnnData, Seurat, SEACell, scVI/scANVI, classifier, and SHAP objects produced by earlier stages.
- Metadata with stable dataset, sample, participant, disease, age, and sex fields.

The scripts now resolve the former server layout relative to `data/`. The historical analysis hierarchy is retained below `data/analysis/` so that path-dependent plotting scripts remain traceable.

## Known missing worker scripts

The supplied stage 4 folder contains orchestration scripts that call the following single-compartment workers, but those worker files were not present in the source folder:

- `17_xgboost_fit_one_celltype.py`
- `18_compute_shap_one_celltype.py`
- `19_xgboost_fit_one_celltype_studyID2.py`
- `20_compute_shap_one_celltype_studyID.py`
- `21_compare_dshap_sshap_one_celltype.py`

Add these files to `04_gene_importance/` or pass their locations through the corresponding command-line options before running those orchestration steps.

## Execution assumptions

- Run commands from the repository root unless a directory README says otherwise.
- Use participant-grouped partitions. Never split cells from the same participant across training and evaluation folds.
- Preserve the random seeds and disease-label mappings recorded in the scripts.
- GPU acceleration is recommended for scVI/scANVI, XGBoost, SHAP, SEACells, and Spectra stages.
- Capture the final Python environment and R `sessionInfo()` used for the paper release.

## Validation completed for this archive

- Python files were checked for syntax without importing optional scientific packages.
- Hard-coded references to the original personal server root were removed.
- The repository was scanned for common credential strings and large data files.
- R syntax could not be parsed in the preparation environment because R was not installed.
- Full numerical reproduction was not possible because the input data and intermediate objects were not supplied.
