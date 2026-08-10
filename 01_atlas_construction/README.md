# 01 - Atlas construction

This stage builds and annotates the pediatric PBMC atlas used in figure 1 and all downstream analyses.

## Suggested order

1. `01_run_aifi_celltypist_annotation.py` - annotate per-dataset AnnData objects and perform doublet/QC handling.
2. `02-merge_and_paper_filter.py` - merge datasets, retain raw counts, and apply the paper gene filters.
3. `03-align_L1_L2_from_aifi.py` - harmonize AIFI labels into the manuscript Level 1 and Level 2 hierarchy.
4. `04-run_scvi_restore_main_only.py` - prepare the main atlas representation used for integration.
5. `05_add_clinical_to_h5ad-binned_age.py` - attach harmonized clinical covariates.
6. `06_prepare_object_for_downstream.py` - construct the downstream analysis object.
7. `07_run_scvi_scanvi.py` - fit scVI/scANVI models and export latent/corrected representations.
8. `run_integration_benchmark_scib.py` - compare PCA, Harmony, scVI, and scANVI representations.
9. `Figure1d_subCluster_circular.R` and `Figure1f_disease_Level_3.R` - generate figure 1 summaries.

The `feature_selection/` subdirectory contains participant-level pseudobulk differential expression, recurrent highly variable gene selection, curated gene mapping, and gene-universe integration steps.

Use `python <script> --help` for scripts with command-line interfaces. The annotation entry script keeps its main configuration near the top of the file.
