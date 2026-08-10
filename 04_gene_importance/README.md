# 04 - Gene importance

This stage trains compartment-specific disease and study classifiers, computes SHAP values, and prioritizes genes with high disease-associated and low study-associated importance.

## Workflow

1. `16_generate_SHAP_input_from_scanvi.py` exports one input object per immune compartment.
2. `17_xgboost_fit_remaining_celltypes.py` orchestrates disease classifier fitting.
3. `18_compute_shap_all_celltypes.py` orchestrates disease SHAP calculation.
4. `19_xgboost_fit_all_celltypes_studyID2.py` fits matched study-identity classifiers.
5. `20_compute_shap_all_celltypes_studyID.py` calculates study SHAP values.
6. `21_compare_dshap_sshap_all_celltypes.py` compares disease and study importance.
7. `22_summarize_xgboost_shap_results.py` produces cross-compartment summaries.
8. `42_export_corrected_pred_test_from_models_v2.py` and `43_fit_unintegrated_xgboost_from_mainobject_same_splits_v4.py` support corrected-versus-uncorrected comparisons.
9. `fig4a.py`, `fig4a.R`, and `fig4b.py` generate figure panels.

The five orchestration scripts in steps 2 to 6 require companion single-compartment worker scripts that were not included in the supplied folder. See the explicit list in the repository-level `REPRODUCIBILITY.md`.
