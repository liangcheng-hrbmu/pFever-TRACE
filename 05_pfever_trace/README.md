# 05 - pFever-TRACE

This stage builds participant-level representations, trains one disease classifier per immune compartment, combines predictions by majority vote, and transfers the framework to whole-blood RNA sequencing.

## Single-cell workflow

1. `single_cell/52_generate_patient_classifier_kfold_splits_keep_all_diseases.py`
2. `single_cell/53_generate_patient_pseudobulks_from_splits_v2.py`
3. `single_cell/54_train_patient_classifiers_and_majority_vote_v2.py`
4. `single_cell/plot_fig5cde_scenario1.R`
5. `single_cell/plot_figS5_scenario1_paper.R`

## Whole-blood workflow

1. `bulk/make_bulk_sample_metadata_with_mapping.py`
2. `bulk/Deconvolution.R`
3. `bulk/52b_generate_bulk_classifier_kfold_splits.py`
4. `bulk/53b_generate_bulk_pseudobulks_from_instaprism_fixed.py`
5. `bulk/54_train_patient_classifiers_and_majority_vote_v2.py`
6. The `Bulk Broad5`, `Bulk fine13`, and supplementary R scripts generate the reported benchmark panels and source data.

Participant identifiers must remain grouped within outer cross-validation folds. Whole-blood analyses require the InstaPrism outputs and the single-cell atlas reference described in the manuscript.
