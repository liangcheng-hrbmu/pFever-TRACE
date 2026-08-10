# External evaluation

## Scripts

1. `plot_patient_disjoint_broad5_withLinearSVC_main.R`
   - Overall weighted F1 comparison for pFever-TRACE, the Med139 baseline, a random baseline, and the 11 Level 1 compartment models.
   - Disease-by-compartment F1 heatmap.
   - Column-normalized confusion matrix.
2. `plot_patient_disjoint_broad5_withLinearSVC_supp.R`
   - Balanced accuracy and Matthews correlation coefficient summaries.
   - Disease-by-compartment precision and recall heatmaps.
3. `plot_patient_disjoint_222_ours_vs_145_withLinearSVC_final.R`
   - Participant-disjoint comparison of pFever-TRACE and the Med139 ridge baseline.

## Run

```bash
bash 06_benchmarking_and_external_validation/02_external_evaluation/run_all_patient_disjoint_plots.sh
```

The scripts expect pFever-TRACE predictions, de-identified evaluation metadata, and Med139 predictions beneath `data/analysis/09_bulk_deconvolution/`. Each script writes source-data TSV files together with PDF and 300 dpi PNG figures.
