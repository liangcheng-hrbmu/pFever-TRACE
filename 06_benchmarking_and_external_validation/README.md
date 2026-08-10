# 06 - Benchmarking and external validation

This stage compares pFever-TRACE with the published 139-transcript multiclass blood RNA signature and evaluates transfer to independent cohorts.

- `01_internal_benchmark/` contains the matched broad five-category and fine 13-disease internal comparisons.
- `02_external_evaluation/` contains participant-disjoint external evaluation, per-compartment metrics, confusion matrices, and the pFever-TRACE versus reference-signature comparison.

Both subdirectories expect outputs produced by `05_pfever_trace/` and local metadata under `data/analysis/09_bulk_deconvolution/`. The shell runner in `02_external_evaluation/` generates the three external-validation figure groups.
