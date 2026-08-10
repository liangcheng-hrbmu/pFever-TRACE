# 03 - Inflammation programs

This stage derives the metacell representation, learns inflammation-related Spectra factors, scores their activity, and generates figure 3.

## Workflow

1. `seacells/split_for_seacells.py`
2. `seacells/SEAcellDefinition.py`
3. `seacells/patientNorm_to_global_tables.py`
4. `spectra/11_run_spectra_from_metacells.py`
5. `spectra/12_select_spectra_factors.py`
6. `ulm_lmem/13_run_ULMsignatures_from_pseudobulk.py`
7. `ulm_lmem/14_run_LMEM_L1_corr2.py`
8. `ulm_lmem/15_run_ULMsignatures_from_pseudobulk_L2_all.py`
9. `ulm_lmem/16_run_LMEM_L2_corr_all.py`
10. `fig3a.R`, `fig3b.R`, and `fig3cde.R`

The Spectra scripts require the curated inflammation gene-set file described in the manuscript. ULM outputs feed the Level 1 and Level 2 mixed-effects analyses and then the figure scripts.
