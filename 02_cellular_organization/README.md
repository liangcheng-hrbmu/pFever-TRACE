# 02 - Cellular organization

This stage reproduces the abundance and cellular-module analyses associated with figure 2.

- `Fig2a.R` computes observed-to-expected enrichment across diseases and broad immune compartments.
- `fig2b.R` summarizes participant-level cell fractions and disease differences.
- `fig2de.R` constructs the 52-state abundance matrix, fits the four-module NMF solution, and prepares participant/module panels.
- `fig2f.R` visualizes cell-state weights within each cellular module.
- `Result2_export_supplementary_tables.R` exports the corresponding source-data and supplementary tables.

Several plotting scripts expect `data/analysis/atlas_objects/integrated_pediatric.rds` and intermediate NMF objects produced by earlier steps. Run from the repository root so relative paths resolve consistently.
