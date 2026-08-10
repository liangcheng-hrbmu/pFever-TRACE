# Local data layout

This directory is intentionally excluded from version control except for this file. Place local datasets and generated intermediate objects here; do not commit controlled participant data.

The scripts currently reference the following broad layout:

```text
data/
|-- processed/
|-- processed_celltypist/
|-- models/
`-- analysis/
    |-- atlas_objects/
    |-- 03_integrated_umap_scvi_v4/
    |-- scANVI/
    |-- inflammation_signatures/
    |-- 08_gene_importance/
    |-- 08_patient_classifier/
    `-- 09_bulk_deconvolution/
```

Use de-identified, stable sample and participant identifiers. The expected metadata field names are defined in each script and commonly include `dataset_id`, `sampleID`, `patient_id`, `disease`, `sex`, and age-group fields.
