# Processing of single cell data

This folder contains 2 workflows of scRNA-seq data of COVID-19 patients ([Liao *et al.*, 2020](https://www.nature.com/articles/s41591-020-0901-9)):
- `sc_pipeline` - **COVID-19 single cell data processing pipeline** that contains the typical Seurat processing and downstream cell type data merging to create aggregated single-cell datasets
- `regulon_pipeline` - **regulon activity inference pipeline** that contains the `VIPER` regulon enrichment pipeline based on single cell data using known (i.e. `DoRothEA`) or inferred (i.e. `pySCENIC`) regulons

Also in this folder one could find `sc_metadata.tsv` tabular file that contains various metadata about raw input files including location and patient type.

The user could obtain the data using [this link](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE145926).