# Processing of single cell data

This folder contains 2 workflows of scRNA-seq data of COVID-19 patients ([Liao *et al.*, 2020](https://www.nature.com/articles/s41591-020-0901-9)):
- `sc_pipeline` - **COVID-19 single cell data processing pipeline** that contains the Seurat processing and cell type data merging
- `regulon_pipeline` - **regulon activity inference pipeline** that contains the `VIPER` regulon enrichment pipeline based on single cell data using known (i.e. `DoRothEA`) or inferred (i.e. `pySCENIC`) regulons

The user could obtain the data using [this link](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE145926).