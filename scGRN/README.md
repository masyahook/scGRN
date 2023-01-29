# Organization of the `scGRN` package

Here the user can find the organization of the whole package:

- `single_cell_processing` - all scripts, pipelines and modules that process single-cell datasets, containing mostly analysis that is done prior to network inference
- `network_inference` - all scripts, pipelines and modules that are required for gene regulatory network inference (GRN) based on `pySCENIC` ([Van de Sande *et al.*, 2020](https://www.nature.com/articles/s41596-020-0336-2)) package
- `network_analysis` - analysis of computed GRNs, includes exploration of network structure, community analysis and gene set enrichment analysis (GSEA)