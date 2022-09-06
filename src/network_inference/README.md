# Network inference pipeline

This folder contains the workflow for computing the gene regulatory network (GRN) based on processed (Seurat) scRNA-seq matrix data. The user can use either `grnboost2` or `genie3` methods for computing GRN. The workflow uses `pySCENIC` ([Van de Sande *et al.*, 2020](https://www.nature.com/articles/s41596-020-0336-2)) package for to parallelize the GRN inference. It is recommended to use `grnboost2` as this method is less memory-heavy and has comparable performance with `geneie3`. The user can compute GRNs in two modes:

- Computing **all possible gene-gene connections** (only scRNA-seq data is needed). In this case the GRN method will detect co-regulatory genes, produce a list of adjacencies and compute pairwise Spearman correlation. As a result the following files are produced:
  - `<data_name>.tsv` - initial list of adjacencies between genes
  - `<data_name>_cor.tsv` - list of adjacencies with computed Spearman correlation
- Computing only **TF-target connections** with additional filtering based on motif enrichment (additional TF and motif data is needed). In this case the GRN method will detect regulons, produce a list of adjacencies, compute Spearman correlation and then filter the obtained connections based on motif enrichment. As a result the following files are produced:
  - `<data_name>_TF.tsv` - initial list of adjacencies between TFs and targets
  - `<data_name>_TF_cor.tsv` - list of adjacencies with computed Spearman correlation
  - `<data_name>_TF_ctx.tsv` - filtered list of adjacencies based on motif enrichment

To run GRN inference for **all gene-gene connections** for a given dataset, please use:

```
  ./infer_GRN.sh <"grnboost2"|"genie3"> <path_to_data> <num_workers> <q_threshold> <logging_folder_path>
```

To run GRN inference for **TF-target only connections** for a given dataset, please use:

```
  ./infer_GRN.sh <"grnboost2"|"genie3"> <path_to_data> <num_workers> <q_threshold> <logging_folder_path> <path_to_tf_file> <path_to_reg_feature_db> <path_to_motif_annotation_data>
```

Alternatively, the user can use custom scripts in `custom_scripts` folder where one can specify the type of data to focus on. For example, to run GRN inference **for a specific patient**, please use:

```
  ./custom_scripts/infer_pat_GRN.sh <"grnboost2"|"genie3"> <patient_ID> <cell_type_ID> <num_workers> <q_threshold> <path_to_tf_file - if in TF-target mode>
```

Similarly, to run GRN inference for **cell type aggregated data**, please use:

```
  ./custom_scripts/infer_agg_GRN.sh <"grnboost2"|"genie3"> <cell_type_ID> <pat_type_ID> <num_workers> <q_threshold> <path_to_tf_file - if run in TF-target mode>
```

In case the user wants to infer **only TF-target enriched** connections, it is **required** to pass **auxiliary data**:

- A list of transcription factors (`path_to_tf_file` parameter) 
- Pre-calculated whole-genome rankings (`path_to_reg_feature_db` parameter)
- Motif annotations (`path_to_motif_annotation_data` parameter)

All these data resources could be obtained from [here](https://resources.aertslab.org/cistarget/). The user could also find additional explanations in the original [`pySCENIC` workflow paper](https://www.nature.com/articles/s41596-020-0336-2#Sec32).

The optimal use of `pySCENIC` package requires careful memory and CPU allocation for each parallelized process. The GRN inference is **memory-heavy** and very high parallelization level could lead to memory overload and stalled processes. The general recommendations regarding choosing `num_workers` parameter (number of parallel processes):

- 8 GB RAM per worker for datasets `num_cells < 10,000`
- 12 GB RAM per worker for datasets `10,000 < num_cells < 20,000`
- 16 GB RAM per worker for datasets `num_cells > 20,000`

For example, the execution on a 96 GB RAM node (48 cores x 2 GB RAM) the recommended `num_workers = 8` by default. For high-memory nodes (384 GB RAM = 48 cores x 8 GB RAM): `num_workers = 32`.

For more information about arbitrary input parameters and output format please look in the corresponding scripts.