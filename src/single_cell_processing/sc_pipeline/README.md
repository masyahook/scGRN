# COVID-19 Single cell processing pipeline

This folder contains the workflow for processing the scRNA-seq matrix data of COVID-19 patients ([Liao *et al.*, 2020](https://www.nature.com/articles/s41591-020-0901-9)). The metadata describing the location of input data is in the `src/single_cell_processing/sc_metadata.csv`. The data is processed using Seurat pipeline with SingleR package to annotate obtained clusters of cells.

To run the Seurat pipeline for **only one** patient, please use in the simplest case:

```
    ./run_seurat_individual.sh <patient_ID> <path_to_metadata> <output_folder>
```

This will process scRNA-seq matrix for the patient `<patient_ID>`. To run the Seurat pipeline for **all** patients **in parallel**, please use in the simplest case:

```
    ./run_seurat.sh <path_to_metadata> <output_folder>
```

After running `single_cell_processing_individual.R` and processing each patient data separately, the next task is to merge patient data *cell type-wise* (e.g. all T cells from all patients into one dataset) to subsequently obtain gene regulatory networks for each *cell type*. Please use in the simplest case:

```
    ./merge_data_by_cell_type.sh <patient_ID> <path_to_metadata> <output_folder>
```

To run the **full** COVID-19 single cell processing pipeline that processes all patients and then merges data, please use in the simplest case:

```
  ./sc_pipeline.sh <path_to_metadata> <output_folder>
```

For more information about arbitrary input parameters and output format please look in the corresponding scripts.
