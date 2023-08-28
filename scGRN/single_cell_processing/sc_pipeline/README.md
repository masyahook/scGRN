# COVID-19 Single cell data processing pipeline

This folder contains the workflow for processing the scRNA-seq matrix data of COVID-19 patients ([Liao *et al.*, 2020](https://www.nature.com/articles/s41591-020-0901-9)). The metadata describing the location of input data is in the `scGRN/single_cell_processing/sc_metadata.csv`. The data is processed using `Seurat` pipeline with `SingleR` package used to annotate obtained clusters of cells.

To run the `Seurat` pipeline for **only one** patient, please use in the simplest case:

```bash
    ./run_seurat_individual.sh <patient_ID> <path_to_metadata> <output_folder>
```

The `run_seurat_individual` could be run on a low-memory cluster (2 GB per core), still at least 20 cores are advised.

This will process scRNA-seq matrix for the patient `<patient_ID>`. To run the Seurat pipeline for **all** patients **in parallel**, please use in the simplest case:

```bash
    ./run_seurat.sh <path_to_metadata> <output_folder>
```

The `run_seurat` could be also run on a low-memory cluster (2 GB per core), 48 cores are advised. The script takes ~2-3 hours to run.

After running `run_seurat.sh` and processing each patient data separately, the next task is to merge patient data *cell type-wise* (e.g. all T cells from all patients into one dataset) to subsequently obtain gene regulatory networks for each *cell type*. Please use in the simplest case:

```bash
    ./run_seurat_merged.sh <path_to_metadata> <output_folder>
```

The `run_seurat_merged` requires high-mem nodes (8 GB per core, a full node with 48 cores is advised) as this pipeline loads into RAM all patient-level datasets, and then separates it into subsets of different cell types. `PRE_MERGED=T` could be used to rerun pipeline from intermediate output (i.e. first integration part of the pipeline was run before) and perform only cell-type specific analysis.

To run the **full** COVID-19 single cell processing pipeline (i.e. `run_seurat + run_seurat_merged`) that processes all patients and then merges data, please use in the simplest case:

```bash
  ./sc_pipeline.sh <path_to_metadata> <output_folder>
```

Also, for quick generation of input commands please use the `notebooks/Generate_sbatch_commands.ipynb` notebook, which will streamline the batch job scheduling on Slurm cluster.

The pipeline will produce several folders in the `<output_folder>` each corresponding to a patient-specific data, or cell type-specific data. Each patient-specific folder will be labeled by `<patient_ID>`, while cell type-specific folder will be stored in `cell_types` folder. The full file structure is described below:

```bash
<output_folder>
├── <patient_ID_1>
│   ├── data
│   │   └── Seurat
│   └── figs
│       └── Seurat
|-- <patient_ID_2>
│   ├── data
│   │   └── Seurat
│   └── figs
│       └── Seurat
|-- cell_types
│   ├── <cell_type_1>
│   │   ├── data
│   │   │   └── Seurat
│   │   └── figs
│   │       └── Seurat
│   └── <cell_type_2>
│       ├── data
│       │   └── Seurat
│       └── figs
│           └── Seurat
└── cell_type_meta.tsv
```

Each `Seurat` folder will contain data to corresponding patient/cell type. For more information about arbitrary input parameters and output format please look in the corresponding scripts.
