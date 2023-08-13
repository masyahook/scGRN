# Regulon activity inference pipeline

This folder contains the workflow for TF-target **regulon activity** using VIPER algorithm ([Alvarez *et al.*, 2016](https://www.nature.com/articles/ng.3593)). This tool uses the scRNA-sequencing data and provided TF-target connections (i.e. regulons) to evaluate the **activity** (or **enrichment**) of transcription factors or rather corresponding TF regulons. The user could use two types of input regulon data:

- `DoRothEA` regulons - database of human and mouse regulons ([article](https://genome.cshlp.org/content/29/8/1363), [package](https://saezlab.github.io/dorothea/))
- `pySCENIC` regulons - regulons inferred from the **same** scRNA-sequencing data (by running `scGRN/network_inference` pipeline using [`pySCENIC`](https://pyscenic.readthedocs.io/en/latest/))

***Please note: this pipeline assumes that the user ran `scGRN/single_cell_processing/sc_pipeline` (i.e. Seurat pipeline) prior to this analysis.***

To run the VIPER pipeline for **only one** patient, please use in the simplest case:

```
    ./run_viper_individual.sh <patient_ID> <path_to_metadata> <output_folder> <pyscenic|dorothea|pyscenic_dorothea>
```

This command will run VIPER for the patient `<patient_ID>` and output a matrix of regulon activities (either based `pyscenic`, `dorothea` or both) with shape (`num_regulons` x `num_cells`). To run the VIPER pipeline for **all** patients, please use in the simplest case:

```
    ./run_viper.sh <path_to_metadata> <output_folder> <pyscenic|dorothea|pyscenic_dorothea>
```

To run VIPER on merged patient data **cell type-wise**, please use in the simplest case:

```
    ./run_viper_merged.sh <path_to_metadata> <output_folder> <pyscenic|dorothea|pyscenic_dorothea>
```

All VIPER pipelines produce `.tsv` files where VIPER scores are stored in a form of matrix for each sample. To speed up data loading we additionally convert `.tsv` files to `.pickle` format. To do that, please run:

```bash
    python all_viper_to_pickle.py -o <output_folder> -r <pyscenic|dorothea|pyscenic_dorothea>
```

And finally, to run the **full** VIPER regulon activity pipeline (first - infer activity based on each patient data, then - on merged data), please use in the simplest case:

```
    ./viper_pipeline.sh <path_to_metadata> <output_folder> <pyscenic|dorothea|pyscenic_dorothea>
```

The pipeline will produce a `regulon` subfolder in the corresponding Seurat-produced data folder:

```
<output_folder>
|-- <patient_ID_1>
|   |-- data
|   |   `-- Seurat
|   |       `-- regulon
|   |           `-- pickle
|   `-- figs
|       `-- Seurat
|           `-- regulon
|-- cell_types
|   |-- <cell_type_1>
|   |   |-- data
|   |   |   `-- Seurat
|   |   |       `-- regulon
|   |   |           `-- pickle
|   |   `-- figs
|   |       `-- Seurat
|   |            `-- regulon
|   `-- <cell_type_2>
|       |-- data
|       |   `-- Seurat
|       |       `-- regulon
|       |           `-- pickle
|       `-- figs
|           `-- Seurat
|              `-- regulon
```

For more information about arbitrary input parameters and output format please look in the corresponding scripts.