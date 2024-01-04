# Reading data

The `scGRN` package provides a set of tools to read and process inferred gene regulatory networks from scRNA-sequencing data. The analysis is based on the [Liao *et al.*, 2020](https://www.nature.com/articles/s41591-020-0901-9) dataset, which contains scRNA-seq data of COVID-19 patients. The uesr can see the sample metadata using:

```python

import scGRN

# Setting file system, also defined in config.py
_PROJ_HOME = (
    "/gpfs/projects/bsc08/shared_projects/scGRN_analysis"
)
_FMETA = (
    f"{_PROJ_HOME}/Data_home/data/GSE145926_RAW/metadata.tsv"
)
_DATA_HOME = f"{_PROJ_HOME}/Data_home/res/covid_19"

# Loading full metadata for all patients
full_meta = scGRN.ana.get_meta(_DATA_HOME, _FMETA)

print(full_meta.shape)
print(full_meta.columns)
```

::: {.table-responsive}

| id   | group   | file                                                                                                                           |   num_cells |   Macrophage |   T_cells |   DC |   Pre-B_cell_CD34- |   Monocyte |   NK_cell |   B_cell |   Epithelial_cells |   BM |   Pro-B_cell_CD34+ |   HSC_-G-CSF |   CMP |   Neutrophils |   GMP |   Erythroblast |   Gametocytes |   Neurons |   Fibroblasts |   Smooth_muscle_cells |   Hepatocytes |   Keratinocytes |   Pro-Myelocyte |
|:-----|:--------|:-------------------------------------------------------------------------------------------------------------------------------|------------:|-------------:|----------:|-----:|-------------------:|-----------:|----------:|---------:|-------------------:|-----:|-------------------:|-------------:|------:|--------------:|------:|---------------:|--------------:|----------:|--------------:|----------------------:|--------------:|----------------:|----------------:|
| C51  | C       | /gpfs/projects/bsc08/shared_projects/scGRN_analysis/Data_home/data/GSE145926_RAW/GSM4475048_C51_filtered_feature_bc_matrix.h5  |        9431 |         8348 |       608 |  215 |                 98 |         70 |        68 |        9 |                  7 |    4 |                  3 |            1 |   nan |           nan |   nan |            nan |           nan |       nan |           nan |                   nan |           nan |             nan |             nan |
| C52  | C       | /gpfs/projects/bsc08/shared_projects/scGRN_analysis/Data_home/data/GSE145926_RAW/GSM4475049_C52_filtered_feature_bc_matrix.h5  |        8696 |         8611 |        13 |   23 |                  3 |         14 |         5 |        2 |                 25 |  nan |                nan |          nan |   nan |           nan |   nan |            nan |           nan |       nan |           nan |                   nan |           nan |             nan |             nan |
| C100 | C       | /gpfs/projects/bsc08/shared_projects/scGRN_analysis/Data_home/data/GSE145926_RAW/GSM4475050_C100_filtered_feature_bc_matrix.h5 |         907 |          338 |       411 |   45 |                  5 |         51 |        20 |       12 |                 18 |  nan |                  2 |          nan |     5 |           nan |   nan |            nan |           nan |       nan |           nan |                   nan |           nan |             nan |             nan |
| C141 | M       | /gpfs/projects/bsc08/shared_projects/scGRN_analysis/Data_home/data/GSE145926_RAW/GSM4339769_C141_filtered_feature_bc_matrix.h5 |        1449 |          197 |       932 |   48 |                  5 |         86 |        96 |       33 |                 37 |  nan |                nan |          nan |     3 |            11 |     1 |            nan |           nan |       nan |           nan |                   nan |           nan |             nan |             nan |
| C142 | M       | /gpfs/projects/bsc08/shared_projects/scGRN_analysis/Data_home/data/GSE145926_RAW/GSM4339770_C142_filtered_feature_bc_matrix.h5 |        1790 |          482 |       996 |   39 |                 13 |         67 |       113 |       20 |                 38 |    1 |                  1 |            2 |     3 |            14 |     1 |            nan |           nan |       nan |           nan |                   nan |           nan |             nan |             nan |
| C144 | M       | /gpfs/projects/bsc08/shared_projects/scGRN_analysis/Data_home/data/GSE145926_RAW/GSM4339772_C144_filtered_feature_bc_matrix.h5 |         452 |           37 |       181 |   41 |                  8 |         73 |        34 |       14 |                 54 |    2 |                  1 |            1 |     1 |             3 |   nan |              2 |           nan |       nan |           nan |                   nan |           nan |             nan |             nan |
| C143 | S       | /gpfs/projects/bsc08/shared_projects/scGRN_analysis/Data_home/data/GSE145926_RAW/GSM4339771_C143_filtered_feature_bc_matrix.h5 |       14933 |         2048 |      1394 |  154 |                 33 |       7489 |       562 |       72 |                145 |  nan |                  1 |           24 |     1 |          3005 |     2 |            nan |             2 |         1 |           nan |                   nan |           nan |             nan |             nan |
| C145 | S       | /gpfs/projects/bsc08/shared_projects/scGRN_analysis/Data_home/data/GSE145926_RAW/GSM4339773_C145_filtered_feature_bc_matrix.h5 |       15550 |         6960 |       719 |  859 |                 46 |       5616 |       421 |       58 |                207 |  nan |                  1 |           26 |   nan |           635 |     2 |            nan |           nan |       nan |           nan |                   nan |           nan |             nan |             nan |
| C146 | S       | /gpfs/projects/bsc08/shared_projects/scGRN_analysis/Data_home/data/GSE145926_RAW/GSM4339774_C146_filtered_feature_bc_matrix.h5 |        2545 |          247 |        61 |   36 |                nan |        127 |        14 |        3 |                417 |  nan |                nan |            2 |   nan |          1632 |   nan |            nan |             1 |       nan |             2 |                     1 |             1 |               1 |             nan |
| C148 | S       | /gpfs/projects/bsc08/shared_projects/scGRN_analysis/Data_home/data/GSE145926_RAW/GSM4475051_C148_filtered_feature_bc_matrix.h5 |        1165 |           98 |       122 |   24 |                nan |        641 |        36 |        8 |                 52 |  nan |                  1 |          nan |     3 |           178 |     1 |            nan |             1 |       nan |           nan |                   nan |           nan |             nan |             nan |
| C149 | S       | /gpfs/projects/bsc08/shared_projects/scGRN_analysis/Data_home/data/GSE145926_RAW/GSM4475052_C149_filtered_feature_bc_matrix.h5 |        1936 |          176 |       681 |   80 |                  1 |        691 |        59 |       38 |                 41 |  nan |                nan |            5 |   nan |           164 |   nan |            nan |           nan |       nan |           nan |                   nan |           nan |             nan |             nan |
| C152 | S       | /gpfs/projects/bsc08/shared_projects/scGRN_analysis/Data_home/data/GSE145926_RAW/GSM4475053_C152_filtered_feature_bc_matrix.h5 |        2557 |          466 |       397 |   41 |                176 |        795 |        74 |      317 |                201 |    6 |                 40 |          nan |   nan |            30 |    13 |            nan |           nan |       nan |           nan |                   nan |           nan |             nan |               1 |

:::

The user should run the [`network_inference`](../scGRN/network_inference/) pipeline to obtain the GRN graphs stored in `pickle` format.
