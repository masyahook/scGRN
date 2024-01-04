# Reading data for analysis

The `scGRN` package provides a set of tools to read and process inferred gene regulatory networks from scRNA-sequencing data. The analysis is based on the [Liao *et al.*, 2020](https://www.nature.com/articles/s41591-020-0901-9) dataset, which contains scRNA-seq data of COVID-19 patients. Thus the API is developed to simplify loading the data from original and derived files.

## Reading metadata

The user can see the sample metadata using [`scGRN.ana.get_meta()`](https://github.com/masyahook/scGRN/blob/04024eb246af2277dee9253b772f0b2141c3fe87/scGRN/network_analysis/_data_processing.py#L18)

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
full_meta = scGRN.ana.get_meta(
    data_home=_DATA_HOME,  # file path to the data home folder
    meta_file=_FMETA  # file path to the patient metadata file
)

# You can also just run without params, as constants are defined in config.py
# full_meta = scGRN.ana.get_meta()

print(full_meta)
```

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

## Getting the list of available data
To get the list of available single-cell datasets for different patients (e.g. C51), patient group (e.g. C / M / S) and cell types (e.g. T cells) use [`scGRN.ana.get_avail_sc_data()`](https://github.com/masyahook/scGRN/blob/04024eb246af2277dee9253b772f0b2141c3fe87/scGRN/network_analysis/_data_processing.py#L46C5-L46C22):

```python

# Getting the list of available cell types for each patient data
avail_sc_data = scGRN.ana.get_avail_sc_data(
    data_home=_DATA_HOME,  # file path to the data home folder
    meta_file=_FMETA  # file path to the patient metadata file
)

# You can also just run without params, as constants are defined in config.py
# avail_sc_data = scGRN.ana.get_avail_sc_data()

# 1 = data sample exists
# 0 = data sample does not exist
# nan = data sample is missing probably due to computation error
print(avail_sc_data)
```

|          |   all_data |   Macrophage |   T_cells |   DC |   Pre-B_cell_CD34- |   Monocyte |   NK_cell |   B_cell |   Epithelial_cells |   BM |   Pro-B_cell_CD34+ |   HSC_-G-CSF |   CMP |   Neutrophils |   GMP |   Erythroblast |   Gametocytes |   Neurons |   Fibroblasts |   Smooth_muscle_cells |   Hepatocytes |   Keratinocytes |   Pro-Myelocyte |
|:---------|-----------:|-------------:|----------:|-----:|-------------------:|-----------:|----------:|---------:|-------------------:|-----:|-------------------:|-------------:|------:|--------------:|------:|---------------:|--------------:|----------:|--------------:|----------------------:|--------------:|----------------:|----------------:|
| all_data |          1 |            1 |         1 |    1 |                  1 |          1 |         1 |        1 |                  1 |    0 |                  0 |            0 |     0 |             1 |     0 |              0 |             0 |         0 |             0 |                     0 |             0 |               0 |               0 |
| C        |          1 |            1 |         1 |    1 |                  1 |          1 |         1 |        1 |                  1 |    0 |                  0 |            0 |     0 |           nan |   nan |            nan |           nan |       nan |           nan |                   nan |           nan |             nan |             nan |
| M        |          1 |            1 |         1 |    1 |                  1 |          1 |         1 |        1 |                  1 |    0 |                  0 |            0 |     0 |             1 |     0 |              0 |           nan |       nan |           nan |                   nan |           nan |             nan |             nan |
| S        |          1 |            1 |         1 |    1 |                  1 |          1 |         1 |        1 |                  1 |    0 |                  0 |            0 |     0 |             1 |     0 |            nan |             0 |         0 |             0 |                     0 |             0 |               0 |               0 |
| C51      |          1 |            1 |         1 |    1 |                  1 |          1 |         1 |        1 |                  1 |    1 |                  1 |            1 |   nan |           nan |   nan |            nan |           nan |       nan |           nan |                   nan |           nan |             nan |             nan |
| C52      |          1 |            1 |         1 |    1 |                  1 |          1 |         1 |        1 |                  1 |  nan |                nan |          nan |   nan |           nan |   nan |            nan |           nan |       nan |           nan |                   nan |           nan |             nan |             nan |
| C100     |          1 |            1 |         1 |    1 |                  1 |          1 |         1 |        1 |                  1 |  nan |                  1 |          nan |     1 |           nan |   nan |            nan |           nan |       nan |           nan |                   nan |           nan |             nan |             nan |
| C141     |          1 |            1 |         1 |    1 |                  1 |          1 |         1 |        1 |                  1 |  nan |                nan |          nan |     1 |             1 |     1 |            nan |           nan |       nan |           nan |                   nan |           nan |             nan |             nan |
| C142     |          1 |            1 |         1 |    1 |                  1 |          1 |         1 |        1 |                  1 |    1 |                  1 |            1 |     1 |             1 |     1 |            nan |           nan |       nan |           nan |                   nan |           nan |             nan |             nan |
| C144     |          1 |            1 |         1 |    1 |                  1 |          1 |         1 |        1 |                  1 |    1 |                  1 |            1 |     1 |             1 |   nan |              1 |           nan |       nan |           nan |                   nan |           nan |             nan |             nan |
| C143     |          1 |            1 |         1 |    1 |                  1 |          1 |         1 |        1 |                  1 |  nan |                  1 |            1 |     1 |             1 |     1 |            nan |             1 |         1 |           nan |                   nan |           nan |             nan |             nan |
| C145     |          1 |            1 |         1 |    1 |                  1 |          1 |         1 |        1 |                  1 |  nan |                  1 |            1 |   nan |             1 |     1 |            nan |           nan |       nan |           nan |                   nan |           nan |             nan |             nan |
| C146     |          1 |            1 |         1 |    1 |                nan |          1 |         1 |        1 |                  1 |  nan |                nan |            1 |   nan |             1 |   nan |            nan |             1 |       nan |             1 |                     1 |             1 |               1 |             nan |
| C148     |          1 |            1 |         1 |    1 |                nan |          1 |         1 |        1 |                  1 |  nan |                  1 |          nan |     1 |             1 |     1 |            nan |             1 |       nan |           nan |                   nan |           nan |             nan |             nan |
| C149     |          1 |            1 |         1 |    1 |                  1 |          1 |         1 |        1 |                  1 |  nan |                nan |            1 |   nan |             1 |   nan |            nan |           nan |       nan |           nan |                   nan |           nan |             nan |             nan |
| C152     |          1 |            1 |         1 |    1 |                  1 |          1 |         1 |        1 |                  1 |    1 |                  1 |          nan |   nan |             1 |     1 |            nan |           nan |       nan |           nan |                   nan |           nan |             nan |               1 |

To get the list of available adjacency lists (i.e. edge lists describing GRN graphs), one can use [`scGRN.ana.get_avail_adj_lists()`](https://github.com/masyahook/scGRN/blob/04024eb246af2277dee9253b772f0b2141c3fe87/scGRN/network_analysis/_data_processing.py#L131)

```python
# Getting the list of available adjacency lists for each patient data
avail_adj_lists = scGRN.ana.get_avail_adj_lists(
    data_home=_DATA_HOME,  # file path to the data home folder
    meta_file=_FMETA,  # file path to the patient metadata file
    filtered=0.95  # quantile threshold for the edge weights, default is 0.95
)

# You can also just run without params, as constants are defined in config.py
# avail_adj_lists = scGRN.ana.get_avail_adj_lists()

# 'all' - gene-gene GRN
# 'TF' - TF-target GRN
# 'ctx' - enriched TF-target GRN
print(avail_adj_lists.keys())
```
```python
> dict_keys(['all', 'TF', 'ctx'])
```
```python
# 1 = data sample exists
# 0 = data sample does not exist
# nan = data sample is missing probably due to computation error
print(avail_adj_lists['all'])
```

|          |   all_data |   Macrophage |   T_cells |   DC |   Pre-B_cell_CD34- |   Monocyte |   NK_cell |   B_cell |   Epithelial_cells |   BM |   Pro-B_cell_CD34+ |   HSC_-G-CSF |   CMP |   Neutrophils |   GMP |   Erythroblast |   Gametocytes |   Neurons |   Fibroblasts |   Smooth_muscle_cells |   Hepatocytes |   Keratinocytes |   Pro-Myelocyte |
|:---------|-----------:|-------------:|----------:|-----:|-------------------:|-----------:|----------:|---------:|-------------------:|-----:|-------------------:|-------------:|------:|--------------:|------:|---------------:|--------------:|----------:|--------------:|----------------------:|--------------:|----------------:|----------------:|
| all_data |          0 |            1 |         1 |    1 |                  1 |          1 |         1 |        1 |                  1 |    0 |                  0 |            0 |     0 |             1 |     0 |              0 |             0 |         0 |             0 |                     0 |             0 |               0 |               0 |
| C        |          1 |            1 |         1 |    1 |                  1 |          1 |         1 |        1 |                  1 |    0 |                  0 |            0 |     0 |           nan |   nan |            nan |           nan |       nan |           nan |                   nan |           nan |             nan |             nan |
| M        |          1 |            1 |         1 |    1 |                  1 |          1 |         1 |        1 |                  1 |    0 |                  0 |            0 |     0 |             1 |     0 |              0 |           nan |       nan |           nan |                   nan |           nan |             nan |             nan |
| S        |          0 |            1 |         1 |    1 |                  1 |          1 |         1 |        1 |                  1 |    0 |                  0 |            0 |     0 |             1 |     0 |            nan |             0 |         0 |             0 |                     0 |             0 |               0 |               0 |
| C51      |          1 |            1 |         1 |    1 |                  1 |          1 |         1 |        1 |                  1 |    1 |                  1 |            1 |   nan |           nan |   nan |            nan |           nan |       nan |           nan |                   nan |           nan |             nan |             nan |
| C52      |          1 |            1 |         1 |    1 |                  1 |          1 |         1 |        1 |                  1 |  nan |                nan |          nan |   nan |           nan |   nan |            nan |           nan |       nan |           nan |                   nan |           nan |             nan |             nan |
| C100     |          1 |            1 |         1 |    1 |                  1 |          1 |         1 |        1 |                  1 |  nan |                  1 |          nan |     1 |           nan |   nan |            nan |           nan |       nan |           nan |                   nan |           nan |             nan |             nan |
| C141     |          1 |            1 |         1 |    1 |                  1 |          1 |         1 |        1 |                  1 |  nan |                nan |          nan |     1 |             1 |     1 |            nan |           nan |       nan |           nan |                   nan |           nan |             nan |             nan |
| C142     |          1 |            1 |         1 |    1 |                  1 |          1 |         1 |        1 |                  1 |    1 |                  1 |            1 |     1 |             1 |     1 |            nan |           nan |       nan |           nan |                   nan |           nan |             nan |             nan |
| C144     |          1 |            1 |         1 |    1 |                  1 |          1 |         1 |        1 |                  1 |    1 |                  1 |            1 |     1 |             1 |   nan |              1 |           nan |       nan |           nan |                   nan |           nan |             nan |             nan |
| C143     |          1 |            1 |         1 |    1 |                  1 |          1 |         1 |        1 |                  1 |  nan |                  1 |            1 |     1 |             1 |     1 |            nan |             1 |         1 |           nan |                   nan |           nan |             nan |             nan |
| C145     |          1 |            1 |         1 |    1 |                  1 |          1 |         1 |        1 |                  1 |  nan |                  1 |            1 |   nan |             1 |     1 |            nan |           nan |       nan |           nan |                   nan |           nan |             nan |             nan |
| C146     |          1 |            1 |         1 |    1 |                nan |          1 |         1 |        1 |                  1 |  nan |                nan |            1 |   nan |             1 |   nan |            nan |             1 |       nan |             1 |                     1 |             1 |               1 |             nan |
| C148     |          1 |            1 |         1 |    1 |                nan |          1 |         1 |        1 |                  1 |  nan |                  1 |          nan |     1 |             1 |     1 |            nan |             1 |       nan |           nan |                   nan |           nan |             nan |             nan |
| C149     |          1 |            1 |         1 |    1 |                  1 |          1 |         1 |        1 |                  1 |  nan |                nan |            1 |   nan |             1 |   nan |            nan |           nan |       nan |           nan |                   nan |           nan |             nan |             nan |
| C152     |          1 |            1 |         1 |    1 |                  1 |          1 |         1 |        1 |                  1 |    1 |                  1 |          nan |   nan |             1 |     1 |            nan |           nan |       nan |           nan |                   nan |           nan |             nan |               1 |

In a similar manner we can get the list of available computed [`NetworkX`](https://networkx.org) graphs using [`scGRN.ana.get_avail_nx_graphs()`](https://github.com/masyahook/scGRN/blob/04024eb246af2277dee9253b772f0b2141c3fe87/scGRN/network_analysis/_data_processing.py#L251)

```python
# Getting the list of available NetworkX graphs for each patient data
avail_nx_graphs = scGRN.ana.get_avail_nx_graphs(
    data_home=_DATA_HOME,  # file path to the data home folder
    meta_file=_FMETA,  # file path to the patient metadata file
    filtered=0.95  # quantile threshold for the edge weights, default is 0.95
)

# You can also just run without params, as constants are defined in config.py
# avail_nx_graphs = scGRN.ana.get_avail_nx_graphs()

# 'all' - gene-gene GRN
# 'TF' - TF-target GRN
# 'ctx' - enriched TF-target GRN
print(avail_nx_graphs['ctx'])
```

|          |   all_data |   Macrophage |   T_cells |   DC |   Pre-B_cell_CD34- |   Monocyte |   NK_cell |   B_cell |   Epithelial_cells |   BM |   Pro-B_cell_CD34+ |   HSC_-G-CSF |   CMP |   Neutrophils |   GMP |   Erythroblast |   Gametocytes |   Neurons |   Fibroblasts |   Smooth_muscle_cells |   Hepatocytes |   Keratinocytes |   Pro-Myelocyte |
|:---------|-----------:|-------------:|----------:|-----:|-------------------:|-----------:|----------:|---------:|-------------------:|-----:|-------------------:|-------------:|------:|--------------:|------:|---------------:|--------------:|----------:|--------------:|----------------------:|--------------:|----------------:|----------------:|
| all_data |          1 |            1 |         1 |    1 |                  1 |          1 |         1 |        1 |                  1 |    0 |                  0 |            0 |     0 |             1 |     0 |              0 |             0 |         0 |             0 |                     0 |             0 |               0 |               0 |
| C        |          1 |            1 |         1 |    1 |                  1 |          1 |         1 |        1 |                  1 |    0 |                  0 |            0 |     0 |           nan |   nan |            nan |           nan |       nan |           nan |                   nan |           nan |             nan |             nan |
| M        |          1 |            1 |         1 |    1 |                  1 |          1 |         1 |        1 |                  1 |    0 |                  0 |            0 |     0 |             1 |     0 |              0 |           nan |       nan |           nan |                   nan |           nan |             nan |             nan |
| S        |          1 |            1 |         1 |    1 |                  1 |          1 |         1 |        1 |                  1 |    0 |                  0 |            0 |     0 |             1 |     0 |            nan |             0 |         0 |             0 |                     0 |             0 |               0 |               0 |
| C51      |          1 |            1 |         1 |    1 |                  1 |          1 |         1 |        1 |                  1 |    1 |                  1 |            0 |   nan |           nan |   nan |            nan |           nan |       nan |           nan |                   nan |           nan |             nan |             nan |
| C52      |          1 |            1 |         1 |    1 |                  1 |          1 |         1 |        0 |                  1 |  nan |                nan |          nan |   nan |           nan |   nan |            nan |           nan |       nan |           nan |                   nan |           nan |             nan |             nan |
| C100     |          1 |            1 |         1 |    1 |                  1 |          1 |         1 |        1 |                  1 |  nan |                  0 |          nan |     1 |           nan |   nan |            nan |           nan |       nan |           nan |                   nan |           nan |             nan |             nan |
| C141     |          1 |            1 |         1 |    1 |                  1 |          1 |         1 |        1 |                  1 |  nan |                nan |          nan |     1 |             1 |     0 |            nan |           nan |       nan |           nan |                   nan |           nan |             nan |             nan |
| C142     |          1 |            1 |         1 |    1 |                  1 |          1 |         1 |        1 |                  1 |    0 |                  0 |            0 |     1 |             1 |     0 |            nan |           nan |       nan |           nan |                   nan |           nan |             nan |             nan |
| C144     |          1 |            1 |         1 |    1 |                  1 |          1 |         1 |        1 |                  1 |    0 |                  0 |            0 |     0 |             1 |   nan |              0 |           nan |       nan |           nan |                   nan |           nan |             nan |             nan |
| C143     |          1 |            1 |         1 |    1 |                  1 |          1 |         1 |        1 |                  1 |  nan |                  0 |            1 |     0 |             1 |     0 |            nan |             0 |         0 |           nan |                   nan |           nan |             nan |             nan |
| C145     |          1 |            1 |         1 |    1 |                  1 |          1 |         1 |        1 |                  1 |  nan |                  0 |            1 |   nan |             1 |     0 |            nan |           nan |       nan |           nan |                   nan |           nan |             nan |             nan |
| C146     |          1 |            1 |         1 |    1 |                nan |          1 |         1 |        1 |                  1 |  nan |                nan |            0 |   nan |             1 |   nan |            nan |             0 |       nan |             0 |                     0 |             0 |               0 |             nan |
| C148     |          1 |            1 |         1 |    1 |                nan |          1 |         1 |        1 |                  1 |  nan |                  0 |          nan |     1 |             1 |     0 |            nan |             0 |       nan |           nan |                   nan |           nan |             nan |             nan |
| C149     |          1 |            1 |         1 |    1 |                  0 |          1 |         1 |        1 |                  1 |  nan |                nan |            1 |   nan |             1 |   nan |            nan |           nan |       nan |           nan |                   nan |           nan |             nan |             nan |
| C152     |          1 |            1 |         1 |    1 |                  1 |          1 |         1 |        1 |                  1 |    1 |                  1 |          nan |   nan |             1 |     1 |            nan |           nan |       nan |           nan |                   nan |           nan |             nan |               0 |

We also implemented two function to extract available patient-specific files, [`scGRN.ana.get_avail_pat_sc`](https://github.com/masyahook/scGRN/blob/04024eb246af2277dee9253b772f0b2141c3fe87/scGRN/network_analysis/_data_processing.py#L107) and [`scGRN.ana.get_avail_pat_nx`](https://github.com/masyahook/scGRN/blob/04024eb246af2277dee9253b772f0b2141c3fe87/scGRN/network_analysis/_data_processing.py#L352) for single-cell data and NetworkX graphs respectively.

```python
# Getting the list of available single-cell data files for patient C51
avail_samples = scGRN.ana.get_avail_pat_sc(pat='C51')
```

```python
['raw_data',
 'raw_data_Macrophage',
 'raw_data_T_cells',
 'raw_data_DC',
 'raw_data_Pre-B_cell_CD34-',
 'raw_data_Monocyte',
 'raw_data_NK_cell',
 'raw_data_B_cell',
 'raw_data_Epithelial_cells',
 'raw_data_BM',
 'raw_data_Pro-B_cell_CD34+',
 'raw_data_HSC_-G-CSF',
 'raw_data_CMP',
 'raw_data_Neutrophils',
 'raw_data_GMP',
 'raw_data_Erythroblast',
 'raw_data_Gametocytes',
 'raw_data_Neurons',
 'raw_data_Fibroblasts',
 'raw_data_Smooth_muscle_cells',
 'raw_data_Hepatocytes',
 'raw_data_Keratinocytes',
 'raw_data_Pro-Myelocyte']
```

```python
# Getting the list of available NetworkX graph files for patient C51
avail_samples = get_avail_pat_nx(
    net_type='ctx',  # could be also 'all' or 'TF'
    pat='C51', # patient ID
)

print(avail_samples)
```

```python
['raw_data',
 'raw_data_Macrophage',
 'raw_data_T_cells',
 'raw_data_DC',
 'raw_data_Pre-B_cell_CD34-',
 'raw_data_Monocyte',
 'raw_data_NK_cell',
 'raw_data_B_cell',
 'raw_data_Epithelial_cells',
 'raw_data_BM',
 'raw_data_Pro-B_cell_CD34+']
```

The user should run the [`network_inference`](../scGRN/network_inference/) pipeline to obtain the GRN graphs stored in `pickle` format.
