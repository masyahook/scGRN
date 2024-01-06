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

Please note, that the user should run the [`network_inference`](../scGRN/network_inference/) pipeline before using the functions listed above. All the successful runs will be visualized in the availability tables.

## Load the data into the memory

To load the single cell RNA-sequencing data matrix into the memory, one can use [`scGRN.ana.get_sc_data()`](https://github.com/masyahook/scGRN/blob/d7eb61b97fdb0a87c74221451c8805d3c46d3280/scGRN/network_analysis/_data_processing.py#L393):

```python
# Loading the single-cell data for patient C51, T cells only
sc_data = scGRN.ana.get_sc_data(
    pat='C51',  # patient ID, look in docstring
    cell_type='T_cells',  # cell type, look in docstring
    # data_home=_DATA_HOME  # optional
)

sc_data
```

```python
>
            AAACCTGCACGGCTAC-1  ...  TTTGCGCCACAGGTTT-1
AL627309.1                   0  ...                   0
AL669831.5                   0  ...                   0
FAM87B                       0  ...                   0
LINC00115                    0  ...                   0
FAM41C                       0  ...                   0
...                        ...  ...                 ...
AC011043.1                   0  ...                   0
AL592183.1                   0  ...                   0
AC007325.4                   0  ...                   0
AL354822.1                   0  ...                   0
AC240274.1                   0  ...                   0

[16848 rows x 608 columns]
```

To load the adjacency list into the memory, one can use [`scGRN.ana.get_adj_list()`](https://github.com/masyahook/scGRN/blob/d7eb61b97fdb0a87c74221451c8805d3c46d3280/scGRN/network_analysis/_data_processing.py#L616):

```python
# Loading the adjacency list for healthy ('C') patients, B cells only, enriched TF-target GRN
adj_list = scGRN.ana.get_adj_list(
    pat='C',  # patient ID
    cell_type='B_cell',  # cell type
    net_type='ctx',  # could be also 'all' or 'TF'
    # data_home=_DATA_HOME  # optional
)

adj_list
```

```python
>
           TF    target  importance  regulation       rho
0      ARID3A    CYFIP2    1.901923         NaN       NaN
1      ARID3A  ZEB2-AS1    1.565187         1.0  1.000000
2      ARID3A     CXCR4    1.664066         1.0  0.852121
3      ARID3A    ATP2A2    2.819717         1.0  1.000000
4      ARID3A    ARID3A    1.000000         NaN       NaN
...       ...       ...         ...         ...       ...
15854   ZNF91    GPR160    2.403347         1.0  0.690849
15855   ZNF91     DMTF1    0.800444         1.0  0.550482
15856   ZNF91     DUOX1    4.431942         1.0  1.000000
15857   ZNF91    UBE2D3    0.320025         1.0  0.433138
15858   ZNF91      BPTF    0.753351         1.0  0.329073

[15859 rows x 5 columns]
```

To load the NetworkX graph into the memory, one can use [`scGRN.ana.get_nx_graph()`](https://github.com/masyahook/scGRN/blob/d7eb61b97fdb0a87c74221451c8805d3c46d3280/scGRN/network_analysis/_data_processing.py#L737):

```python
# Loading the NetworkX graph for all NK cell patient-cells, gene-gene network
nx_graph = scGRN.ana.get_nx_graph(
    pat='all',  # patient ID
    cell_type='NK_cell',  # cell type
    net_type='all',  # could be also 'all' or 'TF'
    # data_home=_DATA_HOME  # optional
)

print('Number of nodes:', nx_graph.number_of_nodes())
print('Number of edges:', nx_graph.number_of_edges())
nx_graph
```

```python
> Number of nodes: 16441
> Number of edges: 3344431
> <networkx.classes.digraph.DiGraph at 0x7f0a982390d0>
```

To load the community analysis information about corresponding GRNs, please use [`scGRN.ana.get_community_info()`](https://github.com/masyahook/scGRN/blob/d7eb61b97fdb0a87c74221451c8805d3c46d3280/scGRN/network_analysis/_data_processing.py#L860):

```python
# Loading the community analysis information for all DC cells
comm_info = scGRN.ana.get_community_info(
    pat='all',  # patient ID
    cell_type='DC',  # cell type
    # data_home=_DATA_HOME  # optional
)

comm_info
```

```python
>
   num_nodes num_edges  ...                 top_links_scores_with_community_17                 top_links_scores_with_community_18
0       2584     13400  ...  MT-CO3 <-> MT-ND1 (score=42.92); MT-ND4 <-> MT...  MT-CO3 <-> MT-ND1 (score=42.92); MT-ND4 <-> MT...
1       1651     10203  ...  CCL8 <-> CCL2 (score=93.95); PSAP <-> LAPTM5 (...  CCL8 <-> CCL2 (score=93.95); PSAP <-> LAPTM5 (...
2       1479     10956  ...  BIRC3 <-> LAMP3 (score=40.53); FSCN1 <-> MARCK...  BIRC3 <-> LAMP3 (score=40.53); FSCN1 <-> MARCK...
3       1229      3445  ...  PDLIM2 <-> TBC1D32 (score=17.13); LINC02345 <-...  PDLIM2 <-> TBC1D32 (score=17.13); LINC02345 <-...
4       1205      8261  ...  DNAH7 <-> EFCAB6 (score=34.36); TUBB4B <-> SMI...  DNAH7 <-> EFCAB6 (score=34.36); TUBB4B <-> SMI...
5       1137     14024  ...  B2M <-> CYBA (score=35.29); HLA-C <-> HLA-B (s...  B2M <-> CYBA (score=35.29); HLA-C <-> HLA-B (s...
6        985      6473  ...  CXCL11 <-> CXCL10 (score=28.07); CXCL10 <-> IF...  CXCL11 <-> CXCL10 (score=28.07); CXCL10 <-> IF...
7        906      2580  ...  MROH1 <-> RECK (score=12.17); VIM-AS1 <-> PPIP...  MROH1 <-> RECK (score=12.17); VIM-AS1 <-> PPIP...
8        885      4447  ...  CCL4 <-> CCL3 (score=69.11); CCL4L2 <-> CCL3L1...  CCL4 <-> CCL3 (score=69.11); CCL4L2 <-> CCL3L1...
9        857     10181  ...  RPS27 <-> RPS29 (score=41.84); CD74 <-> HLA-DR...  RPS27 <-> RPS29 (score=41.84); CD74 <-> HLA-DR...
10       773      2596  ...  LPXN <-> BCL3 (score=14.88); H1F0 <-> EMP2 (sc...  LPXN <-> BCL3 (score=14.88); H1F0 <-> EMP2 (sc...
11       699      2556  ...  IGHV4-34 <-> IGLV3-19 (score=113.10); IGLV3-19...  IGHV4-34 <-> IGLV3-19 (score=113.10); IGLV3-19...
12       500      1637  ...  GATA2 <-> CLU (score=9.55); OFD1 <-> NIPSNAP1 ...  GATA2 <-> CLU (score=9.55); OFD1 <-> NIPSNAP1 ...
13       410      2362  ...  BPIFA1 <-> MSMB (score=13.02); MUC5AC <-> MSMB...  BPIFA1 <-> MSMB (score=13.02); MUC5AC <-> MSMB...
14       397      1240  ...  DCAF4L1 <-> TRIM25 (score=11.89); SND1-IT1 <->...  DCAF4L1 <-> TRIM25 (score=11.89); SND1-IT1 <->...
15       392      1090  ...  FASTKD2 <-> REXO4 (score=8.34); AC022762.2 <->...  FASTKD2 <-> REXO4 (score=8.34); AC022762.2 <->...
16       263       784  ...  TDRD12 <-> FRRS1 (score=6.00); CCL18 <-> MT1X ...  TDRD12 <-> FRRS1 (score=6.00); CCL18 <-> MT1X ...
17       132       439  ...                                                NaN  AP001363.2 <-> NOSTRIN (score=6.88); FOXK2 <->...
18        22        69  ...  RCCD1 <-> EFCAB12 (score=3.47); AC245052.7 <->...                                                NaN

[19 rows x 63 columns]
```

```python
print(comm_info.columns)
```

```python
Index(['num_nodes', 'num_edges', 'main_functions_GO', 'main_functions_KEGG',
       'main_functions_immunological', 'main_functions_hallmark',
       'non_lambert_2018_TF_central_genes', 'non_dorothea_TF_central_genes',
       'new_gene_gene_links_KEGG', 'new_gene_gene_links_hallmark',
       'whole_G_central_genes_scores', 'other_functions_GO',
       'other_functions_KEGG', 'other_functions_immunological',
       'other_functions_hallmark', 'sorted_central_genes_scores',
       'sorted_central_functions_GO', 'sorted_central_functions_KEGG',
       'sorted_central_functions_immunological',
       'sorted_central_functions_hallmark', 'most_frequent_function_words_GO',
       'most_frequent_function_words_KEGG',
       'most_frequent_function_words_immunological',
       'most_frequent_function_words_hallmark', 'all_sorted_genes',
       'top_links_scores_central_genes<->community_0',
       'top_links_scores_central_genes<->community_1',
       'top_links_scores_central_genes<->community_2',
       'top_links_scores_central_genes<->community_3',
       'top_links_scores_central_genes<->community_4',
       'top_links_scores_central_genes<->community_5',
       'top_links_scores_central_genes<->community_6',
       'top_links_scores_central_genes<->community_7',
       'top_links_scores_central_genes<->community_8',
       'top_links_scores_central_genes<->community_9',
       'top_links_scores_central_genes<->community_10',
       'top_links_scores_central_genes<->community_11',
       'top_links_scores_central_genes<->community_12',
       'top_links_scores_central_genes<->community_13',
       'top_links_scores_central_genes<->community_14',
       'top_links_scores_central_genes<->community_15',
       'top_links_scores_central_genes<->community_16',
       'top_links_scores_central_genes<->community_17',
       'top_links_scores_central_genes<->community_18',
       'top_links_scores_with_community_0',
       'top_links_scores_with_community_1',
       'top_links_scores_with_community_2',
       'top_links_scores_with_community_3',
       'top_links_scores_with_community_4',
       'top_links_scores_with_community_5',
       'top_links_scores_with_community_6',
       'top_links_scores_with_community_7',
       'top_links_scores_with_community_8',
       'top_links_scores_with_community_9',
       'top_links_scores_with_community_10',
       'top_links_scores_with_community_11',
       'top_links_scores_with_community_12',
       'top_links_scores_with_community_13',
       'top_links_scores_with_community_14',
       'top_links_scores_with_community_15',
       'top_links_scores_with_community_16',
       'top_links_scores_with_community_17',
       'top_links_scores_with_community_18'],
      dtype='object')
```

For all data-loading functions there is a parameter `tolerate_missing` that handles the missing data. For example:

```python
# Loading the adjacency list for C146 patient, Pre-B_cell_CD34- only, gene-genet network
adj_list = scGRN.ana.get_adj_list(
    pat='C146',  # patient ID
    cell_type='Pre-B_cell_CD34-',  # cell type
    net_type='all',  # could be also 'all' or 'TF'
    tolerate_missing=True,  # if True, the function will return the empty dataframe if the data is missing
    # data_home=_DATA_HOME  # optional
)
```

which raises warning:

```python
/gpfs/home/bsc08/bsc08890/scGRN_analysis/scGRN/network_analysis/_data_processing.py:727: UserWarning: The GRN for pat="C146", cell_type="Pre-B_cell_CD34-", net_type="all" is not found (should be at "/gpfs/projects/bsc08/shared_projects/scGRN_analysis/Data_home/res/covid_19/C146/data/grnboost2/pickle/raw_data_Pre-B_cell_CD34-_cor.pickle"). Returning as `None` instead!
  f'The GRN for pat="{pat}", cell_type="{cell_type}", net_type="{net_type}" is not found '
```
