![scGRN_logo](https://raw.githubusercontent.com/masyahook/scGRN/main/docs/scGRN_logo.png)

# `scGRN` - Gene regulatory network inference and analysis based on scRNA-sequencing data

This repository contains the computational pipeline of the GRN-inference analysis based on expression of single cells. scRNA-sequencing datasets capture rich information about the transciptomic gene levels across thousands of cells and thus could be used to describe the gene interaction. All gene interactions could be depicted as **gene regulatory networks (GRNs)** that summarize genetic communication and regulation as a graph network. The focus of this repository is to provide an end-to-end pipeline to **infer and** <ins>**analyze**</ins> the GRNs that were computed based on scRNA-sequencing data.

As a case study, we focus on the COVID-19 patient dataset ([Liao *et al.*, 2020](https://www.nature.com/articles/s41591-020-0901-9)). The data is available on GEO ([GSE145926](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE145926)) and contains 10xGenomics pipeline samples from the lung immune microenvironment in the bronchoalveolar lavage fluid (BALF) from <span style="color:red">6 severe</span> and <span style="color:#8B8000">3 moderate</span> COVID-19 patients and <span style="color:green">3 healthy</span> controls (look below).

## Key features

- Single cell processing using [`Seurat`](https://satijalab.org/seurat/) and cell type identification using [`SingleR`](https://bioconductor.org/packages/release/bioc/html/SingleR.html)
- Gene regulatory network inference using [`arboreto`](https://arboreto.readthedocs.io) and [`pyscenic`](https://pyscenic.readthedocs.io) (particularly with `GRNBoost2` algorithm)
- Exploratory data analysis (EDA) of the inferred GRNs, description of nodes, edges and other GRN properties
- Rich data visualization of GRNs, on-the-fly comparison with known singaling networks from [NDEx](https://www.ndexbio.org)
- Community detection analysis of inferred GRN using [Louvain](https://python-louvain.readthedocs.io/en/latest/index.html) or [Leiden](https://leidenalg.readthedocs.io/en/stable/index.html) algorithms supported with [wordcloud](https://github.com/amueller/word_cloud) visualization
- Enrichment analysis of gene set clusters using [`EnrichR`](https://maayanlab.cloud/Enrichr/) and [`clusterProfiler`](https://yulab-smu.top/biomedical-knowledge-mining-book/enrichment-overview.html), generation of functional gene networks
- Identification of clinical relevance based on gene-gene interaction, linking tailored gene communcation structure to patient phenotype

## Types of inferred networks

We will work with two types of networks:

- **Gene-gene networks** - GRNs that are generated based on co-expression between **all** genes, i.e. we consider all possible pair-wise gene connections are possible
- **TF regulon networks** - GRNs that are generated based on co-expression and [motif enrichment](https://www.nature.com/articles/nmeth.4463#Abs2), i.e. consider only connections between transcription factors and corresponding targets a.k.a. regulons

![Main pipeline](https://raw.githubusercontent.com/masyahook/scGRN/main/docs/main_pipeline.png)

## Usage

### General structure of the package

The source code is available in the [`scGRN`](https://github.com/masyahook/scGRN/tree/main/scGRN) folder which consists of the following submodules:

- [`single_cell_processing`](scGRN/single_cell_processing/) - processing of single cell RNA-sequencing data
  - [`sc_pipeline`](scGRN/single_cell_processing/sc_pipeline/) - [`Seurat`](https://satijalab.org/seurat/) pipeline for single cell data processing, includes quality control, normalization, dimensionality reduction, clustering, cell type identification and data aggregation
  - [`regulon_pipeline`](scGRN/single_cell_processing/regulon_pipeline/) - regulon enrichment pipeline using [`VIPER`](https://viper.readthedocs.io/en/latest/). Regulon is a set of genes that are regulated by a common regulatory protein, e.g. transcription factor. [`VIPER`](https://viper.readthedocs.io/en/latest/) measures the activity of different regulons based on the co-expression of corresponding genes. Either [`pyscenic`](https://pyscenic.readthedocs.io) or [`DoRothEA`](https://saezlab.github.io/dorothea/) regulons can be used
- [`network_inference`](scGRN/network_inference/) - gene regulatory network inference pipeline using [`pySCENIC`](https://pyscenic.readthedocs.io/en/latest/), look in [Types of inferred networks](#types-of-inferred-networks) section for more details
- [`network_analysis`](scGRN/network_analysis/) - analysis of inferred networks, includes exploratory data analysis (EDA), visualization, community detection and enrichment analysis

### Script-based modules

The work is heavily based on the [Marenostrum](https://www.bsc.es/marenostrum/marenostrum) supercomputer, thus the code is optimized for working in the HPC environment (i.e. [Slurm](https://slurm.schedmd.com)). The [`single_cell_processing`](scGRN/single_cell_processing/) pipeline is implemented in R, while [`network_inference`](scGRN/network_inference/) and [`network_analysis`](scGRN/network_analysis/) are implemented in Python. Both pipelines are wrapped in bash scripts for easy execution. Some scripts are tailored to the Marenostrum file system (file names conventions, Slurm-based workload commands), but could be easily adapted to other HPC environments. The community analysis provided in [`network_analysis`](scGRN/network_analysis/) could be also run at scale using Slurm manager. You can see the examples of the `sbatch` and [`greasy`](https://www.bsc.es/marenostrum/marenostrum/available-software/greasy) commands in the [`notebooks/Generate_sbatch_commands.ipynb`](notebooks/Generate_sbatch_commands.ipynb) notebook. For more details of the pipeline usage, please look in the corresponding `README.md` files in the submodules.

### Jupyter-based modules

The `network_analysis` pipeline is implemented as a Python package that the user can use to analyze the inferred networks. It includes utilities for I/O operations, graph analysis, data aggregation and visualization. As an example, to read the patient sample metadata:

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

```python
(12, 4)
Index(['group', 'file', 'num_cells', 'Macrophage', 'T_cells', 'DC',
       'Pre-B_cell_CD34-', 'Monocyte', 'NK_cell', 'B_cell', 'Epithelial_cells',
       'BM', 'Pro-B_cell_CD34+', 'HSC_-G-CSF', 'CMP', 'Neutrophils', 'GMP',
       'Erythroblast', 'Gametocytes', 'Neurons', 'Fibroblasts',
       'Smooth_muscle_cells', 'Hepatocytes', 'Keratinocytes', 'Pro-Myelocyte'],
      dtype='object')
```

Let's display some patient samples and their cell type composition:

```python
display(full_meta[['group', 'num_cells', 'Macrophage', 'T_cells']])
```

Output:

| id   | group   |   num_cells |   Macrophage |   T_cells |
|:-----|:--------|------------:|-------------:|----------:|
| C51  | C       |        9431 |         8348 |       608 |
| C52  | C       |        8696 |         8611 |        13 |
| C100 | C       |         907 |          338 |       411 |
| C141 | M       |        1449 |          197 |       932 |
| C142 | M       |        1790 |          482 |       996 |
| C144 | M       |         452 |           37 |       181 |
| C143 | S       |       14933 |         2048 |      1394 |
| C145 | S       |       15550 |         6960 |       719 |
| C146 | S       |        2545 |          247 |        61 |
| C148 | S       |        1165 |           98 |       122 |
| C149 | S       |        1936 |          176 |       681 |
| C152 | S       |        2557 |          466 |       397 |

For more examples of use look into [`examples`](examples/) and [`notebooks`](notebooks/).

## Contact

For further information please contact by mail at `mkriukov.job@gmail.com`.
