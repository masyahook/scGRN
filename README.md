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

## Contact

For further information please contact by mail at `mkriukov.job@gmail.com`.