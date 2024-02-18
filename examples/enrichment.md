# Enrichment analysis

**Enrichment analysis** is a statistical method used to identify the biological processes that are over-represented in a set of genes or proteins. It is widely used to understand which underlying biological processes are associated with a given gene set. In our context, we can perform enrichment analysis on **differentially expressed genes** or **community genes**. The former example can help us identify the cell composition in obtained scRNA-seq sample, which biological functions are enriched within each cell cluster. The latter example could hint on the sub-cellular organization of molecular processes.

We use [`EnrichR`](https://maayanlab.cloud/Enrichr/) tool to perform enrichment analysis. It is a comprehensive gene set enrichment analysis web server that includes various gene set libraries, such as Gene Ontology, KEGG, Reactome, and many others. Moreover, it provides the request API to extract the analysis results programmatically that we will use in our work.

## Running enrichment analysis with `run_enrichr.py`

We implemented [`run_enrichr.py`](../scGRN/network_analysis/run_enrichr.py) script to perform enrichment analysis on a given gene set. Let's have a look at the usage of the script:

```bash
> python run_enrichr.py --help
usage: run_enrichr.py [-h] -i IN_PATH -g GENE_COL [-o OUT_PATH] [-c GROUP_COL]
                      [-e ENRICHR_LIBRARY] [-q QUERY] [-n TOP_N]

Running EnrichR.

optional arguments:
  -h, --help            show this help message and exit
  -i IN_PATH, --in_path IN_PATH
                        The path to input data where gene sets are stored
  -g GENE_COL, --gene_col GENE_COL
                        The column name that stores gene names
  -o OUT_PATH, --out_path OUT_PATH
                        The path to output data to store, if not specified the
                        result will be saved in the same folder. Should be a
                        .tsv file
  -c GROUP_COL, --group_col GROUP_COL
                        The column names that store group names
  -e ENRICHR_LIBRARY, --enrichr_library ENRICHR_LIBRARY
                        The EnrichR library to use for enrichment analysis.
  -q QUERY, --query QUERY
                        The The query that can be used to select a subset of
                        original dataset.
  -n TOP_N, --top_n TOP_N
                        Select top_n genes for enrichment analysis, applies
                        for community datasets
```

We can see that script takes the input file with gene sets that are stored as dataframe, the column name of the gene identifier, and the output filename for results. The examples of input file are the following:

- **DE genes**: the output of [`FindAllMarkers`](https://satijalab.org/seurat/reference/findallmarkers) function that we also use internally for [`sc_pipeline`](../scGRN/single_cell_processing/sc_pipeline) and [`regulon_pipeline`](../scGRN/single_cell_processing/regulon_pipeline/) pipelines. The dataframe should contain the following columns: `gene` and `cluster`, and should look something like this:

    ```bash
    > head all_markers.tsv
    "p_val" "avg_log2FC"    "pct.1" "pct.2" "p_val_adj"     "cluster"       "gene"
    "IRF2"  0       8.34662342848101        0.05    0.005   0       "M"     "IRF2"
    "STAT1" 0       5.00663275002127        0.002   0       0       "M"     "STAT1"
    "NR3C1" 0       4.83123705748613        0.092   0.002   0       "M"     "NR3C1"
    "IRF1"  0       4.67667026971982        0.002   0       0       "M"     "IRF1"
    "REST"  0       4.02468677045236        0.215   0.016   0       "M"     "REST"
    "TBX21" 0       3.96278780361005        0.314   0.057   0       "M"     "TBX21"
    "FLI1"  0       3.95119356987923        0.18    0.021   0       "M"     "FLI1"
    "ELF2"  0       3.39472116379728        0.501   0.055   0       "M"     "ELF2"
    "CREB1" 0       2.88222505481395        0.27    0.046   0       "M"     "CREB1"
    ```

    An example of use:

    ```bash
    > python run_enrichr.py -i all_markers.tsv -g gene -c cluster -q "'p_val_adj' < 0.05 & 'avg_log2FC' > 1'"
    ```

- **Community genes**: the output of [`community_ana.py`](../scGRN/network_analysis/community_ana.py) script that we use internally for community detection. The dataframe should contain the the column `all_sorted_genes`, each row of the dataframe represent one cluster of genes. For more details, please look in the [`community.md`](community.md#loading-the-results).

    An example of use:

    ```bash
    > python run_enrichr.py -i raw_data_communities_info.pickle -g all_sorted_genes
    ```

- **General purpose**: the user can also provide any dataframe with `gene_col` column storing gene names and `group_col` column storing group names. 

Look also for more info about `run_enrichr` input parameters [here](https://github.com/masyahook/scGRN/blob/99f1ba91303351cd9948016dfaea7ec78f35c30c/scGRN/network_analysis/_enrichment.py#L15).

## Visual advanced examples

In the [`GRNBoost2_Enrrichment_analysis.ipynb`](../notebooks/GRNBoost2_Enrichment_analysis.ipynb) notebook we gave two detailed examples of running enrichment analysis on **differentially activated TFs** and **community genes in gene-gene regulatory networks**. With powerful [`clusterProfiler`](https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html) visualization capabilities we can easily depict how the biological processes are enriched in the gene sets. Please take a look!
