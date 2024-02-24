# Community analysis

In this example, we will dive into **community analysis** of the inferred gene regulatory networks. Communities, also called *clusters* or *modules*, are groups of vertices which can be easily grouped together (potentially overlapping) and/or are densely connected internally. In the gene-regulatory network the communities should represent groupings of genes according to their mutual function or molecular processes.

In this document we will focus on `T_cells` cell type only for community analysis. The community detection is performed using [Leiden](https://leidenalg.readthedocs.io/en/stable/index.html) algorithm.

## Loading the graph

Let us first read the data.

```python
import scGRN

# Loading the NetworkX graph
nx_graph = scGRN.ana.get_nx_graph(
    pat=None,  # loading all patients
    cell_type='T_cells',  # cell type
    net_type='all',
    # data_home=_DATA_HOME  # optional
)

# Loading the adjacency list of the graph
adj_list = scGRN.ana.get_adj_list(
    pat=None,  ## loading all patients
    cell_type='T_cells',  # cell type
    net_type='all',
    # data_home=_DATA_HOME  # optional
)

# Loading the scRNA-seq data
sc_data = scGRN.ana.get_sc_data(
    pat=None,  # loading all patients
    cell_type='T_cells',  # cell type
    # data_home=_DATA_HOME  # optional
)

print('Number of nodes:', nx_graph.number_of_nodes())
print('Number of edges:', nx_graph.number_of_edges())
print('Number of cells:', sc_data.shape[0])
```

```python
> Number of nodes: 18920
> Number of edges: 4678478
> Number of cells: 23742
```

## Community detection

### Running community detection

The user can run community detection and processing using the [`scGRN.ana.process_community()`](https://github.com/masyahook/scGRN/blob/c86b1219bdb0df6b415d73fb2ff0f0cd4ebf4a1c/scGRN/network_analysis/_community.py#L380) function. This function could take a few minutes/hours to run depending on the size of the graph, so it's recommended to use a script [`community_ana.py`](../scGRN/network_analysis/community_ana.py) to run it on a cluster (also look into [`community_scripts`](../scGRN/network_analysis/community_scripts)). The script will produce a `.pickle` file with the results of community detection and processing.

An example of running `community_ana.py` script is as follows:

```bash
export PYTHONPATH="<PATH_TO_scGRN>:$PYTHONPATH"

# Running for T cells on all patients
python -m scGRN.network_analysis.community_ana --cell_type T_cells --patient all
```

Also we wrote specific wrapper scripts to run community detection and processing on a cluster which are located in [`community_scripts`](../scGRN/network_analysis/community_scripts). The script [`run_community_ana_pat.sh`](../scGRN/network_analysis/community_scripts/run_community_ana_pat.sh) can be used to run community detection on patient-specific data, while [`run_community_ana_agg.sh`](../scGRN/network_analysis/community_scripts/run_community_ana_agg.sh) is used to run community analysis on patient-aggregated data. Below you can see the example of running the script with `sbatch` command:

```bash
sbatch --job-name='T_cells_all_patients_community_ana_leiden' \
    --chdir=/gpfs/home/bsc08/bsc08890/scGRN_analysis/scGRN/network_analysis/community_scripts \
    --ntasks=1 \
    --time='10:00:00' \
    --output=/dev/null \
    --error=/dev/null \
    --cpus-per-task=24 \
    /gpfs/home/bsc08/bsc08890/scGRN_analysis/scGRN/network_analysis/community_scripts/community_ana_agg.sh \
    leiden T_cells all_patients SBATCH 13
```

All scripts will run the underlying [`scGRN.ana.process_community()`](https://github.com/masyahook/scGRN/blob/c86b1219bdb0df6b415d73fb2ff0f0cd4ebf4a1c/scGRN/network_analysis/_community.py#L380) function. Let us have a look at its signature and docstring:

```python
?scGRN.ana.process_community
```

```python
Signature:
scGRN.ana.process_communities(
    cell_type: str,
    pat: Union[str, NoneType] = None,
    algo: str = 'leiden',
    filter_quantile: float = 0.95,
    if_betweenness: bool = True,
    limit_anno_until: int = 50,
    k: int = 5000,
    save_top_intercommunity_links_until: int = 20,
    other_functions_until: int = 20,
    save_top_new_found_cluster_links: int = 20,
    seed: int = 42,
    data_home: str = '/gpfs/projects/bsc08/shared_projects/scGRN_analysis/Data_home/res/covid_19',
)
Docstring:
    Find communities in the passed graph, annotate them by identifying general community properties
    such as number of nodes/edges, as well as biology-related properties, e.g. top important genes
    and functions. At the end we save metadata to .tsv format along with word clouds in the figs
    folder.

    The output has the following format:

    num_nodes num_edges main_functions_GO                                 main_functions_KEGG                                 main_functions_immunological                     main_functions_hallmark                             non_lambert_2018_TF_central_genes                 non_dorothea_TF_central_genes                     new_gene_gene_links_KEGG                         new_gene_gene_links_hallmark                     ... top_links_scores_with_community_6                 top_links_scores_with_community_7                 top_links_scores_with_community_8                 top_links_scores_with_community_9                 top_links_scores_with_community_10                 top_links_scores_with_community_11                 top_links_scores_with_community_12                 top_links_scores_with_community_13                 top_links_scores_with_community_14                 top_links_scores_with_community_15
0 3978     12830     >>> phosphatase activity <<<: PTPRK,DUSP14,PTE... >>> MAPK signaling pathway <<<: JUN,ELK1,GADD4... >>> Genes down-regulated in effector CD8 T cel... >>> Genes regulated by NF-kB in response to TN... DUSP1 (rank=1); IL7R (rank=3); MTRNR2L12 (rank... DUSP1 (rank=1); KLF2 (rank=2); IL7R (rank=3); ... IL7R (KEGG_HEMATOPOIETIC_CELL_LINEAGE & KEGG_J... TXNIP (HALLMARK_APOPTOSIS & HALLMARK_P53_PATHW... ... DUSP1 <-> JUN (score=33.34); IL7R <-> PTGER2 (... DUSP1 <-> JUN (score=33.34); IL7R <-> PTGER2 (... DUSP1 <-> JUN (score=33.34); IL7R <-> PTGER2 (... IL7R <-> PTGER2 (score=36.18); TXNIP <-> IL7R ... DUSP1 <-> JUN (score=33.34); IL7R <-> PTGER2 (... DUSP1 <-> JUN (score=33.34); IL7R <-> PTGER2 (... DUSP1 <-> JUN (score=33.34); IL7R <-> PTGER2 (... DUSP1 <-> JUN (score=33.34); IL7R <-> PTGER2 (... DUSP1 <-> JUN (score=33.34); IL7R <-> PTGER2 (... DUSP1 <-> JUN (score=33.34); IL7R <-> PTGER2 (...
1 1812     11803     >>> regulation of microtubule polymerization o... >>> MAPK signaling pathway <<<: MEF2C,STMN1,MA... >>> Genes down-regulated in naïve CD8 T cells ... >>> Genes involved in the G2/M checkpoint, as ... STMN1 (rank=1); TYMS (rank=2); HMGB2 (rank=3);... STMN1 (rank=1); TYMS (rank=2); HMGB2 (rank=3);... STMN1 (KEGG_MAPK_SIGNALING_PATHWAY) <-> HMGN2 ... HMGN2 (HALLMARK_G2M_CHECKPOINT) <-> H2AFV (); ... ... STMN1 <-> HMGN2 (score=31.23); TYMS <-> PCNA (... STMN1 <-> HMGN2 (score=31.23); TYMS <-> PCNA (... STMN1 <-> HMGN2 (score=31.23); TYMS <-> PCNA (... STMN1 <-> HMGN2 (score=31.23); TYMS <-> PCNA (... STMN1 <-> HMGN2 (score=31.23); TYMS <-> PCNA (... STMN1 <-> HMGN2 (score=31.23); TYMS <-> PCNA (... STMN1 <-> HMGN2 (score=31.23); TYMS <-> PCNA (... STMN1 <-> HMGN2 (score=31.23); TYMS <-> PCNA (... STMN1 <-> HMGN2 (score=31.23); TYMS <-> PCNA (... STMN1 <-> HMGN2 (score=31.23); TYMS <-> PCNA (...
2 1576     6114     >>> immune response <<<: CCL5,IGKV1-12,CST7,IG... >>> Cytokine-cytokine receptor interaction <<<... >>> Genes down-regulated in naïve CD8 T cells ... >>> Genes regulated by NF-kB in response to TN... CCL5 (rank=1); NKG7 (rank=2); LAG3 (rank=3); S... CCL5 (rank=1); NKG7 (rank=2); LAG3 (rank=3); S... NKG7 () <-> GZMA (KEGG_NEUROACTIVE_LIGAND_RECE... NKG7 () <-> GZMA (HALLMARK_ALLOGRAFT_REJECTION... ... NKG7 <-> GZMA (score=41.34); NKG7 <-> CST7 (sc... NKG7 <-> GZMA (score=41.34); NKG7 <-> CST7 (sc... NKG7 <-> GZMA (score=41.34); NKG7 <-> CST7 (sc... NKG7 <-> GZMA (score=41.34); NKG7 <-> CST7 (sc... NKG7 <-> GZMA (score=41.34); NKG7 <-> CST7 (sc... NKG7 <-> GZMA (score=41.34); NKG7 <-> CST7 (sc... NKG7 <-> GZMA (score=41.34); NKG7 <-> CST7 (sc... NKG7 <-> GZMA (score=41.34); NKG7 <-> CST7 (sc... NKG7 <-> GZMA (score=41.34); NKG7 <-> CST7 (sc... NKG7 <-> GZMA (score=41.34); NKG7 <-> CST7 (sc...
3 1494     11988     >>> nucleic acid binding <<<: TRIM27,PARN,SON,... >>> RIG-I-like receptor signaling pathway <<<:... >>> Genes up-regulated in naïve CD8 T cells co... >>> Genes up-regulated in response to low oxyg... ISG20 (rank=1); MX1 (rank=2); ISG15 (rank=3); ... ISG20 (rank=1); MX1 (rank=2); ISG15 (rank=3); ... B2M (KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION)... MYL12A (HALLMARK_ANDROGEN_RESPONSE) <-> MYL12B... ... B2M <-> LY6E (score=51.16); MYL12A <-> MYL12B ... B2M <-> LY6E (score=51.16); LY6E <-> IFITM1 (s... B2M <-> LY6E (score=51.16); MYL12A <-> MYL12B ... B2M <-> LY6E (score=51.16); MYL12A <-> MYL12B ... B2M <-> LY6E (score=51.16); MYL12A <-> MYL12B ... B2M <-> LY6E (score=51.16); MYL12A <-> MYL12B ... B2M <-> LY6E (score=51.16); MYL12A <-> MYL12B ... B2M <-> LY6E (score=51.16); MYL12A <-> MYL12B ... B2M <-> LY6E (score=51.16); MYL12A <-> MYL12B ... B2M <-> LY6E (score=51.16); MYL12A <-> MYL12B ...
4 1481     21413     >>> nucleotide binding <<<: AK6,DDX52,ABCF1,DD... >>> Pathogenic Escherichia coli infection <<<:... >>> Genes down-regulated in naïve CD8 T cells ... >>> Genes encoding components of apical juncti... ACTB (rank=1); GAPDH (rank=2); ACTG1 (rank=3);... ACTB (rank=1); GAPDH (rank=2); ACTG1 (rank=3);... GAPDH (KEGG_GLYCOLYSIS_GLUCONEOGENESIS & KEGG_... GAPDH (HALLMARK_MTORC1_SIGNALING & HALLMARK_HY... ... GAPDH <-> ACTB (score=62.04); PFN1 <-> MYL6 (s... GAPDH <-> ACTB (score=62.04); PFN1 <-> MYL6 (s... GAPDH <-> ACTB (score=62.04); PFN1 <-> MYL6 (s... GAPDH <-> ACTB (score=62.04); PFN1 <-> MYL6 (s... GAPDH <-> ACTB (score=62.04); PFN1 <-> MYL6 (s... GAPDH <-> ACTB (score=62.04); PFN1 <-> MYL6 (s... GAPDH <-> ACTB (score=62.04); PFN1 <-> MYL6 (s... GAPDH <-> ACTB (score=62.04); PFN1 <-> MYL6 (s... GAPDH <-> ACTB (score=62.04); PFN1 <-> MYL6 (s... GAPDH <-> ACTB (score=62.04); PFN1 <-> MYL6 (s...

    The following columns are generated:

    'num_nodes' -> number of nodes in the community
    'num_edges' -> number of edges in the community
    'main_functions_GO' -> main functions based on GO annotation, the following formatting used >>> func_term <<<: gene_1, gene_2, ..., gene_N; >>> func_term <<<: ...
    'main_functions_KEGG' -> same as 'main_functions_GO' but based on KEGG annotation
    'main_functions_immunological' -> same as 'main_functions_GO' but based on MSigDB_immunological annotation
    'main_functions_hallmark' -> same as 'main_functions_GO' but based on MSigDB_hallmark annotation
    'non_lambert_2018_TF_central_genes' -> genes with highest centrality not included in TF Lambert et al, 2018 list, i.e. potential non-TF regulators
    'non_dorothea_TF_central_genes' -> genes with highest centrality not included in TF DoRothEA list, i.e. potential non-TF regulators
    'new_gene_gene_links_KEGG' -> inter-functional connections between genes (if genes are connected, but no common functional terms are found)
    'new_gene_gene_links_hallmark' -> inter-functional connections between genes (if genes are connected, but no common functional terms are found)
    'whole_G_central_genes_scores' -> centrality score for each gene computed on the whole graph, i.e. what's the most central gene in the whole network?
    'other_functions_GO' -> other functions identifed based on non-central genes (based on GO)
    'other_functions_KEGG' -> other functions identifed based on non-central genes (based on KEGG)
    'other_functions_immunological' -> other functions identifed based on non-central genes (based on MSigDB_immunological)
    'other_functions_hallmark' ->  other functions identifed based on non-central genes (based on MSigDB_hallmark)
    'sorted_central_genes_scores' -> top central genes with scores (sorted)
    'sorted_central_functions_GO' -> similar to 'main_functions_GO', but here grouping by gene instead of grouping by function
    'sorted_central_functions_KEGG' -> similar to 'main_functions_KEGG', but here grouping by gene instead of grouping by function
    'sorted_central_functions_immunological' -> similar to 'main_functions_immunological', but here grouping by gene instead of grouping by function
    'sorted_central_functions_hallmark' -> similar to 'main_functions_hallmark', but here grouping by gene instead of grouping by function
    'most_frequent_function_words_GO' -> most frequent functional terms, visualized on the word cloud (based on GO)
    'most_frequent_function_words_KEGG' -> most frequent functional terms, visualized on the word cloud (based on KEGG)
    'most_frequent_function_words_immunological' -> most frequent functional terms, visualized on the word cloud (based on MSigDB_immunological)
    'most_frequent_function_words_hallmark' -> most frequent functional terms, visualized on the word cloud (based on MSigDB_hallmark)
    'all_sorted_genes' -> all genes with corresponding centrality scores (fuller list than 'sorted_central_genes_scores')
    'top_links_scores_central_genes<->community_i' -> gene-gene links with highest importance between the current community and community_i (only central genes)
    'top_links_scores_with_community_0' -> similar to top_links_scores_central_genes.* but any gene is considered here

    :param cell_type: The cell type name of the data - could be either:
        e.g. 'Macrophage', 'T_cells' (the cell type identifier)
        e.g. 'all', 'all_data' (the aggregated data - include all cell types)
        e.g. 'raw_data_Macrophage', 'raw_data_T_cells', 'raw_data' (the data file identifier, i.e. 'raw_data_T_cells'
            corresponds to the same data as 'T_cells')
    :param net_type: The type of data - could be either:
        'all' (all gene-gene connections)
        'TF' (TF-target connections)
        'ctx' (enriched TF-target connections)
    :param pat: The patient identifier - could be either:
        e.g. None (include all patients)
        e.g. 'C51', 'C141' (the patient identifier)
        e.g. 'C', 'M', 'S', 'all_data', 'all' (the identifier of aggregated patient data)
    :param algo: The algorithm used to identify communities, either 'leiden' or 'louvain'
    :param filter_quantile: The quantile threshold
    :param if_betweenness: True if use betweenness centrality as node importance score, False if use
        closeness centrality
    :param limit_anno_until: Number of genes to use to calculate wordcloud
    :param k: Use k nodes to estimate centrality
    :param save_top_intercommunity_links_until: The number of inter-community links to save with highest importance
    :param other_functions_until: Number of other functional terms to save (not included in the main ones)
    :param save_top_new_found_cluster_links: Number of inter-functional connections between genes to save
    :param seed: A random seed
    :param data_home: The path to the data folder

File:      /gpfs/home/bsc08/bsc08890/scGRN_analysis/scGRN/network_analysis/_community.py
Type:      function
```

As we can see the function takes the `cell_type` and the `pat` input parameters indicating the type of data to use to run community detection. Other parameters are optional and configure the community detection algorithm and reporting.

*Please have a look at the format of output dataframe, i.e. which community information is saved.*

### Loading the results

Let's imagine we ran the community detection on the cluster and now we want to process the results. We can load the output using the function [`scGRN.ana.get_community_info()`](https://github.com/masyahook/scGRN/blob/c86b1219bdb0df6b415d73fb2ff0f0cd4ebf4a1c/scGRN/network_analysis/_data_processing.py#L861):

```python
comm_info = scGRN.ana.get_community_info(
    pat=None,  # loading all patients
    cell_type='T_cells',  # cell type
    # data_home=_DATA_HOME  # optional
)
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
       'top_links_scores_with_community_12'],
      dtype='object')
```

```python
print(comm_info)
```

```python
   num_nodes num_edges                                  main_functions_GO  \
0       5888     19001  >>> phosphatase activity <<<: PTPRK,PTEN,MDP1,...   
1       2536     12198  >>> immune response <<<: IGLV1-40,IGLV1-44,IGL...   
2       2011     10242  >>> immune response <<<: HLA-A,HLA-B,HLA-E,HLA...   
3       1628     22856  >>> nucleotide binding <<<: ABCF1,PSMC4,MAP4K1...   
4       1434      6960  >>> immune response <<<: KIR2DL3,FCAR,KIR2DL1,...   
5       1342      8926  >>> regulation of microtubule polymerization o...   
6        920     10911  >>> protein binding <<<: IRF7,AGK,TRIM69,PSMB9...   
7        788     10993  >>> GTP binding <<<: GLUD2,RHOH,GSPT1,EIF2S3,N...   
8        686      2290  >>> cell differentiation <<<: MEIG1,ID4,FOXC1,...   
9        664      1944  >>> RNA binding <<<: MAPT,PPP1R10,HNRNPH1,DDX5...   
10       662      1950  >>> immune response <<<: IGHV3-48,SEMA7A,LIF,B...   
11       172       625  >>> protein binding <<<: DERL3,IGHG1,PARN,MRM1...   
12         2         1  >>> DNA binding <<<: HOXA5; >>> regulation of ...   

                                  main_functions_KEGG  \
0   >>> MAPK signaling pathway <<<: JUN,MEF2C,GADD...   
1   >>> Cytokine-cytokine receptor interaction <<<...   
2   >>> Endocytosis <<<: PIKFYVE,CHMP2B,ARFGAP3,VP...   
3   >>> Pathogenic Escherichia coli infection <<<:...   
4   >>> Cytokine-cytokine receptor interaction <<<...   
5   >>> MAPK signaling pathway <<<: ELK1,STMN1,AKT...   
6   >>> RIG-I-like receptor signaling pathway <<<:...   
7   >>> Cytokine-cytokine receptor interaction <<<...   
8   >>> Arginine and proline metabolism <<<: AZIN2...   
9   >>> Ribosome <<<: RPS10; >>> Oxidative phospho...   
10  >>> Cell adhesion molecules (CAMs) <<<: CD28,J...   
11  >>> Arginine and proline metabolism <<<: ODC1;...   
12                                                      

                         main_functions_immunological  \
0   >>> Genes down-regulated in effector CD8 T cel...   
1   >>> Genes down-regulated in naïve CD8 T cells ...   
2   >>> Genes up-regulated in comparison of Th1 ce...   
3   >>> Genes down-regulated in comparison of heal...   
4   >>> Genes down-regulated in comparison of heal...   
5   >>> Genes down-regulated in naïve CD8 T cells ...   
6   >>> Genes up-regulated in comparison of naive ...   
7   >>> Genes up-regulated in comparison of health...   
8   >>> Genes down-regulated in comparison of CD4 ...   
9   >>> Genes down-regulated in comparison of naiv...   
10  >>> Genes down-regulated in naïve CD8 T cells ...   
11  >>> Genes down-regulated in effector CD8 T cel...   
12  >>> Genes up-regulated in naïve CD8 T cells co...   

                              main_functions_hallmark  \
0   >>> Genes regulated by NF-kB in response to TN...   
1   >>> Genes regulated by NF-kB in response to TN...   
2   >>> Genes up-regulated in response to alpha in...   
3   >>> Genes involved in cholesterol homeostasis....   
4   >>> Genes regulated by NF-kB in response to TN...   
5   >>> Genes involved in the G2/M checkpoint, as ...   
6   >>> Genes up-regulated in response to alpha in...   
7   >>> A subgroup of genes regulated by MYC - ver...   
8   >>> Genes defining early response to estrogen....   
9   >>> A subgroup of genes regulated by MYC - ver...   
10  >>> Genes up-regulated by STAT5 in response to...   
11  >>> Genes up-regulated in response to low oxyg...   
12                                                      

                    non_lambert_2018_TF_central_genes  \
0   DUSP1 (rank=1); IL7R (rank=2); TSC22D3 (rank=4...   
1   CCL5 (rank=1); NKG7 (rank=3); TRDV1 (rank=4); ...   
2   HLA-C (rank=1); HLA-A (rank=2); MALAT1 (rank=3...   
3   ACTG1 (rank=1); ACTB (rank=2); TPM4 (rank=3); ...   
4   CCL2 (rank=1); FTH1 (rank=2); FCER1G (rank=3);...   
5   STMN1 (rank=1); HMGB2 (rank=2); TYMS (rank=3);...   
6   ISG15 (rank=1); SAT1 (rank=2); MX1 (rank=3); T...   
7   EEF1A1 (rank=1); PABPC1 (rank=2); CCR7 (rank=3...   
8   FAM183A (rank=1); C9orf24 (rank=2); DYNLRB2 (r...   
9   RPS10 (rank=1); MT-ATP6 (rank=2); C1orf56 (ran...   
10  CTLA4 (rank=1); SRGN (rank=2); TNFRSF18 (rank=...   
11  IGHG1 (rank=1); IGLV3-19 (rank=2); MT1E (rank=...   
12                                      CIB2 (rank=2)   

                        non_dorothea_TF_central_genes  \
0   DUSP1 (rank=1); IL7R (rank=2); KLF2 (rank=3); ...   
1   CCL5 (rank=1); ZNF683 (rank=2); NKG7 (rank=3);...   
2   HLA-C (rank=1); HLA-A (rank=2); MALAT1 (rank=3...   
3   ACTG1 (rank=1); ACTB (rank=2); TPM4 (rank=3); ...   
4   CCL2 (rank=1); FTH1 (rank=2); FCER1G (rank=3);...   
5   STMN1 (rank=1); HMGB2 (rank=2); TYMS (rank=3);...   
6   ISG15 (rank=1); SAT1 (rank=2); MX1 (rank=3); T...   
7   EEF1A1 (rank=1); PABPC1 (rank=2); CCR7 (rank=3...   
8   FAM183A (rank=1); C9orf24 (rank=2); DYNLRB2 (r...   
9   RPS10 (rank=1); MT-ATP6 (rank=2); C1orf56 (ran...   
10  CTLA4 (rank=1); SRGN (rank=2); TNFRSF18 (rank=...   
11  IGHG1 (rank=1); IGLV3-19 (rank=2); MT1E (rank=...   
12                      HOXA5 (rank=1); CIB2 (rank=2)   

                             new_gene_gene_links_KEGG  \
0   NFKBIA (KEGG_RIG_I_LIKE_RECEPTOR_SIGNALING_PAT...   
1   GZMA (KEGG_NEUROACTIVE_LIGAND_RECEPTOR_INTERAC...   
2   HLA-C (KEGG_ALLOGRAFT_REJECTION & KEGG_ANTIGEN...   
3   ACTB (KEGG_DILATED_CARDIOMYOPATHY & KEGG_ADHER...   
4   FTH1 (KEGG_PORPHYRIN_AND_CHLOROPHYLL_METABOLIS...   
5   PTMA () <-> DUT (KEGG_PYRIMIDINE_METABOLISM); ...   
6   IFITM1 (KEGG_B_CELL_RECEPTOR_SIGNALING_PATHWAY...   
7   EEF1A1 () <-> RPL30 (KEGG_RIBOSOME); RPS12 (KE...   
8   PAPOLG (KEGG_RNA_DEGRADATION) <-> MSMB (); PCD...   
9                                                       
10  SRGN () <-> PDCD1 (KEGG_T_CELL_RECEPTOR_SIGNAL...   
11  IGLV3-19 () <-> HSP90B1 (KEGG_NOD_LIKE_RECEPTO...   
12                                                      

                         new_gene_gene_links_hallmark  ...  \
0   INSIG1 (HALLMARK_MTORC1_SIGNALING & HALLMARK_U...  ...   
1   GZMA (HALLMARK_COMPLEMENT & HALLMARK_INTERFERO...  ...   
2   HLA-B (HALLMARK_INTERFERON_GAMMA_RESPONSE) <->...  ...   
3   ACTB (HALLMARK_APICAL_JUNCTION) <-> CORO1A (HA...  ...   
4   HLA-DRB1 (HALLMARK_INTERFERON_GAMMA_RESPONSE) ...  ...   
5   PTMA () <-> DUT (HALLMARK_E2F_TARGETS & HALLMA...  ...   
6   IFI6 () <-> LY6E (HALLMARK_INTERFERON_GAMMA_RE...  ...   
7   RPL32 () <-> RPL34 (HALLMARK_MYC_TARGETS_V1); ...  ...   
8   PAPOLG () <-> MSMB (HALLMARK_ESTROGEN_RESPONSE...  ...   
9                                                      ...   
10  SRGN (HALLMARK_ALLOGRAFT_REJECTION) <-> PDCD1 ...  ...   
11  MT1G () <-> MT1E (HALLMARK_HYPOXIA & HALLMARK_...  ...   
12                                                     ...   

                    top_links_scores_with_community_3  \
0   FOS <-> JUN (score=36.54); NFKBIA <-> JUNB (sc...   
1   GZMA <-> HCST (score=46.88); XCL1 <-> XCL2 (sc...   
2   HLA-A <-> HLA-B (score=105.29); HLA-B <-> HLA-...   
3                                                 NaN   
4   HLA-DRA <-> CD74 (score=79.94); HLA-DRB1 <-> H...   
5   PTMA <-> DUT (score=35.94); STMN1 <-> HMGN2 (s...   
6   IFI6 <-> LY6E (score=51.59); IFITM1 <-> IFITM2...   
7   RPL13 <-> RPL32 (score=59.80); RPL32 <-> RPL34...   
8   AC097382.2 <-> TBC1D14 (score=19.49); PAPOLG <...   
9   MT-CO3 <-> MT-ND3 (score=67.22); MT-CYB <-> MT...   
10  TNFRSF18 <-> IL2RA (score=24.51); SRGN <-> SNX...   
11  IGHG4 <-> IGHG1 (score=51.20); IGHG1 <-> IGLV3...   
12                        CIB2 <-> HOXA5 (score=2.10)   

                    top_links_scores_with_community_4  \
0   FOS <-> JUN (score=36.54); NFKBIA <-> JUNB (sc...   
1   GZMA <-> HCST (score=46.88); XCL1 <-> XCL2 (sc...   
2   HLA-A <-> HLA-B (score=105.29); HLA-B <-> HLA-...   
3   ACTB <-> CORO1A (score=62.61); B2M <-> TMSB4X ...   
4                                                 NaN   
5   PTMA <-> DUT (score=35.94); STMN1 <-> HMGN2 (s...   
6   IFI6 <-> LY6E (score=51.59); IFITM1 <-> IFITM2...   
7   RPL13 <-> RPL32 (score=59.80); RPL32 <-> RPL34...   
8   AC097382.2 <-> TBC1D14 (score=19.49); PAPOLG <...   
9   MT-CO3 <-> MT-ND3 (score=67.22); MT-CYB <-> MT...   
10  SRGN <-> PDCD1 (score=20.38); TNFRSF18 <-> IL2...   
11  IGHG4 <-> IGHG1 (score=51.20); IGHG1 <-> IGLV3...   
12                        CIB2 <-> HOXA5 (score=2.10)   

                    top_links_scores_with_community_5  \
0   FOS <-> JUN (score=36.54); NFKBIA <-> JUNB (sc...   
1   GZMA <-> HCST (score=46.88); XCL1 <-> XCL2 (sc...   
2   HLA-A <-> HLA-B (score=105.29); HLA-B <-> HLA-...   
3   ACTB <-> ARPC2 (score=51.12); ACTB <-> CORO1A ...   
4   HLA-DRB1 <-> HLA-DRA (score=126.13); S100A12 <...   
5                                                 NaN   
6   IFI6 <-> LY6E (score=51.59); IFITM1 <-> IFITM2...   
7   RPL13 <-> RPL32 (score=59.80); RPL32 <-> RPL34...   
8   AC097382.2 <-> TBC1D14 (score=19.49); PAPOLG <...   
9   MT-CO3 <-> MT-ND3 (score=67.22); MT-CYB <-> MT...   
10  SRGN <-> PDCD1 (score=20.38); TNFRSF18 <-> IL2...   
11  IGHG4 <-> IGHG1 (score=51.20); IGHG1 <-> IGLV3...   
12                        CIB2 <-> HOXA5 (score=2.10)   

                    top_links_scores_with_community_6  \
0   FOS <-> JUN (score=36.54); NFKBIA <-> JUNB (sc...   
1   GZMA <-> HCST (score=46.88); XCL1 <-> XCL2 (sc...   
2   HLA-A <-> HLA-B (score=105.29); HLA-B <-> HLA-...   
3   ACTB <-> ARPC2 (score=51.12); ACTB <-> CORO1A ...   
4   HLA-DRB1 <-> HLA-DRA (score=126.13); S100A12 <...   
5   PTMA <-> DUT (score=35.94); STMN1 <-> HMGN2 (s...   
6                                                 NaN   
7   RPL13 <-> RPL32 (score=59.80); RPL32 <-> RPL34...   
8   AC097382.2 <-> TBC1D14 (score=19.49); PAPOLG <...   
9   MT-CO3 <-> MT-ND3 (score=67.22); MT-CYB <-> MT...   
10  TNFRSF18 <-> IL2RA (score=24.51); SRGN <-> SNX...   
11  IGHG4 <-> IGHG1 (score=51.20); IGHG1 <-> IGLV3...   
12                        CIB2 <-> HOXA5 (score=2.10)   

                    top_links_scores_with_community_7  \
0   FOS <-> JUN (score=36.54); NFKBIA <-> JUNB (sc...   
1   GZMA <-> HCST (score=46.88); XCL1 <-> XCL2 (sc...   
2   HLA-A <-> HLA-B (score=105.29); HLA-B <-> HLA-...   
3   ACTB <-> ARPC2 (score=51.12); ACTB <-> CORO1A ...   
4   HLA-DRB1 <-> HLA-DRA (score=126.13); S100A12 <...   
5   PTMA <-> DUT (score=35.94); STMN1 <-> HMGN2 (s...   
6   IFI6 <-> LY6E (score=51.59); IFITM1 <-> IFITM2...   
7                                                 NaN   
8   AC097382.2 <-> TBC1D14 (score=19.49); PAPOLG <...   
9   MT-CO3 <-> MT-ND3 (score=67.22); MT-CYB <-> MT...   
10  TNFRSF18 <-> IL2RA (score=24.51); SRGN <-> SNX...   
11  IGHG4 <-> IGHG1 (score=51.20); IGHG1 <-> IGLV3...   
12                        CIB2 <-> HOXA5 (score=2.10)   

                    top_links_scores_with_community_8  \
0   FOS <-> JUN (score=36.54); NFKBIA <-> JUNB (sc...   
1   GZMA <-> HCST (score=46.88); XCL1 <-> XCL2 (sc...   
2   HLA-A <-> HLA-B (score=105.29); HLA-B <-> HLA-...   
3   ACTB <-> ARPC2 (score=51.12); ACTB <-> CORO1A ...   
4   HLA-DRB1 <-> HLA-DRA (score=126.13); S100A12 <...   
5   PTMA <-> DUT (score=35.94); STMN1 <-> HMGN2 (s...   
6   IFI6 <-> LY6E (score=51.59); IFITM1 <-> IFITM2...   
7   RPL13 <-> RPL32 (score=59.80); RPL32 <-> RPL34...   
8                                                 NaN   
9   MT-CO3 <-> MT-ND3 (score=67.22); MT-CYB <-> MT...   
10  SRGN <-> PDCD1 (score=20.38); TNFRSF18 <-> IL2...   
11  IGHG4 <-> IGHG1 (score=51.20); IGHG1 <-> IGLV3...   
12                        CIB2 <-> HOXA5 (score=2.10)   

                    top_links_scores_with_community_9  \
0   FOS <-> JUN (score=36.54); NFKBIA <-> JUNB (sc...   
1   GZMA <-> HCST (score=46.88); XCL1 <-> XCL2 (sc...   
2   HLA-A <-> HLA-B (score=105.29); HLA-B <-> HLA-...   
3   ACTB <-> ARPC2 (score=51.12); ACTB <-> CORO1A ...   
4   HLA-DRB1 <-> HLA-DRA (score=126.13); S100A12 <...   
5   PTMA <-> DUT (score=35.94); STMN1 <-> HMGN2 (s...   
6   IFI6 <-> LY6E (score=51.59); IFITM1 <-> IFITM2...   
7   RPL13 <-> RPL32 (score=59.80); RPL32 <-> RPL34...   
8   AC097382.2 <-> TBC1D14 (score=19.49); PAPOLG <...   
9                                                 NaN   
10  SRGN <-> PDCD1 (score=20.38); TNFRSF18 <-> IL2...   
11  IGHG4 <-> IGHG1 (score=51.20); IGHG1 <-> IGLV3...   
12                        CIB2 <-> HOXA5 (score=2.10)   

                   top_links_scores_with_community_10  \
0   FOS <-> JUN (score=36.54); NFKBIA <-> JUNB (sc...   
1   GZMA <-> HCST (score=46.88); XCL1 <-> XCL2 (sc...   
2   HLA-A <-> HLA-B (score=105.29); HLA-B <-> HLA-...   
3   ACTB <-> ARPC2 (score=51.12); ACTB <-> CORO1A ...   
4   HLA-DRB1 <-> HLA-DRA (score=126.13); S100A12 <...   
5   PTMA <-> DUT (score=35.94); STMN1 <-> HMGN2 (s...   
6   IFI6 <-> LY6E (score=51.59); IFITM1 <-> IFITM2...   
7   RPL13 <-> RPL32 (score=59.80); RPL32 <-> RPL34...   
8   AC097382.2 <-> TBC1D14 (score=19.49); PAPOLG <...   
9   MT-CO3 <-> MT-ND3 (score=67.22); MT-CYB <-> MT...   
10                                                NaN   
11  IGHG4 <-> IGHG1 (score=51.20); IGHG1 <-> IGLV3...   
12                        CIB2 <-> HOXA5 (score=2.10)   

                   top_links_scores_with_community_11  \
0   FOS <-> JUN (score=36.54); NFKBIA <-> JUNB (sc...   
1   GZMA <-> HCST (score=46.88); XCL1 <-> XCL2 (sc...   
2   HLA-A <-> HLA-B (score=105.29); HLA-B <-> HLA-...   
3   ACTB <-> ARPC2 (score=51.12); ACTB <-> CORO1A ...   
4   HLA-DRB1 <-> HLA-DRA (score=126.13); S100A12 <...   
5   PTMA <-> DUT (score=35.94); STMN1 <-> HMGN2 (s...   
6   IFI6 <-> LY6E (score=51.59); IFITM1 <-> IFITM2...   
7   RPL13 <-> RPL32 (score=59.80); RPL32 <-> RPL34...   
8   AC097382.2 <-> TBC1D14 (score=19.49); PAPOLG <...   
9   MT-CO3 <-> MT-ND3 (score=67.22); MT-CYB <-> MT...   
10  SRGN <-> PDCD1 (score=20.38); TNFRSF18 <-> IL2...   
11                                                NaN   
12                        CIB2 <-> HOXA5 (score=2.10)   

                   top_links_scores_with_community_12  
0   FOS <-> JUN (score=36.54); NFKBIA <-> JUNB (sc...  
1   GZMA <-> HCST (score=46.88); XCL1 <-> XCL2 (sc...  
2   HLA-A <-> HLA-B (score=105.29); HLA-B <-> HLA-...  
3   ACTB <-> ARPC2 (score=51.12); ACTB <-> CORO1A ...  
4   HLA-DRB1 <-> HLA-DRA (score=126.13); S100A12 <...  
5   PTMA <-> DUT (score=35.94); STMN1 <-> HMGN2 (s...  
6   IFI6 <-> LY6E (score=51.59); IFITM1 <-> IFITM2...  
7   RPL13 <-> RPL32 (score=59.80); RPL32 <-> RPL34...  
8   AC097382.2 <-> TBC1D14 (score=19.49); PAPOLG <...  
9   MT-CO3 <-> MT-ND3 (score=67.22); MT-CYB <-> MT...  
10  SRGN <-> PDCD1 (score=20.38); TNFRSF18 <-> IL2...  
11  IGHG4 <-> IGHG1 (score=51.20); IGHG1 <-> IGLV3...  
12                                                NaN
```

Let us do a quick EDA, we will start with printing the number of nodes and edges in each community:

```python
# Let us see how many nodes and edges we have per community:
for community_id, community_row in comm_info.iterrows():
    print(f'Community {community_id} has {community_row["num_nodes"]} nodes and {community_row["num_edges"]} edges')
```

```python
Community 0 has 5888 nodes and 19001 edges
Community 1 has 2536 nodes and 12198 edges
Community 2 has 2011 nodes and 10242 edges
Community 3 has 1628 nodes and 22856 edges
Community 4 has 1434 nodes and 6960 edges
Community 5 has 1342 nodes and 8926 edges
Community 6 has 920 nodes and 10911 edges
Community 7 has 788 nodes and 10993 edges
Community 8 has 686 nodes and 2290 edges
Community 9 has 664 nodes and 1944 edges
Community 10 has 662 nodes and 1950 edges
Community 11 has 172 nodes and 625 edges
Community 12 has 2 nodes and 1 edges
```

Let us see the top-5 main functions of first 5 communities using [GO](https://geneontology.org) and [KEGG](https://www.kegg.jp) ontologies (with gene sets that were used for functional annotation):

```python
TOP_N = 5
FIRST = 5
for community_id, community_row in comm_info.head(FIRST).iterrows():
    print(f'Community {community_id} has the following top {TOP_N} GO terms:')
    print('\n'.join(community_row['main_functions_GO'].split(';')[:TOP_N]))
    print()
```

```python
Community 0 has the following top 5 GO terms:
>>> phosphatase activity <<<: PTPRK,PTEN,MDP1,PHLPP1,STYX,PPM1A,SYNJ1,PTPN21,PTPRM,PUDP,CDC14B,PTPRO,PTPN3,LHPP,PTPN13,MTMR7,RNGTT,MTMR1,PPM1F,DUSP6,CTDSPL,DUSP2,SSH1,MTMR2,PPM1E,PTPRU,INPP5B,PPM1L,CDC14A,SSH3,PTPRB,DUSP10,DUSP11,PTPN14,PHLPP2,PGAM5,PHOSPHO2,INPPL1,DUSP1,ACP6,DUSP12,PRUNE1,PPM1J
 >>> protein tyrosine phosphatase activity <<<: PTPRK,PTEN,MDP1,EYA2,PTPN21,PTPRM,MTMR3,CDC14B,PTPRO,PTPN3,PTPN13,MTMR7,RNGTT,MTMR1,DUSP6,DUSP2,SSH1,MTMR2,EYA3,PTPRU,CDC14A,SSH3,PTPRB,DUSP10,DUSP11,PTPN14,DUSP1,DUSP12
 >>> protein dephosphorylation <<<: PTPRK,PPP1R3B,MDP1,PPP2R3B,EYA2,PHLPP1,STYX,PPM1A,PPP1R3D,PTPN21,PTPRM,MTMR3,PPP4R1,PPA2,CDC14B,PTPRO,PTPN3,CTDP1,PTPN13,MTMR7,CPPED1,RNGTT,PPM1F,DUSP6,CTDSPL,DUSP2,SSH1,MTMR2,EYA3,PPM1E,PTPRU,CAMK2G,PPM1L,MYH3,PPP2CB,PTEN,CDC14A,SSH3,STYXL1,PTPRB,DUSP10,DUSP11,PTPN14,PHLPP2,DLG1,DUSP1,STK11,DUSP12,PPM1J
 >>> protein tyrosine/serine/threonine phosphatase activity <<<: PTEN,STYX,MTMR3,CDC14B,RNGTT,PPM1F,DUSP6,DUSP2,SSH1,MTMR2,CDC14A,SSH3,STYXL1,DUSP10,DUSP11,DUSP1,DUSP12
 >>> dephosphorylation <<<: PTPRK,PTEN,MDP1,NANP,PHLPP1,STYX,PPM1A,ENTPD6,PTPN21,APTX,PFKFB3,PTPRM,MTMR3,PUDP,THTPA,CDC14B,LPIN3,PTPRO,PTPN3,LHPP,PTPN13,MTMR7,PLPP2,RNGTT,MTMR1,DUSP6,CTDSPL,DUSP2,SSH1,MTMR2,NT5DC3,PTPRU,PNKP,CDC14A,SSH3,STYXL1,PTPRB,DUSP10,DUSP11,PTPN14,PHLPP2,G6PC3,PHOSPHO2,DUSP1,PFKFB4,ACP6,DUSP12,PRUNE1

Community 1 has the following top 5 GO terms:
>>> immune response <<<: IGLV1-40,IGLV1-44,IGLV1-36,HLA-DQB2,CCL5,IGKV3-11,IGKV1-16,CST7,TRDV3,AIRE,CXCR3,IL2RG,TRGV3,IGLV4-69,IGF1R,CD8B,NOTCH1,CX3CL1,GZMA,TRAV8-4,TRAV14DV4,TRAV9-2,TRDV1,CBLB,SPN,ITGAD,CTSW,IL1A,TGFBR3,ATP6V0A2,CD7,SECTM1,CD96,CX3CR1,TNFSF4,TLR3,IL16,FASLG,XCL2,XCL1,CD8A,CCL25,CD244,CXCR6,CCR9,CCR1,CCR2,TLR5,FCGR3A,CD86,CCR5,CTSK,PXDN,MR1
 >>> chemokine activity <<<: CCL5,CX3CL1,CKLF,XCL2,XCL1,CCL25
 >>> exocytosis <<<: CCL5,SLC17A9,CPLX1,STX3,STXBP1,PAK1,SNX19,ARFGEF1,RAB27A,SYTL3,UNC13D,CADPS2,LLGL2,SYTL2,CCR1
 >>> cellular process <<<: TBCD,CCL5,TPST2,ESCO1,DYDC1,NR1H3,SNX13,APH1B,ARFIP1,IGF1R,PAK1,CX3CL1,ARL6,FBXO2,PARVA,HK3,CYFIP2,SMPD3,VCAM1,GSR,SNX14,PRKG2,ELMO3,MLH1,INPP4A,SNX25,CELF2,TAF8,DUSP7,HDAC3,HDAC5,PEX12,OBSCN,MIB2,HMGCL
 >>> protein phosphorylation <<<: PRRT1,RPS6KA1,CCL5,STKLD1,ADCK5,TEC,MOK,CDK8,MAPK1,STK38,TGFBR1,SLK,STYK1,CDKL2,MATK,STK31,TYK2,CDK20,PIK3CG,RPS6KA4,MAP2K6,FYN,IGF1R,PAK1,NEK5,MAP4K5,NPRL2,CTBP1,MTOR,RPS6KA3,MAP3K6,MET,BMPR1B,LIMK1,DAPK2,TGFBR2,GRK4,MAPK14,DYRK2,PKN1,PRKCH,ULK4,PRKG2,MARK2,GRK2,CSNK1G2,CDK13,SCYL3,HCST,TTBK2,MAPK9,PRKCB,TEX14,VRK2,RPS6KC1,GSK3B,CDC7,WNK1,MKNK2,AATK,ABR,PPP1R9B,PPM1D,GRK7,RUNX3,MAPKAPK2,EIF2AK3,MST1R,BIRC6,BMPR2,FGR,CDKL4,MAP4K3,OBSCN,MAP3K19,PIK3CD

Community 2 has the following top 5 GO terms:
>>> immune response <<<: HLA-A,HLA-B,HLA-E,HLA-C,IGKV6-21,IL27RA,MBP,GPR65,CHUK,PNP,TRAV2,IGLV9-49,ADGRE5,CYSLTR2,XBP1,IKBKG,TCF12,TRAV18,TRAV38-2DV8,TRAV40,SEMA4D,TRBV7-9,LCP2,FYB1,CXCR5,SCAP,IGKV1D-33,ZAP70
 >>> antigen processing and presentation <<<: HLA-A,HLA-B,HLA-E,HLA-C,RAB6A,AP3B1,RAB35
 >>> interferon-gamma-mediated signaling pathway <<<: HLA-A,HLA-B,HLA-E,HLA-C,IRF9,JAK1
 >>> antigen processing and presentation of exogenous peptide antigen via MHC class I, TAP-independent <<<: HLA-A,HLA-B,HLA-E,HLA-C,LNPEP
 >>> leukocyte activation <<<: HLA-A,HLA-B,HLA-E,HLA-C,TOLLIP

Community 3 has the following top 5 GO terms:
>>> nucleotide binding <<<: ABCF1,PSMC4,MAP4K1,RUVBL1,RALBP1,CCT8,ATP5F1A,LTB4R,DRG1,RIOK3,CTPS2,PSMC6,MTHFD1,DHX15,ABCG1,PFKL,PSMC1,RAC2,UBE2L3,UBE2K,UBE2A,TUBB6,RAB33A,MYLK4,KIF5B,ITPA,IRAK1,PRKAR1B,XRCC6,UBA6,PRPS2,EIF2AK4,CCT6A,STK17A,BLVRA,MFHAS1,DDX27,HBS1L,ASCC3,ULK3,PGK1,UBE2Z,ARF5,CCT5,MYH9,ABCE1,SRP54,CCT2,MYO1G,HPRT1,UBE2N,FARS2,MAP2K3,UBE2R2,IDH3G,RAP1B,GNAI2,CDK6,UBA1,RIOK2,CDK5,RAB11A,OXSR1,RALA,RAB1B,PAICS,UBE2D2,GK,UBA2,MAPKAPK3,SUCLG1,ACTR3,ARL1,TAOK3,ACTB,PTPA,MARS2,SAR1A,PFKP,TUFM,NUDT2,HK1,PSMC2,MTHFD1L,PSMC5,UBE2E3,RPS6KB2,DDX19A,RAC1,XRCC5,ARF3,DDX31,TTLL7,RAB18,VRK3,RAB5C,NUBP2,GRK6,FARSA,ILK,ARF6,SPHK1,VDAC3,EFTUD2,STK25,HINT1,JAK3,RAB11B,PSMC3,UBE2I,ATP5F1B,DCTPP1,RUVBL2,CSNK1A1,GAPDH,CSK,TIMM44,DLD,MAP2K1,DDX59,DDX46,ARL8B,STK16,DDX49,COASY,MTIF2,PKM,RAB4A,DHX9,PAK2,ACTG1,DDX1,HNRNPU,GALK1,RAB40B,RHOA,SMARCA4,KRAS,UBE2J2,UBE2G1,PRKCD,SLFN11,RAN,RHOF,MYH14,CHFR,PRKAR1A,ACLY,DGUOK,EIF5B,NVL,RAB10,LCK,GUK1,MAP3K2,AK2,GNAI3,ACTR2,UBE2Q1,NAXE,RTCA,CDC42,PMVK,IARS2
 >>> ATP binding <<<: ABCF1,PSMC4,MAP4K1,RUVBL1,RALBP1,CCT8,ATP5F1A,RIOK3,CTPS2,PSMC6,MTHFD1,DHX15,ABCG1,PFKL,PSMC1,UBE2L3,UBE2K,UBE2A,MYLK4,KIF5B,IRAK1,XRCC6,UBA6,PRPS2,EIF2AK4,EIF2B2,CCT6A,STK17A,DDX27,ATP6V1B2,ASCC3,ULK3,PGK1,UBE2Z,CCT5,DNAJA2,MYH9,ABCE1,CCT2,MYO1G,UBE2N,FARS2,MAP2K3,UBE2R2,IDH3G,CDK6,UBA1,RIOK2,CDK5,OXSR1,PAICS,UBE2D2,GK,UBA2,MAPKAPK3,PPP5C,ACTR3,TAOK3,ACTB,PTPA,MARS2,PFKP,HK1,PSMC2,MTHFD1L,PSMC5,UBE2E3,RPS6KB2,DDX19A,XRCC5,DDX31,TTLL7,VRK3,NUBP2,GRK6,FARSA,ILK,SPHK1,STK25,JAK3,PSMC3,UBE2I,ATP5F1B,RUVBL2,CSNK1A1,CSK,TIMM44,MAP2K1,DDX59,DDX46,STK16,DDX49,COASY,PKM,NRBP1,DHX9,TWF2,PAK2,ACTG1,DDX1,HNRNPU,GALK1,ATP5F1D,SMARCA4,UBE2J2,UBE2G1,PRKCD,SLFN11,MYH14,ACLY,DGUOK,NVL,LCK,GUK1,MAP3K2,AK2,ACTR2,UBE2Q1,RTCA,PMVK,IARS2
 >>> protein binding <<<: LSM14A,RNH1,LAIR1,TIMM22,SSBP1,NDUFA3,SMARCB1,MIF,TAF9,NCF4,B2M,CHCHD10,PRPF31,ARHGAP27,PPP1R18,LSM2,CLIC1,YWHAE,FRG1,PUF60,IFI27,PFDN6,ABCF1,CSNK2B,EIF1AY,AATF,NTAN1,PSMB3,HMOX2,PSMC4,BORCS5,MAP4K1,ACTN4,ECH1,HNRNPL,SIRT2,NFKBIB,MRPS12,FLII,DNAJC8,TMEM14C,ATP5IF1,NDUFS3,MTCH2,BOP1,SLC52A2,MFHAS1,CPSF1,VPS28,RUVBL1,RNF187,GTF2H1,LDHA,DCAF11,MCRIP1,MCTS1,PSMG2,PSME1,PSME2,NDUFV2,RALBP1,MRPS6,ACAA2,XRN2,SNRPB2,CHAF1B,NEDD8,CCT8,DSTN,ATP5F1A,TMSB4X,HAUS1,SMAD2,HCCS,ASCC2,DRG1,ELP5,RIOK3,YWHAB,EIF5A,POMP,ACAP1,FKBP3,SLC35E4,CTPS2,PDHA1,PSMC6,TSR2,KPNA3,PPP2R5E,SYAP1,PREX1,CFAP298,DNAJC15,MTHFD1,SUMO3,RRP1B,RALY,DHX15,SAMM50,HDAC8,ATP5ME,PFDN4,PSMA7,PIN4,WDR1,GID8,SLC19A1,EIF2S2,WAS,PSMG1,SH3KBP1,ATXN10,HMGN1,SLC25A5,RCL1,PLGRKT,COMMD8,ITM2A,ABCG1,ARHGEF6,RTRAF,RGS19,ATP5F1C,PARVG,RBM25,RIBC2,ADRM1,PHF20,EXOC5,LGALS1,ANKRD54,FUNDC2,NUDT5,TSPO,PFKL,PDS5A,IPO5,PSMC1,VBP1,ITGB2,RBX1,ATP5F1E,CUBN,TMEM33,GLRX5,ATP5PF,PAX5,YWHAH,RAC2,FAM50A,PSAT1,IL10RB,SNRPD3,TTC7B,UBE2L3,UBE2K,UBE2A,DNAJC1,DNMT3B,OSTF1,ACTR10,SNRPD1,MAPRE1,TUBB6,RAB33A,AIFM1,MYL12B,RSU1,RRP7A,TSSC4,TXN,MYLK4,PPIL1,COX7B,KIF5B,RPS19BP1,ATP6V1E1,CCAR1,DDRGK1,ANXA11,PSMF1,SPAG7,ARHGAP4,ACOT8,CISD2,PCIF1,FKBP1A,PSMB5,HCFC1,CAPRIN1,BID,IRAK1,CRELD2,CCNQ,UQCC2,PRKAR1B,PARVB,CHMP4B,C9orf78,PHF5A,NDUFA4,EIF6,PMM1,DESI1,XRCC6,COX19,PLP2,DYNLRB1,SNX5,EMC6,SLC25A11,SNU13,UBA6,SUSD3,SLIRP,PFN1,SLC52A1,CDC5L,PRPS2,RASSF2,EIF2AK4,ROMO1,COMT,RHNO1,EIF2B2,CAMTA2,TXNDC17,HIRA,LYRM2,UFD1,EXOC6,ARPC5L,PSMB7,FLNA,GSTO1,BLOC1S2,ABTB2,COMTD1,SNW1,SRSF3,SMPD1,TMEM14A,ZSWIM7,BCAT1,NONO,ASL,LAMTOR3,EIF3A,CCT6A,CLEC6A,NIPSNAP2,NDUFS8,SLTM,PGAM1,C22orf39,ARX,TCIRG1,AP2B1,SMAP1,PDZD11,STRAP,BLVRA,THOC1,RTF2,PICALM,VTA1,NAA20,BCCIP,INPP5A,KLHL22,SCAND1,CNOT7,COQ4,EMID1,URM1,KDM4D,LAMP2,EWSR1,COX5A,PSMD3,HBS1L,MAK16,MED4,DPM1,HNRNPK,CCDC107,TLN1,PDAP1,MIEN1,HDAC2,ZMAT5,RNPS1,PRDX3,UCHL3,MRPL41,PSMA4,SRP72,ATP6V1B2,C5,LASP1,EIF4EBP1,HDDC2,VMA21,ESD,PLAA,SNRPC,VPS16,ASCC3,MLST8,ULK3,TOMM22,TCEANC,BUB3,PLRG1,BUD31,VAPA,NUDT16L1,SNX3,TBC1D32,RGS10,GOPC,BECN1,POLR2E,NCLN,STX11,RTL8C,ACP1,ZNF449,TMBIM6,APOO,PSMA3,DKC1,FBXW5,NPDC1,NDUFAF4,ADI1,UBE2V2,CPSF6,AP1B1,MRPL40,PARD3,OAZ1,RNASEH1,SSNA1,SNAP29,PGK1,SPCS3,VAV1,UBE2Z,CLTC,ARF5,CUEDC2,YWHAQ,TYMP,HADHA,LSM8,KISS1R,STMP1,CCT5,GTF2A2,EIF2S1,MRPL54,WDR61,FAM156A,DOHH,CD99,MTMR6,MTIF3,COMMD4,HSD17B10,DNAJA2,HEXB,COA3,TMEM141,C16orf87,CD82,SCAMP2,NET1,PHPT1,TMEM167A,MYH9,MTPN,ABCE1,SUB1,ADGRE2,PGM2,MCF2,CHCHD2,VTI1B,ATP6V1D,SRP54,PRDX5,TPM4,PSMB1,LCMT1,HADHB,AGGF1,GLRX,FIS1,CMKLR1,WASHC3,ANXA4,CPSF3,CCP110,ACYP2,RFXANK,RAB5IF,GOLPH3,SEC11A,CCT2,ASXL1,EXOSC3,MRPS11,FAM32A,TALDO1,FABP5,CEND1,SRP14,SNF8,RASAL3,DPH3,NSMCE1,NDUFA8,C19orf44,PSMD12,NUTF2,YJEFN3,MTDH,EID1,SNRPB,CHERP,POLR2L,HPRT1,OSBPL8,CYB5R3,TOX4,LRRFIP2,AKR1B1,SLCO3A1,UBE2N,TXN2,MORF4L1,PIGT,YWHAZ,PSMD11,SDHD,COX16,CHCHD3,ANXA5,FARS2,STPG4,NDUFS4,ZNHIT1,SERPINB1,MAP2K3,BBS7,METAP2,UBE2R2,DDA1,PTGES2,POLE4,IDH3G,TBC1D7,POLR2K,MPLKIP,COQ9,MPHOSPH6,SAP18,COPS6,MRPL57,TRMT10A,ZDHHC20,TMSB10,MRPL19,DCAF13,SH2D1A,STOML2,SREK1,RDX,TMEM14B,DENR,DCTN3,CAPN1,TRAPPC2L,ACAT2,PSMD9,SLC25A46,RAP1B,GNAI2,CIAPIN1,CDK6,MEA1,COX8A,WDR59,LEO1,LZIC,UBA1,NDFIP2,NDUFB11,TESC,ENO1,ZNF114,RIOK2,ZNF428,LAMTOR1,PHKA2,CDK5,RAB11A,DNPH1,CLTA,MRPL13,GTF3C6,OXSR1,KRCC1,CRADD,CYBA,NDUFB6,MRPL15,NUDT14,TIMM17B,LYPLA1,HAUS5,SFT2D1,MPC1,RALA,PREP,LSM1,TTC7A,CTU2,CCDC85B,RER1,NCKAP1L,NQO1,HDDC3,SMARCC1,WRAP73,ARPC1B,BANF1,MCRS1,GPI,SF3B2,USB1,CDKN2AIPNL,RAB1B,ZCCHC10,PAICS,YIF1A,UBE2D2,COTL1,GK,AKAP13,CMC2,UBA2,UBN1,LYRM4,TMEM258,MAPKAPK3,RDH10,SUCLG1,PPP5C,VDAC1,VPS33B,COA5,ACTR3,RPP30,LAGE3,EMC8,ARL1,ZNF655,WDR73,PARK7,SF3B6,TAOK3,CHURC1,SZRD1,CHMP3,THOC6,PPP2R1A,USP30,ACTB,WASF2,PTPA,NUDT22,CASP3,EMC7,CHCHD1,EMC4,PPP2CA,NOP10,SNRNP40,TIMM23,YIPF3,YIPF5,ARC,PTPN18,SAR1A,PPP6R1,USO1,TMEM160,TMEM147,MAEA,ORAI3,PFKP,PYCARD,ATP5MC1,BLOC1S6,HSPBP1,COPS8,ANXA6,ABI1,SYNCRIP,UQCRFS1,FAAP20,LPXN,ZC3H15,COPB1,OCIAD1,TUFM,SH3RF2,PSMB2,CABCOCO1,UQCRH,AP3M1,PRDM5,SLC1A5,TIMMDC1,VASP,ETF1,MESD,SPCS1,RASGRP4,ITGAV,PPP4C,DELE1,NUDT2,HK1,GNG10,CD2BP2,CTNNBIP1,CORO1A,EMP3,PSMC2,ANAPC16,STAM,POLR2G,TMEM179B,AP1S3,MRPL48,PPIA,MYCBP,DBNL,TRIP4,NOL7,SNX4,SMIM12,PSMC5,SYNJ2,PTPRH,ARAP3,POLD4,PPP1CA,POP4,UBE2E3,CORO1B,HPS1,EIF1B,MCMBP,PPP4R2,SEC23IP,DCPS,SNRPA,TXNDC12,NDUFA5,DMRT2,RNF7,EIF4G2,MOSPD3,MRPS23,GNB2,FERMT3,POP7,BAD,RBM42,CTSC,PDCD5,MAGOH,MPDU1,TRIP6,CAPZB,MRPL11,GYPC,DDX19A,ZDHHC24,RAC1,HACD4,EIF4E,ADAM8,POLDIP2,ZNF511,TIMM50,ECHS1,CALM3,DENND1C,ATP6V1F,INPP5K,NCOA4,TIMM10,GABARAPL2,TERF2IP,UBLCP1,EIF3J,HTATIP2,XRCC5,COX6A1,SRSF9,DYNLL1,NPM3,COQ5,POP5,ANXA2,SFXN1,ARF3,DDX31,HNRNPA2B1,PSMD14,UCP2,RAB18,VRK3,PRPF38A,COPS3,TUBGCP2,AIMP1,GLRX3,SMIM19,ARMC8,TPI1,DEF6,RNASEH2C,TMEM11,CFL1,NUDCD2,ZFAND6,RAB5C,NUP54,ATP1B3,SPG21,AP3D1,ALDOA,JOSD2,ERCC1,FXR1,NUBP2,COPS5,AP2S1,CD72,CCND2,ORMDL2,PIN1,GPX7,ARPP19,TRAF2,STIM1,NDUFAB1,WDR18,SH3GLB1,GPAA1,SCP2,CYC1,SHARPIN,POLD2,NUCB1,TRAPPC1,ZNF202,GRK6,CALM1,FARSA,RAD23A,ENY2,GADD45GIP1,VDR,ILK,ISOC2,RANGRF,MAD2L2,EIF4H,STX8,CCDC32,U2AF2,GOLGA2,ZNF106,CHRAC1,ATXN7L1,PPP1R7,GNG2,BIN1,TM2D2,TERF1,SRI,RWDD3,MRPL4,TMEM9B,ARF6,ICAM3,SPHK1,PPP2R3C,PSMA5,ECHDC1,VPS37C,SRSF2,PSMA6,TATDN1,NDUFB9,HNRNPF,IL12RB1,ERH,EFTUD2,PCBD1,SSBP4,ERP44,STK25,NHP2,ATG4B,SKP1,PCMT1,ARL2BP,RCHY1,ARAP1,AP2A1,HINT1,UBL5,COX6C,AP1M1,CPNE3,E2F4,TMEM208,JAK3,MRPL12,SCNM1,CD320,NDUFA7,RAB11B,PSMC3,HNRNPM,ZMAT2,U2SURP,CRIPT,OGFOD2,HPF1,PFDN1,SRP19,RBPJ,TEN1,REEP5,NDUFAF7,PAFAH1B1,ELOF1,UBE2I,BCL7B,BRK1,PIH1D1,ALDH16A1,SEC13,SLC25A20,GRSF1,NDUFAF3,GSTP1,NDUFV1,THOC7,ATP5F1B,TBCA,DCTPP1,CLINT1,TPM3,VHL,MRPL35,JMJD6,TRMT112,PRR13,ILKAP,CMTM3,RUVBL2,CSNK1A1,MRPL51,GTF2H5,GAPDH,GRAP,LSM7,TIMM13,LSM6,CCS,MZT2A,NDUFB7,AIP,MPC2,CAPZA2,TRAPPC5,CSK,MGST3,ALG3,TIMM44,SPPL2B,PSMD2,ECD,TIPRL,TLCD1,MRPS16,ANXA7,TEX45,CUL2,ZCRB1,DLD,GHITM,MAP2K1,MLF2,WDR83OS,VIM,TPR,PRDX2,ACAA1,LARP7,SHKBP1,TCEA1,CD58,DCTN2,TXNDC11,ZFAND2B,TXNDC15,MYB,ATOX1,TSG101,JPT1,SUMO2,ARL8B,STK16,COPE,LYPLAL1,DDX49,GRB2,DCUN1D5,PRDX6,DPY30,PHB,GPATCH2,PELP1,ARRB2,KEAP1,PSMB6,PRPF40A,COASY,SSR3,PSMD6,TRIP12,RIC8A,CLTB,SEM1,PKM,MDH1,TMC6,ELOC,MRPS14,ORAI2,SNX17,TMEM186,PPM1G,NRBP1,TBCEL,PRMT1,MRPL28,ZNF777,PAFAH1B3,EZR,RPUSD3,JAGN1,MED20,LAMTOR4,COX7A2,RAB4A,POLR2J,DHX9,MIF4GD,KLHL18,SYNGR2,EML2,CCAR2,PDCD10,RTF1,BRD4,MLX,NSL1,NUDT21,STX10,ITGB1BP1,CD79B,DYNC1I2,SIPA1,ARPC3,SLC16A3,ETFA,TMEM101,LSM12,VPS29,FADD,TWF2,ARHGDIA,PAK2,EIPR1,HIKESHI,COPB2,GINS4,ATXN2L,GOLGA7,PA2G4,TTC1,PYCR1,ESYT1,MYL6B,MYL6,ACTG1,IL32,BSG,OTUB1,NABP2,MITD1,TXNDC9,USP4,CLPP,VPS25,RFFL,NLE1,GOLT1B,SSB,METTL5,SRA1,DDX1,YIF1B,NCL,TBCB,HNRNPA3,CAPNS1,COPS9,COX20,COX17,ZDHHC7,DPP3,SNRPG,FAM136A,PSMD1,HNRNPU,SLC12A8,ANAPC15,TIFAB,PDCD6,PPIG,LAT,MUL1,PRELID1,LMAN2,ATP6V0D1,NDUFA12,RNF10,LYPLA2,SLC25A3,CHMP6,TOMM40,GALK1,ANKRD11,RAB40B,ARRDC2,PYM1,FKBP8,KRR1,PDHB,KXD1,PSMD7,SNRPE,REX1BD,ELOB,C1orf131,NDUFB5,BOLA3,MOB1A,MTHFD2,ZNF625,NDUFB10,MRPL27,LRRC59,RHOA,LAMTOR2,ATP5F1D,ATP5MG,CCDC90B,RPL26L1,SMARCA4,BOLA2-SMG1P6,ATP6V0E1,SREK1IP1,CYB5B,PSMD4,PWP1,TAGLN2,DDOST,MPG,ATP5PD,CD74,ATP5MC2,HYAL3,DNAJB11,C1D,C17orf49,UQCRC1,ACADVL,ANKRD13A,HSBP1,GABARAP,ARHGDIB,ETFB,MLEC,KRAS,ZNF8,PTPN11,UBE2J2,DPM3,RBBP4,UBE2G1,UFC1,NAA38,PIAS4,SIRT6,CCR3,ZMPSTE24,ZBTB32,PRKCD,LMAN2L,CNNM4,ANKRD39,UROD,SLFN11,LRRFIP1,PLEKHJ1,STX12,BTBD17,RPA2,RAN,GNG8,SDHB,RXYLT1,TADA3,DENND4B,MYH14,CHFR,TRPV2,PRKAR1A,ACLY,DNAJC7,FDPS,NDUFS7,SNRPF,TRAF4,EBNA1BP2,TET3,MMADHC,CCDC12,CHMP2A,PTRHD1,ECSIT,EXOSC5,PSME3,PPIL3,MEAF6,PPT1,MRPL20,ETHE1,FMNL1,PSMD8,ATP6V0B,FAM104A,RPIA,SNRNP27,EAF2,STAMBP,GNB1,SH3BGRL3,ACAP2,EIF5B,NVL,CNIH4,GPS1,ENSA,UBXN4,RAB10,LCK,COL6A3,HDAC1,EIF3I,LAMTOR5,AP2M1,TERC,TEX264,HDLBP,PFDN2,NDUFS5,SDHC,PARP1,NDUFS2,PRDX1,SLC39A1,GUK1,MAP3K2,SERBP1,EMC3,UCHL5,CAP1,HSPB11,TMCO1,AKR1A1,SUMO1,AK2,NOP58,SUPT7L,SLC4A1AP,GNAI3,ACTR2,CAB39,UBE2Q1,RTN4,ELOVL1,NAXE,DBI,SH2D2A,CALM2,RTCA,KHDRBS1,SRP9,SRM,TMEM50A,TSEN15,COMMD1,MPV17,MTX1,ATRAID,S100A10,COX5B,SFPQ,BZW1,ARPC2,ZBTB8OS,PNKD,CERS2,ANXA9,VAMP8,VAMP5,RNF181,C2orf68,EIF4E2,SSU72,PLEKHB2,CDC42,RASSF5,S100A6,S100A4,CHTOP,SNAPIN,LSM10,EVA1B,RBM45,SVBP,SRSF4,PMVK,POLR3GL,S100A11,NUDC,CSRP1,RBM8A,CAPZA1,TAF12,TMEM183A,RABIF,ARPC5,GNG5,PSMB4,VPS72,ATP5PB,KCNAB2,NFYC,AURKAIP1
 >>> identical protein binding <<<: SSBP1,MIF,B2M,PRPF31,YWHAE,PUF60,IFI27,CSNK2B,ATP5IF1,LDHA,PSME2,XRN2,SMAD2,DRG1,YWHAB,CTPS2,PSMC6,PSMA7,WAS,RTRAF,NUDT5,PFKL,YWHAH,PSAT1,MAPRE1,KIF5B,ITPA,CAT,ARHGAP4,HCFC1,IRAK1,CHMP4B,DESI1,CDC5L,PRPS2,NONO,ASL,PSPH,EWSR1,HNRNPK,PRDX3,ESD,GOPC,ADSL,BECN1,IAH1,SSNA1,YWHAQ,HEXB,MYH9,OXCT1,SUB1,TPM4,FIS1,ANXA4,ACYP2,FABP5,RASAL3,NUTF2,HPRT1,DECR1,YWHAZ,TMEM14B,ZNF114,DNPH1,OXSR1,NUDT14,NQO1,BANF1,PAICS,PPP5C,ACAT1,PARK7,CHMP3,ACTB,PFKP,PYCARD,BLOC1S6,ANXA6,MESD,CORO1A,CORO1B,DCPS,SNRPA,DMRT2,CTSC,EIF3J,ANXA2,HNRNPA2B1,MLYCD,GLRX3,NUP54,ALDOA,TRAF2,STIM1,SH3GLB1,SHARPIN,GOLGA2,BIN1,TERF1,SRI,PCBD1,ADH5,PSMC3,PAFAH1B1,BRK1,SEC13,DCTPP1,JMJD6,RUVBL2,GAPDH,MPC2,CSK,MGST3,VIM,SHKBP1,DCTN2,GRB2,TGFB1,PRDX6,DPY30,ARRB2,KEAP1,PRMT1,PAFAH1B3,EZR,MIF4GD,NUDT21,CD79B,FADD,PAK2,PYCR1,ESYT1,ACTG1,MITD1,USP4,CLPP,NCL,HNRNPU,PDCD6,MUL1,ANKRD11,FKBP8,PSMD4,CD74,C17orf49,ACADVL,HSBP1,KRAS,ZBTB32,TRAF4,PSME3,ETHE1,RPIA,LCK,PARP1,PRDX1,NAXE,DBI,KHDRBS1,SRM,COMMD1,CDC42,RASSF5,S100A4,RBM45
 >>> structural constituent of cytoskeleton <<<: TUBB6,TLN1,ARPC1B,ACTR3,ACTB,VIM,ARPC3,ACTG1,ACTR2,ARPC2,ARPC5

Community 4 has the following top 5 GO terms:
>>> immune response <<<: KIR2DL3,FCAR,KIR2DL1,LILRB2,LST1,HLA-DRA,HLA-DMB,CCL4L2,HLA-DPB1,HLA-DMA,CCL3L1,CCL18,CCL3,CCL4,CCL23,HLA-DRB1,HLA-DPA1,HLA-DQA1,HLA-DQA2,HLA-DQB1,IGKV3-20,DEFB1,RFX1,EDA,CXCL8,THBS1,TNFSF11,CFP,CLEC4E,FTH1,CCL13,CCL7,CCL2,CD36,CCL8,CEBPB,TNFSF13B,CCL22,TNFSF15,CTSL,CEBPG,IRF8,HAMP,TLR2,IL1B,IFNG,PTAFR,IGSF6,CXCL14,C5AR1,CD80,C3,TNFSF14,CXCL2,CXCL9,CXCL10,CXCL11,FCGR1B,PF4V1,CXCL1,CXCL5,HLA-DRB5,MADCAM1,IL1R1,IL13,HRH2,TNFSF13,CCRL2,C1QC,FCGR3B,IL1RN,FCGR1A,SLC11A1,FCGR2B
 >>> chemokine activity <<<: CCL4L2,CCL3L1,CCL18,CCL3,CCL4,CCL23,CXCL8,CCL13,CCL7,CCL2,CCL8,CCL22,CXCL14,CXCL2,CXCL9,CXCL10,CXCL11,PF4V1,CXCL1,CXCL5,CXCL16
 >>> cellular process <<<: SBNO2,CCL4L2,CCL3,RXRA,CCL13,CCL2,MSR1,TPST1,IPO4,CD36,CTNNAL1,DYNC1LI2,HDAC11,MAP1B,PTK2,CXCL14,HACD1,SRD5A3,NAMPT,ADGRG3,GCK,DCLRE1B,SHROOM1,HDAC9,XKR8,MARCO,CAD
 >>> signal transduction <<<: LILRA2,LILRB2,CCL4L2,CCL18,CCL3,CCL4,CCL23,ARHGAP23,FGF11,C1QTNF4,EDA,FAM3B,BEX1,EDN3,IFNA13,TIMP1,IFNA1,SSTR3,DGKQ,IRS2,ARHGAP12,TLR8,NFAM1,BANK1,CXCL8,CD83,GPR146,RASSF4,PGF,GRAPL,ROR2,CCL13,SKAP2,PTGES,CCL7,RIN2,CCL2,TRAF1,TOM1,ASIP,GNG11,CCL8,MPP1,RTKN2,GNA15,LEP,LYN,AKAP12,NRG1,PLAU,GNA11,STARD8,TNFSF13B,OR7C1,CSF2RB,OSGIN2,TAS2R3,SH2B3,NSMAF,CCL22,TNFSF15,ADGRV1,CLEC5A,FGF9,PILRA,ESR2,GPNMB,PLAUR,RIN1,TLR2,PLA2G1B,IL1B,YWHAG,NCK2,PTAFR,IGFBP2,CXCL14,C5AR1,PENK,NPBWR1,FZD2,FAM3D,PLPP3,C3,ZFYVE16,TNFSF14,ANPEP,GAPVD1,RASD1,CXCL9,CXCL10,CXCL11,P2RY1,NAMPT,CSF1R,HTR4,AFDN,SCUBE2,SIGLEC8,FPR2,ADGRG3,SPP1,FPR1,HCAR2,HCAR3,HCAR1,IRAK2,CXCL1,CXCL5,HLA-DRB1,RGS12,RETN,TNFRSF10B,GNB4,TYROBP,ADGRB1,DIRAS1,CXCL16,PDGFA,FFAR2,DOK2,ULK1,VEGFA,ARRB1,LRRK2,MADCAM1,IL1R1,ABL2,IRAK3,IL13,HRH2,GRN,DLG4,TNFSF13,P2RX4,CCRL2,TNFAIP6,RASSF8,ADCY3,IL1RN,RGL1,GPR42,KSR1,FCGR1A,TRIM54,GPBAR1,CSF3R,FCGR2B
 >>> cytokine activity <<<: CCL4L2,CCL3L1,CCL18,CCL3,CCL4,CCL23,C1QTNF4,FAM3B,IFNA16,IFNA2,IFNA13,IFNA7,IFNA17,TIMP1,IFNA1,CXCL8,TNFSF11,IFNA14,IFNA10,CCL13,IFNA8,CCL7,CCL2,CCL8,NRG1,TNFSF13B,CCL22,TNFSF15,IL1B,IFNG,IFNE,CXCL14,INHBB,FAM3D,IFNL1,IFNB1,IFNW1,IFNA21,TNFSF14,CXCL2,CXCL9,CXCL10,CXCL11,NAMPT,SPP1,PF4V1,CXCL1,CXCL5,CXCL16,VEGFA,IL13,GRN,TNFSF13,IL1RN
 ```

```python
TOP_N = 5
FIRST = 5
for community_id, community_row in comm_info.head(FIRST).iterrows():
    print(f'Community {community_id} has the following top {TOP_N} KEGG terms:')
    print('\n'.join(community_row['main_functions_KEGG'].split(';')[:TOP_N]))
    print()
```

```python
 Community 0 has the following top 5 KEGG terms:
>>> MAPK signaling pathway <<<: JUN,MEF2C,GADD45B,MAP3K20,RRAS2,MAP3K1,MAP3K3,MAP3K4,CHP2,MYC,FLNC,TRAF6,TGFB2,DUSP1,DUSP2,DUSP6,PLA2G12A,DUSP10,FGF5,MAPK10,RASGRP2,PLA2G4A,MECOM,MAPK13,MAP2K7,MAP2K5,RASGRF2,MAPK8IP3,MAP3K12,FGFR1,FGF17,FGF18,CACNA1I,FOS,TAOK2,RPS6KA6,PLA2G10,MAP4K4,BRAF,MAP3K11,AKT3,GADD45G,MAPKAPK5,GNA12,GADD45A,MAP3K14,CRK,MAP3K13,CACNB2,SOS2,CACNB3,TNF,MKNK1,PRKCG,RAPGEF2,CACNA1A,CACNA1D,CACNA1F,PPM1A,MAPK7,MAP4K2,MAPK8,MRAS,NR4A1,NTRK2,NTRK1,CHP1
 >>> Primary immunodeficiency <<<: RFXAP,IL7R,CD40LG,AICDA,BTK,RAG1,CD79A,IGLL1,CIITA
 >>> Cytokine-cytokine receptor interaction <<<: CCL26,TNFSF12,TNFSF8,CCL15,FLT4,IL22,CCL17,CCL20,CSF2,CSF2RA,TGFB2,IL23R,CXCL6,CLCF1,CCR10,IFNLR1,CXCR4,LEPR,IL19,LTA,IL17A,IL18,IL15,LTBR,CD40LG,IL18RAP,LIFR,IL24,CNTF,NGFR,ACVR2B,ACVR1,TNF,CXCR2,CCL28,IL7R,TNFRSF12A,CCR6,IL26,EPOR,TNFRSF10A,IFNGR1,TNFRSF10D,AMH,IL2,BMPR1A,CRLF2,CXCL3,TNFRSF21,TSLP,IL4,IL5
 >>> Jak-STAT signaling pathway <<<: PIAS3,SOCS5,SPRY4,IL22,AKT3,MYC,IL24,CSF2,CNTF,CSF2RA,IL23R,BCL2L1,CLCF1,SOS2,IFNLR1,SOCS2,IL7R,LEPR,IL26,IL19,EPOR,IL15,IFNGR1,IL2,CRLF2,SPRY1,SOCS4,PIK3CA,PIK3CB,LIFR,TSLP,IL4,IL5,PIK3R2
 >>> Hematopoietic cell lineage <<<: GP1BA,CD33,TNF,ITGA3,ITGA2B,IL7R,MME,EPOR,CSF2,CSF2RA,ITGA2,CD1D,IL4,IL5

Community 1 has the following top 5 KEGG terms:
>>> Cytokine-cytokine receptor interaction <<<: TNFRSF8,IL21R,CCL5,CD27,CSF1,CX3CL1,XCL1,CCR2,CCL25,TGFBR2,TGFBR1,MET,CXCR3,KITLG,CXCR6,TNFRSF9,TNFSF4,CX3CR1,IL1A,CCR9,TNFRSF1A,CCR5,IL2RG,IL2RB,CCR1,XCL2,BMPR2,BMPR1B,FASLG,IL5RA
 >>> Chemokine signaling pathway <<<: GRK2,CCL5,CX3CR1,GRK7,NRAS,PIK3R5,FGR,CX3CL1,GNGT2,XCL1,CCR9,CCR2,CCL25,PXN,CXCR3,PAK1,CRKL,PLCB2,PRKCB,ADCY7,CXCR6,CCR5,MAPK1,GSK3B,CCR1,XCL2,TIAM2,PIK3CD,GRK4,PIK3CG
 >>> Toll-like receptor signaling pathway <<<: CD86,CCL5,CTSK,PIK3R5,TLR5,MAPK14,TICAM1,MAPK9,MAP2K6,MAPK1,PIK3CD,TLR3,PIK3CG
 >>> NOD-like receptor signaling pathway <<<: CCL5,MAPK9,ERBIN,PSTPIP1,MAPK1,MAPK14
 >>> Cytosolic DNA-sensing pathway <<<: CCL5,POLR3F

Community 2 has the following top 5 KEGG terms:
>>> Endocytosis <<<: PIKFYVE,CHMP2B,ARFGAP3,VPS36,PDCD6IP,HLA-F,HLA-E,RNF41,CHMP1B,USP8,HLA-C,HLA-B,HLA-A,EPS15,FGFR3,RET,VPS37A,ARFGAP2,VPS4B,STAM2,RABEP1
 >>> Cell adhesion molecules (CAMs) <<<: ITGA4,ITGB1,SELPLG,CD6,MAG,HLA-F,HLA-E,CLDN8,HLA-C,HLA-B,ESAM,HLA-A,CLDN15,PTPRC,GLG1
 >>> Antigen processing and presentation <<<: HLA-C,HLA-B,HLA-A,HSPA5,HLA-F,HLA-E,CALR,PDIA3
 >>> Natural killer cell mediated cytotoxicity <<<: ZAP70,PPP3CC,CD247,HLA-E,LCP2,PLCG1,NFATC2,RAF1,HLA-C,HLA-B,HLA-A,SH3BP2,ARAF,IFNAR1,TNFRSF10C,PIK3R1
 >>> Type I diabetes mellitus <<<: HLA-C,HLA-B,HLA-F,HLA-E,HLA-A

Community 3 has the following top 5 KEGG terms:
>>> Pathogenic Escherichia coli infection <<<: WAS,RHOA,ARPC5,ARPC1B,ARPC3,ACTG1,ARPC2,NCL,YWHAQ,ACTB,EZR,TUBB6,ARPC5L,YWHAZ,CDC42
 >>> Hypertrophic cardiomyopathy (HCM) <<<: ITGAV,TPM4,ACTG1,TPM3,LAMA2,TGFB1,ACTB
 >>> Arrhythmogenic right ventricular cardiomyopathy (ARVC) <<<: ACTN4,ITGAV,ACTG1,LAMA2,ACTB
 >>> Dilated cardiomyopathy <<<: ITGAV,TPM4,ACTG1,TPM3,LAMA2,TGFB1,ACTB
 >>> Focal adhesion <<<: PARVG,FLNA,BAD,ACTN4,MAP2K1,ILK,GRB2,ACTG1,MYL12B,ACTB,RAP1B,TLN1,COL6A3,PARVB,VAV1,VASP,ITGAV,CCND2,CDC42,PAK2,RHOA,RAC2,PPP1CA,RAC1,LAMA2

Community 4 has the following top 5 KEGG terms:
>>> Cytokine-cytokine receptor interaction <<<: TNFSF13,CCL3L1,CCL2,CCL3,TNFSF14,CCL4,CCL13,CCL8,CCL7,TNFSF11,CCL18,CCL23,CCL22,CSF1R,IL13RA1,PRLR,IFNE,CSF2RB,CXCL5,CXCL11,CSF3R,CXCL10,INHBB,TNFRSF19,LEP,CXCL9,EDA,CXCL14,TNFSF13B,VEGFA,PDGFA,CCL4L2,PDGFRA,TNFSF15,CXCL1,IFNA2,IFNA1,CXCL8,IL13,IFNA17,IFNA21,IFNA7,IFNA8,IFNA10,IFNA13,IFNA14,IFNA16,CXCL16,IL3RA,IFNG,IFNB1,TNFRSF10B,IL1B,IL1R1,IL6R,CXCL2,PF4V1,IFNL1,IFNW1
 >>> Chemokine signaling pathway <<<: CCL3L1,CCL2,GNB4,CCL3,CCL4,CCL13,CCL8,CCL7,PTK2,CXCL14,SHC4,CCL18,ARRB1,CCL23,CCL22,GNG11,CXCL5,ADCY3,CXCL11,CCL4L2,CXCL10,CXCL1,HCK,CXCL8,WASL,CXCL16,CXCL2,PF4V1,CXCL9,NCF1,LYN
 >>> NOD-like receptor signaling pathway <<<: CXCL1,BIRC3,CCL2,CCL13,CCL8,CCL7,CXCL8,CASP5,NOD2,CARD9,HSP90AA1,NLRC4,IL1B,CXCL2,HSP90AB1
 >>> Porphyrin and chlorophyll metabolism <<<: ALAS1,FTH1
 >>> Natural killer cell mediated cytotoxicity <<<: FCER1G,SHC4,FCGR3B,IFNA2,IFNA1,IFNA17,IFNA21,IFNA7,IFNA8,TYROBP,IFNA10,IFNA13,IFNA14,IFNA16,KIR2DL1,ULBP2,KIR2DL3,IFNG,IFNB1,TNFRSF10B
```

Finally, let us take a look at the top-5 central genes in each community according to Leiden algorithm:

```python
TOP_N = 5

# Looking at central genes in the communities
for community_id, community in comm_info.iterrows():
    print(f"Community {community_id}")
    print(community["all_sorted_genes"].split(';')[:TOP_N])
    print("\n")
```

```python
Community 0
['DUSP1 (score=0.18998955572905765)', ' IL7R (score=0.1638468554537804)', ' KLF2 (score=0.11492169260222583)', ' TSC22D3 (score=0.11486962803428785)', ' FOS (score=0.07371105776759161)']


Community 1
['CCL5 (score=0.33543835396789073)', ' ZNF683 (score=0.15410550633670056)', ' NKG7 (score=0.10819809175100292)', ' TRDV1 (score=0.10617402147363898)', ' CD8A (score=0.06761798903745356)']


Community 2
['HLA-C (score=0.33831786810100817)', ' HLA-A (score=0.2377502730251183)', ' MALAT1 (score=0.2327474127619741)', ' PTPRC (score=0.1756880109160519)', ' MTRNR2L12 (score=0.1676294981043018)']


Community 3
['ACTG1 (score=0.22816198967152546)', ' ACTB (score=0.21815708323032829)', ' TPM4 (score=0.21529789053268528)', ' PFN1 (score=0.2145721303555998)', ' GAPDH (score=0.17829357150363143)']


Community 4
['CCL2 (score=0.20395934613870187)', ' FTH1 (score=0.15592898049565898)', ' FCER1G (score=0.11666055897110021)', ' FTL (score=0.10644153960710624)', ' SOD2 (score=0.07922054758739527)']


Community 5
['STMN1 (score=0.19471379122285665)', ' HMGB2 (score=0.167439647400581)', ' TYMS (score=0.14912295346533552)', ' DUT (score=0.08756385855955125)', ' HIST1H4C (score=0.07516166371720814)']


Community 6
['ISG15 (score=0.2147048155497237)', ' SAT1 (score=0.21059288181479824)', ' MX1 (score=0.18493389376062358)', ' CREM (score=0.15695875738761228)', ' TNFSF10 (score=0.12447934076302508)']


Community 7
['EEF1A1 (score=0.20307412760151444)', ' PABPC1 (score=0.1760801316559486)', ' CCR7 (score=0.10807621301622097)', ' SARAF (score=0.10085970817126913)', ' RPL41 (score=0.08126004313090261)']


Community 8
['FAM183A (score=0.24097835830452044)', ' C9orf24 (score=0.13250096043027276)', ' DYNLRB2 (score=0.1294403892944039)', ' KRT19 (score=0.09845050582661033)', ' CKB (score=0.08643445596960772)']


Community 9
['RPS10 (score=0.15061539372895336)', ' MT-ATP6 (score=0.141663590837218)', ' C1orf56 (score=0.10008749025987342)', ' AC058791.1 (score=0.056722851817929125)', ' JUND (score=0.05522366976072326)']


Community 10
['CTLA4 (score=0.22315591619676342)', ' SRGN (score=0.1451588502269289)', ' TNFRSF18 (score=0.08355109338467887)', ' FOXP3 (score=0.07254618805299592)', ' NEAT1 (score=0.06707238802548939)']


Community 11
['IGHG1 (score=0.1369453044375645)', ' IGLV3-19 (score=0.12242862057103543)', ' MT1E (score=0.11272789817681458)', ' IGLV2-14 (score=0.07567939456484347)', ' IGHV4-34 (score=0.07354661162710698)']


Community 12
['HOXA5 (score=0.0)', ' CIB2 (score=0.0)']
```

## Wordcloud visualization

The best way to depict the communities in gene-regulatory networks is to use [wordcloud](https://github.com/amueller/word_cloud) visualization on top of graph network. We implemented [`scGRN.network_analysis.plot_cloud()`](https://github.com/masyahook/scGRN/blob/99f1ba91303351cd9948016dfaea7ec78f35c30c/scGRN/network_analysis/_plotting.py#L648) function to simplify the process. For example, we can visualize `KEGG` terms in the communities as follows:

```python
def get_central_genes_and_partition(comm_info):
    """Get central genes and partition from community info."""

    # Extract genes sorted by centrality
    all_partition_gene_scores = {
        i: {
            el[: el.find(" ")]: float(el[el.find("=") + 1 : -1])
            for el in comm_info.loc[i, "all_sorted_genes"].split("; ")
        }
        for i in comm_info.index
    }

    # Get partition
    partition = {
        gene: i for i in all_partition_gene_scores for gene in all_partition_gene_scores[i]
    }

    # Pick only top 50 genes for annotation
    central_gene_scores = {
        i: {
            gene_name: gene_score
            for k, (gene_name, gene_score) in enumerate(gene_scores.items())
            if k < limit_anno_until
        }
        for i, gene_scores in all_partition_gene_scores.items()
    }

    # Scale the centrality score
    central_genes = {
        i: dict(
            zip(
                gene_scores.keys(),
                scGRN.util.scale_int(list(gene_scores.values()), 1, 100),
            )
        )
        for i, gene_scores in central_gene_scores.items()
    }

    return central_genes, partition


def extract_genes(gene_str):
    """Extract genes from string in `all_sorted_genes`."""
    return list(map(lambda el: el[: el.find(" ")], gene_str.split("; ")))

limit_anno_until = 50  # use only top 50 genes for annotation

# Load community info
comm_info = scGRN.ana.get_community_info(
    pat=None,  # loading all patients
    cell_type='T_cells',  # cell type
    # data_home=_DATA_HOME  # optional
)

# Load the graph
G = scGRN.ana.get_nx_graph(
    pat=None,  # loading all patients
    cell_type='T_cells',  # cell type
    net_type='all',  # network type
    # data_home=_DATA_HOME  # optional
)

# Select only genes from partition in the graph
G = G.subgraph(list(set(chain(*comm_info["all_sorted_genes"].apply(extract_genes)))))

# Get central genes and partition
central_genes, partition = get_central_genes_and_partition(comm_info)
    
f, ax = plt.subplots(figsize=(25, 45))
cmap = ListedColormap(sns.color_palette(cc.glasbey_bw, n_colors=comm_info.shape[0]).as_hex())

# Getting positions of squeezed graph (i.e. compressed graph with less nodes)
squeezed_G, squeezed_partition = scGRN.ana.squeeze_graph(G, partition)

squeezed_pos = scGRN.ana.netgraph_community_layout(squeezed_G, squeezed_partition, seed=_SEED)
scGRN.ana.plot_cloud(
    G, partition, squeezed_pos, ax=ax, anno_db='KEGG', 
    central_genes=central_genes, limit_anno_until=limit_anno_until, 
    display_func=True
)
nx.draw(
    squeezed_G, 
    squeezed_pos, 
    ax=ax, 
    arrowstyle="->", 
    arrowsize=20, 
    connectionstyle=f'arc3, rad = 0.25', 
    edge_color='k', 
    width=0.4, 
    node_color='k', 
    node_size=50, 
    alpha=0.03
)
nx.draw_networkx_nodes(
    squeezed_G, 
    squeezed_pos, 
    ax=ax, 
    node_size=100, 
    nodelist=list(squeezed_partition.keys()), 
    node_color=list(squeezed_partition.values()), 
    cmap=cmap, 
    alpha=0.01
)

plt.axis('off')
plt.show()
```

![TF-target graph](figs/wordcloud.png)

We see that there are 12 communities with tags such as "Cell adhesion", "Cardiomyopathy", "Receptor signaling", and others.

## Visual advanced examples

Please take a look at the [`GRNBoost2_Community_analysis.ipynb`](../notebooks/GRNBoost2_Community_analysis.ipynb) notebook for detailed examples of running community detection and visualizing the results.
