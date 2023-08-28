# Network analysis

In this folder the user can find different modules and functions targeted for analysis of inferred gene regulatory networks (GRNs). Mainly, we can categorize it into:

- **Exploratory data analysis (EDA)**
    - Loading the data, calculating the general statistics of GRNs
    - The functionalities are implemented in `_data_processing.py`
- **Community analysis**
    - Identifying the substructure in GRN, clustering analysis based on connectivity
    - The functionalities are implemented in `_community.py`
- **Enrichment analysis**
    - Performing functional characterization of obtained gene sets (based on group / community)
    - The functionalities are implemented in `_enrichment.py`
- **Various visualizations**
    - Plotting gene networks with rich information about the structure
    - Plotting communities of genes
    - Plotting enriched terms, functional networks

Each analysis is supported with corresponding notebook in the `notebooks/` folder.

## Community analysis

To simplify computation of different steps, we generated integrated `bash` scripts.

### Community analysis for chosen patient and cell type GRN

To run **community analysis** for a **patient-specific** network, please run:

```bash
cd scGRN_analysis
./scGRN/network_analysis/community_scripts/community_ana_pat.sh <"leiden"|"louvain"> <patient_ID> <cell_type_ID>
```

This will produce a `raw_data_<cell_type>_communities_info.pickle` file with the information summary about the computed communities. 

### Community analysis for aggregated GRN

Moreover, to run **community_analysis** for **aggregated** networks (e.g. cell type networks), please run:

```bash
cd scGRN_analysis
./scGRN/network_analysis/community_scripts/community_ana_agg.sh <"leiden"|"louvain"> <pat_type_ID> <cell_type_ID>
```

This will produce a `raw_data_communities_info.pickle` or `raw_data_<pat_type_ID>_communities_info.pickle` file depending on the input.

### Other details

Eventually, the outputs will be structured, for example, in the following format:

```bash
.
├── <leiden|louvain>_communities
│   ├── raw_data_communities_info.pickle
│   ├── raw_data_C_type_communities_info.pickle
│   ├── raw_data_M_type_communities_info.pickle
│   └── raw_data_S_type_communities_info.pickle
├── nx_graph
│   └── ...
├── pickle
│   └── ...
└── ...
```

During community analysis either [Louvain](https://python-louvain.readthedocs.io/en/latest/index.html) or [Leiden](https://leidenalg.readthedocs.io/en/stable/index.html) algorithms can be used for community detection. Also, for quick generation of input commands please use the `notebooks/Generate_sbatch_commands.ipynb` notebook, which will streamline the batch job scheduling on Slurm cluster.

## Enrichment analysis

The user can perform simple enrichment analysis based on passed gene sets with the help of [`EnrichR`](https://maayanlab.cloud/Enrichr/). As `EnrichR` is a web tool the user should make sure that the instance machine has the internet access when running the pipeline. 

For example, if the data is stored in the following format:

```bash
>>> df.head()
       p_val  avg_log2FC  pct.1  pct.2  p_val_adj cluster   gene
IRF2     0.0    8.346623  0.050  0.005        0.0       M   IRF2
STAT1    0.0    5.006633  0.002  0.000        0.0       M  STAT1
NR3C1    0.0    4.831237  0.092  0.002        0.0       M  NR3C1
IRF1     0.0    4.676670  0.002  0.000        0.0       M   IRF1
REST     0.0    4.024687  0.215  0.016        0.0       M   REST
```

One can run the enrichment analysis with the following command:

```bash
python run_enrichr.py -i <path_to_df> -g gene -c cluster -o <save_path>
```

The result of the command will be stored in `<save_path>` (if `-o` is not specified, the output is saved in the same folder.)

Alternatively, the results of enrichment analysis are stored in the following format:

```bash
>>> df.head()
                                       all_sorted_genes
0     DUSP1 (score=0.18998955572905765); IL7R (score...
1     CCL5 (score=0.33543835396789073); ZNF683 (scor...
2     HLA-C (score=0.33831786810100817); HLA-A (scor...
3     ACTG1 (score=0.22816198967152546); ACTB (score...
4     CCL2 (score=0.20395934613870187); FTH1 (score=...
5     STMN1 (score=0.19471379122285665); HMGB2 (scor...
6     ISG15 (score=0.2147048155497237); SAT1 (score=...
7     EEF1A1 (score=0.20307412760151444); PABPC1 (sc...
8     FAM183A (score=0.24097835830452044); C9orf24 (...
9     RPS10 (score=0.15061539372895336); MT-ATP6 (sc...
10    CTLA4 (score=0.22315591619676342); SRGN (score...
11    IGHG1 (score=0.1369453044375645); IGLV3-19 (sc...
12                  HOXA5 (score=0.0); CIB2 (score=0.0)
```

Then we can run:

```bash
python run_enrichr.py -i <path_to_df> -g all_sorted_genes -o <save_path>
```

## Download NDEx network

We also implemented a helper script to quickly download an [NDEx](https://www.ndexbio.org) network:

```bash
python download_ndex_net.py -i <net_UUID> -f <save_in_folder> -n <save_with_name>
```