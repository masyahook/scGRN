# Graph description

Let's load a gene regulatory network for a certain patient and describe it. We will use the `scGRN` reading data capabilities to do that. We will focus on patient `C51` and cell type `T_cells` as an example.

## Reading the graph and the scRNA-seq data

```python
import scGRN

# Loading the NetworkX graph
nx_graph = scGRN.ana.get_nx_graph(
    pat='C51',  # patient ID
    cell_type='T_cells',  # cell type
    net_type='all',  # could be also 'all' or 'TF'
    # data_home=_DATA_HOME  # optional
)

# Loading the adjacency list of the graph
adj_list = scGRN.ana.get_adj_list(
    pat='C51',  # patient ID
    cell_type='T_cells',  # cell type
    net_type='all',  # could be also 'all' or 'TF'
    # data_home=_DATA_HOME  # optional
)

# Loading the scRNA-seq data
sc_data = scGRN.ana.get_sc_data(
    pat='C51',  # patient ID
    cell_type='T_cells',  # cell type
    # data_home=_DATA_HOME  # optional
)

print('Number of nodes:', nx_graph.number_of_nodes())
print('Number of edges:', nx_graph.number_of_edges())
print('Number of cells:', sc_data.shape[0])
nx_graph
```

```python
> Number of nodes: 11676
> Number of edges: 920050
> Number of cells: 16848
> <networkx.classes.digraph.DiGraph at 0x7f0a96ee93d0>
```

Let us have a look at the adjacency list:

```python
adj_list
```

```python
>
                TF    target    importance  regulation       rho
0           TMSB4X       B2M  4.728644e+01           1  0.611593
1          C1orf56  CDC42SE1  4.410510e+01           1  0.620330
2           MT-CO1    MT-CO3  4.320781e+01           1  0.715067
3           MT-CO1    MT-CO2  4.273359e+01           1  0.682574
4              B2M    TMSB4X  4.261495e+01           1  0.611593
...            ...       ...           ...         ...       ...
920045        DOK2     RPS24  4.892047e-09           1  0.076073
920046     CCNDBP1      RPSA  3.702046e-09           1  0.102179
920047        CMIP     HINT1  3.633490e-09           1  0.034972
920048     SELENOS      FTH1  2.389440e-09           1  0.056361
920049  AC093323.1     RPL36  1.441238e-09           1  0.077012

[920050 rows x 5 columns]
```

Here `TF` is the start of the edge and `target` is the end of the edge (the naming is inherited from `TF-target` networks). `importance` is the importance of the edge or the strength of co-expression, `regulation` is the type of interaction, and `rho` is the correlation between the expression of the TF and the target based on the scRNA-seq data.

## General description

Let us see the top 10 most important edges in the GRN:

```python
adj_list.sort_values(by='importance', ascending=False).head(10)
```

```python
> 
         TF    target  importance  regulation       rho
0    TMSB4X       B2M   47.286441           1  0.611593
1   C1orf56  CDC42SE1   44.105103           1  0.620330
2    MT-CO1    MT-CO3   43.207813           1  0.715067
3    MT-CO1    MT-CO2   42.733589           1  0.682574
4       B2M    TMSB4X   42.614955           1  0.611593
5  CDC42SE1   C1orf56   42.193033           1  0.620330
6     HLA-E      SGO1   41.473206           1  0.365724
7      CUTC  CALCOCO1   38.825905           1  0.399168
8    MT-CO3    MT-CO1   35.109970           1  0.715067
9     AURKB      SGO1   34.115454           1  0.773315
```

We see a lot of mitochondrial genes, which is not surprising as they are highly expressed in the cells. Moreover, we see some immune-related genes, such as `B2M`, `HLA-E`, etc. Interestingly, between genes `C1orf56` and `CDC42SE1` there is a bidirectional edge.

Let us see the genes with the highest out-degree, given the importance is high enough (picking only edges with `importance` > 10):

```python
top_regs = (
    adj_list[adj_list['importance'] > 10]
        .groupby('TF').count()
        .sort_values(by='target', ascending=False)
        .head(10)
)
```

```python
>
        target  importance  regulation  rho
TF                                         
RPL13       29          29          29   29
RPL34       28          28          28   28
RPL41       25          25          25   25
RPS18       24          24          24   24
RPS27       24          24          24   24
RPLP1       23          23          23   23
B2M         23          23          23   23
RPS15A      14          14          14   14
TMSB4X      14          14          14   14
RPL21       12          12          12   12
```

We see that mainly ribosomal genes have the highest out-degree.

## Visualization of gene-gene network

Let us visualize the distribution of importance and correlation values:

```python
import matplotlib.pyplot as plt
import seaborn as sns

pat = 'C51'
cell_type = 'T_cells'
net_type = 'all'
n_bins = 100

f, ax = plt.subplots(1, 2, figsize=(20, 8))
adj_list["importance"].plot(ax=ax[0], bins=n_bins, kind="hist", logy=True)
ax[0].set_title(f"GRNBoost2-inferred gene-gene connection importance")
ax[0].set_xlabel("Gene-gene connection importance")

adj_list["rho"].plot(ax=ax[1], bins=n_bins, density=True, kind="hist", logy=True)
ax[1].set_title(f"Correlation between linked genes")
ax[1].set_xlabel("Spearman correlation")

f.suptitle(
    f"Distributions of GRN edge properties: pat={pat}, net_type={net_type}, cell_type={cell_type}",
    fontsize=20,
)
```

![Importance and correlation distributions](figs/importance_rho_dist.png)

We see that the distribution of importance is highly skewed, which is expected (gene-gene regulation is compartmentalized, genes are regulated by only a few genes). The correlation distribution is also skewed, but not as much as the importance distribution.

Now let us visualize the network itself. We will use the [`scGRN.ana.draw_graph()`](https://github.com/masyahook/scGRN/blob/1ca37b068458b359821b28054d2d304e6b29c165/scGRN/network_analysis/_plotting.py#L240) function:

```python
import networkx as nx

pat = 'C51'
cell_type = 'T_cells'
net_type = 'all'
alpha = 0.5
top_n = 100

out_link_dict = {
    (st, end): info["importance"] for st, end, info in nx_graph.edges(data=True)
}
out_dict_sorted = dict(
    sorted(out_link_dict.items(), key=lambda x: x[1], reverse=True)[:top_n]
)
edges = list(out_dict_sorted.keys())
G = nx_graph.edge_subgraph(edges)

pos = nx.nx_agraph.pygraphviz_layout(
    G, prog="dot"
)  # nx.nx.spring_layout(G, seed=_SEED, k=0.2, iterations=20)
cmap = sns.light_palette("seagreen", as_cmap=True)

f, ax = plt.subplots(figsize=(30, 10))

scGRN.ana.draw_graph(
    G,
    pos=pos,
    ax=ax,
    TF_names=None,  # we will not highlight any special genes, e.g. TFs
    alpha=alpha,
    label_edges=False,
    if_alpha_edges=False,
    cmap=cmap,
)

f.suptitle(
    f"Top {top_n} links with highest importances: pat={pat}, cell_type={cell_type}, net_type={net_type}\n",
    fontsize=20,
)

plt.tight_layout()
```

![Graph with highest importances](figs/top_link_graph.png)

Great! We can see here a cluster of ribosomal genes, mitochondrial genes, and some immune-related genes.

## Visualization of TF-target network

Now in a similar manner, let us load the TF-target network and visualize the TF-target connections on the graphs:

```python
import scGRN

# Loading the NetworkX graph
ctx_nx_graph = scGRN.ana.get_nx_graph(
    pat='C51',  # patient ID
    cell_type='T_cells',  # cell type
    net_type='ctx',
    # data_home=_DATA_HOME  # optional
)

# Loading the adjacency list of the graph
ctx_adj_list = scGRN.ana.get_adj_list(
    pat='C51',  # patient ID
    cell_type='T_cells',  # cell type
    net_type='ctx',
    # data_home=_DATA_HOME  # optional
)

print('Number of nodes in enriched TF-target graph:', ctx_nx_graph.number_of_nodes())
print('Number of edges in enriched TF-target graph:', ctx_nx_graph.number_of_edges())
```

```python
> Number of nodes in enriched TF-target graph: 8639
> Number of edges in enriched TF-target graph: 26205
```

```python
# Transcription factors that we will show
TFs = "STAT1, ETS1, JUND, REL, IRF7, ELF1".split(", ")

pat = 'C51'
cell_type = 'T_cells'
net_type = 'ctx'
alpha = 0.5
top_n = 10

# Making sure that all TFs are in the graph
TFs_to_plot = [tf for tf in TFs if tf in ctx_nx_graph.nodes()]

# Subsetting the initial graph to only chosen TFs and targets
edges = []
for TF in TFs_to_plot:
    out_link_dict = {
        end: info["importance"] for st, end, info in ctx_nx_graph.edges(TF, data=True)
    }
    out_dict_sorted = dict(
        sorted(out_link_dict.items(), key=lambda x: x[1], reverse=True)[:top_n]
    )
    edges += [(TF, target) for target in out_dict_sorted.keys()]
G = ctx_nx_graph.edge_subgraph(edges)

pos = nx.nx.spring_layout(
    G, seed=_SEED, k=0.2, iterations=20
)  # nx.nx_agraph.pygraphviz_layout(G, prog="dot")
cmap = sns.light_palette("seagreen", as_cmap=True)

f, ax = plt.subplots(figsize=(30, 15))

scGRN.ana.draw_graph(
    G,
    pos=pos,
    ax=ax,
    TF_names=TFs_to_plot,
    label_edges=False,
    alpha=alpha,
    if_alpha_edges=False,
    cmap=cmap,
)

ax.set_title(
    f"Graph with TFs (in blue) and their top {top_n} targets (in yellow) with highest importances: "
    f"pat={pat}, cell_type={cell_type}, net_type={net_type}\n"
    f'TFs: {", ".join(TFs_to_plot)}',
    fontsize=20,
)
```

![TF-target graph](figs/tf_target_graph.png)

We can see how TFs (in blue) are connected to their targets (in yellow).

## Calculation of graph properties across all samples

To show aggregated information about the graph properties across all samples, we implemented [`scGRN.ana.get_graph_stats()`](https://github.com/masyahook/scGRN/blob/1ca37b068458b359821b28054d2d304e6b29c165/scGRN/network_analysis/_data_processing.py#L1053) function which computes the following properties in parallel:

- `num_nodes`: The number of nodes in the graph.
- `num_edges`: The number of edges in the graph.
- `num_cells`: The number of cells in the graph.
- `average_degree`: The average degree of the nodes in the graph.
- `radius`: The minimum eccentricity among all nodes in the graph.
- `average_path_length`: The average shortest path length in the graph.
- `diameter`: The maximum eccentricity of any node in the graph.
- `importances`: The importance scores of the nodes in the graph.
- `rhos`: The correlation coefficients of the nodes in the graph.
- `median_importance`: The median importance score of the nodes in the graph.
- `median_rho`: The median correlation coefficient of the nodes in the graph.
- `STD_importance`: The standard deviation of the importance scores of the nodes in the graph.
- `STD_rho`: The standard deviation of the correlation coefficients of the nodes in the graph.
- `tfs`: The transcription factors in the graph (out-degree nodes).
- `num_tfs`: The number of transcription factors in the graph (out-degree nodes).

```python
%%time

_FMETA = (
    f"{_PROJ_HOME}/Data_home/data/GSE145926_RAW/metadata.tsv"
)
_DATA_HOME = f"{_PROJ_HOME}/Data_home/res/covid_19"
_N_JOBS = 48
_Q_THRESH = 0.95

# Loading full metadata for all patients
full_meta = scGRN.ana.get_meta(
    data_home=_DATA_HOME,  # file path to the data home folder
    meta_file=_FMETA  # file path to the patient metadata file
)

graph_stats = scGRN.ana.get_graph_stats(
    meta=full_meta, n_jobs=_N_JOBS, filtered=_Q_THRESH
)
```

```python
CPU times: user 16.5 s, sys: 4.22 s, total: 20.8 s
Wall time: 7min 19s
```

```python
graph_stats.keys()
```

```python
> dict_keys(['all', 'ctx'])
```

```python
graph_stats['all'].keys()
```

```python
> dict_keys(['num_nodes', 'num_edges', 'num_cells', 'average_degree', 'radius', 'average_path_length', 'diameter', 'importances', 'rhos', 'median_importance', 'median_rho', 'STD_importance', 'STD_rho', 'tfs', 'num_tfs'])
```

Now we can show the graph properties:

```python
# Plotting graph properties against number of cells for all GRNs
scGRN.ana.graph_stats_vs_num_cells(graph_stats)
```

![Graph properties vs number of cells](figs/graph_stats_vs_num_cells.png)

We could conclude that the number of network nodes and edges positively correlates with the input size as GRNBoost2 manages to capture bigger signal from larger datasets (left panel). However, the network diameter fluctuates and is intrinsic to each input data (middle panel). For all GRNs every vertex could be reached in less than 10 steps for both types of networks which indicates their high connectivity (right panel).

Finally let us see the edge properties across all samples:

```python
# Plotting graph edge properties
scGRN.ana.graph_edge_stats_vs_num_cells(graph_stats)
```

![Graph edge properties vs number of cells](figs/graph_edge_properties.png)

If we look at the edge attributes, we could observe that GRNBoost2-inferred importances do not co-align monotonically with Spearman correlation hinting at complex nature of genetic interactions (left panel). Interestingly, networks contain low amount of edges with negative regulation (i.e. negative Spearman correlation). The importance median increases with the input data size (middle panel) which again indicates the gain in detected signal in bigger datasets. Meanwhile, the median Spearman correlation decreases with data size due to the fact that we use quantile filtering to obtain final networks (for larger networks the quantile shifts towards smaller values and low-correlation links are captured) (right panel).
