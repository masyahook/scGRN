
"""Community analysis."""

# General
import os
import warnings

# Tools/utils
import multiprocessing
from tqdm.notebook import tqdm
from tqdm import tqdm as tqdm_cli
from itertools import chain  # for aggregate functions

# Data management
import math
import numpy as np
import pandas as pd
import networkx as nx
import igraph as ig
import leidenalg as la
from community import community_louvain

# Visualization
import matplotlib.pyplot as plt
import seaborn as sns
import colorcet as cc
from matplotlib.colors import ListedColormap
from wordcloud import WordCloud, STOPWORDS
from termcolor import colored  # colored text output

from typing import Dict, Tuple, List


stopwords = STOPWORDS.union({
    'regulation', 'activity', 'positive', 'negative', 'catabolic', 'process', 'protein', 'complex', 
    'binding', 'response', 'gene', 'genes', 'encoding', 'defining', 'GeneID', 'regulated',
})


def netgraph_community_layout(
    G: nx.DiGraph, 
    node_to_community: Dict[str, int], 
    community_scale: float = 1., 
    node_scale: float = 2., 
    seed: int = 42
) -> Dict[str, Tuple[float, float]]:
    """
    Compute the node positions of the partitioned graph using `netgraph` layout. It was inspired by
    the following example:

    https://netgraph.readthedocs.io/en/stable/sphinx_gallery_output/plot_10_community_layout.html#sphx-glr-sphinx-gallery-output-plot-10-community-layout-py

    :param G: NetworkX graph
    :param node_to_community: A dictionary where key is the node name and value is the community number
    :param community_scale: A scale factor for positioning communities
    :param node_scale: A scale factor for positioning nodes
    :param seed: A random seed

    :returns: A dictionary where each 
    """

    # assert that there multiple communities in the graph; otherwise abort
    communities = set(node_to_community.values())
    if len(communities) < 2:
        warnings.warn("Graph contains a single community. Unable to compute a community layout. Computing spring layout instead.")
        return nx.spring_layout(G, weight='importance')

    community_size = _get_community_sizes(node_to_community)
    community_centroids = _get_community_positions(G, node_to_community, community_scale, seed=seed)
    relative_node_positions = _get_node_positions(G, node_to_community, node_scale, seed=seed)

    # combine positions
    node_positions = dict()
    for node, community in node_to_community.items():
        xy = community_centroids[node]
        delta = relative_node_positions[node] * community_size[community]
        node_positions[node] = xy + delta

    return node_positions


def _get_community_sizes(node_to_community: Dict[str, int]) -> Dict[int, float]:
    """
    Compute the area of the canvas reserved for each community. Inspired by:

    https://github.com/paulbrodersen/netgraph/blob/e555201737f43f8d3e1027d9e1f8c8c8ab154361/netgraph/_node_layout.py#L1565C9-L1565C9

    :param node_to_community: A dictionary where key is the node name and value is the community number

    :returns: A dictionary where each key is the community number and value is the community size    
    """
    
    def _invert_dict(mydict):
        """Invert a dictionary such that values map to keys."""
        inverse = dict()
        for key, value in mydict.items():
            inverse.setdefault(value, set()).add(key)
        return inverse
    
    scale = (1, 1)
    
    total_nodes = len(node_to_community)
    max_radius = np.linalg.norm(scale) / 2
    scalar = max_radius / total_nodes
    community_to_nodes = _invert_dict(node_to_community)
    community_size = {community : len(nodes) * scalar for community, nodes in community_to_nodes.items()}
    
    return community_size


def _get_community_positions(
    G: nx.DiGraph, 
    node_to_community: Dict[str, int], 
    community_scale: float = 1., 
    seed: int = 42, 
    simple: bool = True
) -> Dict[str, Tuple[float, float]]:
    """
    Compute a centroid position for each community. Inspired by:

    https://github.com/paulbrodersen/netgraph/blob/e555201737f43f8d3e1027d9e1f8c8c8ab154361/netgraph/_node_layout.py#L1575



    """
    
    # create a weighted graph, in which each node corresponds to a community,
    # and each edge weight to the number of edges between communities
    between_community_edges = _find_between_community_edges(G, node_to_community)

    communities = set(node_to_community.values())
    hypergraph = nx.DiGraph()
    hypergraph.add_nodes_from(communities)
    
    if not simple:  
        for (ci, cj), edges in between_community_edges.items():
            hypergraph.add_edge(ci, cj, weight=len(edges))

    # find layout for communities
    pos_communities = nx.spring_layout(hypergraph, scale=community_scale, seed=seed)

    # set node positions to position of community
    pos = dict()
    for node, community in node_to_community.items():
        pos[node] = pos_communities[community]

    return pos


def _find_between_community_edges(
    G: nx.DiGraph, 
    node_to_community: Dict[str, int],
    fixed_community: Union[int, None] = None
) -> Dict[Tuple[int, int], List[Tuple[str, str]]]:
    """
    Calculate all inter-community edges in a graph.

    :param G: NetworkX graph
    :param node_to_community: A dictionary where key is the node name and value is the community number
    :param fixed_community: The number of the community, if specified only edges with this community
        will be extracted

    :returns: A dictionary with the following format:
        {
            (community_i, community_j): [
                (node_1_from_community_i, node_1_from_community_j),
                (node_2_from_community_i, node_2_from_community_j),
                ...
            ],
            ...
        }
    """

    edges = dict()

    for (ni, nj) in G.edges():
        ci = node_to_community[ni]
        cj = node_to_community[nj]
        
        if fixed_community is not None:
            if fixed_community != ci and fixed_community != cj:
                continue

        if ci != cj:
            try:
                edges[(ci, cj)] += [(ni, nj)]
            except KeyError:
                edges[(ci, cj)] = [(ni, nj)]

    return edges

def _get_node_positions(G, node_to_community, node_scale, seed):
    """
    Positions nodes within communities.
    """
    communities = dict()
    for node, community in node_to_community.items():
        try:
            communities[community] += [node]
        except KeyError:
            communities[community] = [node]

    pos = dict()
    for ci, nodes in communities.items():
        subgraph = G.subgraph(nodes)
        pos_subgraph = nx.spring_layout(subgraph, weight='importance', scale=node_scale, seed=seed)
        pos.update(pos_subgraph)

    return pos

def squeeze_graph(G, partition, approximate_size=4000):
    """
    Squeeze graph by picking only top nodes (according to number of connections) in each partition. This
    step is needed to speed up the networkx visualization and show only the general POV on the graph.
    """
    
    #### STEP 1 - filtering nodes
    
    # Getting the number of partitions
    num_partitions = len(set(partition.values()))
    
    # Getting partition parameters
    partition_sizes = {i: len([1 for node, k in partition.items() if k == i]) for i in range(num_partitions)}
    min_partition_size = min(partition_sizes.values())
    
    # Normalizing partition size: divide each partition size by the minimal partition size
    normalized_partition_size = {i: (size // min_partition_size) for i, size in partition_sizes.items()}
    
    # Getting scale factor - to get approximately size of the graph close to approximate_size
    scale_factor = math.ceil(approximate_size / sum(normalized_partition_size.values()))
    squeezed_partition = {i: (size * scale_factor) for i, size in normalized_partition_size.items()}
    
    top_nodes = []
    for i, num_nodes in squeezed_partition.items():
        # Getting partition graph
        partition_i = G.subgraph([node for node, k in partition.items() if k == i])
        
        # Finding inter-community edges
        intercommunity_edges = _find_between_community_edges(G, partition, i)
        
        # Calculating node importance according to number of inter-community edges
        node_importance = {}
        for (part_1, part_2), edges in intercommunity_edges.items():
            for node_1, node_2 in edges:
                curr_node = node_1 if part_1 == i else node_2
                if curr_node in node_importance:
                    node_importance[curr_node] += 1
                else:
                    node_importance[curr_node] = 1
                    
        # Getting top nodes in the partition according to maximum number of inter-community edge (node_importance)
        top_nodes += list(dict(sorted(node_importance.items(), key=lambda x: x[1], reverse=True)[:squeezed_partition[i]]).keys())
    
    filtered_partition = {node: i for node, i in partition.items() if node in top_nodes}
    filtered_G = G.subgraph(top_nodes)
    
    #### STEP 2 - filtering edges
    
    # Setting up the size of the squeezed graph (number of edges)
    keep_num_edges = 20000
    edges_to_keep = \
        list(
            dict(
                sorted(
                    {
                        (st, end): data['importance'] for st, end, data in filtered_G.edges(data=True)
                    }.items(), key=lambda x: x[1], reverse=True)[:keep_num_edges]
            ).keys()
    )
    squeezed_G = filtered_G.edge_subgraph(edges_to_keep)
    squeezed_partition = {node: i for node, i in filtered_partition.items() if node in squeezed_G.nodes()}
    
    return squeezed_G, squeezed_partition


def get_elipsis_mask():
        h, w = 600, 800
        center = (int(w/2), int(h/2))
        radius_x = w // 2
        radius_y = h // 2

        Y, X = np.ogrid[:h, :w]
        mask = ((X - center[0])**2/radius_x**2 + (Y - center[1])**2/radius_y**2 >= 1)*255

        return mask


def plot_cloud(G, partition, squeezed_pos, ax, anno_db, filter_genes=True, 
               limit_anno_until=50, display_func=False, if_betweenness=True, 
               k=3000): 
    """
    Plot word cloud that indicates the function(s) of each gene cluster.
    """
    
    # Loading the gene functional annotation
    gene_func = load_gene_func_db(anno_db, reload=False, as_series=True)
    
    # Reversing partition dict -> {group_1: [gene_1, gene_2, ...], group_2: [gene_3, gene_4, ...], ...}
    partition_genes_ = {}
    for gene, i in partition.items():
        if i not in partition_genes_.keys():
            partition_genes_[i] = [gene]
        else:
            partition_genes_[i] += [gene]
           
    # If display gene function in the word clouds
    if display_func:
            
        # Whether to filter the genes on which we compute the word cloud (most important genes)
        if filter_genes:
            compute_centrality = nx.betweenness_centrality if if_betweenness else nx.closeness_centrality
            distance_metric = {'weight': 'distance'} if if_betweenness else {'distance': 'distance'}
            partition_genes = {}
            t = tqdm(partition_genes_.items())
            for i, genes in t:
                t.set_description(f'Processing cluster {i}, size={G.subgraph(genes).order()}')
                top_len = min(limit_anno_until, len(genes))
                top_gene_scores = dict(
                    sorted(
                        compute_centrality(
                            G.subgraph(genes), k=min(G.subgraph(genes).order(), k), **distance_metric
                        ).items(), 
                        key=lambda x: x[1], reverse=True
                    )[:top_len]
                )
                # Renormalizing centrality scores between 1 and 100, and rounding them to use later when 
                # displaying wordclouds (higher score - higher "frequency" or word size)
                norm_top_gene_scores = dict(
                    zip(
                        top_gene_scores.keys(), list(map(lambda x: int(x), scale(list(top_gene_scores.values()), 1, 100)))
                    )
                )
                partition_genes[i] = norm_top_gene_scores
            print('Filtered genes for generating the function word cloud..')
        else:
            partition_genes = {{gene_: 1 for gene_ in gene_list} for i, gene_list in partition_genes_.items()}
        
        # Computing functional annotation for each cluster as a concatenated list of annotations
        # Each annotation is weighted by its duplication gene_score times (e.g. a gene has score = 2 -> 
        # the functional annotation is duplicated and have bigger font in WordCloud)
        partition_funcs = {
            i: ' '.join(
                chain.from_iterable([
                   gene_func[gene_func.index == gene].to_list()*gene_score 
                        for gene, gene_score in gene_score_list.items()
            ])) for i, gene_score_list in partition_genes.items()
        }

        # Generating word counts from aggregated gene annotation texts -> obtaining main (most frequent) function tokens
        word_counts = {i: WordCloud(max_words=30, min_font_size=15, stopwords=stopwords).process_text(text) for i, text in partition_funcs.items()}
        word_counts = {
            i: (freqs if freqs else {'no found function': 1}) for i, freqs in word_counts.items()
        }  # dealing with no word case
        wordclouds = {
            i: WordCloud(
                max_words=30, min_font_size=15, stopwords=stopwords, background_color='white', mask=get_elipsis_mask()
            ).generate_from_frequencies(freqs) for i, freqs in word_counts.items()
        }
        
    # Display main genes in decreasing order of importance (top `top_len` genes)
    else:
        
        compute_centrality = nx.betweenness_centrality if if_betweenness else nx.closeness_centrality
        distance_metric = {'weight': 'distance'} if if_betweenness else {'distance': 'distance'}
        partition_genes = {}
        t = tqdm(partition_genes_.items())
        for i, genes in t:
            t.set_description(f'Processing cluster {i}, size={G.subgraph(genes).order()}')
            top_len = min(limit_anno_until, len(genes))
            top_gene_scores = dict(
                sorted(
                    compute_centrality(
                        G.subgraph(genes), k=min(G.subgraph(genes).order(), k), **distance_metric
                    ).items(), 
                    key=lambda x: x[1], reverse=True
                )[:top_len]
            )
            # Renormalizing centrality scores between 1 and 100, and rounding them to use later when 
            # displaying wordclouds (higher score - higher "frequency" or word size)
            norm_top_gene_scores = dict(
                zip(
                    top_gene_scores.keys(), list(map(lambda x: int(x), scale(list(top_gene_scores.values()), 1, 100)))
                )
            )
            partition_genes[i] = norm_top_gene_scores
        print('Obtained top genes for generating the gene word cloud..')
        
        wordclouds = {
            i: WordCloud(
                max_words=30, min_font_size=15, background_color='white', mask=get_elipsis_mask()
            ).generate_from_frequencies(gene_score_dict) for i, gene_score_dict in partition_genes.items()
        }
        
    
    # Plotting
    partition_coords = {}
    for gene, coords in squeezed_pos.items():
        if partition[gene] not in partition_coords:
            partition_coords[partition[gene]] = [coords]
        else:
            partition_coords[partition[gene]] += [coords]
    for i, coords in partition_coords.items():
        x, y = zip(*coords)
        min_x, max_x = min(x), max(x)
        min_y, max_y = min(y), max(y)
        ax.imshow(wordclouds[i], interpolation='bilinear', extent=[min_x, max_x, min_y, max_y])
    
    return ax

           
def process_communities(data, pat=None, algo='leiden', filter_quantile=0.95, if_betweenness=True, 
                        limit_anno_until=50, k=5000, save_top_intercommunity_links_until=20, 
                        other_functions_until=20, save_top_new_found_cluster_links=20, seed=42):
    """
    Process graph by finding its communities, annotate its communities, and save everything into .tsv format.
    """
    
    from joblib import Parallel, delayed
    
    def highlight_TFs(word, font_size, position, orientation, font_path, random_state):
        TF_color = (255, 0, 0)  # red
        if word in lambert_TF_names or word in dorothea_TF_names:
            return TF_color
        else:
            r, g, b, alpha = plt.get_cmap('viridis')(font_size / 120)
            return (int(r * 255), int(g * 255), int(b * 255))
    
    print('\nPerforming community analysis..\n\n')
    
    # Setting pathways to files
    _PROJ_PATH = '/gpfs/projects/bsc08/bsc08890'
    _FMETA = os.path.join(_PROJ_PATH, 'data/GSE145926_RAW/metadata.tsv')
    _DATA_HOME = os.path.join(_PROJ_PATH, 'res/covid_19')

    # Loading sample meta data, reordering patients
    full_meta = pd.read_csv(_FMETA, sep='\t', index_col=0)
    
    # Prepare everything to save the figs and dataframe
    if data == 'all_data':
        data = 'raw_data'
    elif 'raw_data_' not in data:
        data = f'raw_data_{data}'
    else:
        pass
    
    if pat is None or pat == 'all_data':
        
        # Cell-type aggregated data
        data_folder = 'all_data' if data == 'raw_data' else data.replace('raw_data_', '')
        
        figs_as = os.path.join(_DATA_HOME, 'cell_types', data_folder, 'figs', 'grnboost2', f'raw_data')
        
        data_to = os.path.join(_DATA_HOME, 'cell_types', data_folder, 'data', 'grnboost2', f'{algo}_communities')
        data_as = os.path.join(data_to, f'raw_data_communities_info.pickle')
        
    elif pat in ['C', 'M', 'S']:
        
        # Patient-type aggregated data
        data_folder = 'all_data' if data == 'raw_data' else data.replace('raw_data_', '')
        
        figs_as = os.path.join(_DATA_HOME, 'cell_types', data_folder, 'figs', 'grnboost2', 
                               f'raw_data_{pat}_type')
        
        data_to = os.path.join(_DATA_HOME, 'cell_types', data_folder, 'data', 'grnboost2', f'{algo}_communities')
        data_as = os.path.join(data_to, f'raw_data_{pat}_type_communities_info.pickle')
        
    else:
        
        # Loading patient-specific data
        figs_as = os.path.join(_DATA_HOME, pat, 'figs', 'grnboost2', f'{data}')
        
        data_to = os.path.join(_DATA_HOME, pat, 'data', 'grnboost2', f'{algo}_communities')
        data_as = os.path.join(data_to, f'{data}_communities_info.pickle')
    
    os.makedirs(data_to, exist_ok=True)
    os.makedirs(os.path.dirname(figs_as), exist_ok=True)
    
    # Loading lists of TFs from Lambert 2018 and DoRothEA, in the latter case we will keep only confident regulons
    lambert_TF_names = pd.read_csv(os.path.join(_PROJ_PATH, 'data/TF_lists/lambert2018.txt'), header=None)[0].to_list()
    dorothea_TF_names = list(
        pd.read_csv(os.path.join(_PROJ_PATH, 'data/TF_lists/dorothea_regulons.tsv'), sep='\t') \
            .loc[lambda x: x['confidence'].isin(['A', 'B', 'C'])]['tf'].unique()
    )
    
    # Loading the graph
    G = get_nx_graph(data=data, data_type='all', pat=pat, get_filtered=filter_quantile)
    print(f"Loaded the graph: {colored('pat', 'green')}='{colored(pat, 'red')}', "
          f"{colored('data', 'green')}='{colored(data, 'red')}', "
          f"{colored('data_type', 'green')}='{colored('all', 'red')}'\n")
    
    
    ###### FINDING COMMUNITIES IN THE GRAPH #######
    
    print('Finding communities in the graph..')
    
    if algo == 'louvain':
        partition = community_louvain.best_partition(G.to_undirected(), weight='importance', random_state=seed)
    elif algo == 'leiden':
        G_igraph = ig.Graph.from_networkx(G.to_undirected())
        la_partition = la.find_partition(G_igraph, la.ModularityVertexPartition, weights='importance', seed=seed)
        partition = {G_igraph.vs[node]['_nx_name']: i for i, cluster_nodes in enumerate(la_partition) for node in cluster_nodes}
    else:
        raise NotImplementedError
        
    num_partitions = len(set(partition.values()))
    print(f'Number of partitions using {algo} algorithm: {colored(num_partitions, "cyan")}\n')
    
    
    ###### FINDING HIGH-CENTRALITY GENES IN THE WHOLE GRAPH
    
    print('Finding high-centrality genes in the whole graph..')
    
    num_workers = max(multiprocessing.cpu_count() // 2, 1)
    whole_G_central_genes = dict(
        sorted(betweenness_centrality_parallel(G, processes=num_workers).items(), key=lambda x: x[1], reverse=True)[:limit_anno_until]
    )
    print(f'Computed the {"betweenness" if if_betweenness else "closeness"} centrality for all genes in the graph\n')
    
    ###### FINDING HIGH-CENTRALITY GENES AND CORRESPONDING FUNCTIONS IN EACH COMMUNITY USING GO ANNOTATION ######
    
    print('Finding high-centrality genes/functions in each cluster..')
    
    # Loading the gene functional annotation
    anno_db_tags = ['GO', 'KEGG', 'immunological', 'hallmark']
    gene_func_dbs = {tag: load_gene_func_db(tag, as_series=True) for tag in anno_db_tags}
    
    # Reversing partition dict -> {group_1: [gene_1, gene_2, ...], group_2: [gene_3, gene_4, ...], ...}
    partition_genes_ = {}
    for gene, i in partition.items():
        if i not in partition_genes_.keys():
            partition_genes_[i] = [gene]
        else:
            partition_genes_[i] += [gene]

    # Whether to filter the genes on which we compute the word cloud (most important genes)
    compute_centrality = nx.betweenness_centrality if if_betweenness else nx.closeness_centrality
    distance_metric = {'weight': 'distance'} if if_betweenness else {'distance': 'distance'}
    all_partition_genes = {}
    norm_partition_genes = {}
    t = tqdm_cli(partition_genes_.items(), ascii=True)
    for i, genes in t:
        t.set_description(f'Processing cluster {i}, size={G.subgraph(genes).order()}')
        gene_scores = dict(
            sorted(
                compute_centrality(
                    G.subgraph(genes), k=min(G.subgraph(genes).order(), k), normalized=True, **distance_metric
                ).items(), 
                key=lambda x: x[1], reverse=True
            )
        )
        all_partition_genes[i] = gene_scores
        central_gene_scores = {gene: gene_scores[gene] for k, gene in enumerate(gene_scores.keys()) if k < limit_anno_until}
        
        # Renormalizing centrality scores between 1 and 100, and rounding them to use later when 
        # displaying wordclouds (higher score - higher "frequency" or word size)
        norm_partition_genes[i] = dict(
            zip(
                central_gene_scores.keys(), 
                list(map(lambda x: int(x), scale(list(central_gene_scores.values()), 1, 100)))
            )
        )
    print('Computed centrality scores for each gene in each community\n')
    
    print('Finding functional annotations for each cluster..')
    
    # Computing functional annotation for each cluster as a concatenated list of annotations
    # Each annotation is weighted by its duplication gene_score times (e.g. a gene has score = 2 -> 
    # the functional annotation is duplicated and have bigger font in WordCloud)
    # We also do it for different functional annotations like GO, KEGG, Hallmark, etc..
    partition_funcs = {
        tag: 
            {
                i: ' '.join(
                    chain.from_iterable([
                       gene_func[gene_func.index == gene].to_list()*gene_score 
                            for gene, gene_score in gene_score_list.items()
                ])) for i, gene_score_list in norm_partition_genes.items()
            } for tag, gene_func in gene_func_dbs.items()
    }
    
    print('Computed functional annotations for each cluster\n')

    
    ###### PLOTTING GENE AND FUNC COMMUNITY CLOUDS ######
    
    print('Plotting clusters..')
    
    # Getting positions of squeezed graph - we do not plot every gene on the figure
    squeezed_G, squeezed_partition = squeeze_graph(G, partition)
    print('Computed a squeezed graph representation..')
    
    squeezed_pos = netgraph_community_layout(squeezed_G, squeezed_partition, seed=seed)  # nx.nx_agraph.pygraphviz_layout(G.to_undirected(), prog="sfdp")  # nx.nx.spring_layout(G, seed=seed, k=0.2, iterations=20)
    partition_coords = {}
    for gene, coords in squeezed_pos.items():
        if partition[gene] not in partition_coords:
            partition_coords[partition[gene]] = [coords]
        else:
            partition_coords[partition[gene]] += [coords]
    print('Computed node positions of the squeezed graph representation..')
    
    cmap = ListedColormap(sns.color_palette(cc.glasbey_bw, n_colors=num_partitions).as_hex())
    
    for plot_type in ['genes'] + list(map(lambda x: f"func_{x}", anno_db_tags)):
    
        if plot_type.startswith('func'):
            # Getting current functional annotation
            curr_partition_funcs = partition_funcs[plot_type[plot_type.find('_') + 1:]]
        
        f, ax = plt.subplots(figsize=(20, 35))
        
        if plot_type == 'genes':
            wordclouds = {
                i: WordCloud(
                    max_words=30, min_font_size=15, background_color='white', mask=get_elipsis_mask()
                ).generate_from_frequencies(gene_score_dict).recolor(color_func=highlight_TFs) 
                    for i, gene_score_dict in norm_partition_genes.items()
            }
        else:
            word_counts = {
                i: WordCloud(max_words=30, min_font_size=15, stopwords=stopwords).process_text(text) for i, text in curr_partition_funcs.items()
            }
            word_counts = {
                i: (freqs if freqs else {'no found function': 1}) for i, freqs in word_counts.items()
            }  # dealing with no word case
            wordclouds = {
                i: WordCloud(
                    max_words=30, min_font_size=15, stopwords=stopwords, background_color='white', mask=get_elipsis_mask()
                ).generate_from_frequencies(freqs) for i, freqs in word_counts.items()
            }
            
        # Plotting clouds
        for i, coords in partition_coords.items():
            x, y = zip(*coords)
            min_x, max_x = min(x), max(x)
            min_y, max_y = min(y), max(y)
            ax.imshow(wordclouds[i], interpolation='bilinear', extent=[min_x, max_x, min_y, max_y])
        print(f'Finished plotting {plot_type} word cloud..')
        
        nx.draw(squeezed_G, squeezed_pos, ax=ax, arrowstyle="->", arrowsize=20, 
                connectionstyle=f'arc3, rad = 0.25', edge_color='gray', width=0.4, 
                node_color='k', node_size=50, alpha=0.02)
        nx.draw_networkx_nodes(squeezed_G, squeezed_pos, ax=ax, node_size=100, 
                               nodelist=list(squeezed_partition.keys()), 
                               node_color=list(squeezed_partition.values()), 
                               cmap=cmap, alpha=0.005)
        print(f'Finished plotting {plot_type} nodes..')

        ax.set_title(f'Found communities ({pat}, "all", {data}), '
                     f'annotation - {plot_type}', 
                     fontsize=30)
        plt.axis('off')

        plt.savefig(f'{figs_as}_{plot_type}.png', bbox_inches='tight', dpi=400)
            
    print('Finished plotting..\n')
            
    
    ###### SAVING DATAFRAME CONTAINING INFORMATION ABOUT EACH COMMUNITY ######

    def compute_community_info(i):
        """
        Parallel saving of the dataframe.
        """

        # Getting information for each community
        genes = list(all_partition_genes[i].keys())
        community_subgraph = G.subgraph(genes)

        communities_i = pd.Series(dtype='object')

        # Setting tqdm logs
        # t.set_description(f'Saving info about {i} cluster, size={community_subgraph.order()}')

        # Getting information about cluster genes
        central_genes_and_scores = {
            gene: all_partition_genes[i][gene] for k, gene in enumerate(genes) if k < limit_anno_until
        }

        non_lambert_TFs = [
            f'{gene} (rank={k})' for k, gene in enumerate(central_genes_and_scores.keys(), start=1) if gene not in lambert_TF_names
        ]
        non_dorothea_TFs = [
            f'{gene} (rank={k})' for k, gene in enumerate(central_genes_and_scores.keys(), start=1) if gene not in dorothea_TF_names
        ]

        # Filling dataframe with the information
        communities_i['num_nodes'] = community_subgraph.number_of_nodes()
        communities_i['num_edges'] = community_subgraph.number_of_edges()
        communities_i['all_sorted_genes'] = '; '.join(
            f'{gene} (score={score})' for gene, score in all_partition_genes[i].items()
        )
        communities_i['sorted_central_genes_scores'] = '; '.join(
            f'{gene} (score={score:.2f})' for gene, score in central_genes_and_scores.items()
        )
        communities_i['non_lambert_2018_TF_central_genes'] = '; '.join(non_lambert_TFs)
        communities_i['non_dorothea_TF_central_genes'] = '; '.join(non_dorothea_TFs)
        communities_i['whole_G_central_genes_scores'] = '; '.join(
            f'{gene} (score={score:.2f})' for gene, score in whole_G_central_genes.items()
        )

        # Filling information about newly found gene-gene links (based on absence in KEGG and Hallmark)
        top_cluster_links = set()

        iter_i = 0

        for st, end, edge_info in sorted(community_subgraph.edges(data=True), 
                                         key=lambda t: t[2]['importance'], 
                                         reverse=True):

            # If the current (reverse directed) link was not encountered previously..
            if (end, st) not in [(uniq_st, uniq_end) for uniq_st, uniq_end, _ in top_cluster_links]:
                top_cluster_links.add((st, end, edge_info['importance']))
                iter_i += 1
            if iter_i == save_top_new_found_cluster_links:
                break

        for anno_tag in ['KEGG', 'hallmark']:

            curr_db = load_gene_func_db(anno_tag)
            tmp_list = []

            # if `st` gene and `end` gene have non-overlapping annotations..
            for st, end, imp in top_cluster_links:
                st_anno_IDs = set(curr_db[curr_db.index == st]['ID'])
                end_anno_IDs = set(curr_db[curr_db.index == end]['ID'])
                if len(st_anno_IDs.intersection(end_anno_IDs)) == 0 and \
                        (len(st_anno_IDs) != 0 or len(end_anno_IDs) != 0):
                    tmp_list.append(f"{st} ({' & '.join(st_anno_IDs)}) <-> {end} ({' & '.join(end_anno_IDs)})")

            communities_i[f'new_gene_gene_links_{anno_tag}'] = '; '.join(tmp_list)

        # Filling information about cluster functions
        for tag, gene_func in gene_func_dbs.items():

            curr_partition_funcs = partition_funcs[tag]

            # Filling main functions - non duplicates at the top 
            main_functions = list(dict.fromkeys([  # dropping duplicates, but preserving order
                func for gene in central_genes_and_scores.keys() 
                    for func in gene_func[gene_func.index == gene].to_list()
            ]))
            gene_with_main_functions = [
                ','.join(
                    gene_func[gene_func == func].loc[lambda x: x.index.isin(genes)].index.to_list()
                ) for func in main_functions
            ]
            main_functions = [
                f'>>> {func} <<<: {gene}' for gene, func in zip(gene_with_main_functions, main_functions)
            ]
            communities_i[f'main_functions_{tag}'] = '; '.join(main_functions)  # saving..

            # Saving functions corresponding to each gene
            central_functions_per_gene = [
                f">>> {gene} <<<: {' & '.join(gene_func[gene_func.index == gene].to_list())}" for gene in central_genes_and_scores.keys()
            ]
            communities_i[f'sorted_central_functions_{tag}'] = '; '.join(central_functions_per_gene)  # saving..

            # Saving most frequent function words
            freq_words = WordCloud(
                max_words=30, min_font_size=15, stopwords=stopwords
            ).process_text(curr_partition_funcs[i])
            freq_words = dict(
                sorted(freq_words.items(), key=lambda x: x[1], reverse=True)
            ) if freq_words else {'no found function': 1}  # dealing with no word case
            communities_i[f'most_frequent_function_words_{tag}'] = '; '.join(freq_words.keys())  # saving

            # Saving other functions present in this cluster
            other_functions = list(dict.fromkeys([  # dropping duplicates, but preserving order
                func for gene in genes if gene not in central_genes_and_scores.keys() 
                    for func in gene_func[gene_func.index == gene].to_list() if func not in main_functions
            ]))[:other_functions_until]
            genes_with_other_functions = [
                ','.join(
                    gene_func[gene_func == func].loc[lambda x: x.index.isin(genes)].index.to_list()
                ) for func in other_functions
            ]
            other_functions = [
                f'>>> {func} <<<: {gene}' for gene, func in zip(genes_with_other_functions, other_functions)
            ]
            communities_i[f'other_functions_{tag}'] = '; '.join(other_functions)  # saving

        # Filling information about top inter-community links
        # t_sub = tqdm(range(num_partitions), ascii=True, leave=False)
        for k in range(num_partitions):  # t_sub:
            # t_sub.set_description(f'Extracting top inter-community links with {k}')

            if i != k:
                genes_in_k = list(all_partition_genes[k].keys())

                # Getting the subgraph that contains central genes in community_i and all genes in comunity_k
                G_central_i_k = G.subgraph(list(central_genes_and_scores.keys()) + genes_in_k)
                # Getting the subgraph that contains all genes from community_i and community_k
                G_i_k = G.subgraph(genes + genes_in_k)

                # Creating two helper sets that allow us to keep only unique links
                links_central_i_k = set()
                links_i_k = set()

                iter_i = 0

                # Getting out top links from the second subgraph
                for st, end, edge_info in sorted(G_central_i_k.edges(data=True), 
                                                 key=lambda t: t[2]['importance'], 
                                                 reverse=True):
                    # If the current (reverse directed) link was not encountered previously..
                    if (end, st) not in [(uniq_st, uniq_end) for uniq_st, uniq_end, _ in links_central_i_k] and \
                            ((st in genes and end not in genes) or (end in genes and st in genes)):
                        links_central_i_k.add((st, end, edge_info['importance']))
                        iter_i += 1
                    if iter_i == save_top_intercommunity_links_until:
                        break

                iter_i = 0

                # Getting out top links from the second subgraph
                for st, end, edge_info in sorted(G_i_k.edges(data=True), 
                                                 key=lambda t: t[2]['importance'], 
                                                 reverse=True):
                    # If the current (reverse directed) link was not encountered previously..
                    if (end, st) not in [(uniq_st, uniq_end) for uniq_st, uniq_end, _ in links_i_k] and \
                            ((st in genes and end not in genes) or (end in genes and st in genes)):
                        links_i_k.add((st, end, edge_info['importance']))
                        iter_i += 1
                    if iter_i == save_top_intercommunity_links_until:
                        break

                # Adding top links to the dataframe
                communities_i[f'top_links_scores_central_genes<->community_{k}'] = \
                    '; '.join(f'{st} <-> {end} (score={score:.2f})' for st, end, score in links_central_i_k)
                communities_i[f'top_links_scores_with_community_{k}'] = \
                    '; '.join([f'{st} <-> {end} (score={score:.2f})' for st, end, score in links_i_k])

        return communities_i
    
    print('Saving info dataframe..')
    
    t = tqdm_cli(range(num_partitions), ascii=True)
    
    # Getting dataframe
    result = Parallel(n_jobs=num_workers)(delayed(compute_community_info)(i) for i in t)
    communities_df = pd.concat(result, axis=1).T.reindex(
        columns=[
            'num_nodes', 'num_edges',
            'main_functions_GO', 'main_functions_KEGG', 'main_functions_immunological', 'main_functions_hallmark', 
            'non_lambert_2018_TF_central_genes', 'non_dorothea_TF_central_genes', 
            'new_gene_gene_links_KEGG', 'new_gene_gene_links_hallmark',
            'whole_G_central_genes_scores',
            'other_functions_GO', 'other_functions_KEGG', 'other_functions_immunological', 'other_functions_hallmark',
            'sorted_central_genes_scores',
            'sorted_central_functions_GO', 'sorted_central_functions_KEGG', 'sorted_central_functions_immunological', 'sorted_central_functions_hallmark', 
            'most_frequent_function_words_GO', 'most_frequent_function_words_KEGG', 'most_frequent_function_words_immunological', 'most_frequent_function_words_hallmark',
            'all_sorted_genes'] + 
            [f'top_links_scores_central_genes<->community_{i}' for i in range(num_partitions)] + 
            [f'top_links_scores_with_community_{i}' for i in range(num_partitions)
            ]
    )
    
    # Saving dataframe
    communities_df.to_pickle(data_as)
    print(f"Saved the data to {data_as}!\n")