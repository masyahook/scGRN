# General
import os

# Tools/utils
import itertools
import multiprocessing
from tqdm.notebook import tqdm
from tqdm import tqdm as tqdm_cli
from functools import reduce  # for aggregate functions
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
import pygraphviz as pgv
import colorcet as cc
from matplotlib.colors import ListedColormap
from wordcloud import WordCloud, STOPWORDS
from termcolor import colored  # colored text output

from sklearn.preprocessing import MinMaxScaler


scale = lambda x, min_y, max_y: list(MinMaxScaler(feature_range=(min_y, max_y)).fit_transform(np.expand_dims(np.array(x), axis=1))[:, 0])
stopwords = STOPWORDS.union({'regulation', 'activity', 'positive', 'negative', 
                                 'catabolic', 'process', 'protein', 'complex', 
                                 'binding', 'response'})

def save_pickle(f, fn):
    """
    Save object as a pickle file (usually used for dicts).
    
    f: file object
    fn: file name
    """
    import pickle
    
    with open(fn, 'wb') as fo:
        pickle.dump(f, fo)
        

def load_pickle(fn):
    """
    Load object from pickle file (usually used for dicts).
    
    fn: file name
    
    return: loaded object
    """
    import pickle
    
    with open(fn, 'rb') as fo:
        f = pickle.load(fo)
        
    return f


def load_gene_func_db_mapping():
    """
    Load downloaded gene function databases like MSigDB, GO, SIGNOR, etc mappings.
    """
    
    path_to_dbs = '/gpfs/projects/bsc08/bsc08890/data/Gene_func_associations'    

    file_mapping = {
        'GO': os.path.join(path_to_dbs, 'GO_annotation.tsv'),
        'DoRothEA': os.path.join(path_to_dbs, 'dorothea_regulons.tsv'),
        'MSigDB_hallmark_gene_sets_h': os.path.join(path_to_dbs, 'MSigDB', 'h.all.v7.5.1.symbols.gmt.txt'),
        'MSigDB_curated_gene_sets_c2_all': os.path.join(path_to_dbs, 'MSigDB', 'c2.all.v7.5.1.symbols.gmt.txt'), 
        'MSigDB_curated_gene_sets_c2_cp': os.path.join(path_to_dbs, 'MSigDB', 'c2.cp.v7.5.1.symbols.gmt.txt'),
        'MSigDB_curated_gene_sets_c2_cp_kegg': os.path.join(path_to_dbs, 'MSigDB', 'c2.cp.kegg.v7.5.1.symbols.gmt.txt'),
        'MSigDB_curated_gene_sets_c2_cp_reactome': os.path.join(path_to_dbs, 'MSigDB', 'c2.cp.reactome.v7.5.1.symbols.gmt.txt'),
        'MSigDB_curated_gene_sets_c2_cp_wikipathways': os.path.join(path_to_dbs, 'MSigDB', 'c2.cp.wikipathways.v7.5.1.symbols.gmt.txt'),
        'MSigDB_regulatory_target_gene_sets_c3_all': os.path.join(path_to_dbs, 'MSigDB', 'c3.all.v7.5.1.symbols.gmt.txt'),
        'MSigDB_immunologic_signature_gene_sets_c7_all': os.path.join(path_to_dbs, 'MSigDB', 'c7.all.v7.5.1.symbols.gmt.txt')
    }

    return file_mapping


def load_gene_func_db(db, reload=False, as_series=False):
    """
    Load data from gene function database like MSigDB, GO, DoRothEA, or etc. The output will be in a Pandas Dataframe format.
    """
    
    def short_to_long_tag(db):
        """
        Transforming the short annotation db tag to long format.
        """
        if db == 'GO':
            return 'GO'
        elif db == 'DoRothEA':
            return 'DoRothEA'
        elif db == 'KEGG':
            return 'MSigDB_curated_gene_sets_c2_cp_kegg'
        elif db == 'hallmark':
            return 'MSigDB_hallmark_gene_sets_h'
        elif db == 'reactome':
            return 'MSigDB_curated_gene_sets_c2_cp_reactome'
        elif db == 'wikipathways':
            return 'MSigDB_curated_gene_sets_c2_cp_wikipathways'
        elif db == 'immunological':
            return 'MSigDB_immunologic_signature_gene_sets_c7_all'
        elif db == 'curated':
            return 'MSigDB_curated_gene_sets_c2_all'
        elif db == 'canonical':
            return 'MSigDB_curated_gene_sets_c2_cp'
        else:
            raise NotImplementedError(f"The tag '{db}' not found in the downloaded databases..")
       
    # Getting the db tag
    db = short_to_long_tag(db)            
    
    # Dealing with path names
    db_path = load_gene_func_db_mapping()[db]
    db_folder = db_path[:db_path.rfind('/')]
    
    if db.startswith('MSigDB'):
        
        if reload:
    
            # Loading saved data from https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp
            with open(db_path, 'r') as f:
                entries = [f_line.strip().split('\t') for f_line in f.readlines()]

            # Converting data to pandas format
            entries_df = pd.DataFrame(columns=['ID', 'link', 'genes'])

            for i, entry in enumerate(entries):
                entries_df.loc[i, 'ID'] = entry[0]
                entries_df.loc[i, 'link'] = entry[1]
                entries_df.loc[i, 'genes'] = entry[2:]
                
            # Adding some information about each functional group by scraping from the web
            for i, row in tqdm_cli(entries_df.iterrows(), leave=False, total=entries_df.shape[0], ascii=True):
                brief, full = pd.read_html(
                    row['link']
                )[1].loc[
                    lambda x: x[0].isin(['Brief description', 'Full description or abstract'])
                ][1].values
                
                entries_df.loc[i, 'brief_desc'] = brief
                entries_df.loc[i, 'full_desc'] = full
                
            # Transforming data to format when gene is the index
            per_gene_entries_df = pd.DataFrame(columns=['gene_name', 'ID', 'link', 'brief_desc', 'full_desc'])
            for i, row in entries_df.iterrows():
                num_genes = len(row['genes'])
                row_df = pd.DataFrame(
                    dict(
                        gene_name=row['genes'], ID=[row['ID']]*num_genes, link=[row['link']]*num_genes,
                        brief_desc=[row['brief_desc']]*num_genes, full_desc=[row['full_desc']]*num_genes
                    )
                )
                per_gene_entries_df = pd.concat([per_gene_entries_df, row_df], axis=0)

            out = per_gene_entries_df.set_index('gene_name')
            out.to_csv(os.path.join(db_folder, f'{db}.tsv'), sep='\t')
            
        else:
            
            # Reading previously saved data
            out = pd.read_csv(os.path.join(db_folder, f'{db}.tsv'), sep='\t', index_col=0)
            
        if as_series:
            out = out['brief_desc']
        
    elif db == 'GO':
        
        # Loading gene GO description - it was downloaded using ENSEMBL BioMart
        entries_df = pd.read_csv(db_path, sep='\t', index_col=0) \
            .drop(columns=['GO term evidence code']) \
            .dropna(subset=['Gene name', 'GO term name', 'GO domain']) \
            .drop_duplicates()
        entries_df = entries_df.rename(
            columns={
                'Gene name': 'gene_name', 'GO term accession': 'ID', 'GO term name': 'brief_desc',
                'GO term definition': 'full_desc', 'GO domain': 'domain'
            }
        )
        
        out = entries_df.set_index('gene_name')
        
        if as_series:
            out = out[
                (out['domain'] == 'molecular_function') | (out['domain'] == 'biological_process')
            ]['brief_desc']
        
    elif db == 'DoRothEA':
        
        # Loading DoRothEA database of TF-target gene associations
        out = pd.read_csv(db_path, sep='\t')
        
    else:
        
        raise NotImplementedError
        
    return out


def fancy_draw_network_edge_labels(
    G,
    pos,
    edge_labels=None,
    label_pos=0.5,
    font_size=10,
    font_color="k",
    font_family="sans-serif",
    font_weight="normal",
    alpha=None,
    bbox=None,
    horizontalalignment="center",
    verticalalignment="center",
    ax=None,
    rotate=True,
    clip_on=True,
    rad=0
):
    """Draw edge labels.

    Parameters
    ----------
    G : graph
        A networkx graph

    pos : dictionary
        A dictionary with nodes as keys and positions as values.
        Positions should be sequences of length 2.

    edge_labels : dictionary (default={})
        Edge labels in a dictionary of labels keyed by edge two-tuple.
        Only labels for the keys in the dictionary are drawn.

    label_pos : float (default=0.5)
        Position of edge label along edge (0=head, 0.5=center, 1=tail)

    font_size : int (default=10)
        Font size for text labels

    font_color : string (default='k' black)
        Font color string

    font_weight : string (default='normal')
        Font weight

    font_family : string (default='sans-serif')
        Font family

    alpha : float or None (default=None)
        The text transparency

    bbox : Matplotlib bbox, optional
        Specify text box properties (e.g. shape, color etc.) for edge labels.
        Default is {boxstyle='round', ec=(1.0, 1.0, 1.0), fc=(1.0, 1.0, 1.0)}.

    horizontalalignment : string (default='center')
        Horizontal alignment {'center', 'right', 'left'}

    verticalalignment : string (default='center')
        Vertical alignment {'center', 'top', 'bottom', 'baseline', 'center_baseline'}

    ax : Matplotlib Axes object, optional
        Draw the graph in the specified Matplotlib axes.

    rotate : bool (deafult=True)
        Rotate edge labels to lie parallel to edges

    clip_on : bool (default=True)
        Turn on clipping of edge labels at axis boundaries

    Returns
    -------
    dict
        `dict` of labels keyed by edge

    Examples
    --------
    >>> G = nx.dodecahedral_graph()
    >>> edge_labels = nx.draw_networkx_edge_labels(G, pos=nx.spring_layout(G))

    Also see the NetworkX drawing examples at
    https://networkx.org/documentation/latest/auto_examples/index.html

    See Also
    --------
    draw
    draw_networkx
    draw_networkx_nodes
    draw_networkx_edges
    draw_networkx_labels
    """
    import matplotlib.pyplot as plt
    import numpy as np

    if ax is None:
        ax = plt.gca()
    if edge_labels is None:
        labels = {(u, v): d for u, v, d in G.edges(data=True)}
    else:
        labels = edge_labels
    text_items = {}
    for (n1, n2), label in labels.items():
        (x1, y1) = pos[n1]
        (x2, y2) = pos[n2]
        (x, y) = (
            x1 * label_pos + x2 * (1.0 - label_pos),
            y1 * label_pos + y2 * (1.0 - label_pos),
        )
        pos_1 = ax.transData.transform(np.array(pos[n1]))
        pos_2 = ax.transData.transform(np.array(pos[n2]))
        linear_mid = 0.5*pos_1 + 0.5*pos_2
        d_pos = pos_2 - pos_1
        rotation_matrix = np.array([(0,1), (-1,0)])
        ctrl_1 = linear_mid + rad*rotation_matrix@d_pos
        ctrl_mid_1 = 0.5*pos_1 + 0.5*ctrl_1
        ctrl_mid_2 = 0.5*pos_2 + 0.5*ctrl_1
        bezier_mid = 0.5*ctrl_mid_1 + 0.5*ctrl_mid_2
        (x, y) = ax.transData.inverted().transform(bezier_mid)

        if rotate:
            # in degrees
            angle = np.arctan2(y2 - y1, x2 - x1) / (2.0 * np.pi) * 360
            # make label orientation "right-side-up"
            if angle > 90:
                angle -= 180
            if angle < -90:
                angle += 180
            # transform data coordinate angle to screen coordinate angle
            xy = np.array((x, y))
            trans_angle = ax.transData.transform_angles(
                np.array((angle,)), xy.reshape((1, 2))
            )[0]
        else:
            trans_angle = 0.0
        # use default box of white with white border
        if bbox is None:
            bbox = dict(boxstyle="round", ec=(1.0, 1.0, 1.0), fc=(1.0, 1.0, 1.0))
        if not isinstance(label, str):
            label = str(label)  # this makes "1" and 1 labeled the same

        t = ax.text(
            x,
            y,
            label,
            size=font_size,
            color=font_color,
            family=font_family,
            weight=font_weight,
            alpha=alpha,
            horizontalalignment=horizontalalignment,
            verticalalignment=verticalalignment,
            rotation=trans_angle,
            transform=ax.transData,
            bbox=bbox,
            zorder=1,
            clip_on=clip_on,
        )
        text_items[(n1, n2)] = t

    ax.tick_params(
        axis="both",
        which="both",
        bottom=False,
        left=False,
        labelbottom=False,
        labelleft=False,
    )

    return text_items


def draw_graph(G, pos, ax, TF_names=None, label_edges=True, node_size=1200,
               if_alpha_edges=False, plot_cmap=True, cmap=plt.cm.plasma):
    """
    Draw GRN using NetworkX.
    
    """
    import networkx as nx
    import matplotlib as mpl
    
    seed = 42
    alpha = 0.5  
    
    nodes = nx.draw_networkx_nodes(G, pos, node_color='pink', ax=ax, node_size=node_size)
    if TF_names is not None:
        nx.draw_networkx_nodes(G.subgraph(TF_names), pos=pos, node_color='limegreen', ax=ax, node_size=node_size)
    nx.draw_networkx_labels(G, pos, ax=ax)
    
    if G.edges():
    
        edges, importances = zip(*nx.get_edge_attributes(G, 'importance').items())
        edges, rhos = zip(*nx.get_edge_attributes(G, 'rho').items())
        widths = scale(importances, 1, 10)
        
        curved_mask = [reversed(edge) in G.edges() for edge in G.edges()]
        
        edge_meta = {
            'curved_edge': [e for m, e in zip(curved_mask, edges) if m],
            'curved_importance': [e for m, e in zip(curved_mask, importances) if m],
            'curved_rho': [e for m, e in zip(curved_mask, rhos) if m],
            'curved_width': [e for m, e in zip(curved_mask, widths) if m],
            'straight_edge': [e for m, e in zip(curved_mask, edges) if not m],
            'straight_importance': [e for m, e in zip(curved_mask, importances) if not m],
            'straight_rho': [e for m, e in zip(curved_mask, rhos) if not m],
            'straight_width': [e for m, e in zip(curved_mask, widths) if not m]
        }

        mpl_straight_edges = nx.draw_networkx_edges(
            G, pos, ax=ax, edgelist=edge_meta['straight_edge'], arrowstyle="->", arrowsize=30,
            edge_color=edge_meta['straight_rho'], edge_cmap=cmap, width=edge_meta['straight_width'], node_size=node_size
        )
        mpl_curved_edges = nx.draw_networkx_edges(
            G, pos, ax=ax, edgelist=edge_meta['curved_edge'], connectionstyle=f'arc3, rad = 0.25', arrowstyle="->", 
            arrowsize=30, edge_color=edge_meta['curved_rho'], edge_cmap=cmap, width=edge_meta['curved_width'], node_size=node_size
        )
        
        if mpl_curved_edges is None:
            mpl_curved_edges = []
        if mpl_straight_edges is None:
            mpl_straight_edges = []
                
        if plot_cmap:
            pc = mpl.collections.PatchCollection(mpl_straight_edges + mpl_curved_edges, cmap=cmap)
            pc.set_array(rhos)
            cbar = plt.colorbar(pc)
            cbar.ax.set_ylabel('Spearman correlation', rotation=270, fontsize=17, labelpad=17)
    
        if label_edges:
            edge_weights = nx.get_edge_attributes(G, 'importance')
            curved_edge_labels = {edge: f'{edge_weights[edge]:.1f}' for edge in edge_meta['curved_edge']}
            straight_edge_labels = {edge: f'{edge_weights[edge]:.1f}' for edge in edge_meta['straight_edge']}
            fancy_draw_network_edge_labels(G, pos, ax=ax, edge_labels=curved_edge_labels, rotate=False, rad = 0.25)
            fancy_draw_network_edge_labels(G, pos, ax=ax, edge_labels=straight_edge_labels, rotate=False)

        if if_alpha_edges:
            edge_importances = [w_dict['importance'] for st, end, w_dict in G.edges(data=True)]
            alpha_importances = [(w - min(edge_importances))/(max(edge_importances) - min(edge_importances))*alpha + alpha for w in edge_importances]
            # set alpha value for each edge
            for i in range(len(alpha_importances)):
                edges[i].set_alpha(alpha_importances[i])
        
    plt.axis('off')
    
    return ax


def process_ndex_net(raw_G, G_name, remove_single_nodes=True):
    """
    Process ndex net
    """
    
    import networkx as nx
    import re
    
    G = nx.DiGraph()
        
    # Remove empty nodes
    raw_G.remove_nodes_from([node[0] for node in raw_G.nodes(data=True) if not node[1]])
    
    # Define the network type
    node_0 = list(raw_G.nodes(data=True))[0]
    if 'bel_function_type' in node_0[1].keys():
        source = 'causalbionet'
    elif '__gpml:textlabel' in node_0[1].keys():
        source = 'wikipathways'
    else:
        source = 'signor' 
        
    if source == 'causalbionet':
        tmp_nodes = {}
        node_duplicates = {}
        for node_i, node_info in raw_G.nodes(data=True):
            if ':' in node_info['name']:
                name = '+'.join(re.findall(r'(?<=\:).+?(?=\))', node_info['name'], re.IGNORECASE))
            else:
                name = re.findall(r'(?<=\().+?(?=\))', node_info['name'], re.IGNORECASE)[0]
            name = name[:name.find(',')] if '(' in name else name
            if name not in node_duplicates:
                node_duplicates[name] = 0
            else:
                node_duplicates[name] += 1
            new_name = node_duplicates[name]*'_' + name + node_duplicates[name]*'_'
            tmp_nodes[node_i] = {'old_name': name, 'new_name': new_name}
            attrs = {'type': node_info['bel_function_type'], 'full_name': node_info['name']}
            G.add_node(new_name, **attrs)
        for st_i, end_i, edge_info in raw_G.edges(data=True):
            st_name, end_name = tmp_nodes[st_i]['new_name'], tmp_nodes[end_i]['new_name']
            inter = edge_info['interaction']
            attrs = {'regulation': 'up' if 'increases' in inter.lower() else 'down' if 'decreases' in inter.lower() else None, 'info': inter}
            G.add_edge(st_name, end_name, **attrs)
    elif source == 'wikipathways':
        # Removing empty nodes
        raw_G.remove_nodes_from([node_i for node_i, node_info in raw_G.nodes(data=True) if 'name' not in node_info])
        # Removing phosphorus
        raw_G.remove_nodes_from([node_i for node_i, node_info in raw_G.nodes(data=True) if node_info['name'] == 'P'])
        tmp_nodes = {}
        node_duplicates = {}
        for node_i, node_info in raw_G.nodes(data=True):
            name = node_info['name']
            if name not in node_duplicates:
                node_duplicates[name] = 0
            else:
                node_duplicates[name] += 1
            new_name = node_duplicates[name]*'_' + name + node_duplicates[name]*'_'
            tmp_nodes[node_i] = {'old_name': name, 'new_name': new_name}
            attrs = {'type': 'None', 'full_name': 'None'}
            G.add_node(new_name, **attrs)
        for st_i, end_i, edge_info in raw_G.edges(data=True):
            st_name, end_name = tmp_nodes[st_i]['new_name'], tmp_nodes[end_i]['new_name']
            attrs = {'regulation': None, 'info': None}
            G.add_edge(st_name, end_name, **attrs)
    elif source == 'signor':
        tmp_nodes = {}
        node_duplicates = {}
        for node_i, node_info in raw_G.nodes(data=True):
            name = node_info['name']
            if name not in node_duplicates:
                node_duplicates[name] = 0
            else:
                node_duplicates[name] += 1
            new_name = node_duplicates[name]*'_' + name + node_duplicates[name]*'_'
            tmp_nodes[node_i] = {'old_name': name, 'new_name': new_name}
            attrs = {'node_type': node_info['type'], 'full_name': None}
            G.add_node(new_name, **attrs)
        for st_i, end_i, edge_info in raw_G.edges(data=True):
            inter = edge_info['interaction']
            attrs = {'regulation': 'up' if 'up' in inter.lower() else 'down' if 'down' in inter.lower() else None, 'info': edge_info['mechanism'] if 'mechanism' in edge_info else None}
            st_name, end_name = tmp_nodes[st_i]['new_name'], tmp_nodes[end_i]['new_name']
            G.add_edge(st_name, end_name, **attrs)
    else:
        raise NameError(f'No identified source for {G_name}')
        
    if remove_single_nodes:
        single_node_sets = list(filter(lambda x: len(x) == 1, nx.weakly_connected_components(G)))
        single_nodes = reduce(lambda x, y: x.union(y), single_node_sets) if len(single_node_sets) > 0 else set()
        G.remove_nodes_from(single_nodes)

    return G
    
    
def get_tf_targ_ctx(df):
    tf_target_dict = {'TF': [], 'target': [], 'importance': []}
    tf_target_info = (
        df.droplevel(axis=0, level=1).droplevel(axis=1, level=0)['TargetGenes']
          .map(set)  # transform each list into set
          .groupby('TF').agg(lambda x: reduce(lambda a, b: a.union(b), x))  # combine all targets per TF
    )
    for tf, target_info in tf_target_info.iteritems():
        tf_target_dict['TF'] += [tf for target_name, score in target_info]
        tf_target_dict['target'] += [target_name for target_name, score in target_info]
        tf_target_dict['importance'] += [score for target_name, score in target_info]
    return pd.DataFrame(tf_target_dict)


def style_data_availability(df):
    return df.style.apply(lambda x: ["background-color: green" if v == '+' else "background-color: red" if v == '-' else 'background: white' for v in x], axis=1)


def get_adj_list(pat, data, data_type, method='grnboost2', is_filter=False):
    _DATA_HOME = '/gpfs/projects/bsc08/bsc08890/res/covid_19'
    data_suffix = 'TF_cor' if data_type == 'TF' else 'TF_ctx' if data_type == 'ctx' else 'cor'
    filter_suffix = '' if not is_filter else '_filtered'
    return pd.read_pickle(os.path.join(_DATA_HOME, pat, 'data', method, 'pickle', f'{data}_{data_suffix}{filter_suffix}.pickle'))


def get_nx_graph(pat, data, data_type, method='grnboost2', is_filter=False):
    _DATA_HOME = '/gpfs/projects/bsc08/bsc08890/res/covid_19'
    data_suffix = 'TF_cor' if data_type == 'TF' else 'TF_ctx' if data_type == 'ctx' else 'cor'
    filter_suffix = '' if not is_filter else '_filtered'
    return nx.read_gpickle(os.path.join(_DATA_HOME, pat, 'data', method, 'nx_graph', f'{data}_{data_suffix}{filter_suffix}.gpickle'))


def netgraph_community_layout(G, node_to_community, community_scale=1., node_scale=2., seed=42):
    """
    Compute the node positions for a modular graph.
    """

    # assert that there multiple communities in the graph; otherwise abort
    communities = set(node_to_community.values())
    if len(communities) < 2:
        warnings.warn("Graph contains a single community. Unable to compute a community layout. Computing spring layout instead.")
        return nx.spring_layout(G, weight='importance', **kwargs)

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


def _get_community_sizes(node_to_community):
    """
    Compute the area of the canvas reserved for each community.
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


def _get_community_positions(G, node_to_community, community_scale, seed, simple=True):
    """
    Compute a centroid position for each community.
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

def _find_between_community_edges(G, node_to_community, fixed_community=None):
    """Convert the graph into a weighted network of communities."""
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

def squeeze_graph(G, partition):
    """
    Squeeze graph by picking only top nodes (according to number of connections) in each partition. This
    step is needed to speed up the networkx visualization and show only the general POV on the graph.
    """
    
    #### STEP 1 - filtering nodes
    
    # Setting up the size of the squeezed graph (number of nodes)
    approximate_size = 4000
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
        word_counts = {i: WordCloud(stopwords=stopwords).process_text(text) for i, text in partition_funcs.items()}
        word_counts = {
            i: (freqs if freqs else {'no found function': 1}) for i, freqs in word_counts.items()
        }  # dealing with no word case
        wordclouds = {
            i: WordCloud(
                max_font_size=40, stopwords=stopwords, background_color='white', mask=get_elipsis_mask()
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
                max_font_size=80, background_color='white', mask=get_elipsis_mask()
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

           
def process_communities(pat, data, algo='leiden', if_betweenness=True, limit_anno_until=50, 
                        k=5000, save_top_intercommunity_links_until=20, other_functions_until=20, 
                        save_top_new_found_cluster_links=20, seed=42):
    """
    Process graph by finding its communities, annotate its communities, and save everything into .tsv format.
    """
    
    print('\nPerforming community analysis..\n\n')
    
    # Setting pathways to files
    _PROJ_PATH = '/gpfs/projects/bsc08/bsc08890'
    _FMETA = os.path.join(_PROJ_PATH, 'data/GSE145926_RAW/metadata.tsv')
    _DATA_HOME = os.path.join(_PROJ_PATH, 'res/covid_19')

    # Loading sample meta data, reordering patients
    full_meta = pd.read_csv(_FMETA, sep='\t', index_col=0)

    # Getting information about all patients
    _ALL_PATIENTS = full_meta.index.to_list()
    _ALL_FIG_DIRS = {
        pat: os.path.join(_DATA_HOME, pat, 'figs/grnboost2') for pat in _ALL_PATIENTS
    }
    
    # Getting plot titles
    dtype_title = 'gene-gene links'
    data_title = 'all data' if data == 'raw_data' else 'all data, only HVGs' if 'HVG' in data else data.replace('raw_data_', '').replace('_', ' ')
    
    # Loading lists of TFs from Lambert 2018 and DoRothEA, in the latter case we will keep only confident regulons
    lambert_TF_names = pd.read_csv(os.path.join(_PROJ_PATH, 'data/TF_lists/lambert2018.txt'), header=None)[0].to_list()
    dorothea_TF_names = list(
        pd.read_csv(os.path.join(_PROJ_PATH, 'data/TF_lists/dorothea_regulons.tsv'), sep='\t') \
            .loc[lambda x: x['confidence'].isin(['A', 'B', 'C'])]['tf'].unique()
    )
    
    # Loading the graph
    G = get_nx_graph(pat, data, 'all', is_filter=True)
    print(f"Loaded the graph: {colored('pat', 'green')}='{colored(pat, 'red')}', "
          f"{colored('data', 'green')}='{colored(data, 'red')}', "
          f"{colored('data_type', 'green')}='{colored('all', 'red')}'\n")
    
    
    ###### FINDING COMMUNITIES IN THE GRAPH #######
    
    print('Finding communities in the graph..')
    
    if algo == 'louvain':
        partition = community_louvain.best_partition(G.to_undirected(), weight='importance', random_state=seed)
    else:
        G_igraph = ig.Graph.from_networkx(G.to_undirected())
        la_partition = la.find_partition(G_igraph, la.ModularityVertexPartition, weights='importance', seed=seed)
        partition = {G_igraph.vs[node]['_nx_name']: i for i, cluster_nodes in enumerate(la_partition) for node in cluster_nodes}
        
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
        top_len = min(limit_anno_until, len(genes))
        gene_scores = dict(
            sorted(
                compute_centrality(
                    G.subgraph(genes), k=min(G.subgraph(genes).order(), k), **distance_metric
                ).items(), 
                key=lambda x: x[1], reverse=True
            )
        )
        all_partition_genes[i] = gene_scores
        central_gene_scores = {gene: gene_scores[gene] for k, gene in enumerate(gene_scores.keys()) if k < top_len}
        
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
        
        f, ax = plt.subplots(figsize=(25, 45))
        
        if plot_type == 'genes':
            wordclouds = {
                i: WordCloud(
                    max_font_size=80, background_color='white', mask=get_elipsis_mask()
                ).generate_from_frequencies(gene_score_dict) for i, gene_score_dict in norm_partition_genes.items()
            }
        else:
            word_counts = {
                i: WordCloud(stopwords=stopwords).process_text(text) for i, text in curr_partition_funcs.items()
            }
            word_counts = {
                i: (freqs if freqs else {'no found function': 1}) for i, freqs in word_counts.items()
            }  # dealing with no word case
            wordclouds = {
                i: WordCloud(
                    max_font_size=40, stopwords=stopwords, background_color='white', mask=get_elipsis_mask()
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
                connectionstyle=f'arc3, rad = 0.25', edge_color='k', width=0.4, 
                node_color='k', node_size=50, alpha=0.01)
        nx.draw_networkx_nodes(squeezed_G, squeezed_pos, ax=ax, node_size=100, 
                               nodelist=list(squeezed_partition.keys()), 
                               node_color=list(squeezed_partition.values()), 
                               cmap=cmap, alpha=0.01)
        print(f'Finished plotting {plot_type} nodes..')

        ax.set_title(f'Found communities ({pat}, {dtype_title}, {data_title}), '
                     f'annotation - {plot_type}', 
                     fontsize=30)
        plt.axis('off')

        plt.savefig(os.path.join(_ALL_FIG_DIRS[pat], f"{pat}_{data}_all_communities_{plot_type}.png"), bbox_inches='tight', dpi=400)
            
    print('Finished plotting..\n')
            
    
    ###### SAVING DATAFRAME CONTAINING INFORMATION ABOUT EACH COMMUNITY ######
    
    print('Saving info dataframe..')
    
    # Create a directory where the dataframe will be stored
    save_to_folder = os.path.join(_DATA_HOME, pat, 'data', 'grnboost2', f'{algo}_communities')
    os.makedirs(save_to_folder, exist_ok=True)
    
    communities_df = pd.DataFrame(
        index=pd.Series(range(num_partitions), name='community_i'), 
        columns=[
            'num_nodes', 'num_edges',
            'main_functions_GO', 'main_functions_KEGG', 'main_functions_immunological', 'main_functions_hallmark', 
            'sorted_central_genes', 'sorted_central_gene_scores'
            'sorted_central_functions_GO', 'sorted_central_functions_KEGG', 'sorted_central_functions_immunological', 'sorted_central_functions_hallmark', 
            'most_frequent_function_words_GO', 'most_frequent_function_words_KEGG', 'most_frequent_function_words_immunological', 'most_frequent_function_words_hallmark',
            'non_lambert_2018_TF_central_genes', 'non_dorothea_TF_central_genes', 
            'other_functions_GO', 'genes_with_other_functions_GO',
            'other_functions_KEGG', 'genes_with_other_functions_KEGG',
            'other_functions_immunological', 'genes_with_other_functions_immunological',
            'other_functions_hallmark', 'genes_with_other_functions_hallmark',
            'all_sorted_genes',
            'new_gene_gene_links_KEGG', 'new_gene_gene_links_hallmark',
            'whole_G_central_genes', 'whole_G_central_gene_scores'] + 
            [f'top_{n}_between_central_genes_and_community_{i}' for i in range(num_partitions) for n in ['links', 'link_scores']] + 
            [f'top_{n}_with_community_{i}' for i in range(num_partitions) for n in ['links', 'link_scores']]
    )
    
    t = tqdm_cli(range(num_partitions), ascii=True)

    for i in t:
        # Getting information for each community
        genes = list(all_partition_genes[i].keys())
        community_subgraph = G.subgraph(genes)
        
        # Setting tqdm logs
        t.set_description(f'Saving info about {i} cluster, size={community_subgraph.order()}')
        
        # Getting information about cluster genes
        central_genes_and_scores = {
            gene: all_partition_genes[i][gene] for k, gene in enumerate(genes) if k < top_len
        }
                
        non_lambert_TFs = [gene for gene in central_genes_and_scores.keys() if gene not in lambert_TF_names]
        non_dorothea_TFs = [gene for gene in central_genes_and_scores.keys() if gene not in dorothea_TF_names]

        # Filling dataframe with the information
        communities_df.loc[i, 'num_nodes'] = community_subgraph.number_of_nodes()
        communities_df.loc[i, 'num_edges'] = community_subgraph.number_of_edges()
        communities_df.loc[i, 'all_sorted_genes'] = '; '.join(genes)
        communities_df.loc[i, 'sorted_central_genes'] = '; '.join(central_genes_and_scores.keys())
        communities_df.loc[i, 'sorted_central_gene_scores'] = '; '.join(
            [f'{score:.2f}' for score in central_genes_and_scores.values()]
        )
        communities_df.loc[i, 'non_lambert_2018_TF_central_genes'] = '; '.join(non_lambert_TFs)
        communities_df.loc[i, 'non_dorothea_TF_central_genes'] = '; '.join(non_dorothea_TFs)
        communities_df.loc[i, 'whole_G_central_genes'] = '; '.join(whole_G_central_genes.keys())
        communities_df.loc[i, 'whole_G_central_gene_scores'] = '; '.join(
            [f'{score:.2f}' for score in whole_G_central_genes.values()]
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
                if len(st_anno_IDs.intersection(end_anno_IDs)) == 0:
                    tmp_list.append(f"{st} ({' & '.join(st_anno_IDs)}) <-> {end} ({' & '.join(end_anno_IDs)})")
                    
            communities_df.loc[i, f'new_gene_gene_links_{anno_tag}'] = '; '.join(tmp_list)
        
        # Filling information about cluster functions
        for tag, gene_func in gene_func_dbs.items():
        
            curr_partition_funcs = partition_funcs[tag]
            
            central_functions = list(dict.fromkeys([  # dropping duplicates, but preserving order
                func for gene in central_genes_and_scores.keys() 
                    for func in gene_func[gene_func.index == gene].to_list()
            ]))
            central_functions_per_gene = [
                ' & '.join(gene_func[gene_func.index == gene].to_list()) for gene in central_genes_and_scores.keys()
            ]

            freq_words = WordCloud(stopwords=stopwords).process_text(curr_partition_funcs[i])
            freq_words = dict(
                sorted(freq_words.items(), key=lambda x: x[1], reverse=True)
            ) if freq_words else {'no found function': 1}  # dealing with no word case

            other_functions = list(dict.fromkeys([  # dropping duplicates, but preserving order
                func for gene in genes if gene not in central_genes_and_scores.keys() 
                    for func in gene_func[gene_func.index == gene].to_list() if func not in central_functions
            ]))[:other_functions_until]
            genes_with_other_functions = [
                ' & '.join(gene_func[gene_func == func].index.to_list()) for func in other_functions
            ]
            
            # Filling dataframe with the information
            communities_df.loc[i, f'main_functions_{tag}'] = '; '.join(central_functions)
            communities_df.loc[i, f'sorted_central_functions_{tag}'] = '; '.join(central_functions_per_gene)
            communities_df.loc[i, f'most_frequent_function_words_{tag}'] = '; '.join(freq_words.keys())
            communities_df.loc[i, f'other_functions_{tag}'] = '; '.join(other_functions)
            communities_df.loc[i, f'genes_with_other_functions_{tag}'] = '; '.join(genes_with_other_functions)

        # Filling information about top inter-community links
        t_sub = tqdm_cli(range(num_partitions), ascii=True, leave=False)
        for k in t_sub:
            t_sub.set_description(f'Extracting top inter-community links with {k}')
            
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
                    if (end, st) not in [(uniq_st, uniq_end) for uniq_st, uniq_end, _ in links_central_i_k]:
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
                    if (end, st) not in [(uniq_st, uniq_end) for uniq_st, uniq_end, _ in links_central_i_k]:
                        links_central_i_k.add((st, end, edge_info['importance']))
                        iter_i += 1
                    if iter_i == save_top_intercommunity_links_until:
                        break
                        
                # Adding top links to the dataframe
                communities_df.loc[i, f'top_links_between_central_genes_and_community_{k}'] = \
                    '; '.join([f'{st}->{end}' for st, end, _ in links_central_i_k])
                communities_df.loc[i, f'top_link_scores_between_central_genes_and_community_{k}'] = \
                    '; '.join([f'{score:.2f}' for _, _, score in links_central_i_k])
                
                communities_df.loc[i, f'top_links_with_community_{k}'] = \
                    '; '.join([f'{st}->{end}' for st, end, _ in links_i_k])
                communities_df.loc[i, f'top_link_scores_with_community_{k}'] = \
                    '; '.join([f'{score:.2f}' for _, _, score in links_i_k])
    
    
    for k, arr in time_dict.items():
        print(f"Computing '{k}' takes around {np.mean(arr):.4f} seconds..")
    
    # Saving dataframe
    communities_df.to_pickle(os.path.join(save_to_folder, f'{pat}_{data}_all_communities_info.pickle'))
    print(f"Saved the data to {os.path.join(save_to_folder, f'{pat}_{data}_all_communities_info.pickle')}!\n")
        

def betweenness_centrality_parallel(G, processes=None):
    """Parallel betweenness centrality  function"""
    from multiprocessing import Pool
    
    def chunks(l, n):
        """Divide a list of nodes `l` in `n` chunks"""
        l_c = iter(l)
        while 1:
            x = tuple(itertools.islice(l_c, n))
            if not x:
                return
            yield x
    
    p = Pool(processes=processes)
    node_divisor = len(p._pool) * 4
    node_chunks = list(chunks(G.nodes(), int(G.order() / node_divisor)))
    num_chunks = len(node_chunks)
    bt_sc = p.starmap(
        nx.betweenness_centrality_subset,
        zip(
            [G] * num_chunks,
            node_chunks,
            [list(G)] * num_chunks,
            [True] * num_chunks,
            ['distance'] * num_chunks
        ),
    )

    # Reduce the partial solutions
    bt_c = bt_sc[0]
    for bt in bt_sc[1:]:
        for n in bt:
            bt_c[n] += bt[n]
    return bt_c
