import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx

from functools import reduce  # for aggregate functions

from sklearn.preprocessing import MinMaxScaler

scale = lambda x, min_y, max_y: list(MinMaxScaler(feature_range=(min_y, max_y)).fit_transform(np.expand_dims(np.array(x), axis=1))[:, 0])


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
    from functools import reduce
    
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


def netgraph_community_layout(G, node_to_community, community_scale=3., node_scale=1., seed=42):
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


def _get_community_positions(G, node_to_community, community_scale, seed):
    """
    Compute a centroid position for each community.
    """
    
    # create a weighted graph, in which each node corresponds to a community,
    # and each edge weight to the number of edges between communities
    between_community_edges = _find_between_community_edges(G, node_to_community)

    communities = set(node_to_community.values())
    hypergraph = nx.DiGraph()
    hypergraph.add_nodes_from(communities)
    for (ci, cj), edges in between_community_edges.items():
        hypergraph.add_edge(ci, cj, weight=len(edges))

    # find layout for communities
    pos_communities = nx.spring_layout(hypergraph, scale=community_scale, seed=seed)

    # set node positions to position of community
    pos = dict()
    for node, community in node_to_community.items():
        pos[node] = pos_communities[community]

    return pos

def _find_between_community_edges(G, node_to_community):
    """Convert the graph into a weighted network of communities."""
    edges = dict()

    for (ni, nj) in G.edges():
        ci = node_to_community[ni]
        cj = node_to_community[nj]

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