import os
import typing

import numpy as np
import pandas as pd
import networkx as nx
from pandas.io.formats.style import Styler

from sklearn.preprocessing import MinMaxScaler


def scale(
        x: typing.Iterable[float],
        min_x: float,
        max_x: float
) -> list:
    """
    Scale the input array `x` between `min_x` and `max_x`.

    :param x: Iterable array to scale
    :param min_x: The left boundary of scaled array
    :param max_x: The right boundary of scaled array

    :return: The scaled array `x`
    """

    return list(
        MinMaxScaler(
            feature_range=(min_x, max_x)
        ).fit_transform(
            np.expand_dims(np.array(x), axis=1)
        )[:, 0]
    )


def scale_int(
        x: typing.Iterable[float],
        min_x: float,
        max_x: float
) -> list:
    """
    Scale the input array `x` between `min_x` and `max_x` with additional rounding of output values.

    :param x: Iterable array to scale
    :param min_x: The left boundary of scaled array
    :param max_x: The right boundary of scaled array

    :return: The scaled array `x` with only integer values
    """

    return [
        int(el) for el in list(
            MinMaxScaler(
                feature_range=(min_x, max_x)
            ).fit_transform(
                np.expand_dims(np.array(x), axis=1)
            )[:, 0]
        )
    ]


def is_non_empty(fn: str) -> bool:
    """
    Check if non-empty file exists.

    :param fn: The filepath to the file

    :return: True of non-empty file exists, False otherwise
    """

    return os.path.exists(fn) and (os.stat(fn).st_size != 0)


def style_bool_df(df: pd.DataFrame) -> Styler:
    """
    Style (HTML coloring) the Pandas dataframe consisting of True (in green), False (in red) and NaN (in yellow) values.

    :param df: Pandas dataframe to style

    :return Styled dataframe
    """

    check, missing, cross = u'\u2713', '?', u'\u2715'  # green, red, yellow

    return df.applymap(  # replace boolean values with check, cross and question marks
        lambda x: check if x is True else missing if x is False else cross
    ).style.apply(
        lambda x: [
            "background-color: green"
            if v == check else "background-color: red" if v == missing
            else "background-color: yellow" for v in x
        ], axis=1
    )


def save_pickle(f: typing.Any, fn: str):
    """
    Save object as a pickle file (usually used for dicts).

    :param f: File object
    :param fn: File name
    """

    import pickle

    with open(fn, 'wb') as fo:
        pickle.dump(f, fo)


def load_pickle(fn: str) -> typing.Any:
    """
    Load object from pickle file (usually used for dicts).

    :param fn: File name

    return: Loaded object
    """

    import pickle

    with open(fn, 'rb') as fo:
        f = pickle.load(fo)

    return f


def _process_ndex_net(
        raw_G: nx.MultiGraph,
        net_fn: str,
        remove_single_nodes: bool = True
) -> nx.DiGraph:
    """
    Process the loaded `.cx` file, convert it to NetworkX object and adapt it to the pipeline. Each raw network is
        scanned for type identifiers and subsequently processed.

    :param raw_G: The raw NetworkX graph
    :param net_fn: The `.cx` file path
    :param remove_single_nodes: True if remove isolated single nodes in the graph

    :return Processed NetworkX graph
    """

    import re
    from functools import reduce  # for aggregate functions

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
            new_name = name + node_duplicates[name] * '*'
            tmp_nodes[node_i] = {'old_name': name, 'new_name': new_name}
            attrs = {'type': node_info['bel_function_type'], 'full_name': node_info['name']}
            G.add_node(new_name, **attrs)
        for st_i, end_i, edge_info in raw_G.edges(data=True):
            st_name, end_name = tmp_nodes[st_i]['new_name'], tmp_nodes[end_i]['new_name']
            inter = edge_info['interaction']
            attrs = {
                'regulation': 'up' if 'increases' in inter.lower() else 'down' if 'decreases' in inter.lower() else None,
                'info': inter}
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
            new_name = name + node_duplicates[name] * '*'
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
            new_name = name + node_duplicates[name] * '*'
            tmp_nodes[node_i] = {'old_name': name, 'new_name': new_name}
            attrs = {'node_type': node_info['type'], 'full_name': None}
            G.add_node(new_name, **attrs)
        for st_i, end_i, edge_info in raw_G.edges(data=True):
            inter = edge_info['interaction']
            attrs = {'regulation': 'up' if 'up' in inter.lower() else 'down' if 'down' in inter.lower() else None,
                     'info': edge_info['mechanism'] if 'mechanism' in edge_info else None}
            st_name, end_name = tmp_nodes[st_i]['new_name'], tmp_nodes[end_i]['new_name']
            G.add_edge(st_name, end_name, **attrs)

    else:

        raise NameError(f'No identified source for "{net_fn}"')

    if remove_single_nodes:
        single_node_sets = list(filter(lambda x: len(x) == 1, nx.weakly_connected_components(G)))
        single_nodes = reduce(lambda x, y: x.union(y), single_node_sets) if len(single_node_sets) > 0 else set()
        G.remove_nodes_from(single_nodes)

    return G


def load_ndex_net(fn: str, remove_single_nodes: bool = True):
    """
    Load NDEx network from file

    :param fn: The file path to `.cx` file
    :param remove_single_nodes: True if remove isolated single nodes in the graph

    :return The processed NetworkX graph of NDEx network
    """

    import ndex2

    raw_G = ndex2.create_nice_cx_from_file(fn).to_networkx(mode='default')

    return _process_ndex_net(raw_G, fn, remove_single_nodes)
