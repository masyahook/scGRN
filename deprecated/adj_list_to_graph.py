import argparse
import os
import re

from functools import reduce  # for aggregate functions

import pandas as pd
import networkx as nx

import warnings


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Create a save a GRN graph based in pickle format.')
    parser.add_argument('-f', '--fn', type=str, help='The full path to input pickle file with list of adjacencies (either ending with _cor.pickle, _TF_cor.pickle, _TF_ctx.pickle)', required=True)
    args = parser.parse_args()
    
    # Creaing a nx_graph folder within the patient-method directory
    curr_dir = os.path.dirname(os.path.dirname(args.fn))
    os.chdir(curr_dir)
    os.makedirs('nx_graph', exist_ok=True)
    short_fn = args.fn[args.fn.rfind('/') + 1:].replace('.pickle', '')
    
    # Reading the list of adjacencies
    try:
        q_thresh_suffix = re.findall(r'(?<=\_)0\_\d+(?=\.)', args.fn)[0]  # getting q_thresh suffix
    except:
        warnings.warn(f'Cannot infer q_thresh from the input filename: {args.fn}. Saving files without q_thresh..')
        q_thresh_suffix = ''
    adj_list = pd.read_pickle(args.fn)
    filtered_adj_list = pd.read_pickle(f'pickle/{short_fn}_filtered_{q_thresh_suffix}.pickle')
    
    # Adding distance metric as 1 / importance to each graph
    adj_list['distance'] = 1 / adj_list['importance']
    filtered_adj_list['distance'] = 1 / filtered_adj_list['importance']
    
    # Creating graphs
    graph = nx.from_pandas_edgelist(adj_list, 'TF', 'target', ['importance', 'rho', 'distance'], create_using=nx.DiGraph)
    filtered_graph = nx.from_pandas_edgelist(filtered_adj_list, 'TF', 'target', ['importance', 'rho', 'distance'], create_using=nx.DiGraph)
    
    # Saving
    nx.write_gpickle(graph, f'nx_graph/{short_fn}.gpickle')
    nx.write_gpickle(filtered_graph, f'nx_graph/{short_fn}_filtered_{q_thresh_suffix}.gpickle')
    