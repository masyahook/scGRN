import argparse
import os

from functools import reduce  # for aggregate functions

import pandas as pd
import networkx as nx


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Process the list of adjacencies starting from saving it to pickle format, filtering the unimportant/false connections, and creating/saving a Networkx graph.')
    parser.add_argument('-f', '--fn', type=str, help='The full path to filename', required=True)
    parser.add_argument('-q', '--q_thresh', type=float, help='The quantile threshold used to filter out unimportant/false connections.', default=0.8)
    args = parser.parse_args()
    
    # Defining file and path names 
    curr_dir = os.path.dirname(args.fn)
    fn = args.fn[args.fn.rfind('/') + 1:]
    short_fn = fn.replace('.tsv', '')
    os.chdir(curr_dir)
    os.makedirs('pickle', exist_ok=True)
    os.makedirs('nx_graph', exist_ok=True)
    
    ####### 1st step #######
    # Process ctx-format file
    if 'ctx' in fn:
        df = pd.read_csv(args.fn, sep='\t', index_col=[0, 1], header=[0, 1], skipinitialspace=True)
        all_df = pd.read_csv(os.path.join(curr_dir, fn.replace('ctx', 'cor')), sep='\t')
        
        df[('Enrichment', 'Context')] = df[('Enrichment', 'Context')].apply(lambda s: eval(s))
        df[('Enrichment', 'TargetGenes')] = df[('Enrichment', 'TargetGenes')].apply(lambda s: eval(s))

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

        out_df = pd.DataFrame(tf_target_dict).merge(all_df, how='left')
    else:
        out_df = pd.read_csv(args.fn, sep='\t')
    
    # Saving
    out_df.to_pickle(f'pickle/{short_fn}.pickle')
    
    ####### 2nd step #######
    adj_list = out_df
    
    # Filtering
    thresh = adj_list['importance'].quantile(args.q_thresh)
    filtered_adj_list = adj_list.loc[lambda x: x.importance > thresh]
    
    # Saving
    filtered_adj_list.to_pickle(f'pickle/{short_fn}_filtered.pickle')
    
    ####### 3rd step #######
    # Creating graphs
    filtered_graph = nx.from_pandas_edgelist(filtered_adj_list, 'TF', 'target', ['importance', 'rho'], create_using=nx.DiGraph)
    graph = nx.from_pandas_edgelist(adj_list, 'TF', 'target', ['importance', 'rho'], create_using=nx.DiGraph)
    
    # Saving
    nx.write_gpickle(graph, f'nx_graph/{short_fn}.gpickle')
    nx.write_gpickle(filtered_graph, f'nx_graph/{short_fn}_filtered.gpickle')
    