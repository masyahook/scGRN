import argparse
import sys
import os

import pandas as pd

# Setting working directory as home
home_dir = os.path.expanduser('~')
os.chdir(os.path.expanduser('~'))

# Adding home directory to the PYTHONPATH
sys.path.insert(0, home_dir)

from src.func import process_communities

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Run community analysis on a specific GRN.')
    parser.add_argument('-d', '--data', default='raw_data', type=str, help="The name of the data - starts with 'raw_data' and could end with cell type", required=True)
    parser.add_argument('-p', '--patient', type=str, default=None, help='Patient ID')
    parser.add_argument('-a', '--algo', type=str, default='leiden', help="The name of the clustering algorithm - either 'leiden' or 'louvain'")
    parser.add_argument('-ib', '--if_betweenness', type=bool, default=True, help='True if use betweeness metric to compute centrality, False if use closeness centrality')
    parser.add_argument('-la', '--limit_anno_until', type=int, default=50, help='The upper limit of genes that we will base the cluster annotation on')
    parser.add_argument('-k', '--k_centrality', type=int, default=5000, help='The upper number of nodes to compute the centrality on - approximates the centrality and speeds up the process')
    parser.add_argument('-il', '--top_intercom_links', type=int, default=20, help='The number of inter-community links to store in the output info dataframe')
    parser.add_argument('-of', '--other_functions_until', type=int, default=20, help='The number of other functions to include in the final dataframe.')
    parser.add_argument('-nl', '--save_top_new_found_cluster_links', type=int, default=20, help='The number of newly found gene interactions in each functional cluster to save to final dataframe.')
    parser.add_argument('-s', '--seed', type=int, default=42, help='The random seed to use')
    args = parser.parse_args()
    
    # Running process communities script
    process_communities(
        pat=args.patient, data=args.data, algo=args.algo, if_betweenness=args.if_betweenness, 
        limit_anno_until=args.limit_anno_until, k=args.k_centrality, 
        save_top_intercommunity_links_until=args.top_intercom_links, 
        other_functions_until=args.other_functions_until,
        save_top_new_found_cluster_links=args.save_top_new_found_cluster_links,
        seed=args.seed
    )
    
    print('Success!')