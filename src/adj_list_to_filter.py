import argparse
import os

import pandas as pd


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Filter the list of adjacencies depending on the -q parameter (quantile level).')
    parser.add_argument('-f', '--fn', type=str, help='The full path to the pickle-formatted list of adjacencies', required=True)
    parser.add_argument('-q', '--q_thresh', type=float, help='The quantile threshold used to filter out unimportant/false connections.', default=0.8)
    args = parser.parse_args()
    
    # Defining path and file names
    curr_dir = os.path.dirname(args.fn)
    os.chdir(curr_dir)
    short_fn = args.fn[args.fn.rfind('/') + 1:].replace('.pickle', '')
    
    # Reading the list of adjacencies
    adj_list = pd.read_pickle(args.fn)
    
    # Filtering
    thresh = adj_list['importance'].quantile(args.q_thresh)
    filtered_adj_list = adj_list.loc[lambda x: x.importance > thresh]
    
    # Saving
    filtered_adj_list.to_pickle(f'{short_fn}_filtered.pickle')
    
    
    