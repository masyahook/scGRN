import argparse
import os

from functools import reduce  # for aggregate functions

import pandas as pd
import networkx as nx
from tqdm import tqdm
import multiprocessing
from joblib import Parallel, delayed

from func import post_process_adj_list


def run_post_process_adj_list(fn, q):
    """
    Helper function to deal with empty data.
    """
    
    try:
        # Run post-processing
        post_process_adj_list(fn, q)
    except pd.errors.EmptyDataError:
        print(f'The file {fn} is empty!')


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Process the list of adjacencies starting from saving it to pickle format, filtering the unimportant/false connections, and creating/saving a Networkx graph.')
    parser.add_argument('-q', '--q_thresh', type=float, help='The quantile threshold used to filter out unimportant/false connections.', default=0.95)
    args = parser.parse_args()
    
    # Setting up the path to all the data
    _PROJ_PATH = '/gpfs/projects/bsc08/bsc08890/res/covid_19'
    
    pat_spec_data_folders = [folder for folder in os.listdir(_PROJ_PATH) if folder not in ['cell_types', '.ipynb_checkpoints']]
    cell_agg_data_folders = os.listdir(os.path.join(_PROJ_PATH, 'cell_types'))
    
    data_files = [
        os.path.join(_PROJ_PATH, folder, 'data', 'grnboost2', f) for folder in pat_spec_data_folders for f in os.listdir(os.path.join(_PROJ_PATH, folder, 'data', 'grnboost2')) if f.endswith('_cor.tsv') or f.endswith('_TF_cor.tsv') or f.endswith('_TF_ctx.tsv')
    ] + [
        os.path.join(_PROJ_PATH, 'cell_types', folder, 'data', 'grnboost2', f) for folder in cell_agg_data_folders for f in os.listdir(os.path.join(_PROJ_PATH, 'cell_types', folder, 'data', 'grnboost2')) if f.endswith('_cor.tsv') or f.endswith('_TF_cor.tsv') or f.endswith('_TF_ctx.tsv')
    ]
    
    num_cores = max(1, multiprocessing.cpu_count() - 10)
    
    Parallel(n_jobs=num_cores)(delayed(run_post_process_adj_list)(fn, args.q_thresh) for fn in tqdm(data_files))
    
    print('Success!')
    