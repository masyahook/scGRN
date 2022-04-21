import argparse
import os

from functools import reduce  # for aggregate functions

import pandas as pd
from tqdm import tqdm
import multiprocessing
from joblib import Parallel, delayed

def to_pickle(fn):
    """
    Helper function to run the script
    """
    
    # Dealing with path names
    curr_dir = os.path.dirname(fn)
    short_pickle_fn = fn[fn.rfind('/') + 1:].replace('.tsv', '.pickle')
    full_pickle_fn = os.path.join(curr_dir, 'pickle', short_pickle_fn)
    
    # Running
    df = pd.read_csv(fn, sep='\t')
    df.to_pickle(full_pickle_fn)


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Save the output of DoRothEA to pickle format, to subsequently speed up loading this file in the network analysis.')
    parser.add_argument('-a', '--aggregated', type=bool, help='The type of data, cell type aggregated or not', required=True)
    args = parser.parse_args()
    
    if not args.aggregated:
    
        # Setting up the path to all the data
        _DATA_PATH = '/gpfs/projects/bsc08/bsc08890/res/covid_19'

        # Defining path names, creating pickle folders
        pat_spec_data_folders = [folder for folder in os.listdir(_DATA_PATH) if folder not in ['cell_types', '.ipynb_checkpoints']]
        _ = [os.makedirs(os.path.join(_DATA_PATH, pat, 'data', 'Seurat', 'pickle'), exist_ok=True) for pat in pat_spec_data_folders]
        data_files = [
            os.path.join(_DATA_PATH, pat, 'data', 'Seurat', f) for pat in pat_spec_data_folders for f in os.listdir(os.path.join(_DATA_PATH, pat, 'data', 'Seurat')) if f.startswith('dorothea') and f.endswith('.tsv')
        ]

        # Running
        num_cores = max(1, multiprocessing.cpu_count() - 10)

        Parallel(n_jobs=num_cores)(delayed(to_pickle)(fn) for fn in tqdm(data_files))
        
    else:
        
        # Setting up the path to all the data
        _DATA_PATH = '/gpfs/projects/bsc08/bsc08890/res/covid_19/cell_types'
        
        # Defining path names, creating pickle folders
        type_spec_data_folders = [folder for folder in os.listdir(_DATA_PATH) if folder not in ['.ipynb_checkpoints']]
        _ = [os.makedirs(os.path.join(_DATA_PATH, t, 'data', 'Seurat', 'pickle'), exist_ok=True) for t in type_spec_data_folders]
        data_files = [
            os.path.join(_DATA_PATH, t, 'data', 'Seurat', f) for t in type_spec_data_folders for f in os.listdir(os.path.join(_DATA_PATH, t, 'data', 'Seurat')) if f.startswith('dorothea') and f.endswith('.tsv')
        ]

        # Running
        num_cores = max(1, multiprocessing.cpu_count() - 10)

        Parallel(n_jobs=num_cores)(delayed(to_pickle)(fn) for fn in tqdm(data_files))

    print('Success!')