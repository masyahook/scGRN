import argparse
import os

import pandas as pd
from tqdm import tqdm
import multiprocessing
from joblib import Parallel, delayed


""" Convert all saved .tsv VIPER files to .pickle format, for faster IO operations later """


def to_pickle(fn):
    """
    Helper function to run the script
    """
    
    # Dealing with path names
    curr_dir = os.path.dirname(fn)
    short_pickle_fn = fn[fn.rfind('/') + 1:].replace('.tsv', '.pickle')
    full_pickle_fn = os.path.join(curr_dir, 'pickle', short_pickle_fn)
    
    # Running
    try:
        df = pd.read_csv(fn, sep='\t')
        df.to_pickle(full_pickle_fn)
    except BaseException as e:
        print(f'Caught an error for {fn}:\n{e}')


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Save the output of VIPER to pickle format, to subsequently speed up loading this file in the network analysis.')
    parser.add_argument('-o', '--outdir', type=str, help='The output directory folder', default='/gpfs/projects/bsc08/shared_projects/scGRN_analysis/Data_home/res/covid_19')
    parser.add_argument('-r', '--regulon', type=str, help='The regulon used', required=True, default='pyscenic')
    args = parser.parse_args()
    
    # Saving pat-specific VIPER data
    DATA_PATH = '/gpfs/projects/bsc08/bsc08890/res/covid_19'

    # Defining path names, creating pickle folders
    pat_spec_data_folders = [folder for folder in os.listdir(DATA_PATH) if folder not in ['cell_types', '.ipynb_checkpoints']]
    _ = [os.makedirs(os.path.join(DATA_PATH, pat, 'data', 'Seurat', 'regulon', 'pickle'), exist_ok=True) for pat in pat_spec_data_folders]
    data_files = [
        os.path.join(DATA_PATH, pat, 'data', 'Seurat', f) for pat in pat_spec_data_folders for f in os.listdir(os.path.join(DATA_PATH, pat, 'data', 'Seurat')) if f.startswith(args.regulon) and f.endswith('.tsv')
    ]

    # Running
    num_cores = max(1, multiprocessing.cpu_count() - 5)

    Parallel(n_jobs=num_cores)(delayed(to_pickle)(fn) for fn in tqdm(data_files))
        
    # Saving cell type-specific VIPER data
    DATA_PATH = os.path.join(args.outdir, 'cell_types')
    
    # Defining path names, creating pickle folders
    type_spec_data_folders = [folder for folder in os.listdir(DATA_PATH) if folder not in ['.ipynb_checkpoints']]
    _ = [os.makedirs(os.path.join(DATA_PATH, t, 'data', 'Seurat', 'pickle'), exist_ok=True) for t in type_spec_data_folders]
    data_files = [
        os.path.join(DATA_PATH, t, 'data', 'Seurat', f) for t in type_spec_data_folders for f in os.listdir(os.path.join(DATA_PATH, t, 'data', 'Seurat')) if f.startswith(args.regulon) and f.endswith('.tsv')
    ]

    # Running
    num_cores = max(1, multiprocessing.cpu_count() - 5)

    Parallel(n_jobs=num_cores)(delayed(to_pickle)(fn) for fn in tqdm(data_files))

    print('Success!')
