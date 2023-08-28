import argparse
import sys
import os

import pandas as pd

from distributed import Client, LocalCluster
from arboreto.utils import load_tf_names
from arboreto.algo import grnboost2, genie3

if __name__ == '__main__':
    
    # Setting working directory as home
    home_dir = os.path.expanduser('~')
    os.chdir(os.path.expanduser('~'))
    
    parser = argparse.ArgumentParser(description='Run GRN inference by using GRNBoost2.')
    parser.add_argument('-i', '--in_fn', default='res/PilotWorkflow/data/norm_data.tsv', type=str, 
                        help='The filepath to the input expression matrix', required=True)
    parser.add_argument('-o', '--out_fn', type=str, help='The filepath to the output list of link importances', required=True)
    parser.add_argument('-n', '--num_workers', type=int, default=18, help='The number of workers for parallel execution')
    args = parser.parse_args()

    # Reading the DataFrame with expression values
    ex_df = pd.read_csv(args.in_fn, sep='\t').T
        
    # instantiate a custom Dask distributed Client
    local_cluster = LocalCluster(n_workers=args.num_workers, threads_per_worker=1)
    client = Client(local_cluster)

    # compute the GRN
    network = grnboost2(expression_data=ex_df,
                        verbose=True,
                        client_or_address=client)
    
    client.close()
    local_cluster.close()

    # write the GRN to file
    network.to_csv(args.out_fn, sep='\t', index=False)