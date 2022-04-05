import argparse
import sys
import os

import pandas as pd

from distributed import Client, LocalCluster
from arboreto.utils import load_tf_names
from arboreto.algo import grnboost2, genie3

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Run GRN inference by using GRNBoost2/GENIE3.')
    parser.add_argument('-m', '--method', default='grnboost2', type=str, help='The method used for GRN inference', required=True)
    parser.add_argument('-i', '--in_fn', default='res/PilotWorkflow/data/norm_data.tsv', type=str, help='The filepath to the input expression matrix', required=True)
    parser.add_argument('-o', '--out_fn', type=str, help='The filepath to the output list of link importances', required=True)
    parser.add_argument('-t', '--tf_fn', type=str, default=None, help='The filepath to the list of transcription factors')
    parser.add_argument('-n', '--num_workers', type=int, default=18, help='The number of workers for parallel execution')
    args = parser.parse_args()
    
    method = eval(args.method)

    # Reading the DataFrame with expression values
    ex_df = pd.read_csv(args.in_fn, sep='\t').T
             
    # Reading the list of TFs from the file
    if args.tf_fn is not None:
        tf_names = pd.read_csv(args.tf_fn, header=None)[0].to_list()
    else:
        tf_names = None
        
    # instantiate a custom Dask distributed Client
    local_cluster = LocalCluster(n_workers=args.num_workers, threads_per_worker=1)
    client = Client(local_cluster)

    # compute the GRN
    network = method(expression_data=ex_df,
                     tf_names=tf_names,
                     verbose=True,
                     client_or_address=client)
    
    client.close()
    local_cluster.close()

    # write the GRN to file
    network.to_csv(args.out_fn, sep='\t', index=False)