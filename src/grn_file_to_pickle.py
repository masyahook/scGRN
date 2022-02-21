import argparse
import os

from functools import reduce  # for aggregate functions

import pandas as pd


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Save the output of pyScenic ctx step to pickle format, to subsequently speed up loading this file in the network anlaysis.')
    parser.add_argument('-f', '--fn', type=str, help='The full path to filename', required=True)
    args = parser.parse_args()
    
    # Creaing a pickle folder within the patient-method directory
    curr_dir = os.path.dirname(args.fn)
    fn = args.fn[args.fn.rfind('/') + 1:]
    os.chdir(curr_dir)
    os.makedirs('pickle', exist_ok=True)
    
    if 'ctx' in fn:
        df = pd.read_csv(args.fn, sep='\t', index_col=[0, 1], header=[0, 1], skipinitialspace=True)
        TF_df = pd.read_csv(os.path.join(curr_dir, fn.replace('ctx', 'cor')), sep='\t')
        
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

        out_df = pd.DataFrame(tf_target_dict).merge(TF_df, how='left')
    else:
        out_df = pd.read_csv(args.fn, sep='\t')
    
    # Saving pickle file in the 'pickle' subdirectory
    out_df.to_pickle('pickle/' + args.fn[args.fn.rfind('/') + 1:].replace('.tsv', '.pickle'))