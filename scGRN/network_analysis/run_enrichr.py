import argparse
import os

# Setting working directory as home
home_dir = os.path.expanduser('~')
os.chdir(os.path.expanduser('~'))

from ._enrichment import run_enrichr

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Running EnrichR.')
    parser.add_argument('-d', '--data', type=str, help='The data to compute on', required=True)
    parser.add_argument('--is_community', help='Whether to compute community enrichments', action='store_true')
    parser.add_argument('--is_pat_group', dest='is_community', action='store_false')
    parser.add_argument('--on_targets', help='Whether to compute community enrichments', action='store_true')
    parser.add_argument('--on_tfs', dest='on_targets', action='store_false')
    parser.add_argument('--fixed_tf', type=str, help='The fixed TF to fix the analysis on while in on_targets mode', default=None)
    parser.add_argument('-p', '--positive', type=bool, help='Whether to use positive markers only', default=True)
    parser.add_argument('-gt', '--group_types', type=str, help='The group type to compute on: either "all", "M_C", "S_C", or "S_M"', default='all')
    parser.add_argument('-dt', '--data_type', type=str, help='The data type to compute on: either "all" or "ctx"', default='all')
    parser.add_argument('-n', '--n', type=int, help='The top number of community genes', default=50)
    parser.add_argument('-l', '--lib', type=str, help='The library to use in EnrichR', default='MSigDB_Hallmark_2020')  # 'BioPlanet_2019', # 'Reactome_2016'
    args = parser.parse_args()
    
    run_enrichr(
        args.data, is_communities=args.is_community, group_types=args.group_types,
        on_targets=args.on_targets, is_positive_markers=args.positive, choose_fixed_tf=args.fixed_tf,
        data_type=args.data_type, top_n=args.n, algo='leiden', enrichr_library=args.lib
    )