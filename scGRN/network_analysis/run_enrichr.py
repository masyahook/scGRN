"""Script to run EnrichR analysis."""
import argparse
import os

# Setting working directory as home
home_dir = os.path.expanduser('~')
os.chdir(os.path.expanduser('~'))

from _enrichment import run_enrichr

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Running EnrichR.')
    parser.add_argument('-i', '--in_path', type=str, help='The path to input data where gene sets are stored', required=True)
    parser.add_argument('-g', '--gene_col', type=str, help='The column name that stores gene names', required=True)
    parser.add_argument('-o', '--out_path', type=str, help='The path to output data to store, if not specified the result will be save in the same folder')
    parser.add_argument('-c', '--group_col', type=str, help='The column names that store group names')
    parser.add_argument('-e', '--enrichr_library', type=str, help='The EnrichR library to use for enrichment analysis.', default='Reactome_2016')
    parser.add_argument('-q', '--query', type=str, help='The The query that can be used to select a subset of original dataset.', default=None)
    parser.add_argument('-n', '--top_n', type=str, help='Select top_n genes for enrichment analysis, applies for community datasets', default=50)
    args = parser.parse_args()

    if args.out_path is None:
        args.out_path = os.path.join(
            os.path.dirname(args.in_path), 
            f'enrichr_{os.path.basename(args.in_path)}'
        )
    
    run_enrichr(
        in_path=args.in_path, gene_col=args.gene_col, out_path=args.out_path,
        group_col=args.group_col, enrichr_library=args.enrichr_library,
        query=args.query, top_n=args.top_n
    )
