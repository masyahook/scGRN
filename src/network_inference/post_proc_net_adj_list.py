import argparse

from ._post_process import post_process_adj_list


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Process the list of adjacencies starting from saving it to pickle format, filtering the unimportant/false connections, and creating/saving a Networkx graph.')
    parser.add_argument('-f', '--fn', type=str, help='The full path to filename', required=True)
    parser.add_argument('-q', '--q_thresh', type=float, help='The quantile threshold used to filter out unimportant/false connections.', default=0.9)
    args = parser.parse_args()
    
    # Run post-processing
    post_process_adj_list(args.fn, args.q_thresh)
