"""Post process all adjacency lists by filtering, saving to pickle and saving NetworkX graph."""

import argparse
import multiprocessing
import os

import pandas as pd
from joblib import Parallel, delayed
from tqdm import tqdm

from _post_process import post_process_adj_list


def run_post_process_adj_list(fn, q):
    """
    Helper function to deal with empty data.

    :param fn: filename of pyscenic output
    :param q: the quantile threshold
    """

    try:
        # Run post-processing
        post_process_adj_list(fn, q)
    except pd.errors.EmptyDataError:
        print(f"The file {fn} is empty!")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Process all graph adjacency lists computed in the specified data folder. "
        "Useful, when encountered errors during GRN inference and want to just run "
        "the last part of the pipeline."
    )
    parser.add_argument(
        "-f",
        "--fn",
        type=str,
        help="The path to folder containing the data.",
        required=True,
    )
    parser.add_argument(
        "-q",
        "--q_thresh",
        type=float,
        help="The quantile threshold used to filter out unimportant/false connections.",
        default=0.95,
    )
    args = parser.parse_args()

    PAT_FOLDERS = [
        folder
        for folder in os.listdir(args.fn)
        if folder not in ["cell_types", ".ipynb_checkpoints"]
    ]
    CELL_TYPES_FOLDERS = os.listdir(os.path.join(args.fn, "cell_types"))

    DATA_FILES = [
        os.path.join(args.fn, folder, "data", "grnboost2", f)
        for folder in PAT_FOLDERS
        for f in os.listdir(os.path.join(args.fn, folder, "data", "grnboost2"))
        if f.endswith("_cor.tsv")
        or f.endswith("_TF_cor.tsv")
        or f.endswith("_TF_ctx.tsv")
    ] + [
        os.path.join(args.fn, "cell_types", folder, "data", "grnboost2", f)
        for folder in CELL_TYPES_FOLDERS
        for f in os.listdir(
            os.path.join(args.fn, "cell_types", folder, "data", "grnboost2")
        )
        if f.endswith("_cor.tsv")
        or f.endswith("_TF_cor.tsv")
        or f.endswith("_TF_ctx.tsv")
    ]

    NUM_CORES = max(1, multiprocessing.cpu_count() - 10)

    Parallel(n_jobs=NUM_CORES)(
        delayed(run_post_process_adj_list)(fn, args.q_thresh) for fn in tqdm(DATA_FILES)
    )

    print("Success!")
