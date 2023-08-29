""" Convert all saved .tsv VIPER files to .pickle format, for faster IO operations later """

import argparse
import multiprocessing
import os
import re

import pandas as pd
from joblib import Parallel, delayed
from tqdm import tqdm


def _safe_listdir(path):
    """Safe way to list directory, in case it is not existent."""

    try:
        return os.listdir(path)
    except FileNotFoundError as e:
        print(
            f"The path {path} does not exist. Did you run the `regulon_pipeline.sh` for that data?"
            f"Error:\n{e}"
        )
        return []


def to_pickle(fn):
    """
    Helper function to run the script
    """

    # Dealing with path names
    curr_dir = os.path.dirname(fn)
    short_pickle_fn = fn[fn.rfind("/") + 1 :].replace(".tsv", ".pickle")
    full_pickle_fn = os.path.join(curr_dir, "pickle", short_pickle_fn)

    # Running
    try:
        df = pd.read_csv(fn, sep="\t")
        df.to_pickle(full_pickle_fn)
    except BaseException as e:
        print(f"Caught an error for {fn}:\n{e}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Save the output of VIPER to pickle format, to subsequently speed up loading this file in the network analysis."
    )
    parser.add_argument(
        "-o",
        "--outdir",
        type=str,
        help="The output directory folder",
        default="/gpfs/projects/bsc08/shared_projects/scGRN_analysis/Data_home/res/covid_19",
    )
    parser.add_argument(
        "-r",
        "--regulon",
        type=str,
        help="The regulon used",
        required=True,
        default="pyscenic",
    )
    args = parser.parse_args()

    # Saving pat-specific VIPER data
    DATA_PATH = args.outdir

    # Defining path names, creating pickle folders
    pattern = re.compile("^C\d{2,3}$")  # either C51, C52, or C141, C142, etc
    pat_spec_data_folders = [
        folder for folder in os.listdir(DATA_PATH) if pattern.match(folder)
    ]
    _ = [
        os.makedirs(
            os.path.join(DATA_PATH, pat, "data", "Seurat", "regulon", "pickle"),
            exist_ok=True,
        )
        for pat in pat_spec_data_folders
    ]
    data_files = [
        os.path.join(DATA_PATH, pat, "data", "Seurat", "regulon", f)
        for pat in pat_spec_data_folders
        for f in _safe_listdir(
            os.path.join(DATA_PATH, pat, "data", "Seurat", "regulon")
        )
        if f.startswith(args.regulon) and f.endswith(".tsv")
    ]

    # Running
    num_cores = max(1, multiprocessing.cpu_count() - 5)

    Parallel(n_jobs=num_cores)(delayed(to_pickle)(fn) for fn in tqdm(data_files))

    # Saving cell type specific VIPER data
    DATA_PATH = os.path.join(args.outdir, "cell_types")

    # Defining path names, creating pickle folders
    type_spec_data_folders = [
        folder
        for folder in next(os.walk(DATA_PATH))[1]
        if folder not in [".ipynb_checkpoints"]
    ]
    _ = [
        os.makedirs(
            os.path.join(DATA_PATH, t, "data", "Seurat", "regulon", "pickle"),
            exist_ok=True,
        )
        for t in type_spec_data_folders
    ]
    data_files = [
        os.path.join(DATA_PATH, t, "data", "Seurat", "regulon", f)
        for t in type_spec_data_folders
        for f in _safe_listdir(os.path.join(DATA_PATH, t, "data", "Seurat", "regulon"))
        if f.startswith(args.regulon) and f.endswith(".tsv")
    ]

    # Running
    num_cores = max(1, multiprocessing.cpu_count() - 5)

    Parallel(n_jobs=num_cores)(delayed(to_pickle)(fn) for fn in tqdm(data_files))

    print("Success!")
