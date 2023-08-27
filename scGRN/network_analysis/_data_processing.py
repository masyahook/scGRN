"""Data loading / processing for network analysis."""

import itertools
import os
import warnings
from typing import Union

import igraph as ig
import networkx as nx
import numpy as np
import pandas as pd
from joblib import Parallel, delayed

from ..config import _DATA_HOME, _META_FILE
from ..utils import is_non_empty


def get_meta(
    data_home: str = _DATA_HOME,
    meta_file: str = _META_FILE
) -> pd.DataFrame:
    """
    Obtain metadata information about the patients.

    :param data_home: The filepath to the data home folder
    :param meta_file: The filepath to metadata about patients

    :return: The metadata dataframe
    """

    # Getting metadata with file locations
    full_meta = pd.read_csv(meta_file, sep='\t', index_col=0)

    # Filling metadata with the information about present cell types
    for patient_id in full_meta.index:
        cells_metadata = pd.read_csv(
            os.path.join(data_home, patient_id, 'data', 'Seurat', 'cells_metadata.tsv'),
            sep='\t'
        )
        full_meta.loc[patient_id, 'num_cells'] = len(cells_metadata)
        for cell_type, cell_count in cells_metadata['cell_type'].value_counts().iteritems():
            full_meta.loc[patient_id, cell_type] = cell_count

    return full_meta


def get_avail_sc_data(
    data_home: str = _DATA_HOME,
    meta_file: str = _META_FILE
) -> pd.DataFrame:
    """
    Obtain dataframe that contains the information about the scRNA-seq matrix data availability.

    :param data_home: The filepath to the data home folder
    :param meta_file: The filepath to metadata about patients

    :return: Boolean-value dataframe indicating the [presence (True) / absence of file due to failed computation (False)
        / absence of data (NaN)] for the corresponding scRNA-seq matrix data
    """

    # Patient types
    pat_types = ['C', 'M', 'S']

    # Getting metadata
    full_meta = get_meta(data_home, meta_file)

    avail_sc_data = pd.DataFrame(
        columns=['all_data'] + full_meta.columns[3:].to_list(),  # cell-types, e.g. 'T_cells', 'Macrophage', 'all_data'
        index=['all_data'] + pat_types + full_meta.index.to_list()  # patients, e.g. 'C51', 'C141', 'S', 'all_data'
    )

    # Checking the presence of patient data files: e.g. 'C141', 'C51', 'C152'
    for pat in avail_sc_data.index[4:]:
        for d in ['all_data'] + full_meta.loc[pat][3:].dropna().index.to_list():
            d_fn = f'raw_data_{d}' if d != 'all_data' else 'raw_data'
            fn = os.path.join(data_home, pat, 'data', 'Seurat', f'{d_fn}.tsv')
            if is_non_empty(fn):
                avail_sc_data.loc[pat, d] = True
            else:
                avail_sc_data.loc[pat, d] = False

    # Checking the presence of cell type aggregated data files: 'all_data', 'C', 'M', 'S'
    for pat in avail_sc_data.index[:4]:
        # Getting cell types present in the corresponding data
        if pat == 'all_data':
            cell_types = avail_sc_data.columns
        else:
            cell_types = ['all_data'] + \
                         full_meta[full_meta['group'] == pat].iloc[:, 3:].loc[:, lambda x: ~x.isna().all()] \
                         .columns.to_list()
        for d in cell_types:
            d_fn = 'raw_data' if pat == 'all_data' else f'raw_data_{pat}_type'
            fn = os.path.join(data_home, 'cell_types', d, 'data', 'Seurat', f'{d_fn}.tsv')
            if is_non_empty(fn):
                avail_sc_data.loc[pat, d] = True
            else:
                avail_sc_data.loc[pat, d] = False

    return avail_sc_data


def get_avail_pat_sc(
    pat: str,
    as_data_fn: bool = True
) -> list:
    """
    Get a list of available scRNA-seq datasets corresponding to `pat`.

    :param pat: The patient identifier - could be either:
        e.g. 'C51', 'C141' (the patient identifier)
        e.g. 'C', 'M', 'S', 'all_data' (the identifier of aggregated patient data)
    :param as_data_fn: Return list of data file names or just cell type names:
        True - ['raw_data', 'raw_data_Macrophage', 'raw_data_T_cells', ... ]
        False - ['all_data', 'Macrophage', 'T_cells', ... ]

    :return: A list of available data files corresponding to `net_type` and `pat`
    """

    avail_sc = get_avail_sc_data()
    cell_types = avail_sc.loc[pat].dropna().loc[lambda x: x].index
    if as_data_fn:
        return cell_types.map(lambda x: f'raw_data_{x}' if x != 'all_data' else 'raw_data').to_list()
    else:
        return cell_types.to_list()


def get_avail_adj_lists(
    data_home: str = _DATA_HOME,
    meta_file: str = _META_FILE,
    method: str = 'grnboost2',
    filtered: float = 0.95,
) -> dict:
    """
    Obtain dictionary containing the information about GRN adjacency list availability.

    :param data_home: The filepath to the data home folder
    :param meta_file: The filepath to metadata about patients
    :param method: The GRN inference method, could be either 'genie3' or 'grnboost2'
    :param filtered: The quantile threshold

    :return: A dictionary with the following structure:
        {
            'all': pd.DataFrame,  # gene-gene GRN adjacency list availability
            'TF': pd.DataFrame,  # TF-target GRN adjacency list availability
            'ctx': pd.DataFrame  # enriched TF-target GRN adjacency list availability
        }
        Each dataframe has boolean values indicating the [presence (True) / absence of file due to failed computation
            (False)/ absence of data (NaN)] for each GRN adjacency list
    """

    # Formatting filtered parameter
    filtered_suffix = str(filtered).replace('.', '_')

    # Getting info about available scRNA-seq data
    avail_sc_data = get_avail_sc_data(data_home, meta_file)

    avail_adj_lists = dict(
        zip(
            ['all', 'TF', 'ctx'],
            [  # creating 3 separate dataframes
                pd.DataFrame(columns=avail_sc_data.columns, index=avail_sc_data.index),
                pd.DataFrame(columns=avail_sc_data.columns, index=avail_sc_data.index),
                pd.DataFrame(columns=avail_sc_data.columns, index=avail_sc_data.index)
             ]
        )
    )

    for net_type, data in avail_adj_lists.items():
        net_type_suffix = 'cor' if net_type == 'all' else 'TF_cor' if net_type == 'TF' else 'TF_ctx'

        # Checking the presence of patient data files: e.g. 'C141', 'C51', 'C152'
        for pat in data.index[4:]:
            for d in avail_sc_data.loc[pat].dropna().index.to_list():
                d_fn = f'raw_data_{d}' if d != 'all_data' else 'raw_data'
                fn = os.path.join(data_home, pat, 'data', method, f'{d_fn}_{net_type_suffix}.tsv')
                fn_pickled = os.path.join(data_home, pat, 'data', method, 'pickle',
                                          f'{d_fn}_{net_type_suffix}.pickle')
                filtered_fn_pickled = os.path.join(data_home, pat, 'data', method, 'pickle',
                                                   f'{d_fn}_{net_type_suffix}_filtered_{filtered_suffix}.pickle')
                if is_non_empty(fn) and is_non_empty(fn_pickled) and is_non_empty(filtered_fn_pickled):
                    data.loc[pat, d] = True
                else:
                    data.loc[pat, d] = False

        # Checking the presence of cell type aggregated data files: 'all_data', 'C', 'M', 'S'
        for pat in data.index[:4]:
            for d in avail_sc_data.loc[pat].dropna().index.to_list():
                d_fn = 'raw_data' if pat == 'all_data' else f'raw_data_{pat}_type'
                fn = os.path.join(data_home, 'cell_types', d, 'data', method, f'{d_fn}_{net_type_suffix}.tsv')
                fn_pickled = os.path.join(data_home, 'cell_types', d, 'data', method, 'pickle',
                                          f'{d_fn}_{net_type_suffix}.pickle')
                filtered_fn_pickled = os.path.join(data_home, 'cell_types', d, 'data', method, 'pickle',
                                                   f'{d_fn}_{net_type_suffix}_filtered_{filtered_suffix}.pickle')
                if is_non_empty(fn) and is_non_empty(fn_pickled) and is_non_empty(filtered_fn_pickled):
                    data.loc[pat, d] = True
                else:
                    data.loc[pat, d] = False

    return avail_adj_lists


def get_avail_nx_graphs(
    data_home: str = _DATA_HOME,
    meta_file: str = _META_FILE,
    method: str = 'grnboost2',
    filtered: float = 0.95
) -> dict:
    """
    Obtain dictionary containing the information about GRN NetworkX graph availability.

    :param data_home: The filepath to the data home folder
    :param meta_file: The filepath to metadata about patients
    :param method: The GRN inference method, could be either 'genie3' or 'grnboost2'
    :param filtered: The quantile threshold

    :return: A dictionary with the following structure:
        {
            'all': pd.DataFrame,  # gene-gene GRN NetworkX graph availability
            'TF': pd.DataFrame,  # TF-target GRN NetworkX graph availability
            'ctx': pd.DataFrame  # enriched TF-target NetworkX graph availability
        }
        Each dataframe has boolean values indicating the [presence (True) / absence of file due to failed computation
            (False)/ absence of data (NaN)] for each GRN NetworkX graph
    """

    # Formatting filtered parameter
    filtered_suffix = str(filtered).replace('.', '_')

    # Getting info about available scRNA-seq data
    avail_sc_data = get_avail_sc_data(data_home, meta_file)

    avail_nx_graphs = dict(
        zip(
            ['all', 'TF', 'ctx'],
            [  # creating 3 separate dataframes
                pd.DataFrame(columns=avail_sc_data.columns, index=avail_sc_data.index),
                pd.DataFrame(columns=avail_sc_data.columns, index=avail_sc_data.index),
                pd.DataFrame(columns=avail_sc_data.columns, index=avail_sc_data.index)
             ]
        )
    )

    for net_type, data in avail_nx_graphs.items():
        net_type_suffix = 'cor' if net_type == 'all' else 'TF_cor' if net_type == 'TF' else 'TF_ctx'

        # Checking the presence of patient data files: e.g. 'C141', 'C51', 'C152'
        for pat in data.index[4:]:
            for d in avail_sc_data.loc[pat].dropna().index.to_list():
                d_fn = f'raw_data_{d}' if d != 'all_data' else 'raw_data'
                fn = os.path.join(data_home, pat, 'data', method, 'nx_graph', f'{d_fn}_{net_type_suffix}.gpickle')
                fn_filtered = os.path.join(data_home, pat, 'data', method, 'nx_graph',
                                           f'{d_fn}_{net_type_suffix}_filtered_{filtered_suffix}.gpickle')
                if is_non_empty(fn) and is_non_empty(fn_filtered):
                    data.loc[pat, d] = True
                else:
                    data.loc[pat, d] = False

        # Checking the presence of cell type aggregated data files: 'all_data', 'C', 'M', 'S'
        for pat in data.index[:4]:
            for d in avail_sc_data.loc[pat].dropna().index.to_list():
                d_fn = 'raw_data' if pat == 'all_data' else f'raw_data_{pat}_type'
                fn = os.path.join(data_home, 'cell_types', d, 'data', method, 'nx_graph',
                                  f'{d_fn}_{net_type_suffix}.gpickle')
                fn_filtered = os.path.join(data_home, 'cell_types', d, 'data', method, 'nx_graph',
                                           f'{d_fn}_{net_type_suffix}_filtered_{filtered_suffix}.gpickle')
                if is_non_empty(fn) and is_non_empty(fn_filtered):
                    data.loc[pat, d] = True
                else:
                    data.loc[pat, d] = False

    return avail_nx_graphs


def get_avail_pat_nx(
    net_type: str,
    pat: str,
    as_data_fn: bool = True
) -> list:
    """
    Get a list of available GRN NetworkX graphs corresponding to `net_type` and `pat`.

    :param net_type: The type of data - could be either:
        'all' (all gene-gene connections)
        'TF' (TF-target connections)
        'ctx' (enriched TF-target connections)
    :param pat: The patient identifier - could be either:
        e.g. 'C51', 'C141' (the patient identifier)
        e.g. 'C', 'M', 'S', 'all_data' (the identifier of aggregated patient data)
    :param as_data_fn: Return list of data file names or just cell type names:
        True - ['raw_data', 'raw_data_Macrophage', 'raw_data_T_cells', ... ]
        False - ['all_data', 'Macrophage', 'T_cells', ... ]

    :return: A list of available data files corresponding to `net_type` and `pat`
    """

    avail_Gs = get_avail_nx_graphs()
    cell_types = avail_Gs[net_type].loc[pat].dropna().loc[lambda x: x].index
    if as_data_fn:
        return cell_types.map(lambda x: f'raw_data_{x}' if x != 'all_data' else 'raw_data').to_list()
    else:
        return cell_types.to_list()


def get_sc_data(
    cell_type: str,
    pat: str = None,
    data_home: str = _DATA_HOME,
    tolerate_missing: bool = True
) -> Union[pd.DataFrame, None]:
    """
    Load scRNA-seq matrix data.

    :param cell_type: The cell type name of the data - could be either:
        e.g. 'Macrophage', 'T_cells' (the cell type identifier)
        e.g. 'all', 'all_data' (the aggregated data - include all cell types)
        e.g. 'raw_data_Macrophage', 'raw_data_T_cells', 'raw_data' (the data file identifier, i.e. 'raw_data_T_cells'
            corresponds to the same data as 'T_cells')
    :param pat: The patient identifier - could be either:
        e.g. None (include all patients)
        e.g. 'C51', 'C141' (the patient identifier)
        e.g. 'C', 'M', 'S', 'all_patients', 'all_data', 'all' (the identifier of aggregated patient data)
    :param data_home: The filepath to the data home folder
    :param tolerate_missing: True if tolerate missing data file (and output None in this case), False otherwise

    :return Pandas dataframe of corresponding VIPER score matrix
    """

    # Loading data that includes all patients
    if pat is None or 'all' in pat:

        data_home = os.path.join(data_home, 'cell_types')
        data_folder = 'all_data' if cell_type in ['all', 'raw_data'] else cell_type.replace('raw_data_', '')

        fn = os.path.join(
            data_home, data_folder, cell_type, 'Seurat', f'raw_data.tsv'
        )

    # Loading data that includes a specific patient type: 'control', 'moderate' or 'severe'
    elif pat in ['C', 'M', 'S']:

        data_home = os.path.join(data_home, 'cell_types')
        data_folder = 'all_data' if cell_type in ['all', 'raw_data'] else cell_type.replace('raw_data_', '')

        fn = os.path.join(
            data_home, data_folder, cell_type, 'Seurat', f'raw_data_{pat}_type.tsv'
        )

    # Loading patient-specific data
    else:

        # Include all cell types
        if cell_type == 'all_data' or cell_type == 'raw_data' or cell_type == 'all':
            tag = 'raw_data'
        # Specific cell type
        else:
            # Formatting input params for accessing the correct file
            if 'raw_data_' not in cell_type:
                tag = f'raw_data_{cell_type}'
            else:
                tag = cell_type

        fn = os.path.join(
            data_home, pat, 'data', 'Seurat', f'{tag}.tsv'
        )

    try:
        data = pd.read_csv(fn, sep='\t')
        return data
    except FileNotFoundError:
        if tolerate_missing:
            warnings.warn(f'The single-cell data not found for pat="{pat}", cell_type="{cell_type}" '
                          f'(should be at "{fn}"). Returning as `None` instead!')
            return None
        else:
            raise FileNotFoundError(
                f'The single-cell data not found for pat="{pat}", cell_type="{cell_type}" '
                f'(should be at "{fn}")..'
            )


def get_num_cells(pat: str, cell_type: str, meta: pd.DataFrame) -> int:
    """
    Get number of cells present in the corresponding scRNA-seq dataset.
    :param cell_type: The cell type name of the data - could be either:
        e.g. 'Macrophage', 'T_cells' (the cell type identifier)
        e.g. 'all', 'all_data' (the aggregated data - include all cell types)
        e.g. 'raw_data_Macrophage', 'raw_data_T_cells', 'raw_data' (the data file identifier, i.e. 'raw_data_T_cells'
            corresponds to the same data as 'T_cells')
    :param pat: The patient identifier - could be either:
        e.g. None (include all patients)
        e.g. 'C51', 'C141' (the patient identifier)
        e.g. 'C', 'M', 'S', 'all_patients', 'all_data', 'all' (the identifier of aggregated patient data)
    :param meta: The dataframe containing metadata about patients (cell count, patient types, file paths)
    :return: The number of cells
    """
    if pat is None or 'all' in pat:
        curr_meta = meta
    elif pat in ['C', 'M', 'S']:
        curr_meta = meta[meta['group'] == pat]
    else:
        curr_meta = meta.loc[pat]
    if cell_type == 'all' or cell_type == 'all_data' or cell_type == 'raw_data':
        col_key = 'num_cells'
    else:
        col_key = cell_type.replace("raw_data_", "")
    return int(curr_meta[col_key]) if isinstance(curr_meta[col_key], (int, float)) else int(curr_meta[col_key].sum())


def get_viper_mat(
    cell_type: str,
    pat: str = None,
    regulon: str = 'pyscenic',
    data_home: str = _DATA_HOME,
    tolerate_missing: bool = True
) -> Union[pd.DataFrame, None]:
    """
    Load VIPER score TF activity matrix.

    :param cell_type: The cell type name of the data - could be either:
        e.g. 'Macrophage', 'T_cells' (the cell type identifier)
        e.g. 'all', 'all_data' (the aggregated data - include all cell types)
        e.g. 'raw_data_Macrophage', 'raw_data_T_cells', 'raw_data' (the data file identifier, i.e. 'raw_data_T_cells'
            corresponds to the same data as 'T_cells')
    :param pat: The patient identifier - could be either:
        e.g. None (include all patients)
        e.g. 'C51', 'C141' (the patient identifier)
        e.g. 'C', 'M', 'S', 'all_patients', 'all_data', 'all' (the identifier of aggregated patient data)
    :param regulon: The regulon type, could be either 'pyscenic' (inferred by `pyscenic`) or 'dorothea' (obtained from
        `DoRothEA` database)
    :param data_home: The filepath to the data home folder
    :param tolerate_missing: True if tolerate missing data file (and output None in this case), False otherwise

    :return Pandas dataframe of corresponding VIPER score matrix
    """

    # Loading data that includes all patients
    if pat is None or 'all' in pat:

        data_home = os.path.join(data_home, 'cell_types')
        data_folder = 'all_data' if cell_type in ['all', 'raw_data'] else cell_type.replace('raw_data_', '')

        fn = os.path.join(
            data_home, data_folder, 'data', 'Seurat', 'regulon', 'pickle', f'{regulon}_raw_data.pickle'
        )

    # Loading data that includes a specific patient type: 'control', 'moderate' or 'severe'
    elif pat in ['C', 'M', 'S']:

        data_home = os.path.join(data_home, 'cell_types')
        data_folder = 'all_data' if cell_type in ['all', 'raw_data'] else cell_type.replace('raw_data_', '')

        fn = os.path.join(
            data_home, data_folder, 'data', 'Seurat', 'regulon', 'pickle', f'{regulon}_raw_data_{pat}_type.pickle'
        )

    # Loading patient-specific data
    else:

        # Include all cell types
        if cell_type == 'all_data' or cell_type == 'raw_data' or cell_type == 'all':
            tag = 'raw_data'
        # Specific cell type
        else:
            # Formatting input params for accessing the correct file
            if 'raw_data_' not in cell_type:
                tag = f'raw_data_{cell_type}'
            else:
                tag = cell_type

        fn = os.path.join(
            data_home, pat, 'data', 'Seurat', 'regulon', 'pickle', f'{regulon}_{tag}.pickle'
        )

    try:
        mat = pd.read_pickle(fn)
        return mat
    except FileNotFoundError:
        if tolerate_missing:
            warnings.warn(f'The VIPER matrix not found for pat="{pat}", cell_type="{cell_type}" '
                          f'(should be at "{fn}"). Returning as `None` instead!')
            return None
        else:
            raise FileNotFoundError(
                f'The VIPER matrix not found for pat="{pat}", cell_type="{cell_type}" '
                f'(should be at "{fn}")..'
            )


def get_adj_list(
    cell_type: str,
    net_type: str,
    pat: Union[str, None] = None,
    method: str = 'grnboost2',
    filtered: Union[float, None] = None,
    data_home: str = _DATA_HOME,
    tolerate_missing: bool = True
) -> Union[pd.DataFrame, None]:
    """
    Load adjacency list from patient-specific data, or from cell type aggregated data.

    :param cell_type: The cell type name of the data - could be either:
        e.g. 'Macrophage', 'T_cells' (the cell type identifier)
        e.g. 'all', 'all_data' (the aggregated data - include all cell types)
        e.g. 'raw_data_Macrophage', 'raw_data_T_cells', 'raw_data' (the data file identifier, i.e. 'raw_data_T_cells'
            corresponds to the same data as 'T_cells')
    :param net_type: The type of data - could be either:
        'all' (all gene-gene connections)
        'TF' (TF-target connections)
        'ctx' (enriched TF-target connections)
    :param pat: The patient identifier - could be either:
        e.g. None (include all patients)
        e.g. 'C51', 'C141' (the patient identifier)
        e.g. 'C', 'M', 'S', 'all_patients', 'all_data', 'all' (the identifier of aggregated patient data)
    :param method: The GRN inference method, could be either 'genie3' or 'grnboost2'
    :param filtered: The quantile threshold
    :param data_home: The filepath to the data home folder
    :param tolerate_missing: True if tolerate missing data file (and output None in this case), False otherwise

    :return Pandas dataframe which represents GRN adjacency list
    """

    # Formatting input params for accessing the correct file
    net_suffix = 'TF_cor' if net_type == 'TF' else 'TF_ctx' if net_type == 'ctx' else 'cor'
    filter_suffix = f"_filtered_{str(filtered).replace('.', '_')}" if filtered is not None else ''

    # Loading data that includes all patients
    if pat is None or 'all' in pat:

        data_home = os.path.join(data_home, 'cell_types')
        data_folder = 'all_data' if cell_type in ['all', 'raw_data'] else cell_type.replace('raw_data_', '')

        fn = os.path.join(
            data_home, data_folder, 'data', method, 'pickle', f'raw_data_{net_suffix}{filter_suffix}.pickle'
        )

    # Loading data that includes a specific patient type: 'control', 'moderate' or 'severe'
    elif pat in ['C', 'M', 'S']:

        data_home = os.path.join(data_home, 'cell_types')
        data_folder = 'all_data' if cell_type in ['all', 'raw_data'] else cell_type.replace('raw_data_', '')

        fn = os.path.join(
            data_home, data_folder, 'data', method, 'pickle',
            f'raw_data_{pat}_type_{net_suffix}{filter_suffix}.pickle'
        )

    # Loading patient-specific data
    else:

        # Include all cell types
        if cell_type == 'all_data' or cell_type == 'raw_data' or cell_type == 'all':
            tag = 'raw_data'
        # Specific cell type
        else:
            # Formatting input params for accessing the correct file
            if 'raw_data_' not in cell_type:
                tag = f'raw_data_{cell_type}'
            else:
                tag = cell_type

        fn = os.path.join(
            data_home, pat, 'data', method, 'pickle', f'{tag}_{net_suffix}{filter_suffix}.pickle'
        )

    try:
        adj_list = pd.read_pickle(fn)
        return adj_list
    except FileNotFoundError:
        if tolerate_missing:
            warnings.warn(f'The GRN for pat="{pat}", cell_type="{cell_type}", net_type="{net_type}" is not found '
                          f'(should be at "{fn}"). Returning as `None` instead!')
            return None
        else:
            raise FileNotFoundError(
                f'The GRN for pat={pat}, cell_type={cell_type}, net_type={net_type} is not found '
                f'(should be at "{fn}")..'
            )


def get_nx_graph(
    cell_type: str,
    net_type: str,
    pat: str = None,
    method: str = 'grnboost2',
    filtered: float = None,
    data_home: str = _DATA_HOME,
    tolerate_missing: bool = True
) -> Union[nx.DiGraph, None]:
    """
    Load NetworkX graph from patient-specific data, or from cell type aggregated data.

    :param cell_type: The cell type name of the data - could be either:
        e.g. 'Macrophage', 'T_cells' (the cell type identifier)
        e.g. 'all', 'all_data' (the aggregated data - include all cell types)
        e.g. 'raw_data_Macrophage', 'raw_data_T_cells', 'raw_data' (the data file identifier, e.g. 'raw_data_T_cells'
            corresponds to the same data as 'T_cells')
    :param net_type: The type of data - could be either:
        'all' (all gene-gene connections)
        'TF' (TF-target connections)
        'ctx' (enriched TF-target connections)
    :param pat: The patient identifier - could be either:
        e.g. None (include all patients)
        e.g. 'C51', 'C141' (the patient identifier)
        e.g. 'C', 'M', 'S', 'all_patients', 'all_data', 'all' (the identifier of aggregated patient data)
    :param method: The GRN inference method, could be either 'genie3' or 'grnboost2'
    :param filtered: The quantile threshold
    :param data_home: The filepath to the data home folder
    :param tolerate_missing: True if tolerate missing data file (and output None in this case), False otherwise

    :return NetworkX object of the corresponding GRN
    """

    # Formatting input params for accessing the correct file
    net_suffix = 'TF_cor' if net_type == 'TF' else 'TF_ctx' if net_type == 'ctx' else 'cor'
    filter_suffix = f"_filtered_{str(filtered).replace('.', '_')}" if filtered is not None else ''

    # Loading data that includes all patients
    if pat is None or 'all' in pat:

        data_home = os.path.join(data_home, 'cell_types')
        data_folder = 'all_data' if cell_type in ['all', 'raw_data'] else cell_type.replace('raw_data_', '')

        fn = os.path.join(
            data_home, data_folder, 'data', method, 'nx_graph', f'raw_data_{net_suffix}{filter_suffix}.gpickle'
        )

    # Loading data that includes a specific patient type: 'control', 'moderate' or 'severe'
    elif pat in ['C', 'M', 'S']:

        # Loading patient-type aggregated data
        data_home = os.path.join(data_home, 'cell_types')
        data_folder = 'all_data' if cell_type in ['all', 'raw_data'] else cell_type.replace('raw_data_', '')

        fn = os.path.join(
            data_home, data_folder, 'data', method, 'nx_graph',
            f'raw_data_{pat}_type_{net_suffix}{filter_suffix}.gpickle'
        )

    # Loading patient-specific data
    else:

        # Include all cell types
        if cell_type == 'all_data' or cell_type == 'raw_data' or cell_type == 'all':
            tag = 'raw_data'
        # Specific cell type
        else:
            # Formatting input params for accessing the correct file
            if 'raw_data_' not in cell_type:
                tag = f'raw_data_{cell_type}'
            else:
                tag = cell_type

        fn = os.path.join(
            data_home, pat, 'data', method, 'nx_graph', f'{tag}_{net_suffix}{filter_suffix}.gpickle'
        )

    try:
        G = nx.read_gpickle(fn)
        return G
    except FileNotFoundError:
        if tolerate_missing:
            warnings.warn(f'The GRN for pat="{pat}", cell_type="{cell_type}", net_type="{net_type}" is not found '
                          f'(should be at "{fn}"). Returning as `None` instead!')
            return None
        else:
            raise FileNotFoundError(
                f'The GRN for pat={pat}, cell_type={cell_type}, net_type={net_type} is not found '
                f'(should be at "{fn}")..'
            )


def get_community_info(
    cell_type: str,
    pat: str = None,
    method: str = 'grnboost2',
    algo: str = 'leiden',
    data_home: str = _DATA_HOME,
    tolerate_missing: bool = True
) -> pd.DataFrame:
    """
    Load community summary information that is produced by running `community_ana.py`.

    :param cell_type: The cell type name of the data - could be either:
        e.g. 'Macrophage', 'T_cells' (the cell type identifier)
        e.g. 'all', 'all_data' (the aggregated data - include all cell types)
        e.g. 'raw_data_Macrophage', 'raw_data_T_cells', 'raw_data' (the data file identifier, e.g. 'raw_data_T_cells'
            corresponds to the same data as 'T_cells')
    :param net_type: The type of data - could be either:
        'all' (all gene-gene connections)
        'TF' (TF-target connections)
        'ctx' (enriched TF-target connections)
    :param pat: The patient identifier - could be either:
        e.g. None (include all patients)
        e.g. 'C51', 'C141' (the patient identifier)
        e.g. 'C', 'M', 'S', 'all_patients', 'all_data', 'all' (the identifier of aggregated patient data)
    :param method: The GRN inference method, could be either 'genie3' or 'grnboost2'
    :param algo: The community detection algorithm used
    :param data_home: The filepath to the data home folder
    :param tolerate_missing: True if tolerate missing data file (and output None in this case), False otherwise

    :returns: A dataframe with all information about communities in queried graph network. The 
        format of the dataframe could be found in `_community.py/process_communities`
    """

    # Loading data that includes all patients
    if pat is None or 'all' in pat:

        data_home = os.path.join(data_home, 'cell_types')
        data_folder = 'all_data' if cell_type in ['all', 'raw_data'] else cell_type.replace('raw_data_', '')

        fn = os.path.join(
            data_home, data_folder, 'data', method, f'{algo}_communities', 
            f'raw_data_communities_info.pickle'
        )

    # Loading data that includes a specific patient type: 'control', 'moderate' or 'severe'
    elif pat in ['C', 'M', 'S']:

        # Loading patient-type aggregated data
        data_home = os.path.join(data_home, 'cell_types')
        data_folder = 'all_data' if cell_type in ['all', 'raw_data'] else cell_type.replace('raw_data_', '')

        fn = os.path.join(
            data_home, data_folder, 'data', method, f'{algo}_communities', 
            f'raw_data_{pat}_type_communities_info.pickle'
        )

    # Loading patient-specific data
    else:

        # Include all cell types
        if cell_type == 'all_data' or cell_type == 'raw_data' or cell_type == 'all':
            tag = 'raw_data'
        # Specific cell type
        else:
            # Formatting input params for accessing the correct file
            if 'raw_data_' not in cell_type:
                tag = f'raw_data_{cell_type}'
            else:
                tag = cell_type

        fn = os.path.join(
            data_home, data_folder, 'data', method, f'{algo}_communities', 
            f'{tag}_type_communities_info.pickle'
        )

    try:
        df = pd.read_pickle(fn)
        return df
    except FileNotFoundError:
        if tolerate_missing:
            warnings.warn(f'The community info df for pat="{pat}", cell_type="{cell_type}" is not found '
                          f'(should be at "{fn}"). Returning as `None` instead!')
            return None
        else:
            raise FileNotFoundError(
                f'The community info df for pat={pat}, cell_type={cell_type} is not found '
                f'(should be at "{fn}")..'
            )


def _compute_graph_stats(
    curr_pat_: str,
    curr_ctype_: str,
    curr_ntype_: str,
    _meta: pd.DataFrame,
    _q_thresh: float
) -> dict:
    """
    Compute graph statistics for given patient, cell type and net type.

    :param curr_ctype_: The cell type name of the data - could be either:
        e.g. 'Macrophage', 'T_cells' (the cell type identifier)
        e.g. 'all', 'all_data' (the aggregated data - include all cell types)
        e.g. 'raw_data_Macrophage', 'raw_data_T_cells', 'raw_data' (the data file identifier, e.g. 'raw_data_T_cells'
            corresponds to the same data as 'T_cells')
    :param curr_ntype_: The type of data - could be either:
        'all' (all gene-gene connections)
        'TF' (TF-target connections)
        'ctx' (enriched TF-target connections)
    :param curr_pat_: The patient identifier - could be either:
        e.g. None (include all patients)
        e.g. 'C51', 'C141' (the patient identifier)
        e.g. 'C', 'M', 'S', 'all_data', 'all' (the identifier of aggregated patient data)
    :param _meta: The dataframe containing metadata about patients (cell count, patient types, file paths)
    :param _q_thresh: The quantile threshold

    :return: A dictionary containing computed graph properties
    """

    out = {}
    try:
        full_G = get_nx_graph(curr_ctype_, curr_ntype_, pat=curr_pat_,
                              filtered=_q_thresh, tolerate_missing=False)
        if not nx.is_connected(full_G.to_undirected()):
            Gcc = sorted(nx.connected_components(full_G.to_undirected()), key=len, reverse=True)
            G_biggest = full_G.subgraph(Gcc[0])
        else:
            G_biggest = full_G
        G_igraph = ig.Graph.from_networkx(G_biggest.to_undirected())
        out = {
            'num_nodes': full_G.number_of_nodes(),
            'num_edges': full_G.number_of_edges(),
            'num_cells': get_num_cells(curr_pat_, curr_ctype_, _meta),
            'average_degree': np.mean([degree for node, degree in full_G.degree(full_G.nodes())]),
            'radius': G_igraph.radius(),
            'average_path_length': G_igraph.average_path_length(),
            'diameter': G_igraph.diameter(),
            'importances': [info['importance'] for st, end, info in full_G.edges(data=True)],
            'rhos': [info['rho'] for st, end, info in full_G.edges(data=True)],
            'median_importance': np.median([info['importance'] for st, end, info in full_G.edges(data=True)]),
            'median_rho': np.median([info['rho'] for st, end, info in full_G.edges(data=True)]),
            'STD_importance': np.std([info['importance'] for st, end, info in full_G.edges(data=True)]),
            'STD_rho': np.std([info['rho'] for st, end, info in full_G.edges(data=True)]),
            'tfs': set([st for st, _, _ in full_G.edges(data=True)]),
            'num_tfs': len(set([st for st, _, _ in full_G.edges(data=True)]))
        }
    except BaseException as e:
        warnings.warn(f'Caught an error with the message:\n{e}')

    return out


def get_graph_stats(meta: pd.DataFrame, n_jobs: int, filtered: float = None) -> dict:
    """
    Compute statistics of graph properties with available GRNs.

    :param meta: The dataframe containing metadata about patients (cell count, patient types, file paths)
    :param n_jobs: The number of parallel executions
    :param filtered: The quantile threshold or None, if not to use filtered GRN

    :return: A dictionary containing computed graph properties for all available GRNs
    """

    n_types = ['all', 'ctx']  # computing only for these types of networks

    def chunks(lst, n):
        """Yield successive n-sized chunks from lst."""
        for i in range(0, len(lst), n):
            yield lst[i:i + n]

    stat_types = [
        'num_nodes', 'num_edges', 'num_cells', 'average_degree', 'radius',
        'average_path_length', 'diameter', 'importances', 'rhos',
        'median_importance', 'median_rho', 'STD_importance', 'STD_rho',
        'tfs', 'num_tfs'
    ]
    num_all_sets = get_avail_nx_graphs()['all'].sum().sum()

    # Computing statistics in parallel
    all_out = Parallel(n_jobs=n_jobs)(
        delayed(_compute_graph_stats)(_pat, _ctype, _ntype, meta, filtered) for _ntype in n_types
        for _pat in get_avail_nx_graphs()[_ntype].index
        for _ctype in get_avail_pat_nx(_ntype, _pat)
    )
    graph_stats = {
        _ntype: {
            curr_stat_type: [
                tmp_dict[curr_stat_type] for tmp_dict in el if len(tmp_dict) != 0
            ] for curr_stat_type in stat_types
        } for _ntype, el in zip(n_types, chunks(all_out, num_all_sets))
    }

    # Concatenating lists together
    for _n_type in n_types:
        graph_stats[_n_type]['importances'] = list(itertools.chain(*graph_stats[_n_type]['importances']))
        graph_stats[_n_type]['rhos'] = list(itertools.chain(*graph_stats[_n_type]['rhos']))

    return graph_stats
