import os
import pandas as pd

from ..utils import is_non_empty


def get_meta(
        data_home: str = '/gpfs/projects/bsc08/bsc08890/res/covid_19',
        meta_file: str = '/gpfs/projects/bsc08/bsc08890/data/GSE145926_RAW/metadata.tsv'
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
        data_home: str = '/gpfs/projects/bsc08/bsc08890/res/covid_19',
        meta_file: str = '/gpfs/projects/bsc08/bsc08890/data/GSE145926_RAW/metadata.tsv'
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


def get_avail_adj_lists(
        data_home: str = '/gpfs/projects/bsc08/bsc08890/res/covid_19',
        meta_file: str = '/gpfs/projects/bsc08/bsc08890/data/GSE145926_RAW/metadata.tsv',
        filtered: float = 0.95
) -> dict:
    """
    Obtain dictionary containing the information about GRN adjacency list availability.

    :param data_home: The filepath to the data home folder
    :param meta_file: The filepath to metadata about patients
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

    for d_type, data in avail_adj_lists.items():
        d_type_suffix = 'cor' if d_type == 'all' else 'TF_cor' if d_type == 'TF' else 'TF_ctx'

        # Checking the presence of patient data files: e.g. 'C141', 'C51', 'C152'
        for pat in data.index[4:]:
            for d in avail_sc_data.loc[pat].dropna().index.to_list():
                d_fn = f'raw_data_{d}' if d != 'all_data' else 'raw_data'
                fn = os.path.join(data_home, pat, 'data', 'grnboost2', f'{d_fn}_{d_type_suffix}.tsv')
                fn_pickled = os.path.join(data_home, pat, 'data', 'grnboost2', 'pickle',
                                          f'{d_fn}_{d_type_suffix}.pickle')
                filtered_fn_pickled = os.path.join(data_home, pat, 'data', 'grnboost2', 'pickle',
                                                   f'{d_fn}_{d_type_suffix}_filtered_{filtered_suffix}.pickle')
                if is_non_empty(fn) and is_non_empty(fn_pickled) and is_non_empty(filtered_fn_pickled):
                    data.loc[pat, d] = True
                else:
                    data.loc[pat, d] = False

        # Checking the presence of cell type aggregated data files: 'all_data', 'C', 'M', 'S'
        for pat in data.index[:4]:
            for d in avail_sc_data.loc[pat].dropna().index.to_list():
                d_fn = 'raw_data' if pat == 'all_data' else f'raw_data_{pat}_type'
                fn = os.path.join(data_home, 'cell_types', d, 'data', 'grnboost2', f'{d_fn}_{d_type_suffix}.tsv')
                fn_pickled = os.path.join(data_home, 'cell_types', d, 'data', 'grnboost2', 'pickle',
                                          f'{d_fn}_{d_type_suffix}.pickle')
                filtered_fn_pickled = os.path.join(data_home, 'cell_types', d, 'data', 'grnboost2', 'pickle',
                                                   f'{d_fn}_{d_type_suffix}_filtered_{filtered_suffix}.pickle')
                if is_non_empty(fn) and is_non_empty(fn_pickled) and is_non_empty(filtered_fn_pickled):
                    data.loc[pat, d] = True
                else:
                    data.loc[pat, d] = False

    return avail_adj_lists


def get_avail_nx_graphs(
        data_home: str = '/gpfs/projects/bsc08/bsc08890/res/covid_19',
        meta_file: str = '/gpfs/projects/bsc08/bsc08890/data/GSE145926_RAW/metadata.tsv',
        filtered: float = 0.95
) -> dict:
    """
    Obtain dictionary containing the information about GRN NetworkX graph availability.

    :param data_home: The filepath to the data home folder
    :param meta_file: The filepath to metadata about patients
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

    for d_type, data in avail_nx_graphs.items():
        d_type_suffix = 'cor' if d_type == 'all' else 'TF_cor' if d_type == 'TF' else 'TF_ctx'

        # Checking the presence of patient data files: e.g. 'C141', 'C51', 'C152'
        for pat in data.index[4:]:
            for d in avail_sc_data.loc[pat].dropna().index.to_list():
                d_fn = f'raw_data_{d}' if d != 'all_data' else 'raw_data'
                fn = os.path.join(data_home, pat, 'data', 'grnboost2', 'nx_graph', f'{d_fn}_{d_type_suffix}.gpickle')
                fn_filtered = os.path.join(data_home, pat, 'data', 'grnboost2', 'nx_graph',
                                           f'{d_fn}_{d_type_suffix}_filtered_{filtered_suffix}.gpickle')
                if is_non_empty(fn) and is_non_empty(fn_filtered):
                    data.loc[pat, d] = True
                else:
                    data.loc[pat, d] = False

        # Checking the presence of cell type aggregated data files: 'all_data', 'C', 'M', 'S'
        for pat in data.index[:4]:
            for d in avail_sc_data.loc[pat].dropna().index.to_list():
                d_fn = 'raw_data' if pat == 'all_data' else f'raw_data_{pat}_type'
                fn = os.path.join(data_home, 'cell_types', d, 'data', 'grnboost2', 'nx_graph',
                                  f'{d_fn}_{d_type_suffix}.gpickle')
                fn_filtered = os.path.join(data_home, 'cell_types', d, 'data', 'grnboost2', 'nx_graph',
                                           f'{d_fn}_{d_type_suffix}_filtered_{filtered_suffix}.gpickle')
                if is_non_empty(fn) and is_non_empty(fn_filtered):
                    data.loc[pat, d] = True
                else:
                    data.loc[pat, d] = False

    return avail_nx_graphs


def get_avail_nx(
        d_type: str,
        pat: str
) -> pd.Series:
    """
    Get a list of available GRN NetworkX graphs corresponding to `d_type` and `pat`.

    :param d_type: The type of data - could be either:
        'all' (all gene-gene connections)
        'TF' (TF-target connections)
        'ctx' (enriched TF-target connections)
    :param pat: The patient identifier - could be either:
        e.g. 'C51', 'C141' (the patient identifier)
        e.g. 'C', 'M', 'S', 'all_data' (the identifier of aggregated patient data)

    :return: A list of available data files corresponding to `d_type` and `pat`
    """

    avail_Gs = get_avail_nx_graphs()
    return avail_Gs[d_type].loc[pat].dropna().loc[lambda x: x].index.map(
        lambda x: f'raw_data_{x}' if x != 'all_data' else 'raw_data')
