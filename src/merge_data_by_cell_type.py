import argparse
import os
import shutil

import pandas as pd


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Merge datasets by cell type from COVID-19 data.')
    args = parser.parse_args()
    
    # Setting up pathnames
    _PROJ_PATH = '/gpfs/projects/bsc08/bsc08890'
    _FMETA = os.path.join(_PROJ_PATH, 'data/GSE145926_RAW/metadata.tsv')
    _DATA_HOME = os.path.join(_PROJ_PATH, 'res/covid_19')
    
    # Reading metadata
    full_meta = pd.read_csv(_FMETA, sep='\t', index_col=0)
    
    cell_type_metadata = {}
    cell_type_data = {}
    for patient_id in full_meta.index:
        print(f'Processing {patient_id}..')
        # Getting cell type information for the current patient
        pat_metadata = pd.read_csv(os.path.join(_DATA_HOME, patient_id, 'data', 'Seurat', 'cells_metadata.tsv'), sep='\t')
        pat_metadata['patient'] = patient_id
        cell_types = list(pat_metadata['cell_type'].unique())
        
        for cell_type in ['all_data'] + cell_types:            
            # Creating pathnames for the current cell type
            final_filedir = os.path.join(_DATA_HOME, 'cell_types', cell_type, 'data', 'Seurat')
            os.makedirs(final_filedir, exist_ok=True)
            curr_fn = 'raw_data.tsv' if cell_type == 'all_data' else f'raw_data_{cell_type}.tsv'
            
            curr_data = pd.read_csv(os.path.join(_DATA_HOME, patient_id, 'data', 'Seurat', curr_fn), sep='\t')
            
            # Selecting only cell-type specific cells (or all data)
            if cell_type in cell_type_metadata.keys():
                if cell_type == 'all_data':
                    cell_type_metadata[cell_type] = pd.concat([
                        cell_type_metadata[cell_type], pat_metadata
                    ], axis=0)
                else:
                    cell_type_metadata[cell_type] = pd.concat([
                        cell_type_metadata[cell_type], pat_metadata[pat_metadata['cell_type'] == cell_type]
                    ], axis=0)
                cell_type_data[cell_type] = pd.concat([
                    cell_type_data[cell_type], curr_data
                ], axis=1)
            else:
                if cell_type == 'all_data':
                    cell_type_metadata[cell_type] = pat_metadata
                else:
                    cell_type_metadata[cell_type] = pat_metadata[pat_metadata['cell_type'] == cell_type] 
                cell_type_data[cell_type] = curr_data
                
    print('\n\nSaving..')
            
    # Saving metadata for each cell-type specific set of cells
    for cell_type in cell_type_metadata.keys():
        final_filedir = os.path.join(_DATA_HOME, 'cell_types', cell_type, 'data', 'Seurat')
        cell_type_metadata[cell_type].dropna().to_csv(os.path.join(final_filedir, f'all_cells_metadata.tsv'), sep='\t', index_label=False)
        cell_type_data[cell_type].dropna().to_csv(os.path.join(final_filedir, f'all_raw_data.tsv'), sep='\t', index_label=False)
        
