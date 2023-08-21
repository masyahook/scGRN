"""Enrichment analysis."""

import json
import requests
import sys
import io 


def run_enrichr(data, is_communities=False, is_positive_markers=True, group_types = 'all', on_targets=False, choose_fixed_tf=None,
                data_type='all', top_n=50, algo='leiden', enrichr_library='MSigDB_Hallmark_2020'):
    """
    Run enrichment analysis with Enrichr.
    """
    
    out_folder = 'community_ana' if is_communities else 'cohort_ana'
    
    if is_communities == True:
        
        print('Running EnrichR on communities..')
    
        algo = 'leiden'
        _DATA_HOME = '/gpfs/projects/bsc08/bsc08890/res/covid_19'

        if data_type == 'all':
            community_data = pd.read_pickle(os.path.join(
                _DATA_HOME, 'cell_types', data, 'data', 'grnboost2', f'{algo}_communities', 
                f'raw_data_communities_info.pickle'
            ))
        else:
            community_data = pd.read_pickle(os.path.join(
                _DATA_HOME, 'cell_types', data, 'data', 'grnboost2', f'{algo}_communities', 
                f'raw_data_{data_type}_type_communities_info.pickle'
            ))

        df = pd.concat([
            pd.DataFrame({
                'cluster': f'cluster_{i}',
                'gene': [el[: el.find(' ')] for el in vals.split('; ')][:top_n]
            }) for i, vals in community_data['all_sorted_genes'].iteritems()
        ], axis=0).reset_index(drop=True)
        
    else:
        
        if on_targets:
            
            print('Running EnrichR on targets between 3 group types..')
            
            types = ['C', 'M', 'S']
            
            df = pd.concat([
                pd.read_csv(
                    f'/gpfs/home/bsc08/bsc08890/tmp/cohort_ana/tmp_enrichr_{data}_{t}_{choose_fixed_tf}_target_list.tsv', 
                    header=None, names=['gene']
                ).assign(cluster=t) for t in types
            ], axis=0)
            
        else:
        
            if group_types == 'all':
                print('Running EnrichR on TFs between 3 group types..')
                df = pd.read_csv(f'/gpfs/home/bsc08/bsc08890/tmp/tf_markers_df_{data}.tsv', sep='\t')
            else:
                print('Running EnrichR on 2 group types..')
                if group_types == 'M_S':
                    group_types = 'S_M'
                if group_types == 'C_M':
                    group_types = 'M_C'
                if group_types == 'C_S':
                    group_types = 'S_C'
                df_1 = pd.read_csv(f'/gpfs/home/bsc08/bsc08890/tmp/tf_markers_df_{group_types}_{data}.tsv', sep='\t')
                df_1['gene'] = df_1.index
                df_2 = df_1.copy()
                df_2['avg_log2FC'] = - df_2['avg_log2FC']
                df_1['cluster'], df_2['cluster'] = group_types.split('_')

                df = pd.concat([df_1, df_2], axis=0)


            if is_positive_markers:
                df = df[(df['p_val_adj'] < 0.05) & (df['avg_log2FC'] > 1)]
            else:
                df = df[(df['p_val_adj'] < 0.05) & (df['avg_log2FC'] < -1)]

    cluster_dfs = {}
    for cl in df['cluster'].unique():
        
        print(f'Processing {cl}..')

        ENRICHR_URL = 'http://amp.pharm.mssm.edu/Enrichr/addList'
        genes_str = '\n'.join(df[df['cluster'] == cl]['gene'])
        description = f"{data}_{data_type}_{cl}"
        
        if is_communities == True:
            filename = f'tmp/{out_folder}/tmp_enrichr_{data}_{data_type}_{cl}.tsv'
        elif on_targets:
            filename = f'tmp/{out_folder}/tmp_enrichr_{data}_{data_type}_{choose_fixed_tf}_target_{cl}.tsv'
        elif group_types == 'all':
            filename = f'tmp/{out_folder}/tmp_enrichr_{data}_{data_type}_{cl}.tsv'
        else:
            filename = f'tmp/{out_folder}/tmp_enrichr_{data}_2_groups_{cl}.tsv'
            

        payload = {
          'list': (None, genes_str),
          'description': (None, description)
        }
        response = requests.post(ENRICHR_URL, files=payload)

        if not response.ok:
            raise Exception('Error analyzing gene list')

        job_id = json.loads(response.text)

        ################################################################################
        # Get enrichment results
        #
        ENRICHR_URL = 'http://amp.pharm.mssm.edu/Enrichr/export'
        query_string = '?userListId=%s&filename=%s&backgroundType=%s'
        user_list_id = str(job_id['userListId'])
        gene_set_library = str(enrichr_library)
        url = ENRICHR_URL + query_string % (user_list_id, filename, gene_set_library)

        response = requests.get(url, stream=True)

        print('        Enrichr API : Downloading file of enrichment results: Job Id:', job_id)
        with open(filename, 'wb') as f:
            for chunk in response.iter_content(chunk_size=1024):
                if chunk:
                    f.write(chunk)
                    
        print(f'        Saved to {filename}')
                    
        cluster_dfs[cl] = pd.read_csv(filename, sep='\t')

    return cluster_dfs
