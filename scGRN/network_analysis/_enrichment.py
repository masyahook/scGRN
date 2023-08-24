"""Enrichment analysis."""

import json
import sys
from io import StringIO
from pathlib import Path

import pandas as pd
import requests

ENRICHR_POST = 'http://amp.pharm.mssm.edu/Enrichr/addList'
ENRICHR_GET = 'http://amp.pharm.mssm.edu/Enrichr/export'


def run_enrichr(
    in_path: str,
    gene_col: str,
    out_path: str,
    group_col: str = None,
    enrichr_library: str = 'Reactome_2016',
    query: str = None,
    top_n: int = 50,
) -> pd.DataFrame:
    """
    Run enrichment analysis with Enrichr on the passed dataset `in_path` with genes in `gene_col` 
    and groups in `group_col`. This function will run Enrichr on the gene set for each group and 
    save results in the passed file. A user can run the function for the following cases:

        1. Finding functional differences between two or more groups of cells based on sets of genes
            In this case the dataset has, for example, the following format:
                            p_val	avg_log2FC	pct.1	pct.2	p_val_adj	cluster	gene
            IRF2	0.000000e+00	8.346623	0.050	0.005	0.000000e+00	M	IRF2
            STAT1	0.000000e+00	5.006633	0.002	0.000	0.000000e+00	M	STAT1
            NR3C1	0.000000e+00	4.831237	0.092	0.002	0.000000e+00	M	NR3C1
            IRF1	0.000000e+00	4.676670	0.002	0.000	0.000000e+00	M	IRF1
            REST	0.000000e+00	4.024687	0.215	0.016	0.000000e+00	M	REST

            Here `gene_col="gene"` and `group_col="cluster"

        2. Getting functional annotation for communities detected in gene regulatory networks
            In this case the dataset has, for example, the following format:
                num_nodes	num_edges	main_functions_GO	                                all_sorted_genes
            0	3978	    12830	    >>> phosphatase activity <<<: PTPRK,DUSP14,PTE...	DUSP1 (score=0.19707084599626928); KLF2 (score...
            1	1812	    11803	    >>> regulation of microtubule polymerization o...	STMN1 (score=0.19654841042005425); TYMS (score...
            2	1576	    6114	    >>> immune response <<<: CCL5,IGKV1-12,CST7,IG...	CCL5 (score=0.2095641475565237); NKG7 (score=0...
            3	1494	    11988	    >>> nucleic acid binding <<<: TRIM27,PARN,SON,...	ISG20 (score=0.2054951704917856); MX1 (score=0...
            4	1481	    21413	    >>> nucleotide binding <<<: AK6,DDX52,ABCF1,DD...	ACTB (score=0.31799517570308644); GAPDH (score...
            5	1398	    6500	    >>> ferric iron binding <<<: FTH1,FTL; >>> iro...	FTH1 (score=0.17875800169417477); CCL2 (score=...

            Here `gene_col="all_sorted_genes"` and `group_col=None` (all genes stored within one cell)

    The output of the function is stored in `out_path` in the following format:

            Term	                        Overlap	P-value	        Adjusted P-value	Old P-value	Old Adjusted P-value	Odds Ratio	Combined Score	Genes
        0	TNF-alpha Signaling via NF-kB	8/200	5.628165e-10	1.407041e-08	    0	        0	                    35.827899	763.065002	    CEBPB;BCL6;MAFF;FOSB;FOS;KLF2;ATF3;RELB
        1	Hypoxia	                        4/200	2.471030e-04	1.544394e-03	    0	        0	                    14.945578	124.133566	    MAFF;FOS;ETS1;ATF3
        2	G2-M Checkpoint	                4/200	2.471030e-04	1.544394e-03	    0	        0	                    14.945578	124.133566	    TFDP1;E2F1;E2F2;MYBL2
        3	KRAS Signaling Up	            4/200	2.471030e-04	1.544394e-03	    0	        0	                    14.945578	124.133566	    MAFB;IRF8;ETV4;ETS1
        4	UV Response Up	                3/158	1.848501e-03	9.242505e-03	    0	        0	                    13.696313	86.196108	    FOSB;FOS;ATF3

    :param in_path: The path to input data where gene sets are stored
    :param gene_col: The column name that stores gene names
    :param out_path: The path to output data to store
    :param group_col: The column names that store group names
    :param enrichr_library: The EnrichR library to use for enrichment analysis. The list of libraries
        is available here: https://maayanlab.cloud/Enrichr/#libraries
    :param query: The query that can be used to select a subset of original dataset. For example, it
        is useful when a user wants to run Enrichr only on positively expressed (and statistically 
        significant) genes within each group, i.e. in this case:
            `query="'p_val_adj' < 0.05 & 'avg_log2FC' > 1"`
    :param top_n: If a user wants to run EnrichR only on the first `top_n` genes, it can use this 
        parameter. Useful in the 2nd case above (running only on the most central genes within each
        community)

    :returns: The results from EnrichR that was performed on each group separately
    """

    if in_path.endswith('.tsv'):
        in_df = pd.read_csv(in_path, sep='\t')
    elif in_path.endswith('.csv'):
        in_df = pd.read_csv(in_path, sep=',')
    elif in_path.endswith('.pickle'):
        in_df = pd.read_pickle(in_path)
    else:
        raise NotImplementedError(f'Unsupported file format {in_df}. Currently .csv, .tsv and .pickle are supported')
    
    # Selecting a subset of records if needed
    if query is not None:
        in_df = in_df.query(query)

    # Creating a group -> gene list mapping
    if group_col is not None:
        group_to_genes = in_df.groupby(group_col)[gene_col].agg(lambda x: x.to_list())
    else:
        group_col = 'group_to_col'
        group_to_genes = (
            in_df.assign(**{
                group_col: [f'group_{i}' for i in range(in_df.shape[0])],  # creating a group column
                gene_col: in_df[gene_col].map(lambda lst: [  # selecting the first `top_n` genes from the full list
                    el[: el.find(' ')] for el in lst.split('; ')
                ][:top_n]
                )
            })
            ).set_index(group_col)[gene_col] 

    print(f'\nRunning EnrichR analysis on dataframe: {in_path}.\n\n')
    print(f'Overall, detected {len(group_to_genes)} clusters: {group_to_genes.index.to_list()}\n')

    out_df = pd.DataFrame()
    for group, gene_list in group_to_genes.items():
        
        print(f'    Processing {group}..')

        # Uploading the gene set to EnrichR
        genes_str = '\n'.join(gene_list)
        description = Path(in_path).stem           

        payload = {
          'list': (None, genes_str),
          'description': (None, description)
        }
        response = requests.post(ENRICHR_POST, files=payload)

        if not response.ok:
            raise Exception('Error analyzing gene list')

        job_id = json.loads(response.text)

        # Getting the functional set from EnrichR
        query_string = '?userListId=%s&filename=%s&backgroundType=%s'
        user_list_id = str(job_id['userListId'])
        url = ENRICHR_GET + query_string % (user_list_id, out_path, enrichr_library)

        response = requests.get(url, stream=True).text

        # Concatenating the result for the current group to others
        out_df = pd.concat([
            out_df,
            pd.read_csv(StringIO(response), sep='\t').assign(**{group_col: group})
        ])

    # Saving the results to file
    out_df.to_csv(out_path, index=False)

    print(f'\nSuccess! Saved at {out_path}')

    return out_df
