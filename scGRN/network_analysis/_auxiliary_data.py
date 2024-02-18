"""Load external data for network analysis."""

import os
from typing import Union

import pandas as pd
from tqdm import tqdm as tqdm_cli

from ..config import _INPUT_DATA_HOME


def load_gene_func_db_mapping(
    path_to_dbs: str = f"{_INPUT_DATA_HOME}/Gene_func_associations",
) -> dict:
    """
    Load downloaded gene function databases like MSigDB, GO, SIGNOR, etc mappings.

    The user can download MSigDB datasets of gene set mappings here:
    https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp.
    DoRothEA regulons could be accessed here: https://saezlab.github.io/dorothea/

    The GO annotation was obtained from BioMart: https://www.ensembl.org/biomart/martview/

    return: The file paths to the local database files
    """

    file_mapping = {
        "GO": os.path.join(path_to_dbs, "GO_annotation.tsv"),
        "DoRothEA": os.path.join(path_to_dbs, "dorothea_regulons.tsv"),
        "MSigDB_hallmark_gene_sets_h": os.path.join(
            path_to_dbs, "MSigDB", "h.all.v7.5.1.symbols.gmt.txt"
        ),
        "MSigDB_curated_gene_sets_c2_all": os.path.join(
            path_to_dbs, "MSigDB", "c2.all.v7.5.1.symbols.gmt.txt"
        ),
        "MSigDB_curated_gene_sets_c2_cp": os.path.join(
            path_to_dbs, "MSigDB", "c2.cp.v7.5.1.symbols.gmt.txt"
        ),
        "MSigDB_curated_gene_sets_c2_cp_kegg": os.path.join(
            path_to_dbs, "MSigDB", "c2.cp.kegg.v7.5.1.symbols.gmt.txt"
        ),
        "MSigDB_curated_gene_sets_c2_cp_reactome": os.path.join(
            path_to_dbs, "MSigDB", "c2.cp.reactome.v7.5.1.symbols.gmt.txt"
        ),
        "MSigDB_curated_gene_sets_c2_cp_wikipathways": os.path.join(
            path_to_dbs, "MSigDB", "c2.cp.wikipathways.v7.5.1.symbols.gmt.txt"
        ),
        "MSigDB_regulatory_target_gene_sets_c3_all": os.path.join(
            path_to_dbs, "MSigDB", "c3.all.v7.5.1.symbols.gmt.txt"
        ),
        "MSigDB_immunologic_signature_gene_sets_c7_all": os.path.join(
            path_to_dbs, "MSigDB", "c7.all.v7.5.1.symbols.gmt.txt"
        ),
    }

    return file_mapping


def load_gene_func_db(
    db: str, reload: bool = False, as_series: bool = False
) -> Union[pd.DataFrame, pd.Series]:
    """
    Load data from gene function database like MSigDB, GO, DoRothEA, or etc. The output will be in a
        Pandas Dataframe
    format.

    :param db: The db tag - either a short version, or long version
    :param reload: Only relevant for MSigDB - whether to reload the data from the web or not
    :param as_series: Only relevant for MSigDB - return as pd.Series with index denoting gene names,
        and values denoting the functional annotation

    return: Function mapping of the gene set
    """

    def short_to_long_tag(_db):
        """
        Transforming the short annotation db tag to long format.
        """
        if _db == "GO":
            return "GO"
        elif _db == "DoRothEA":
            return "DoRothEA"
        elif _db == "KEGG":
            return "MSigDB_curated_gene_sets_c2_cp_kegg"
        elif _db == "hallmark":
            return "MSigDB_hallmark_gene_sets_h"
        elif _db == "reactome":
            return "MSigDB_curated_gene_sets_c2_cp_reactome"
        elif _db == "wikipathways":
            return "MSigDB_curated_gene_sets_c2_cp_wikipathways"
        elif _db == "immunological":
            return "MSigDB_immunologic_signature_gene_sets_c7_all"
        elif _db == "curated":
            return "MSigDB_curated_gene_sets_c2_all"
        elif _db == "canonical":
            return "MSigDB_curated_gene_sets_c2_cp"
        else:
            raise NotImplementedError(
                f"The tag '{_db}' not found in the downloaded databases.."
            )

    # Getting the db tag
    db = short_to_long_tag(db)

    # Dealing with path names
    db_path = load_gene_func_db_mapping()[db]
    db_folder = db_path[: db_path.rfind("/")]

    if db.startswith("MSigDB"):
        # if reload == True -> download the data from the web, use the local saved version otherwise
        if reload:
            # Loading saved data from https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp
            with open(db_path, "r") as f:
                entries = [f_line.strip().split("\t") for f_line in f.readlines()]

            # Converting data to pandas format
            entries_df = pd.DataFrame(columns=["ID", "link", "genes"])

            for i, entry in enumerate(entries):
                entries_df.loc[i, "ID"] = entry[0]
                entries_df.loc[i, "link"] = entry[1]
                entries_df.loc[i, "genes"] = entry[2:]

            # Adding some information about each functional group by scraping from the web
            for i, row in tqdm_cli(
                entries_df.iterrows(),
                leave=False,
                total=entries_df.shape[0],
                ascii=True,
            ):
                brief, full = (
                    pd.read_html(row["link"])[1]
                    .loc[
                        lambda x: x[0].isin(
                            ["Brief description", "Full description or abstract"]
                        )
                    ][1]
                    .values
                )

                entries_df.loc[i, "brief_desc"] = brief
                entries_df.loc[i, "full_desc"] = full

            # Transforming data to format when gene is the index
            per_gene_entries_df = pd.DataFrame(
                columns=["gene_name", "ID", "link", "brief_desc", "full_desc"]
            )
            for i, row in entries_df.iterrows():
                num_genes = len(row["genes"])
                row_df = pd.DataFrame(
                    dict(
                        gene_name=row["genes"],
                        ID=[row["ID"]] * num_genes,
                        link=[row["link"]] * num_genes,
                        brief_desc=[row["brief_desc"]] * num_genes,
                        full_desc=[row["full_desc"]] * num_genes,
                    )
                )
                per_gene_entries_df = pd.concat([per_gene_entries_df, row_df], axis=0)

            out = per_gene_entries_df.set_index("gene_name")
            out.to_csv(os.path.join(db_folder, f"{db}.tsv"), sep="\t")

        else:
            # Reading previously saved data
            out = pd.read_csv(
                os.path.join(db_folder, f"{db}.tsv"), sep="\t", index_col=0
            )

        if as_series:
            out = out["brief_desc"]

    elif db == "GO":
        # Loading gene GO description - it was downloaded using ENSEMBL BioMart
        entries_df = (
            pd.read_csv(db_path, sep="\t", index_col=0)
            .drop(columns=["GO term evidence code"])
            .dropna(subset=["Gene name", "GO term name", "GO domain"])
            .drop_duplicates()
        )
        entries_df = entries_df.rename(
            columns={
                "Gene name": "gene_name",
                "GO term accession": "ID",
                "GO term name": "brief_desc",
                "GO term definition": "full_desc",
                "GO domain": "domain",
            }
        )

        out = entries_df.set_index("gene_name")

        if as_series:
            out = out[
                (out["domain"] == "molecular_function")
                | (out["domain"] == "biological_process")
            ]["brief_desc"]

    elif db == "DoRothEA":
        # Loading DoRothEA database of TF-target gene associations
        out = pd.read_csv(db_path, sep="\t")

    else:
        raise NotImplementedError

    return out
