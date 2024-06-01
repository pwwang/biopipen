"""Do gene name conversion"""
from __future__ import annotations

import re
import contextlib
import pandas as pd
from mygene import MyGeneInfo

mygene = MyGeneInfo()


class QueryGenesNotFound(ValueError):
    """When genes cannot be found"""


def gene_name_conversion(
    genes: list[str],
    infmt: str | list[str],
    outfmt: str,
    dup: str = "first",
    species: str = "human",
    notfound: str = "na",
    suppress_messages: bool = False,
):
    """Convert gene names using MyGeneInfo

    Args:
        genes: A character/integer vector of gene names/ids
        species: A character vector of species names
        infmt: A character vector of input gene name formats
            See the available scopes at
            https://docs.mygene.info/en/latest/doc/data.html#available-fields
            You can use ensg as a shortcut for ensembl.gene
        outfmt: A character vector of output gene name formats
        dup: How to deal with duplicate gene names found.
            first: keep the first one (default), sorted by score descendingly
            last: keep the last one, sorted by score descendingly
            all: keep all of them, each will be a separate row
            <X>: combine them into a single string, separated by X
        notfound: How to deal with gene names that are not found
            error: stop with an error message
            use-query: use the query gene name as the converted gene name
            skip: skip the gene names that are not found
            ignore: Same as "skip"
            na: use NA as the converted gene name (default)
        suppress_messages: Suppress the messages while querying

    Returns:
        A dataframe with the query gene names and the converted gene names
        When a gene name is not found, the converted name will be "NA"
        When duplicate gene names are found, the one with the highest score will be kept
    """
    notfound = notfound.lower()
    if notfound not in ("error", "use-query", "skip", "ignore", "na"):
        raise ValueError(
            "`notfound` of `gene_name_conversion` must be one of "
            "'error', 'use-query', 'skip', 'ignore', 'na'"
        )

    if infmt in ["ensg", "ensmusg"]:
        infmt = "ensembl.gene"
    if outfmt in ["ensg", "ensmusg"]:
        outfmt = "ensembl.gene"

    orig_genes = genes[:]
    if infmt == "ensembl.gene":
        # Remove version numbers from ensembl gene ids
        genes = [re.sub("\\..*", "", gene) for gene in genes]

    query_df = pd.DataFrame({"query": genes, "orig": orig_genes})

    if suppress_messages:
        with contextlib.redirect_stdout(None):
            out = mygene.querymany(
                genes,
                scopes=infmt,
                fields=outfmt,
                species=species,
                as_dataframe=True,
                df_index=False,
            )
    else:
        out = mygene.querymany(
            genes,
            scopes=infmt,
            fields=outfmt,
            species=species,
            as_dataframe=True,
            df_index=False,
        )

    if out.shape[0] == 0:
        return pd.DataFrame({"query": genes, "converted": ["NA"] * len(genes)})

    if dup == "first":
        out = (
            out
            .sort_values("_score", ascending=False)
            .groupby("query")
            .head(1)
            .reset_index(drop=True)
        )
    elif dup == "last":
        out = (
            out
            .sort_values("_score", ascending=False)
            .groupby("query")
            .tail(1)
            .reset_index(drop=True)
        )
    elif dup != "all":
        out = (
            out
            .sort_values("_score", ascending=False)
            .groupby("query")
            .agg({outfmt: lambda x: f"{dup}".join([str(x) for x in x.unique()])})
            .reset_index()
        )

    out = pd.merge(query_df, out, on="query", how="left")
    out = out.drop(columns=["query"]).rename(columns={"orig": "query"})

    if notfound == "error":
        if out[outfmt].isnull().any():
            nagenes = out[out[outfmt].isnull()]["query"].tolist()
            raise QueryGenesNotFound(f"Query genes not found: {','.join(nagenes)}")
    elif notfound == "use-query":
        out[outfmt] = out[outfmt].combine_first(out["query"])
    elif notfound in ["skip", "ignore"]:
        out = out.dropna(subset=[outfmt])
    else:  # notfound == "na"
        out[outfmt] = out[outfmt].fillna("NA")

    return out
