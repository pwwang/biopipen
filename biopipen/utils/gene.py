"""Do gene name conversion"""
from mygene import MyGeneInfo
from datar.all import (
    c,
    f,
    group_by,
    desc,
    arrange,
    slice_head,
    tibble,
    left_join,
    mutate,
    is_na,
    across,
    if_else,
    filter,
    pull,
    select,
)

mygene = MyGeneInfo()


class QueryGenesNotFound(Exception):
    """When genes cannot be found"""


def gene_name_conversion(
    genes,
    species,
    infmt,
    outfmt,
    notfound,
):
    """Convert gene names using MyGeneInfo

    Args:
        genes: A sequence of genes
        species: The species to limit the query
            Supported: human, mouse, rat, fruitfly, nematode, zebrafish,
            thale-cress, frog and pig

        infmt: What's the original gene name format
            Available fields
            https://docs.mygene.info/en/latest/doc/query_service.html#available-fields
        outfmt: What's the target gene name format
        notfound: What to do if a conversion cannot be done.
            use-query: Ignore the conversion and use the original name
            skip: Ignore the conversion and skip the entire row in input file
            error: Report error

    Returns:
        A dataframe with two columns, query and `outfmt`.
    """
    out = (
        mygene.querymany(
            genes,
            scopes=infmt,
            fields=outfmt,
            as_dataframe=True,
            df_index=False,
            species=species,
        )
        >> group_by(f.query)
        >> arrange(desc(f._score))
        >> slice_head(1)
        >> select(~c(f._id, f._score, f.notfound))
    )
    if isinstance(outfmt, str):
        outfmt = [of.strip() for of in outfmt.split(",")]
    out = tibble(query=genes) >> left_join(out, by=f.query)
    if notfound == "use-query":
        out = out >> mutate(
            across(
                outfmt,
                lambda col, query: if_else(is_na(col), query, col),
                query=f.query,
            )
        )
    elif notfound == "error" and any(is_na(out[outfmt[0]])):
        nagenes = out >> filter(is_na(f[outfmt[0]])) >> pull(f.query)
        raise QueryGenesNotFound(nagenes)
    elif notfound == "skip":
        out = out >> filter(~is_na(f[outfmt[0]]))

    return out
