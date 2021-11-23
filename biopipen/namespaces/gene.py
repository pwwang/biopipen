"""Gene related processes"""

from ..core.proc import Proc
from ..core.config import config

class GeneNameConversion(Proc):
    """Convert gene names back and forth using MyGeneInfo

    Input:
        infile: The input file with original gene names

    Output:
        outfile: The output file with converted gene names

    Envs:
        inopts: Options to read `in.infile` for `pandas.read_csv()`
            See https://pandas.pydata.org/docs/reference/api/pandas.read_csv.html
        outopts: Options to write `out.outfile` for `pandas.to_csv()`
            See https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.to_csv.html
        notfound: What to do if a conversion cannot be done.
            use-query: Ignore the conversion and use the original name
            skip: Ignore the conversion and skip the entire row in input file
            error: Report error
        genecol: The index (0-based) or name of the column where
            genes are present
        output: How to output
            keep: Keep the original name column and add new converted columns
            drop: Drop the original name column, and add the converted names
            replace: Drop the original name column, and insert
                the converted names at the original position
            only: Only keep the query and the converted name columns
        infmt: What's the original gene name format
            Available fields
            https://docs.mygene.info/en/latest/doc/query_service.html#available-fields
        outfmt: What's the target gene name format
        species: Limit gene query to certain species.
            Supported: human, mouse, rat, fruitfly, nematode, zebrafish,
            thale-cress, frog and pig
    """
    input = "infile:file"
    output = "outfile:file:{{in.infile | basename}}"
    lang = config.lang.python
    envs = {
        "inopts": {"sep": "\t", "index_col": False},
        "outopts": {"sep": "\t", "index": False},
        "notfound": "error",
        "genecol": 0,
        "output": "keep",
        "infmt": ["symbol", "alias"],
        "outfmt": "symbol",
        "species": "human",
    }
    script = "file://../scripts/gene/GeneNameConversion.py"
