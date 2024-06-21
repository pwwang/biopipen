"""Gene related processes"""

from ..core.proc import Proc
from ..core.config import config


class GeneNameConversion(Proc):
    """Convert gene names back and forth using MyGeneInfo

    Input:
        infile: The input file with original gene names
            It should be a tab-separated file with header

    Output:
        outfile: The output file with converted gene names

    Envs:
        notfound (choice): What to do if a conversion cannot be done.
            - use-query: Ignore the conversion and use the original name
            - skip: Ignore the conversion and skip the entire row in input file
            - ignore: Same as skip
            - error: Report error
            - na: Use NA
        dup (choice): What to do if a conversion results in multiple names.
            - first: Use the first name, sorted by matching score descendingly (default)
            - last: Use the last name, sorted by matching score descendingly
            - combine: Combine all names using `;` as separator
        genecol: The index (1-based) or name of the column where genes are present
        output (choice): How to output.
            - append: Add the converted names as new columns at the end using `envs.outfmt`
                as the column name.
            - replace: Drop the original name column, and insert
                the converted names at the original position.
            - converted: Only keep the converted names.
            - with-query: Output 2 columns with original and converted names.
        infmt: What's the original gene name format
            Available fields
            https://docs.mygene.info/en/latest/doc/query_service.html#available-fields
        outfmt: What's the target gene name format. Currently only a single format
            is supported.
        species: Limit gene query to certain species.
            Supported: human, mouse, rat, fruitfly, nematode, zebrafish,
            thale-cress, frog and pig
    """  # noqa: E501
    input = "infile:file"
    output = "outfile:file:{{in.infile | basename}}"
    lang = config.lang.rscript
    envs = {
        "notfound": "error",
        "genecol": 1,
        "dup": "first",
        "output": "append",
        "infmt": ["symbol", "alias"],
        "outfmt": "symbol",
        "species": "human",
    }
    script = "file://../scripts/gene/GeneNameConversion.R"


class GenePromoters(Proc):
    """Get gene promoter regions by specifying the flanking regions of TSS

    Input:
        infile: The input file with gene ids/names

    Output:
        outfile: The output file with promoter regions in BED format

    Envs:
        up (type=int): The upstream distance from TSS
        down (type=int): The downstream distance from TSS
            If not specified, the default is `envs.up`
        notfound (choice): What to do if a gene is not found.
            - skip: Skip the gene
            - error: Report error
        refgene: The reference gene annotation file in GTF format
        header (flag): Whether the input file has a header
        genecol (type=int): The index (1-based) of the gene column
        match_id (flag): Should we match the genes in `in.infile` by `gene_id`
            instead of `gene_name` in `envs.refgene`
        sort (flag): Sort the output by chromosome and start position
        chrsize: The chromosome size file, from which the chromosome order is
            used to sort the output
    """
    input = "infile:file"
    output = "outfile:file:{{in.infile | stem}}-promoters.bed"
    lang = config.lang.rscript
    envs = {
        "up": 2000,
        "down": None,
        "notfound": "error",
        "refgene": config.ref.refgene,
        "header": True,
        "genecol": 1,
        "match_id": False,
        "sort": False,
        "chrsize": config.ref.chrsize,
    }
    script = "file://../scripts/gene/GenePromoters.R"
