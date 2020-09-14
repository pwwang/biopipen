"""Rank calculations"""
from diot import Diot
from . import opts, proc_factory

# pylint: disable=invalid-name

pRank = proc_factory(
    desc='Convert values to ranks',
    config=Diot(annotate="""
    @name:
        pRank
    @description:
        Convert values to ranks.
    @input:
        `infile:file`: The input file
    @output:
        `outfile:file`: The output file with ranks.
    @args:
        `na`: Where to put the `NA` values.
            - `"first"` : Put `NA` first
            - `"last"`  : Put `NA` last (default)
            - `"remove"`: Remove `NA` values
            - `"keep"`  : keep `NA` values
        `tie`: How to deal with ties
            - `"average"` : Use average ranks (default)
            - `"first"`   : Use the ranks come first
            - `"last"`    : Use the ranks come last
            - `"random"`  : Use the random ranks
            - `"max"`     : Use the max ranks
            - `"min"`     : Use the min ranks
        `byrow`: Calculate ranks by row (instead of by column)? Default: `True`
        `reverse`: Take the reverse rank? Default: `True`
            - Large number gets higher rank (smaller rank index)
            - `args.na` remains the same.
        `inopts`: The input options:
            - `cnames`: Whether the input file has header. Default: `True`
            - `rnames`: Whether the input file has row names. Default: `True`
            - `delimit`: The separator of columns. Default: `\t`
    """),
    input='infile:file',
    output='outfile:file:{{i.infile | fn}}.rank.txt',
    lang=opts.Rscript,
    args=Diot(
        na='last',  # keep,         first,   remove
        tie='average',  # "average", "first", "last", "random", "max", "min"
        byrow=True,  # else by column
        reverse=True,  # large number ranks higher
        inopts=Diot(cnames=True, rnames=True, delimit="\t")
    )
)

pRankProduct = proc_factory(
    desc='Calculate the rank product of a set of ranks.',
    config=Diot(annotate="""
    @name:
        pRankProduct
    @description:
        Calculate the rank product of a set of ranks. Refer to [here](https://en.wikipedia.org/wiki/Rank_product)
    @input:
        `infile:file`: The input file
            - Format:
            ```
                        Case1	Case2	...
            Feature1	8.2  	10.1 	...
            Feature2	2.3  	8.0  	...
            ...
            ```
            - Or instead of values, you can also have ranks in the input file:
            ```
                        Rank1	Rank2	...
            Feature1	2    	1    	...
            Feature2	3    	2    	...
            ...
            ```
    @output:
        `outfile:file`: The output file with original ranks, rank products and p-value if required
    @args:
        `informat`: The input format of the values. Whether they are real values (value) or ranks (rank). Default: value
            - Records will be ordered descendingly by value (Larger value has higher rank (lower rank index)).
        `pval`:     Whether to calculate the p-value or not. Default: True
        `plot`:     Number of rows to plot. Default: 0 (Don't plot)
        `cex`:      Font size for plotting. Default: 0.9
        `cnheight`: Colname height. Default: 80
        `rnwidth`:  Rowname width. Default: 50
        `devpars`:  device parameters for the plot. Default: `Diot(res=300, width=2000, height=2000)`
        `inopts`:   Options for reading the input file. Default: `Diot(cnames=True, rnames=True, delimit="\t")`
    """),
    input="infile:file",
    output="outdir:dir:{{i.infile | fn}}.rp",
    lang=opts.Rscript,
    args=Diot(
        informat="value",
        pval=True,
        cex=.9,
        cnheight=80,
        rnwidth=50,
        inopts=Diot(cnames=True, rnames=True, delimit='\t'),
        devpars=Diot(res=300, width=2000, height=2000),
    )
)
