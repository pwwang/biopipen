"""Microarray data analysis"""
from diot import Diot
#from .utils import plot, txt, dirnamePattern
from .rnaseq import pBatchEffect, pCoexp, pExprStats # pylint: disable=unused-import
from .utils import dirpat2name
from . import opts, proc_factory

# pylint: disable=invalid-name

pCELDir2Matrix = proc_factory(
    desc='Merge expression files to a matrix.',
    config=Diot(annotate="""
    @name:
        pCELDir2Matrix
    @description:
        Convert CEL files to expression matrix
        File names will be used as sample names (colnames)
    @input:
        `indir:file`:  the directory containing the CEL files, could be gzipped
            - If you have files, then use `pFiles2Dir` first
        `sifile:File`: the sample infor file, extensions for samples are not necessary.
            - So that it can be also used by `pMArrayDEG`
    @output:
        `outfile:file`: the expression matrix file
        `outdir:dir`:   the directory containing expr file and plots
    @args:
        `pattern`  : The pattern to filter files. Default `'*'`
        `norm`     : The normalization method. Default: rma (mas5)
        `transfm`  : The extra tranformer for the expression values after the nomralization.
            - `Note the expression values have been done with log`
        `cdffile`  : The cdffile. Default: ''
        `fn2sample`: The transformer to transform file name
    """),
    input="indir:file, sifile:file",
    output="outfile:file:{{i.indir, args.pattern | *dirpat2name}}.expr.txt",
    lang=opts.Rscript,
    args=Diot(
        fn2sample='function(fn) unlist(strsplit(fn, ".", fixed = T))[1]',
        pattern='*',
        norm='rma',  # mas5,
        transfm=None,  # mas5,
        cdffile='',
    ),
    envs=Diot(dirpat2name=dirpat2name),
)

pMArrayDEG = proc_factory(
    desc='Detect DEGs from microarray data.',
    config=Diot(annotate="""
    @name:
        pMArrayDEG
    @description:
        Detect DEGs from microarray data.
    @input:
        `efile:file`: The expression matrix file
        `gfile:file`: The sample information file
    @output:
        `outfile:file`: The file with DEGs
        `outdir:dir`  : The directory containing files and plots
    @args:
        `tool`    : The tool to use. Default: `limma`
        `annofile`: The annotation file for the probes. Default: ``
            - If not provided, raw probe name will be used.
        `filter`: The filter of the data. Default: `[0, 0]`
            - 1st element: the `rowSums` of the expression matrix
            - 2nd element: how many samples(columns) have to reach the `rowSums` given by the 1st element
        `pval`  : The pval cutoff for DEGs. Default: `0.05`
            - You may specify which kind of measurement to use
            - `p:0.05`: using cutoff 0.05 for pvalue
            - `q:0.05`: using cutoff 0.05 for qvalue(fdr)
            - If no prefix is used, then it defaults to pvalue
        `plot`  : What kind of plots to generate.
            - `mdsplot` : The MDS plot
            - `volplot` : The volcano plot
            - `maplot`  : The MA plot
            - `heatmap` : The heatmap. (use `args.hmrows` to determine how many genes to plot)
        `ggs`   : The extra ggplot element for each plot (should be able to concatenate by `+`).
            - `maplot`  : for MA plot
            - `heatmap` : for heatmap. Default: `Diot(theme = {'axis.text.y': 'r:element_blank()'})`
            - `volplot` : for volcano plot
        `devpars`: The parameters for plotting device. Default: `Diot(res = 300, width = 2000, height = 2000)`
    @requires:
        `r-limma`
    """),
    lang=opts.Rscript,
    input="efile:file, gfile:file",
    output=[
        "outfile:file:{{i.efile | fn2}}-{{i.gfile | fn2}}-DEGs/"
        "{{i.efile | fn2}}-{{i.gfile | fn2}}.degs.txt",
        "outdir:dir:{{i.efile | fn2}}-{{i.gfile | fn2}}-DEGs"
    ],
    args=Diot(
        tool='limma',
        annofile='',
        filter=[0, 0],
        pval=0.05,
        hmrows=100,
        plot=Diot(mdsplot=True,
                  volplot=Diot(fccut=2, pcut=0.05),
                  maplot=False,
                  heatmap=False),
        ggs=Diot(
            maplot=Diot(),
            heatmap=Diot(theme={'axis.text.y': 'r:element_blank()'}),
            volplot=Diot(ylab={0: '-log10(p-value)'})
        ),
        devpars=Diot(res=300, width=2000, height=2000),
    )
)
