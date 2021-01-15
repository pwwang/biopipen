"""Analysis of expression data from RNA-seq"""

from diot import Diot
from .utils import dirpat2name, fs2name
from . import opts, proc_factory

# pylint: disable=invalid-name

pExprDir2Matrix = proc_factory(
    desc='Merge expression files to a matrix.',
    config=Diot(annotate="""
    @name:
        pExprDir2Matrix
    @description:
        Convert expression files to expression matrix
        File names will be used as sample names (colnames)
        Each gene and its expression per line.
        Suppose each expression file has the same rownames and in the same order.
    @input:
        `indir:file`:  the directory containing the expression files, could be gzipped
    @output:
        `outfile:file`: the expression matrix file
        `outdir:dir`:   the directory containing expr file and plots
    @args:
        `pattern` : The pattern to filter files. Default `'*'`
        `namefunc`: Transform filename (no extension) as column name. Default: "function(fn) fn"
        `header`  : Whether each expression file contains header. Default: `False`
        `exrows`  : Rows to be excluded, regular expression applied. Default: `["^Sample", "^Composite", "^__"]`
        `boxplot` : Whether to plot a boxplot. Default: False
        `heatmap` : Whether to plot a heatmap. Default: False
        `histplot`: Whether to plot a histgram. Default: False
        `devpars` : Parameters for png. Default: `{'res': 300, 'width': 2000, 'height': 2000}`
        `boxplotggs`: The ggplot parameters for boxplot. Default: `['r:ylab("Expression")']`
            - See ggplot2 documentation.
        `heatmapggs`: The ggplot parameters for heatmap. Default: `['r:theme(axis.text.y = element_blank())']`
        `histplotggs`: The ggplot parameters for histgram. Default: `['r:labs(x = "Expression", y = "# Samples")']`
    """),
    lang=opts.Rscript,
    input="indir:file",
    output=[
        "outfile:file:{{i.indir, args.pattern | dirpat2name}}.dir/"
        "{{i.indir, args.pattern | dirpat2name}}.expr.txt",
        "outdir:dir:{{i.indir, args.pattern | dirpat2name}}.dir"
    ],
    args=Diot(
        pattern='*',
        fn2sample='function(fn) fn',
        exrows=["^Sample", "^Composite", "^__"],
        hmrows=500,
        plot=Diot(boxplot=False, heatmap=False, histogram=False),
        devpars=Diot(res=300, width=2000, height=2000),
        ggs=Diot(
            boxplot=Diot(ylab={0: "Expression"}),
            heatmap=Diot(theme={'axis.text.y': 'r:element_blank()'}),
            histogram=Diot(labs={
                'x': 'Expression',
                'y': '# Genes'
            })
        )
    ),
    envs=Diot(dirpat2name=dirpat2name)
)

pExprFiles2Mat = proc_factory(
    desc='Merge expression to a matrix from single samples.',
    config=Diot(annotate="""
    @name:
        pExprFiles2Mat
    @description:
        Merge expression to a matrix from single samples.
    @input:
        `infiles:files`: The expression file from single samples, typically with 2 columns: gene and expression value
    @output:
        `outfile:file`: the expression matrix file
    @args:
        `fn2sample`: Transform filename (no extension) as column name. Default: "function(fn) fn"
        `inopts`   : Options to read input files. Default: `Diot(rname = True, cnames = True)`
    """),
    input='infiles:files',
    output='outfile:file:{{i.infiles | fs2name}}.expr.txt',
    lang=opts.Rscript,
    args=Diot(
        inopts=Diot(cnames=True, rnames=True),
        fn2sample='function(fn) unlist(strsplit(fn, ".", fixed=T))[1]',
    ),
    envs=Diot(fs2name=fs2name)
)

pExprStats = proc_factory(
    desc='Expression profile statistics',
    config=Diot(annotate="""
    @input:
        infile: The expression matrix (rows are genes and columns are samples).
        gfile: The sample information file. Determines whether or not to do subgroup or batch stats.
            - If not provided, will do for all samples
    @output:
        outdir: The directory containing the plots
            - If `args.filter` is given, a filtered expression matrix will be generated in `outdir`.
    @args:
        inopts (Diot): Options to read `infile`.
        tsform (string): An R function in string to transform the expression matrix (i.e take log).
        annocols (list|str): Annotation columns, either column names or indexes (1-based)
            - Will be excluded in the stats.
        filter (str): An R function in string to filter the expression data before stats.
        plot (Diot): Which plot to do?
        ggs (Diot): The ggs for each plot.
        params (Diot): The params for each ggplot function.
        devpars (Diot): Parameters for png. Default: `{'res': 300, 'width': 2000, 'height': 2000}`
    """),
    input='infile:file, gfile:file',
    output='outdir:dir:{{i.infile | fn2}}.stats',
    lang=opts.Rscript,
    args=Diot(
        inopts=Diot(cnames=True, rnames=True, dup="ignore"),
        annocols=None,
        tsform=None,
        filter=None,
        plot=Diot(boxplot=True,
                  histogram=True,
                  qqplot=True,
                  violin=True,
                  pca=True),
        ggs=Diot(boxplot=Diot(ylab={0: "Expression"}),
                 violin=Diot(ylab={0: "Expression"}),
                 histogram=Diot(labs={
                     'x': 'Expression',
                     'y': '# Genes'
                 }),
                 pca=Diot(),
                 qqplot=Diot()),
        params=Diot(boxplot=Diot(),
                    histogram=Diot(),
                    qqplot=Diot(),
                    pca=Diot(),
                    violin=Diot()),
        devpars=Diot(res=300, width=2000, height=2000),
    )
)

pBatchEffect = proc_factory(
    desc='Try to remove batch effect of expression data.',
    config=Diot(annotate="""
    @input:
        expre: The expression file, generated by pExprDir2Matrix
        batche: The batch file defines samples and batches.
    @output:
        outfilee: the expression matrix file
        outdir: the directory containing expr file and plots
    @args:
        tool (str)   : The tool used to remove batch effect.
        inopts (Diot): Options to read the input file
        inlog (bool): Whether the input values has already been in log scale.
    """),
    input="expr:file, batch:file",
    output="outfile:file:{{i.expr | fn2}}-{{i.batch | fn}}/{{i.expr | fn2}}.expr.txt, \
            outdir:dir:{{i.expr | fn2}}-{{i.batch | fn}}",
    lang=opts.Rscript,
    args=Diot(
        tool='combat',
        inlog=False,
        inopts=Diot(cnames=True, rnames=True, delimit="\t", dup="drop"),
        params=Diot({'mod': None, 'par.prior': True, 'prior.plots': False})
    )
)

pUnitConversion = proc_factory(
    desc='Convert RNAseq data in different units back and forth',
    config=Diot(annotate="""
    @description:
        Convert expression between units
        See here: https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/ and
        https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/#fpkm
        Available converstions:
        - `count -> cpm, fpkm/rpkm, fpkm-uq/rpkm-rq, tpm, tmm`
        - `fpkm/rpkm -> count, tpm, cpm`
        - `tpm -> count, fpkm/rpkm, cpm`
        - `cpm -> count, fpkm/rpkm, tpm`
        NOTE that during some conversions, `sum(counts/effLen)` is approximated to `sum(counts)/sum(effLen) * length(effLen))`
    @input:
        `infile`: the expression matrix
            - rows are genes, columns are samples
    @output:
        `outfile`: the converted expression matrix
    @args:
        `inunit` : The unit of input expression values.
        `outunit`: The unit of output expression values.
        `meanfl` : A file containing the mean fragment length for each sample by rows, without header.
            - Or a fixed universal estimated number (1 used by TCGA).
        `nreads`:    The estimatied total number of reads for each sample.
            - Or you can pass a file with the number for each sample by rows, without header.
            - In converting `fpkm/rpkm -> count`: it should be total reads of that sample
            - In converting `cpm -> count`: it should be total reads of that sample
            - In converting `tpm -> count`: it should be total reads of that sample
            - In converting `tpm -> cpm`: it should be total reads of that sample
            - In converting `tpm -> fpkm/rpkm`: it should be `sum(fpkm)` of that sample
        `inform`:    Transform the input to the unit specified. (sometimes the values are log transformed)
            - For example, if the `inunit` is `tpm`, but it's actually `log2(expr+1)` transformed,
            - Then this should be `function(expr) 2^(expr - 1)`.
        `outform`: Transform to be done on the output expression values.
        datacols (list): The data columns of the `i.infile`.
            - If `None`, all columns will be used
            - Can be columns names or indexes (1-based, without rownames)
            - Note that the annotation columns (columns other than data) will be put before the transformed data in output file.
        `refexon`: the exome gff file, for RPKM/FPKM
            - `gene_id` is required for gene names
    @requires:
        [edgeR](https://bioconductor.org/packages/release/bioc/html/edger.html) if cpm or rpkm is chosen
        [coseq](https://rdrr.io/rforge/coseq/man/transform_RNAseq.html) if tmm is chosen
    """),
    lang=opts.Rscript,
    input='infile:file',
    output="outfile:file:{{i.infile | fn2}}.{{args.outunit}}.txt",
    args=Diot(inunit='count',
              outunit='tpm',
              inopts=Diot(rnames=True, cnames=True),
              meanfl=1,
              nreads=50000000,
              datacols=None,
              refexon=opts.refexon,
              inform=None,
              outform=None)
)

pRNASeqDEG = proc_factory(
    desc='Detect DEGs from RNA-seq data.',
    config=Diot(annotate="""
    @input:
        efile: The expression matrix
            - Columns other than samples in gfile will be used as annotations
            - See `args.mapping`
        gfile: The group information
            - For example:
                ```
                Sample	[Patient	]Group
                sample1	[patient1	]group1
                sample2	[patient1	]group1
                sample3	[patient2	]group1
                sample4	[patient2	]group2
                sample5	[patient3	]group2
                sample6	[patient3	]group2
                ```
            - If it is not paired comparison, you should not include Patient.
        meta: The meta information for genes
            - If provided, `args.meta` will be ignored
            - See `args.meta`
    @output:
        outfile:file: The DEG list
        outdir:file:  The output directory containing deg list and plots
    @args:
        tool  : The tool used to detect DEGs. Default: 'deseq2' (edger is also available).
        inopts: Options to read `infile`. Default: `Diot(cnames = True, rnames = True)`
        cutoff: The cutoff used to filter the results. Default: `0.05`
            - `0.05` implies `{"by": "p", "value": "0.05", "sign": "<"}`
        ggs   : The ggs for each plot. Default:
            - For heatmap: should be the `draw` argument from `plot.heatmap2` in `plot.r`
            - Not available for `mdsplot`.
            - Others are empty `Diot()`s
            - To disable a plot: `ggs.heatmap = FALSE`
        params: Parameters for each plot. Default:
            - `volplot`: `Diot(pcut = 0.05, logfccut = 2)`
            - `maplot` : `Diot(pcut = 0.05)`
            - `heatmap`: `Diot(ngenes = None, <other arguments for plot.heatmap2's params>)`, all genes in heatmap or a number for up/down genes in heatmap
        devpars: Parameters for png. Default: `{'res': 300, 'width': 2000, 'height': 2000}`
        meta: Meta information for genes, assume genes are used as rownames. This could be:
            - A column names or indexs (1-based, not including rownames) to specify the meta information in the `i.efile` matrix
            - A file with the same row names as `i.efile`
        gscol: The gene symbol column, used to show in some plots.
            - A column name or a index relative to meta.
            - If not provided, assuming rownames of `i.efile`
    """),
    lang=opts.Rscript,
    input="efile:file, gfile:file, meta",
    output=["outfile:file:{{i.efile | stem | stem}}-{{i.gfile | stem}}.DEGs"
            "/{{i.efile | stem | stem}}-{{i.gfile | stem}}.degs.xls",
            "outdir:dir:{{i.efile | stem | stem}}-{{i.gfile | stem}}.DEGs"],
    args=Diot(tool='deseq2',
              inopts=Diot(cnames=True, rnames=True, dup='ignore'),
              gscol=None,
              meta=None,
              cutoff=0.05,
              ggs=Diot(mdsplot=Diot(),
                       volplot=Diot(),
                       maplot=Diot(),
                       heatmap=Diot(),
                       qqplot=Diot(labs={
                           'x': 'Expected',
                           'y': 'Observed -log10(PValue)'
                       })),
              params=Diot(volplot=Diot(logfccut=2),
                          maplot=Diot(),
                          heatmap=Diot(ngenes=None, show_row_names=False)),
              devpars=Diot(res=300, width=2000, height=2000))
)

pCoexp = proc_factory(
    desc="Get co-expression of gene pairs in the expression matrix.",
    config=Diot(annotate="""
    @name:
        pCoexp
    @description:
        Get co-expression of gene pairs in the expression matrix.
    """),
    input="infile:file",
    output="outfile:file:{{i.infile | fn}}.coexp, \
            outpval:file:{{i.infile | fn}}.pval",
    lang=opts.Rscript,
    args=Diot(
        method='pearson',
        pval=False,
    )
)

pExprSimulate = proc_factory(
    desc="Simulate expression values",
    config=Diot(annotate="""
    @input:
        seed: The seed for random initialization. Needs an integer.
    @output:
        outfile: The simulated gene expression matrix
    @args:
        nsamples (int): The number of samples
        ngenes (int): The number of genes
        slabel (str): The prefix for Sample Label
        glabel (str): The prefix for Gene Label
    """),
    input='seed:var',
    output='outfile:file:exprsim.{{i.seed | ?isinstance:int | !:"noseed"}}.txt',
    lang=opts.Rscript,
    args=Diot(
        nsamples=100,
        ngenes=1000,
        slabel='Sample',
        glabel='Gene',
        params=Diot(),
    )
)
