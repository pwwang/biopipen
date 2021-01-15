"""Some statistic processes"""
from diot import Diot
from . import opts, proc_factory

# pylint: disable=invalid-name

pStats = proc_factory(
    desc='Data statistics',
    config=Diot(annotate="""
    @name:
        pStats
    """),
    lang=opts.Rscript,
    input='infile:file',
    output='outdir:dir:{{i.infile | stem}}.stats',
    args=Diot(inopts=Diot(rnames=True,
                          cnames=True,
                          delimit="\t"),
              types=Diot(),
              groups=[],
              ignore=[],
              devpars=Diot(res=300,
                           height=2000,
                           width=2000))
)

pAdjust = proc_factory(
    desc='Calcuate adjusted pvalues',
    config=Diot(annotate="""
                @input:
                    infiles: The input files with pvalues
                        - They must have pvalues at the same column
                @output:
                    outfile: The output file with all input files concatenated and last column the adjusted pvalue
                @args:
                    inopts (Diot): The options to read the input
                    pcol (int|str): The column of the pvalue, could be the column index(1-based)
                        - The index is used to subset the input to get the pvalues in `R`
                        - Note if `args.inopts.rnames` is `True`
                    method (str): The method used to adjust the pvalues
                        - See `?p.adjust` in `R`
                """),
    input='infiles:files',
    output='outfile:file:{{i.infiles | ?commonprefix | ![0] | stem2 | =commprefix}}.padj.txt',
    lang=opts.Rscript,
    args=Diot(method='BH', pcol=0, inopts=Diot(rnames=False, cnames=True))
)

pDeCov = proc_factory(
    desc='Remove covariate effects using linear regression',
    config=Diot(annotate="""
                @description:
                    Adjust covariate effects using linear regression.
                    Residues will be returned.
                @input:
                    infile: The input data, features as rows and samples as columns
                    covfile: The covariate data, covariate as columns and samples as rows
                @output:
                    outfile: The output file with the residues
                @args:
                    tool (str): Tool used to do the job.
                        - linear: using linear model
                        - peer: using peer, see: https://github.com/PMBio/peer/wiki/Tutorial
                    nthread (int): Number of threads to use for regression
                        - Only available for tool `linear`
                    nk (int): Number of hidden confounders, only for tool `peer`
                """),
    input='infile:file, covfile:file',
    output='outfile:file:{{i.infile | stem2}}.decov.txt',
    lang=opts.Rscript,
    args=Diot(nthread=1, tool='linear', nk=10)
)

pMetaPval = proc_factory(
    desc='Calculate meta p-values.',
    config=Diot(annotate="""
    @name:
        pMetaPval
    @description:
        Calculate a meta-pvalue using different methods
    @input:
        `infile:file`: The infile containing multiple pvalues for each entries.
            - Could be two types (see `args.intype`)
            - `matrix`: A matrix with rows as entries (rownames are optional), columns as cases
                ```
                        Case1   Case2   ...  CaseN
                Entry1  1.33e-2 NA      ...  1.77e-10
                Entry2  2.66e-2 4.22e-5 ...  1.71e-3
                ... ...
                EntryM  NA      0.00013 ...  4.11e-3
                ```
            - `melt`: Rows are entries from each case
            `args.rnames` should be `False`, but entry names are required in the first column,
            column names are optional (have to be 3 columns)
                ```
                Entry   Pvalue   Case
                Entry1  1.33e-2  Case1
                Entry1  1.77e-10 CaseN
                Entry2  2.66e-2  Case1
                Entry2  4.22e-5  Case2
                Entry2  1.71e-3  CaseN
                ... ...
                EntryM  0.00013  Case2
                EntryM  4.11e-3  CaseN
                ```
    @output:
        `outfile:file`: The output file containing the meta-pvalues. Default: `{{i.infile | fn}}.meta{{i.infile | ext}}`
    @args:
        `intype`: The type of the input file. Default: `matrix` (see `i.infile`)
        `inopts`: The input options to read the input file. Default: `Diot(rnames = True, cnames = True)`
        `method`: The method used to calculate the meta-pvalue. Default: sumlog (Fisher's method, aka `fisher`)
            - Other available methods: logitp, sumz, votep, sump, meanp and wilkinsonp
            - See: https://www.rdocumentation.org/packages/metap/versions/0.8
        `na`    : How to deal with `NA` p-values. Default: `skip` (just don't count it)
            - Or a numeric value to replace it with (e.g.: `1`).
    @requires:
        [`r-matep`](https://www.rdocumentation.org/packages/metap/)
    """),
    input='infile:file',
    output='outfile:file:{{i.infile | fn}}.meta{{i.infile | ext}}',
    lang=opts.Rscript,
    args=Diot(
        intype='matrix',
        inopts=Diot(rnames=True, cnames=True),
        method='sumlog',
        na='skip',#0, 1,
    )
)

pSurvival = proc_factory(
    desc='Correlation Coefficients between variables',
    config=Diot(annotate="""
    @input:
        infile: The input file (header is required).
            - col1: rownames if args.inopts.rnames = True
            - col2: the survival time
            - col3: the status. 0/1 for alive/dead or 1/2 for alive dead
            - col4: var1.
            - ... other variables
            - Note that covariales should be put in covfile, each variable in this file will be analyzed separately.
        covfile: The covariate file. If not provided, use `args.covfile`
    @output:
        outfile: The outfile containing the pvalues and other details.
        outdir : The output directory containing the pval files and plots, including:
            - `outfile`
            - cutplot: `<feature>.cut.png`, only available when `args.cutting = 'maxstat'`
            - forest plot: `<feature>.forest.png`, only available using cox
            - data with subpopulations: `<feature>.data.txt`
    @args:
        `inunit`    : The time unit in input file. Default: days
        `outunit`   : The output unit for plots. Default: days
        `inopts`    : The options for input file
            - `rnames`: Whether input file has row names.
        `combine`   : Whether combine variables in the same plot.
            - `nrow`: The number of rows. Default: 1
            - `ncol`: The number of cols. Default: 1
        `devpars`   : The device parameters for png.
            - The height and width are for each survival plot.
            - If args.combine is True, the width and height will be multiplied by `max(combine.ncol, combine.nrow)`
        `covfile`   : The covariant file. Require rownames in both this file and input file.
        `ggs`       : Extra ggplot2 elements for main plot. `ggs$table` is for the risk table.
        `params`    : The params for `ggsurvplot`. Extra arguments:
            - cut: Method to cut the continous value into subpopulations. One of mean, median, q25, q75, asis and maxstat (default).
            - cutlabels: Labels for the subpopulations (from small to large).
        `method`    : The method to do survival analysis. km (Kaplan Merier) or cox.
    @requires:
        [`r-survival`](https://rdrr.io/cran/survival/)
        [`r-survminer`](https://rdrr.io/cran/survminer/)
    """),
    lang=opts.Rscript,
    input='infile:file, covfile:file',
    output=[
        'outfile:file:{{i.infile | fn2}}.dir/{{i.infile | fn2}}.survival.txt',
        'outdir:dir:{{i.infile | fn2}}.dir'
    ],
    args=Diot(
        inunit='days',  # months, weeks, years
        outunit='days',
        method='cox',  # or km
        covfile=None,
        inopts=Diot(rnames=True),
        # params for arrange_ggsurvplots.
        # Typically nrow or ncol is set.
        # If combine.ncol = 3, that means {ncol: 3, nrow: 1}.
        # If ncol is not set, then it defaults to 1.
        # If empty, the figures will not be combined
        combine=Diot(ncol=2),
        devpars=Diot(res=300, height=2000, width=2000),
        params=Diot(),
        ggs=Diot(table=Diot())
    )
)

pPostSurvival = proc_factory(
    desc="Post survival analysis: statistics on variables in different groups",
    config=Diot(annotate="""
    @name:
        pPostSurvival
    @description:
        Statistic comparison between groups after survival analysis.
    @input:
        `infile:file`: The result file from `pSurvival`
        `survfile:file`: The survival data. See format of infile of `pSurvival`
    @output:
        `outfile:file`: The output excel file.
    @args:
        `covfile`: The covariant file. Require rownames in both this file and input file.
        `methods`: A list of testing methods
            - `wilcox`: Wilcox rank sum test
            - `t`: t-test
            - `chisq`: chisquare-test
        `inopts`: The input options for `i.survfile`.
            - `rnames`: whether the file has row names. This has to be True if `args.covfile` provided.
    """),
    input='infile:file, survfile:file',
    output='outfile:file:{{i.infile | fn2}}.stats.xlsx',
    lang=opts.Rscript,
    args=Diot(
        chi2n=10,
        inopts=Diot(rnames=True),
        covfile=None,
    )
)

pBin = proc_factory(
    desc="Bin the data",
    config=Diot(annotate="""
    @name:
        pBin
    @description:
        Bin the data in columns.
    @input:
        `infile:file`: The input file
    @output:
        `outfile:file`: The output file. Default: `{{i.infile | stem}}.binned{{i.infile | ext}}`
    @args:
        `inopts`: The input options.
            - `delimit`: The delimiter. Default: `\t`
            - `rnames`: Whether input file has row names. Default: `False`
            - `cnames`: Whether input file has column names. Default: `True`
            - Other arguments available for `read.table`
        `binopts`: The default bin options:
            - `nbin`: Number of bins.
            - `step`: The step of binning.
            - `nan`:  What to do if the value is not a number. Default: `skip`
                - `skip/keep`: Keep it
                - `as0`: Treat it as 0
            - `out`: The out value. Default: `step`
                - `step`: Use the step breaks
                - `lower/min`: Use the min value of the records in the bin
                - `upper/max`: Use the max value of the records in the bin
                - `mean`: Use the mean value of the records in the bin
                - `median`: Use the median value of the records in the bin
                - `binno`: Use the bin number (empty bins will be skipped).
        `cols`: The detailed bin options for each column.
            - If not provided (`None`), all columns will use `binopts`.
            - If column specified, only the specified column will be binned.
            - Column indices can be used. It's 1-based.
    """),
    input='infile:file',
    output='outfile:file:{{i.infile | stem}}.binned{{i.infile | ext}}',
    lang=opts.Rscript,
    args=Diot(
        inopts=Diot(delimit="\t", rnames=False, cnames=True),
        binopts=Diot(
            nbin=None, step=None, nan='skip',
            out='step',  # lower/min, upper/max, mean, median, binno
        ),
        cols=None
    )
)

pQuantileNorm = proc_factory(
    desc='Do quantile normalization',
    config=Diot(annotate="""
    @name:
        pQuantileNorm
    @description:
        Do quantile normalization
    @input:
        `infile:file`: The input matrix
    @output:
        `outfile:file`: The output matrix. Default: `{{i.infile | bn}}`
    """),
    input='infile:file',
    output='outfile:file:{{i.infile | bn}}',
    lang=opts.Rscript,
    args=Diot(inopts=Diot(rnames=True,
                          cnames=True,
                          delimit="\t",
                          skip=0))
)
# pQuantileNorm.envs.rimport = rimport

pChiSquare = proc_factory(
    desc="Do chi-square test.",
    config=Diot(annotate="""
    @name:
        pChiSquare
    @description:
        Do chi-square test.
    @input:
        `infile:file`: The input file.
    @output:
        `outfile:file` : The output file containing Xsquare, df, pval and method
        `obsvfile:file`: The observation matrix
        `exptfile:file`: The expectation matrix
    @args:
        `intype`: The type of the input file:
            - `count` (default): The contingency table
            ```
            #         | Disease | Healthy |
            # --------+---------+---------+
            #   mut   |   40    |   12    |
            # non-mut |   23    |   98    |
            # --------+---------+---------+
            ```
            - `raw`: The raw values:
            ```
            # Contingency table rows: Mut, Non
            # Contingency table cols: Disease, Healthy
            #
            #         | S1 | S2 | ... | Sn |
            # --------+----+----+-----+----+
            # Disease | 1  | 0  | ... | 1  |
            # Healthy | 0  | 1  | ... | 0  |
            # --------+----+----+-----+----+
            # Mut     | 1  | 0  | ... | 1  |
            # Non     | 0  | 1  | ... | 0  |
            ```
        `ctcols`: The colnames of contingency table if input file is raw values
            - You may also specify them in the head of the input file
    """),
    input="infile:file",
    output="outfile:file:{{i.infile | fn2}}.chi2.txt, \
            obsvfile:file:{{i.infile | fn2}}.obsv.txt, \
            exptfile:file:{{i.infile | fn2}}.expt.txt",
    lang=opts.Rscript,
    args=Diot(
        intype='cont',  # raw
        ctcols=''
    )
)

pFisherExact = proc_factory(
    desc="Do fisher exact test.",
    config=Diot(annotate="""
    @name:
        pFisherExact
    @description:
        Do fisher exact test.
    @input:
        `infile:file`: The input file.
    @output:
        `outfile:file` : The output file containing confInt1, confInt2, oddsRatio, pval, alternative and method.
    @args:
        `intype`: The type of the input file:
            - `count` (default): The contingency table
            ```
            #         | Disease | Healthy |
            # --------+---------+---------+
            #   mut   |   40    |   12    |
            # non-mut |   23    |   98    |
            # --------+---------+---------+
            ```
            - `raw`: The raw values:
            ```
            # Contingency table rows: Disease, Healthy
            # Contingency table cols: Mut, Non
            #
            #    | Disease Healthy | Mut  Non  |
            # ---+--------+--------+-----+-----+
            # S1 |    1   |    0   |  0  |  1  |
            # S2 |    0   |    1   |  1  |  0  |
            # .. |   ...  |   ...  | ... | ... |
            # Sn |    0   |    1   |  0  |  1  |
            #
            ```
        `ctcols`: The colnames of contingency table if input file is raw values
            - You may also specify them in the head of the input file
    """),
    input="infile:file",
    output="outfile:file:{{i.infile | fn2}}.fexact.txt",
    lang=opts.Rscript,
    args=Diot(
        intype='cont',  # raw
        ctcols=[]
    )
)

pPWFisherExact = proc_factory(
    desc="Do pair-wise fisher exact test.",
    config=Diot(annotate="""
    @name:
        pPWFisherExact
    @description:
        Do pair-wise fisher exact test.
        Commonly used for co-occurrence/mutual-exclusivity analysis.
        P-value indicates if the pairs are significantly co-occurred or mutually exclusive.
        Co-occurrence: Odds ratio > 1
        Mutual-exclusivity: Odds ratio < 1
    @input:
        `infile:file`: The input file.
    @output:
        `outfile:file` : The output file containing confInt1, confInt2, oddsRatio, pval, qval, alternative and method.
    @args:
        `intype`: The type of the input file:
            - `pairs`: The contingency table
            ```
            #
            # A+	B+	4
            # A-	B-	175
            # A+	B-	12
            # A-	B+	1
            #
            ```
            - `raw` (default): The raw values:
            ```
            #
            #    | A | B | ... | X |
            # ---+---+---+-----+---+
            # S1 | 1 | 0 | ... | 1 |
            # S2 | 0 | 1 | ... | 0 |
            # .. | 0 | 0 | ... | 1 |
            # Sn | 0 | 1 | ... | 1 |
            #
            ```
        `padj`: The p-value adjustment method, see `p.adjust.methods` in R. Default: `BH`
        `rnames`: If the input file has rownames for `raw` input type.
    """),
    input="infile:file",
    output="outfile:file:{{i.infile | fn2}}.pwfexact.txt",
    lang=opts.Rscript,
    args=Diot(
        rnames=True,
        intype='raw',  # pairs,
        padj='BH',
    )
)

pMediation = proc_factory(
    desc="Do mediation analysis.",
    config=Diot(annotate="""
    @input:
        `infile:file`: The input file (a matrix or data.frame). Example:
            ```
                V1   V2   V3
            S1   1    2    3
            S2   4    1    8
            ... ...
            Sn   3    3    1
            ```
        `casefile:file`: The mediation cases. Example:
            ```
            Case1   V3~V2|V1[   lm]
            Case2   V3~V1|V2[   glm]
            ```
            - No column names, but implies `Case`, `Formua`, and `Model`.
            - `\t` as delimiter.
            - This file is optional. If it is not provided, `args.case` is required.
            - If this file is provided, `args.case` is ignored
            - For different models, model for Mediator comes last. For Case2, `V2` is the mediator
    @output:
        `outfile:file`: The result file.
        `outdir:dir`  : The output directory containing output file and plots.
    @args:
        `inopts`: The options for input file.
            - `cnames`: Whether the input file has column names
            - `rnames`: Whether the input file has row names
        `medopts`: The options for mediation analysis.
            - `boot`: Use bootstrap?
            - `sims`: How many time simulations?
        `cov`: The covariate file. Default: ``
        `pval`: The pvalue cutoff. Default: `0.05`
        `fdr` : Method to calculate fdr. Use `False` to disable. Default: `True` (`BH`)
        `plot`: Parameters for `plot.mediate`? Use `False` to disable plotting. Default: `Diot()`
            - Only case with pvalue < `args.pval` will be plotted.
            - To plot all cases, use `args.pval = 1`
        `nthread`: Number of threads to use for different cases.
        `devpars`: device parameters for the plot.
        `case`   : Define cases, each case should have `model` and `fmula`.
            - If you only have one case, then it could be: `Diot(model = 'lm', fmula = 'Y~X|M')`
            In this case, `{{i.infile | fn2}}` will be used as case name
            - For multiple cases, this should be a dict of cases:
            `Diot(Case1 = Diot(model='lm', fmula='Y~X|M'), Case2 = ...)`
    """),
    lang=opts.Rscript,
    input='infile:file, casefile:file',
    output=[
        'outfile:file:{{i.infile | fn2}}.mediation/'
        '{{i.infile | fn2}}.mediation.txt',
        'outdir:dir:{{i.infile | fn2}}.mediation'
    ],
    args=Diot(
        inopts=Diot(cnames=True, rnames=True),
        medopts=Diot(boot=True, sims=500),
        cov='',
        pval=0.05,
        fdr=True,  # BH, or other methods for p.adjust,
        plot=Diot(),
        case=Diot(model='lm', fmula='Y~X|M'),
        nthread=1,
        devpars=Diot(res=300, width=2000, height=2000),
    )
)

pModeration = proc_factory(
    desc="Do moderation analysis.",
    config=Diot(annotate="""
    @input:
        `infile:file`: The input file (a matrix or data.frame). Example:
            ```
                V1   V2   V3
            S1   1    2    3
            S2   4    1    8
            ... ...
            Sn   3    3    1
            ```
        `casefile:file`: The moderation cases. Example:
            ```
            Case1   V3~V1|V2[   lm]
            Case2   V3~V1|V2[   glm]
            ```
            - No column names, but implies `Case`, `Formua`, and `Model`.
            - `\t` as delimiter.
            - This file is optional. If it is not provided, `args.case` is required.
            - If this file is provided, `args.case` is ignored
            - Moderator comes last (V2 in the above examples)
    @output:
        `outfile:file`: The result file.
        `outdir:dir`  : The output directory containing output file and plots.
    @args:
        `inopts`: The options for input file.
            - `cnames`: Whether the input file has column names
            - `rnames`: Whether the input file has row names
        `pval`: The pvalue cutoff.
        `fdr` : Method to calculate fdr. Use `False` to disable.
        `plot`: Parameters for `plot.mediate`? Use `False` to disable plotting.
            - Only case with pvalue < `args.pval` will be plotted.
            - To plot all cases, use `args.pval = 1`
        `nthread`: Number of threads to use for different cases.
        `devpars`: device parameters for the plot.
        `case`   : Define cases, each case should have `model` and `fmula`.
            - If you only have one case, then it could be: `Diot(model = 'lm', fmula = 'Y~X|M')`
            In this case, `{{i.infile | fn2}}` will be used as case name
            - For multiple cases, this should be a dict of cases:
            `Diot(Case1 = Diot(model='lm', fmula='Y~X|M'), Case2 = ...)`
    """),
    lang=opts.Rscript,
    input='infile:file, casefile:file',
    output=[
        'outfile:file:{{i.infile | fn2}}.moderation/'
        '{{i.infile | fn2}}.moderation.txt',
        'outdir:dir:{{i.infile | fn2}}.moderation'
    ],
    args=Diot(
        inopts=Diot(cnames=True, rnames=True),
        pval=0.05,
        fdr=True,  # BH, or other methods for p.adjust,
        plot=Diot(),
        case=Diot(model='lm', fmula='Y~X|M'),
        nthread=1,
        devpars=Diot(res=300, width=2000, height=2000),
    )
)

pLiquidAssoc = proc_factory(
    desc='Do liquid association analysis',
    config=Diot(annotate="""
    @name:
        pLiquidAssoc
    @description:
        Do liquid association analysis
    @input:
        `infile:file`: The input file with input data, where LA will be done on rows.
            ```
                S1   S2 ... ... Sn
            G1   1    2  ... ... 9
            G2   3    1  ... ... 1
            ... ...
            Gm   9    2  ... ... 3
            ```
        `casefile:file`: Defining the groups (X, Y, Z) and the cases. If case (3rd col) is not provided, all will be treated as one case.
            - Group "Z" is required. You can also specify group "X", then the rest will be group "Y"
            ```
            G1   X   Case1
            G2   X   Case1
            Gx   Z   Case1
            Gy   Z   Case1
            ```
    @output:
        `outfile:file`: The results of the analysis
        `outdir:dir`  : The output directory containing the result file and plots.
    @args:
        `inopts` : The options for reading input file
        `zcat`   : Whether the group "Z" is categorical. Default: `False`
            - If it is, then `stein lemma` is not suitable, we will calculate LA manually (`E(g'(z)`)
        `pval`   : The pval cutoff. Default: `0.05`
        `fdr`    : The method to calculate FDR. Use `False` to disable. Default: `True` (BH)
        `nthread`: The number of threads to use. Default: `1` (WCGNA requires)
        `plot`   : Whether do plotting or not. Default: `False`
        `devpars`: device parameters for the plot. Default: `Diot(res = 300, width = 2000, height = 2000)`
    @requires:
        [r-fastLiquidAssociation](https://github.com/pwwang/fastLiquidAssociation)
    """),
    lang=opts.Rscript,
    input='infile:file, casefile:file',
    output=[
        'outfile:file:{{i.infile | fn2}}.la/{{i.infile | fn2}}.la.txt',
        'outdir:dir:{{i.infile | fn2}}.la'
    ],
    args=Diot(
        inopts=Diot(cnames=True, rnames=True),
        zcat=False,
        pval=0.05,
        fdr=True,  # BH, or other methods for p.adjust,
        fdrfor='case',  # all,
        nthread=1,
        plot=False,
        ggs=Diot(),
        devpars=Diot(res=300, width=2000, height=2000),
    )
)

pHypergeom = proc_factory(
    desc="Do hypergeometric test.",
    config=Diot(annotate="""
    @input:
        `infile:file`: The input file, could be raw data (presence (1) and absence (0) of elements) or number of overlapped elements and elements in each category.
            - Set `args.intype` as `raw` if it is raw data. The population size `args.N` is required
            - Set `args.intype` as `numbers` (or any string except `raw`) if it is numbers. You can specified explicit header: `k` = overlapped elements, `m` = size of set 1, `n` = size of set 2 and `N` = the population size. If `N` not included, then `args.N` is required
    @output:
        `outfile:file`: The output file
    @args:
        `intype`: the type of input file. Default: `raw`. See `infile:file`
        `inopts`: The options for input file.
            - `cnames`: Whether the input file has column names
            - `rnames`: Whether the input file has row names
        `N`: The population size. Default: `None`
    """),
    input='infile:file',
    output='outfile:file:{{i.infile | fn2}}.hypergeom.txt',
    lang=opts.Rscript,
    args=Diot(
        intype='raw',  # numbers,
        inopts=Diot(cnames=True, rnames=True),
        N=None,
    )
)

pChow = proc_factory(
    desc="Do Chow-Test",
    config=Diot(annotate="""
    @description:
        Do Chow-Test
    @input:
        infile: The entire dataset used to calculate correlations. Rownames and colnames are required.
            ```
                S1  S2  S3  S4 ... Sn
            G1  1   2   1   4  ... 9
            G2  2   3   1   1  ... 3
            ... ...
            Gm  3   9   1   7  ... 8
            ```
        `colfile:file`: The sample groups, between which you want to compare the correlations.
            - It can be stacked, for example:
              ```
              S1    SG1-Group1
              S2    SG1-Group2
              S2    SG1-Group3
              S3    SG2-Group1
              S3    SG2-Group2
              ... ...
              Sn    SG2-Group3
              ```
            - Or unstacked, for example:
              ```
                    SG1     SG2
              S1    Group1  NA
              S2    Group2  Group1
              S3    Group1  Group2
              ... ...
              Sn    Group3  Group1
              ```
            - See `args.stacked`
        `casefile:file`: Assign the cases to compare. If not provided, it will do for every possible combination.
            - If `i.colfile` and `i.rowfile` are stacked, this should be a 3-column file:
              ```
              Case_1    SG1    GeneGroup1
              Case_2    SG2    GeneGroup2
              ```
            - Otherwise, this should be a 4-column file:
              ```
              Case1     SG1-Group1:SG1-Group2:SG1-Group3     Kinase:TF
              Case2     SG2-Group1:SG2-Group2   LncRNA:Gene
              ```
            - The GeneGroups can be omitted, then it will use all groups defined in row file or all possible row pairs
            - If this file is not provided, will exhaust all possible cases
        `rowfile:file`: Specify groups for rows, then the correlation will be only done within the pairs, each of which is from different groups (only 2 allowed). If not provided, it will investigate for every possible row pairs.
            - If `args.stacked` is true, it should be like:
              ```
              G1	Kinase
              G2	Kinase
              G2    LncRNA
              G3    Gene
              ... ...
              Gm	TF
              ```
            - Otherwise, it should be like:
              ```
                    GeneGroup1  GeneGroup2
              G1    Kinase      NA
              G2    Kinase      LncRNA
              G3    NA          Gene
              ... ...
              Gm    TF          NA
              ```
            - If this file is not provided, will exhaust all possible pairs
    @output:
        `outfile:file`: Significant trios (Case/ColumnGroup ~ RowGroup1 ~ RowGroup2)
        `outdir:dir`: The output directory containing plots and other output files.
    @args:
        stacked: Whether `i.rowfile` and `i.colfile` are stacked or not.
            - If they are stacked, no column names (header) should be in the file
            - Otherwise, rownames and column names are required
            - If one is stacked and the other is not, you can use a dict: `{col: False, row:True}`
        nthread: Number of threads to use.
        `inopts`: The input options for `infile`. `cnames` and `rnames` have to be `True`
        `pval`  : The pvalue cutoff to define correlation change significance. Default: `0.05`
        `fdr`   : Calculate FDR or not. Use `False` to disable. If `True` will use `BH` method, otherwise, specify the method (see `R`'s `p.adjust`).
        `plot`  : Plot the correlation for cases? Default: `False`
        `ggs`   : `ggs` items for the plot.
        `devpars`: The device parameters for the plot.
    """),
    input='infile:file, colfile:file, casefile:file, rowfile:file',
    output='outfile:file:{{i.infile | fn}}.chow/{{i.infile | fn}}.chow.txt, \
            outdir:dir:{{i.infile | fn}}.chow',
    lang=opts.Rscript,
    args=Diot(
        stacked=False,
        inopts=Diot(cnames=True, rnames=True),
        pval=0.05,
        fdr=True,
        plot=True,
        nthread=1,
        devpars=Diot(res=300, width=2000, height=2000),
        ggs=Diot(),
    )
)

pAnovaModel = proc_factory(
    desc="Run anova test on model to test the significance with/without terms",
    config=Diot(annotate="""
    @name:
        pAnovaModel
    """),

    input='infile:file, casefile:file',
    output='outfile:file:{{i.infile | fn}}.anova/{{i.infile | fn}}.anova.txt, \
            outdir:dir:{{i.infile | fn}}.anova',
    lang=opts.Rscript,
    args=Diot(
        model='lm',
        fmula=None,
        inopts=Diot(cnames=True, rnames=True),
        cov=''  # requires inopts.rnames,
    )
)

pCorr = proc_factory(
    desc='Calculate correlation coefficients.',
    config=Diot(annotate="""
    @input:
        infile: The input file of data to calculate correlations.
            - The columns are instances and rows variables when `args.byrow = True`.
            - Transposed otherwise.
        groupfile: The group file defining the groups of variables, then the correlation will be only calculated between variables in different groups.
            - First column is the variable names, and second is the groups. Only two groups allowed.
            - Overlaps allowed in different groups.
            - If this is not provided, `args.groupfile` will be checked.
            - If none of them provided, correlation will be calculated for each pair of variables.
            - We can also specify variables directly here by 'var1,var2;var3,var4'
    @output:
        outfile: The output file containing the correlation coefficients
        outdir : The output directory containing the outfile, pvalue file and the plot
    @args:
        outfmt: The output format. Could be `matrix` or `pairs`.
        method: The method used to calculate the correlation coefficient.
            - Could be `pearson`, `spearman` or `kendall`
        byrow: Whether the rows are variables or not.
        pval:  Whether output the pvalues as well.
            - Will be generated at `<outdir>/<i.infile | stem>.<args.method>.pval.txt`
        groupfile: The group file. See `in.groupfile`
        inopts: The input options:
            - `cnames`: Required if `args.byrow = False`. Whether the input file has header.
            - `rnames`: Required if `args.byrow = True`. Whether the input file has row names.
            - `delimit`: The separator of columns.
        plot   : Whether output a correlation plot.
        params : The params for `plot.heatmap` in `utils/plot.r`
        ggs    : The extra ggplot2 statements.
        devpars: The parameters for the plot device.
    @requires:
        R packages: `ggplot2` and `reshape`
    """),
    lang=opts.Rscript,
    input='infile:file, groupfile:var',
    output=[
        'outfile:file:{{i.infile | stem}}.{{args.method}}/'
        '{{i.infile | stem}}.{{args.method}}.txt',
        'outdir:dir:{{i.infile | stem}}.{{args.method}}'
    ],
    args=Diot(
        outfmt='pairs',
        method='pearson',
        byrow=True,
        pval=False,
        groupfile=None,
        inopts=Diot(cnames=True, rnames=True, delimit="\t"),
        plot=False,
        params=Diot(),  # the parameters for plot.heatmap2
        devpars=Diot(height=2000, width=2000, res=300)
    )
)

pCorr2 = proc_factory(
    desc='Calculate correlation coefficient between instances of two files',
    config=Diot(annotate="""
    @name:
        pCorr2
    @description:
        Calculate correlation coefficient between instances of two files
        Don't do it between instances within the same file.
    @input:
        `infile1:file`: The first file. See input of `pCorr`
        `infile2:file`: The second file.
            - must have same number of columns with `infile1`
    @output:
        `outfile:file`: The output file.
        `outdir:dir`  : The output directory containing output file and other files:
            - pvalues/fdr file and plots
    @args:
        `pval`  : Whether output pvalue. Default: `False`
        `fdr`   : Whether output qvalue. Default: `False`
        `outfmt`: The output format. `pairs` (default) or `matrix`
        `plot`  : Whether plot a heatmap or not. Default: `False`
        `params`: The params for `plot.heatmap` in `utils/plot.r`
        `ggs`:    The extra ggplot2 statements.
        `devpars`:The parameters for the plot device. Default: `Diot(height = 2000, width = 2000, res = 300)`
    """),
    lang=opts.Rscript,
    input='infile1:file, infile2:file',
    output=[
        'outfile:file:{{i.infile1 | fn2}}-{{i.infile2 | fn2}}.corr/'
        '{{i.infile1 | fn2}}-{{i.infile2 | fn2}}.corr.txt',
        'outdir:dir:{{i.infile1 | fn2}}-{{i.infile2 | fn2}}.corr'
    ],
    args=Diot(
        inopts1=Diot(),
        inopts2=Diot(),
        method='pearson',  # spearman, kendall,
        pval=False,
        fdr=False,  # method: 'BH',
        outfmt='pairs',  # matrix,
        plot=False,
        params=Diot(),  # the parameters for plot.heatmap,
        ggs=Diot(),  # extra ggplot statements,
        devpars=Diot(height=2000, width=2000, res=300),
    )
)

pDiffCorr = proc_factory(
    desc='Test correlation differences using Fisher Z method.',
    config=Diot(annotate="""
    @input:
        infile: The entire dataset used to calculate correlations. Rownames and colnames are required.
            ```
                S1  S2  S3  S4 ... Sn
            G1  1   2   1   4  ... 9
            G2  2   3   1   1  ... 3
            ... ...
            Gm  3   9   1   7  ... 8
            ```
        `colfile:file`: The sample groups, between which you want to compare the correlations.
            - It can be stacked, for example:
              ```
              S1    Healthy
              S2    Healthy
              S2    Healthy-1
              S3    Disease
              S3    Disease-1
              ... ...
              Sn    Disease
              ```
            - Or unstacked, for example:
              ```
                    SampleGroup1    SampleGroup2
              S1    Healthy         NA
              S2    Healthy         Healthy
              S3    Disease         Disease
              ... ...
              Sn    Disease         NA
              ```
            - See `args.stacked`
        `casefile:file`: Assign the cases to compare. If not provided, it will do for every possible combination.
            - If `i.colfile` and `i.rowfile` are stacked, this should be a 2-column file:
              ```
              SampleGroup1    GeneGroup1
              SampleGroup2    GeneGroup2
              ```
            - otherwise, it should be:
              ```
              Healthy:Disease     Kinase:TF
              Healthy-1:Disease-1   LncRNA:Gene
              ```
            - The GeneGroups can be omitted, then it will use all groups defined in row file or all possible row pairs
            - If this file is not provided, will exhaust all possible cases
        `rowfile:file`: Specify groups for rows, then the correlation will be only done within the pairs, each of which is from different groups (only 2 allowed). If not provided, it will investigate for every possible row pairs.
            - If `args.stacked` is true, it should be like:
              ```
              G1	Kinase
              G2	Kinase
              G2    LncRNA
              G3    Gene
              ... ...
              Gm	TF
              ```
            - Otherwise, it should be like:
              ```
                    GeneGroup1  GeneGroup2
              G1    Kinase      NA
              G2    Kinase      LncRNA
              G3    NA          Gene
              ... ...
              Gm    TF          NA
              ```
            - If this file is not provided, will exhaust all possible pairs
    @output:
        `outfile:file`: The pairs under different cases that their correlations have been changed significantly
        `outdir:dir`: The output directory containing plots and other output files.
    @args:
        stacked: Whether `i.rowfile` and `i.colfile` are stacked or not.
            - If they are stacked, no column names (header) should be in the file
            - Otherwise, rownames and column names are required
            - If one is stacked and the other is not, you can use a dict: `{col: False, row:True}`
        nthread: Number of threads to use.
        `inopts`: The input options for `infile`. `cnames` and `rnames` have to be `True`
        `method`: The method to calculate the correlation. Default: `pearson`
        `pval`  : The pvalue cutoff to define correlation change significance. Default: `0.05`
        `fdr`   : Calculate FDR or not. Use `False` to disable. If `True` will use `BH` method, otherwise, specify the method (see `R`'s `p.adjust`).
        `plot`  : Plot the correlation for cases? Default: `False`
        `ggs`   : `ggs` items for the plot.
        `devpars`: The device parameters for the plot.
    """),
    lang=opts.Rscript,
    input='infile:file, colfile:file, casefile:file, rowfile:file',
    output=[
        'outfile:file:{{i.infile | fn2}}.diffcorr/'
        '{{i.infile | fn2}}.diffcorr.txt',
        'outdir:dir:{{i.infile | fn2}}.diffcorr'
    ],
    args=Diot(
        stacked=False,
        inopts=Diot(cnames=True, rnames=True),
        method='pearson',  # spearman,
        pval=0.05,
        fdr=True,  # BH,
        plot=False,
        nthread=1,
        ggs=Diot(),  # extra ggplot statements,
        devpars=Diot(height=2000, width=2000, res=300),
    )
)

pBootstrap = proc_factory(
    desc='Do bootstrapping',
    config=Diot(annotate="""
    @name:
        pBootstrap
    @description:
        Do bootstrapping resampling
    @input:
        `infile:file`: The input data file
    @output:
        `outfile:file`: The output file with the bootstrapped statistics values
            - depends on the `args.stats` function
        `outdir:dir`: The directory to save the outfile and figures.
    @args:
        `inopts`:  The options to read the input file. Default: `Diot(cnames = True, rnames = True)`
        `params`:  Other parameters for `boot` function from R `boot` package
        `nthread`: # of threads(cores) to use. Default: `1`
        `n`: Sampling how many times? Default: `1000`
        `stats`: The function to generate statistics for output. Default: `function(x) x`
            - Default to use all data
            - This function can return a multiple statistics in a vector
            - The argument `x` is the data generate for each sampling.
            - Unlink the `statistic` argument from `boot`, to make it convenient, we don't put the `index` here.
        `plot`: Plot the statistics? Default: `all` (plot all statistics)
            - You may also specify indices. For example: `[1, 2]` to plot the 1st and 2nd statistics
            - Use `False` or `None` to disable plotting
        `devpars`: The device parameters for the plot.
    """),
    lang=opts.Rscript,
    input='infile:file',
    output=[
        'outfile:file:{{i.infile | fn2}}.boot/{{i.infile | fn2}}.boot.txt',
        'outdir:dir:{{i.infile | fn2}}.boot'
    ],
    args=Diot(
        inopts=Diot(cnames=True, rnames=True),
        params=Diot(),
        nthread=1,
        n=1000,
        stats='function(x) x',
        plot='all',
        devpars=Diot(height=2000, width=2000, res=300),
    )
)

pPCA = proc_factory(
    desc='Perform PCA analysis using PCAtools',
    config=Diot(annotate="""
    @description:
        Perform PCA analysis using PCAtools.
        See: https://bioconductor.org/packages/release/bioc/vignettes/PCAtools/inst/doc/PCAtools.html
    @input:
        infile: The matrix to do the analysis
            - Columns are features, rows are samples
            - If not, you may use `args.params.transposed = True`
        metafile: The metadata file
            - Columns are features, rows are samples
            - rnames and cnames are required
    @output:
        outfile: The file with the components
        oudir  : The directory containing the output file and plots.
    @args:
        devpars (Diot)        : The parameters for plotting device.
        inopts  (Diot)        : The options to read the input files.
        na      (bool|number): How to deal with `NA` values.
            - A logistic/boolean value will remove them (use `complete.cases`)
            - Otherwise, it will be replaced by the given value.
        plots (Diot): Plotting parameters supported by `PCAtools`. Set `False` to disable a plot.
            - scree: A scree plot. See `?screeplot`
                - Since `PCAtools::screeplot` returns a ggplot object, we have `args.plots.scree.ggs` to extend the plot.
            - bi: A bi-plot. See `?biplot`
            - pairs: A pairs plot. See `?pairsplot`
            - loadings: A loadings plot. See `?plotloadings`
            - eigencor: An eigencor plot. See `?eigencorplot`
        params (Diot): Other parameters for `PCAtools::pca`
        npc (int|str): The number of PCs to write to the output file.
            - A fixed number of PCs; or
            - one of `horn` or `elbow` to determine the optimal number of PCs.
    """),
    input='infile:file, metafile:file',
    output=[
        'outfile:file:{{i.infile | stem | stem}}.pca/'
        '{{i.infile | stem | stem}}.pcs.txt',
        'outdir:dir:{{i.infile | stem | stem}}.pca'
    ],
    lang=opts.Rscript,
    args=Diot(devpars=Diot(height=2000, width=2000, res=300),
              inopts=Diot(cnames=True, rnames=True),
              params=Diot(),
              na=0,
              plots=Diot(scree=True,
                         bi=True,
                         pairs=True,
                         loadings=True,
                         eigencor=True),
              npc='elbow')
)
