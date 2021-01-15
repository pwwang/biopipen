"""A set of algorithms or models"""
from diot import Diot
from pyppl import Proc
from . import opts, proc_factory

# pylint: disable=invalid-name
pRWR = proc_factory(
    desc='Do random walk with restart (RWR).',
    config=Diot(annotate="""
    @input:
        `Wfile:file`: The adjecent matrix
        `Efile:file`: The start vector
    @output:
        `outfile:file`: The output of final probabilities
    @args:
        `c`:       The restart probability. Default: 0.1
        `eps`:     The convergent cutoff || R(i+1) - R(i) ||. Default: 1e-5
        `niter`:   Max iterations to stop. Default: 10000
        `normW`:   Weather to normalize W or not, default True.
            - Laplacian normalization is used (more to add).
        `normE`:   Weather to normalize E or not, default True.
            - E will be normalized as: E = E/sum(E)
    @requires:
        [NetPreProc]
        desc = "Package for the pre-processing and normalization of graphs."
        url = "https://cran.r-project.org/web/packages/NetPreProc/index.html"
        when = "[[ {{args.normW}} == True ]]"
        validate = '[[ $({{proc.lang}} --vanilla -e 'library(NetPreProc)' 2>&1) == \
            *"Network Pre-Processing package"* ]]'
        install = '{{proc.lang}} -e \'install.packages("NetPreProc", \
            repos="https://cran.rstudio.com")\''
    """),
    input="Wfile:file, Efile:file",
    output="outfile:file:{{i.Wfile | fn2}}.rwr.txt",
    lang=opts.Rscript,
    args=Diot(
        c=0.1,
        eps=1e-5,
        niter=10000,
        normW=True,
        normE=True,
    ))

pAR = proc_factory(
    desc='Affinity Regression.',
    config=Diot(annotate="""
    @name:
        pAR
    @description:
        Affinity Regression.
        Ref: https://www.nature.com/articles/nbt.3343
        ```
                b           c        d          d
            _________    _______    ____       ____
            |       |    |  W  |    |  |       |  |
          a |   D   |  b |_____|  c |Pt|  =  a |Y |   <=>
            |_______|               |__|       |  |
                                               |__|
        kronecker(P, YtD)*vec(W) = vec(YtY)             <=>
        X*vec(W) = vec(YtY)
        WPt:
            c           d              d
            _______    ____          _____
            |  W  |    |  |          |   |
          b |_____|  c |Pt|    =>  b |___|
                       |__|

        YtDW:
        WtDtY:
              b             a        d               d
            _______    _________   ____           _____
            |  Wt |    |       |   |  |           |   |
          c |_____|  b |   Dt  | a |Y |     =>  c |___|
                       |_______|   |  |
                                   |__|
        ```
    @input:
        `D:file` : The D matrix, with `b` as the cnames and `a` as rnames
        `Pt:file`: The Pt matrix, with `d` as the cnames and `c` as rnames
        `Y:file`:  The Y matrix, with `d` as the cnames and `a` as rnames
            - All input files could be gzipped
    @output:
        `W:file`:  The interaction matrix
        `outdir:dir`: The output directory
    @args:
        seed (int):  The seed for sampling the training set.
        tfrac (float): The fraction of samples used for training.
        nthread (int): The number of threads to use
        method (str): The solver
        nfold (int): N-fold cross validation to use
        svdP (int): Dimension to reduce using SVD on P
    """),
    input='D:file, Pt:file, Y:file',
    output=[
        'W:file:{{i.D, i.Pt, i.Y | *commonprefix '
        '                        | ? | =_ | !:stem2(i.D) }}.AR/W.txt',
        'outdir:dir:{{i.D, i.Pt, i.Y | *commonprefix '
        '                            | ? | =_ | !:stem2(i.D) }}.AR'
    ],
    lang=opts.Rscript,
    args=Diot(
        seed=None,
        tfrac=.5,
        inopts=Diot(cnames=True, rnames=True),
        svdP=0,
        predY=True,
        WPt=True,
        WtDtY=True,
        nfold=3,
        nthread=1,
        method='glmnet',  # admm,
    )
)

pColoc = proc_factory(
    desc="Bayes Factor colocalisation analyses using R `coloc` package.",
    config=Diot(annotate="""
    @description:
        Bayes Factor colocalisation analyses using R `coloc` package.
        `coloc` package can accept multiple formats of input. Here we adopt the one using pvalues.
        `coloc.abf(dataset1=list(pvalues=p1,N=nrow(X1),type="quant"), dataset2=list(pvalues=p2,N=nrow(X2),type="quant"), MAF=maf)`
    @input:
        `infile:file`: The input file including the MAF, pvalues of 1st and 2nd phenotypes
            - The first 6 columns are in BED6 format.
            - 7th : MAF
            - 8th : Pvalues for the 1st phenotype
            - 9th : Pvalues for the 2nd phenotype
            - This file could have a header with the names for phenotypes
            - Snps have to be on the same chromosome, and sorted by positions.
    @output:
        `outfile:file`: The output file including:
            - # snps, PP.H0.abf, PP.H1.abf, PP.H2.abf, PP.H3.abf and PP.H4.abf
        `outdir:dir`  : The output directory containing the output file and plots.
    @args:
        `plot`: Do manhattan plot? Default: `True`
    """),
    input='infile:file',
    output=[
        'outfile:file:{{i.infile | fn2}}.coloc/{{i.infile | fn2}}.coloc.txt',
        'outdir:dir:{{i.infile | fn2}}.coloc'
    ],
    args=Diot(
        inopts=Diot(cnames=True, rnames=False),
        plot=True,
        ggs=Diot(),
        params=Diot(),
        devpars=Diot(res=300, height=2000, width=2000),
        hifile='',
        hilabel=True,
    ),
    lang=opts.Rscript,
)
