"""EQTL analysis"""
from diot import Diot
from . import params, proc_factory

# pylint: disable=invalid-name

pMatrixeQTL = proc_factory(
    desc='Use matrix eQTL to call eQTLs',
    config=Diot(annotate="""
    @name:
        pMatrixeQTL
    @description:
        Call eQTLs using Matrix eQTL
    @input:
        `snpfile:file`: The genotype file, rows are snps and columns are samples
        `expfile:file`: The expression file, rows are genes
        `covfile:file`: The covariant file, columns are covariants
    @output:
        `outfile:file`: The matrix eqtl output file
    @args:
        `model`: The model to use, either modelLINEAR(default) or modelANOVA
        `pval` : The pvalue cutoff (if `cisopts.dist` > 0, will be used as pval for trans-eQTL)
            - Set this to 0 and `cisopts.cispv` to do ciseqtl analysis only.
        `fdr`  : Calculate FDR or not (default: True)
        `cisopts`: Options for calling cis-, trans-eQTL
            - `snppos` : The snp position file (columns are: snp, chr, pos)
            - `genepos`: The gene position file (columns are: gene, chr, start, end)
            - `dist`   : The distance to define cis-eQTL. (default: 0 (don't do cis-, trans- calling)
            - `cispv`  : The pvalue cutoff for cis-eQTL (`pval` will not work). Default: `1e-3`
    @requires:
        [`Matrix-eQTL (R)`](http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/)
    """),
    input='snpfile:file, expfile:file, covfile:file',
    output='''outfile:file:{{i.snpfile | fn}}-{{i.expfile | fn}}.eqtl.txt,
              cisfile:file:{{i.snpfile | fn}}-{{i.expfile | fn}}.ciseqtl.txt''',
    lang=params.Rscript.value,
    args=Diot(
        model='modelLINEAR',  # or modelANOVA
        pval=1e-5,
        fdr=True,
        cisopts=Diot(
            snppos='',
            genepos=params.refgene.value,
            dist=0,  # 0 don't do cis-, trans- calls
            cispv=1e-3
        )
    )
)
