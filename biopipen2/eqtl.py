"""EQTL analysis"""
from diot import Diot
from . import opts, proc_factory

# pylint: disable=invalid-name

pMatrixeQTL = proc_factory(
    desc='Use matrix eQTL to call eQTLs',
    config=Diot(annotate="""
    @description:
        Call eQTLs using Matrix eQTL.
        By default, this is only calling cis-eQTLs.
        See: http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/
    @input:
        snpfile: The genotype file, rows are snps and columns are samples
        expfile: The expression file, rows are genes
        covfile: The covariant file, columns are covariants
    @output:
        outfile: The cis-eQTL file if `snppos` and `genepos` are provided. Otherwise it'll be empty.
        alleqtl: All eQTL file
    @args:
        model (str): The model to use, either modelLINEAR or modelANOVA
        pval (float): The pvalue cutoff for the results
        transp (float): The pvalue for trans-eQTLs. If this is `False`, then trans-eQTL will not be calculated.
            - Either this is set, or snppos, genepos and dist are required
            - If this is `True`, `1e-5` will be used as cutoff
        fdr (bool): Do FDR calculation or not (save memory if not).
        snppos (str): The path of the SNP positions. Could be BED, GFF, VCF or
            - custom file with `snp`, `chr`, `pos` as first 3 columns.
        genepos (str): The path of Gene positions. Could be BED or GFF file.
        dist (int): The distance to decide the SNP-gene pairs for cis-eQTL calculation.
    @requires:
        [`Matrix-eQTL (R)`](http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/)
    """),
    input='snpfile:file, expfile:file, covfile:file',
    output='''outfile:file:{{i.snpfile | fn}}-{{i.expfile | fn}}.ciseqtl.txt,
              alleqtl:file:{{i.snpfile | fn}}-{{i.expfile | fn}}.eqtl.txt''',
    lang=opts.Rscript,
    args=Diot(
        model='modelLINEAR',  # or modelANOVA
        pval=1e-3,
        transp=False,
        fdr=True,
        snppos='',
        genepos=opts.refgene,
        dist=250000,
    )
)
