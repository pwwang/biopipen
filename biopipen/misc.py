"""Some misc processes"""
from diot import Diot
from . import opts, proc_factory

# pylint: disable=invalid-name

pGEP70 = proc_factory(
    desc='Do GEP70 plot for multiple mylenoma 70-gene-signatures',
    config=Diot(annotate="""
    @name:
        pGEP70
    @description:
        Calculate GEP70 scores for multiple mylenoma 70-gene-signatures and do survival plot.
        Or add one gene to it to see the survival plot.
        see: https://www.ncbi.nlm.nih.gov/pubmed/17105813
    @input:
        `exprfile:file`: The input file with expression matrix
            - Columns are samples, rows are genes
            - make sure the expression values are log2-scale normalized
        `survfile:file`: The survival data file (header is required).
            - col1: rownames
            - col2: the survival time
            - col3: the status. 0/1 for alive/dead or 1/2 for alive dead
        `gene`: An extra gene to be added to the plot.
            - If not provided, only do the GEP70 plot.
    @output:
        `outdir:file`: The output directory, containing:
            - The survival data file.
            - The GEP70 plot and results
            - The GEP70 + gene plot and results if `i.gene` is provided.
    @args:
        `gep70`: The GEP70 genes.
            - Column 1: genes
            - Column 2: how is the regulated (up or down)
        `inunit`  : The time unit in input file. Default: days
        `outunit` : The output unit for plots. Default: days
        `params`  : The params for `ggsurvplot`. Default: `Diot({'risk.table': True, 'conf.int': True, 'font.legend': 13, 'pval': 'Log-rank p = {pval}'})`
            - You may do `ylim.min` to set the min ylim. Or you can set it as 'auto'. Default: 0.
        `ggs`     : Extra ggplot2 elements for main plot. `ggs.table` is for the risk table.
        `devpars` : The device parameters for png. Default: `{res:300, height:2000, width:2000}`
    """),
    input='exprfile:file, survfile:file, gene',
    output='outdir:dir:{{i.survfile | fn2}}.{{args.name}}{{i.gene}}',
    lang='Rscript',
    args=Diot(
        name='GEP70',
        gep70=opts.gep70,
        inunit='days',  # months, weeks, years,
        outunit='days',
        params=Diot({
            'font.legend': 13,
            'pval': 'Log-rank p = {pval}',
            'risk.table': True
        }),
        devpars=Diot(res=300, height=2000, width=2000),
        ggs=Diot(table=Diot()),
    )
)

pNCBI = proc_factory(
    desc='The NCBI E-Utils',
    config=Diot(annotate="""
    @name:
        pNCBI
    @description:
        The NCBI E-Utils
    @input:
        `term`: The term or the id argument for esearch or efetch
    @output:
        `outfile:file`: The output file
    @args:
        `prog`   : The program to use, esearch (Default) or efetch
        `apikey` : The api key for E-utils
            - Without API key, we can only query 3 time in a second
            - With it, we can do 10.
        `sleep`  : Sleep sometime after job done. Default: `0.15`
            - Because of the limit of # queries/sec, we need to sleep sometime after the job is done
            - At the same time, we also have to limit # jobs to run at the same time. typically: `pNCBI.forks = 10`
        `db`     : The database to query. Default: `pubmed`. Available databases:
            - annotinfo, assembly, biocollections, bioproject, biosample, biosystems, blastdbinfo, books,
            - cdd, clinvar, clone, dbvar, gap, gapplus, gds, gencoll, gene, genome, geoprofiles, grasp, gtr,
            - homologene, ipg, medgen, mesh, ncbisearch, nlmcatalog, nuccore, nucest, nucgss, nucleotide,
            - omim, orgtrack, pcassay, pccompound, pcsubstance, pmc, popset, probe, protein, proteinclusters,
            - pubmed, pubmedhealth, seqannot, snp, sparcle, sra, structure, taxonomy, unigene
        `joiner` : The delimit to use if the field is a list
        `record` : A function to transform the record.
    @requires:
        [python-eutils](https://github.com/biocommons/eutils)
    """),
    input='term',
    output='outfile:file:{% import re %}{{i.term | re.sub: r"[^\\w_]", "_", _ \
                                                 | re.sub: r"_+", "_", _ \
                                                 | [:255]}}.{{args.prog}}.txt',
    lang=opts.python,
    args=Diot(
        prog='esearch',
        apikey=opts.ncbikey,
        sleep=.15,
        # annotinfo, assembly, biocollections, bioproject, biosample,
        # biosystems, blastdbinfo, books,
        # cdd, clinvar, clone, dbvar, gap, gapplus, gds, gencoll, gene,
        # genome, geoprofiles, grasp, gtr,
        # homologene, ipg, medgen, mesh, ncbisearch, nlmcatalog, nuccore,
        # nucest, nucgss, nucleotide,
        # omim, orgtrack, pcassay, pccompound, pcsubstance, pmc, popset,
        # probe, protein, proteinclusters,
        # pubmed, pubmedhealth, seqannot, snp, sparcle, sra, structure,
        # taxonomy, unigene
        db='pubmed',
        joiner='|',
        record=None,
    )
)
pNCBI.errhow = 'retry'
