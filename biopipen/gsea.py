"""Gene set enrichment analysis"""

from diot import Diot
from . import opts, proc_factory

# pylint: disable=invalid-name

pGMT2Mat = proc_factory(
    desc='Convert a GMT file to a matrix.',
    config=Diot(annotate="""
    @name:
        pGMT2Mat
    @description:
        Convert a GMT file to a matrix.
        Rownames of GMT file will be the column names of output matrix.
    @input:
        `infile:file`: The input file in GMT format.
    @output:
        `outfile:file`: output matrix file
    """),
    input="infile:file",
    output="outfile:file:{{i.infile | fn}}.gmat",
    lang=opts.python,
)

pExprMat2GCT = proc_factory(
    desc='Convert expression matrix to GCT file.',
    config=Diot(annotate="""
    @name:
        pExprMat2GCT
    @description:
        Convert expression matrix to GCT file.
        Refer to http://software.broadinstitute.org/cancer/software/genepattern/file-formats-guide#GCT for file format
    @input:
        `expfile:file`: the input expression matrix file. Samples as columns, genes as rows.
    @output:
        `outfile:file`: the gct file
    """),
    input='expfile:file',
    output='outfile:file:{{ i.expfile | fn }}.gct',
    lang=opts.python,
)

pSampleinfo2CLS = proc_factory(
    desc='Convert sample infomation to cls file.',
    config=Diot(annotate="""
    @name:
        pSampleinfo2CLS
    @description:
        Convert sample infomation to cls file.
        Refer to http://software.broadinstitute.org/cancer/software/genepattern/file-formats-guide#CLS for file format
        NOTE that the order of samples must be the same as in GMT file in further analysis.
    @input:
        `sifile:file`: the sample information file.
            - Headers are: [Sample, ]Patient, Group, Batch
            - Rows are samples
    @output:
        `outfile:file`: the cls file
    """),
    input='sifile:file',
    output='outfile:file:{{ i.sifile | fn }}.cls',
    #pSampleinfo2CLS.envs.txtSampleinfo = txt.sampleinfo.py
    lang=opts.python,
)

pSSGSEA = proc_factory(
    desc='Do single-sample GSEA.',
    config=Diot(annotate="""
    @name:
        pSSGSEA
    @description:
        Single sample GSEA
        Refer to http://software.broadinstitute.org/cancer/software/genepattern/file-formats-guide#GCT for GCT file format
        Refer to http://software.broadinstitute.org/cancer/software/genepattern/file-formats-guide#GMT for GMT file format
    @input:
        `gctfile:file`: the expression file
        `gmtfile:file`: the gmtfile for gene sets
    @output:
        `outdir:file`: the output directory
            - `report.txt`: the enrichment report for each Gene set.
            - `RES_<GeneSet>.png`: the running ES plot for <GeneSet>
            - `normP_<GeneSet>.png`: the norminal P value plot for <GeneSet>
    @args:
        `weightexp`: Exponential weight employed in calculation of enrichment scores. Default: 0.75
        `nperm`:     Number of permutations. Default: 10000
    """),
    input="gctfile:file, gmtfile:file",
    output="outdir:file:{{i.gctfile | fn}}-{{i.gmtfile | fn}}-ssGSEA",
    args=Diot(
        weightexp=1,
        nperm=1000,
        seed=-1,
    )
)

pGSEA = proc_factory(
    desc='Multiple-sample GSEA',
    config=Diot(annotate="""
    @description:
        Multiple-sample GSEA with 3+ samples using fgsea
    @input:
        exprfile: The expression matrix.
            - Must be normalized and log2 transformed in most cases
            - See `args.ranking` to select the right ranking method for your data
        saminfo: The sample information file
        gmtfile: the gmtfile for gene sets
            - Refer to http://software.broadinstitute.org/cancer/software/genepattern/file-formats-guide#GMT for GMT file format
    @output:
        outfile: The output table for top pathways
        outdir: The output directory containing the other outputs
    @args:
        inopts (Diot): Options to read the expression file
            - `cnames` and `rnames` will be forced to be `True`
        ranking (str): The method to rank the genes. See: https://www.gsea-msigdb.org/gsea/doc/GSEAUserGuideTEXT.htm#_Metrics_for_Ranking Could be one of
            - `s2n`: Signal2Noise
            - `ttest`: tTest
            - `roc`: Ratio_of_Classes
            - `doc`: Diff_of_Classes
            - `log2roc`: log2_Ratio_of_Classes
        minsize (int): Minimal size of a gene set to test. All pathways below the threshold are excluded.
        top (int): Top N pathways to plot
        nthread (int): Number of threads to use
        seed (int): Seed for random number generator
        devpars (Diot): Device parameters for plotting
    """),
    input="exprfile:file, saminfo:file, gmtfile:file",
    output=("outfile:file:{{i.exprfile | stem2}}-{{i.gmtfile | stem2}}.GSEA/"
            "{{i.exprfile | stem2}}-{{i.gmtfile | stem2}}.gsea.txt, "
            "outdir:dir:{{i.exprfile | stem2}}-{{i.gmtfile | stem2}}.GSEA"),
    args=Diot(
        ranking="s2n",
        nthread=1,
        seed=8525,
        top=20,
        minsize=10,
        inopts=Diot(dup='ignore'),
        devpars=Diot(res=300, height=2000, width=2000)
    ),
    lang=opts.Rscript,
)

pEnrichr = proc_factory(
    desc='Gene set enrichment analysis using Enrichr APIs',
    config=Diot(annotate="""
    @input:
        `infile:file`: The gene list, each per line
    @output:
        `outdir:dir`:  The output directory, containing the tables and figures.
    @args:
        `top`:     Top N pathways used to plot
        `genecol`: The columns index containing the genes
        `inopts`:  The input options.
            - `delimit`: The delimit of input file.
            - `skip`:    Skip first N lines.
            - `comment`: Line comment mark.
            - Other parameters fit `bioprocs.utils.tsvio.TsvReader`
        `libs`:  The databases to do enrichment against.
            - A full list can be found here: http://amp.pharm.mssm.edu/Enrichr/#stats
            - Multiple dbs separated by comma (,)
        `plot`: Whether to plot the result.
        `devpars`: Parameters for png.
        include: A lambda function to include the records(genes)
            - argument is `bioprocs.utils.tsvio2.TsvRecord`
    """),
    input='infile:file',
    output='outdir:dir:{{i.infile | stem}}.enrichr',
    lang=opts.python,
    errhow='retry',
    args=Diot(
        inopts=Diot(delimit='\t', skip=0, comment='#'),
        top=20,
        cutoff=1.1,
        genecol=0,
        include=None,
        nthread=1,
        Rscript=opts.Rscript,
        pathview=Diot(),  # Diot(fccol = 2)
        libs="KEGG_2019_Human",
        devpars=Diot(res=300, width=2000, height=2000),
        plot=True
    )
)

pGene2Pathway = proc_factory(
    desc='Find pathways that genes are present.',
    config=Diot(annotate="""
    @name:
        pGene2Pathway
    @description:
        Find pathways where genes are present
    @input:
        `infile:file`: The input file
    @output:
        `outfile:file`: The output file, Default: `{{i.infile | fn}}-pw{{i.infile | ext}}`
    @args:
        `inopts`: Reading options for input file, Default: `Diot(cnames = True)`
        `genecol`: Index or name of the gene column, Default: `0`
        `libs`: Libraries of the pathways, Default: `KEGG_2019_Human`
    """),
    input='infile:file',
    output='outfile:file:{{i.infile | fn}}-pw{{i.infile | ext}}',
    lang=opts.python,
    args=Diot(
        inopts=Diot(cnames=True),
        genecol=0,
        libs="KEGG_2019_Human",
    )
)

pTargetEnrichr = proc_factory(
    desc='Do gene set enrichment analysis for target genes.',
    config=Diot(annotate="""
    @name:
        pTargetEnrichr
    @description:
        Use APIs from http://amp.pharm.mssm.edu/Enrichr/help#api&q=1 to analyze a gene list
    @input:
        `infile:file`: The target genes with regulators
            - Format (RegulatorStatus and TargetStatus are optional):
            ```
            Regulator	Target	RegulatorStatus	TargetStatus	Relation
            has-mir-22	Gene	+	+	+
            ```
    @output:
        `outdir:dir`:  The output directory, containing the tables and figures.
    @args:
        `dbs`       : The databases to do enrichment against. Default: KEGG_2016
            - A full list can be found here: http://amp.pharm.mssm.edu/Enrichr/#stats
            - Multiple dbs separated by comma (,)
        `rmtags`    : Remove pathway tags in the plot. Default: True
            - For example: change "Lysine degradation_Homo sapiens_hsa00310" to "Lysine degradation".
        `enrplot`   : Whether to plot the result. Default: True
        `enrn`      : Top N pathways used to plot. Default: 10
        `netplot`   : Whether to plot the network. Default: True
        `netn`      : Top N pathways used to plot the network. Default: 5
            - Must <= `enrn`. If `netn` >= `enrn`, `netn` = `enrn`
        `title`     : The title for the plot. Default: "Gene enrichment: {db}"
    @requires:
        [`python-mygene`](https://pypi.python.org/pypi/mygene/3.0.0)
        [`graphviz`](https://pypi.python.org/pypi/graphviz)
    """),
    input="infile:file",
    output="outdir:dir:{{i.infile | fn}}.tenrichr",
    lang=opts.python,
    args=Diot(
        inopts=Diot(delimit='\t',
                    skip=0,
                    comment='#',
                    ftype='head'),
        genecol="COL2",
        dbs="KEGG_2016",
        norm=False,
        rmtags=True,
        enrplot=True,
        enrn=10,
        netplot=True,
        netn=5,
        title="Gene enrichment: {db}",
        cachedir=opts.cachedir,
    ),
    errhow='retry',
)
