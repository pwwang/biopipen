"""Gene set enrichment analysis"""
from os import path
from pyppl import Proc
from diot import Diot
from . import params, proc_factory

pGMT2Mat = proc_factory(
	desc = 'Convert a GMT file to a matrix.',
	config = Diot(annotate = """
	@name:
		pGMT2Mat
	@description:
		Convert a GMT file to a matrix.
		Rownames of GMT file will be the column names of output matrix.
	@input:
		`infile:file`: The input file in GMT format.
	@output:
		`outfile:file`: output matrix file
	"""))
pGMT2Mat.input  = "infile:file"
pGMT2Mat.output = "outfile:file:{{i.infile | fn}}.gmat"
pGMT2Mat.lang   = params.python.value

pExprMat2GCT = proc_factory(
	desc = 'Convert expression matrix to GCT file.',
	config = Diot(annotate = """
	@name:
		pExprMat2GCT
	@description:
		Convert expression matrix to GCT file.
		Refer to http://software.broadinstitute.org/cancer/software/genepattern/file-formats-guide#GCT for file format
	@input:
		`expfile:file`: the input expression matrix file. Samples as columns, genes as rows.
	@output:
		`outfile:file`: the gct file
	"""))
pExprMat2GCT.input  = 'expfile:file'
pExprMat2GCT.output = 'outfile:file:{{ i.expfile | fn }}.gct'
pExprMat2GCT.lang   = params.python.value

pSampleinfo2CLS = proc_factory(
	desc = 'Convert sample infomation to cls file.',
	config = Diot(annotate = """
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
	"""))
pSampleinfo2CLS.input              = 'sifile:file'
pSampleinfo2CLS.output             = 'outfile:file:{{ i.sifile | fn }}.cls'
#pSampleinfo2CLS.envs.txtSampleinfo = txt.sampleinfo.py
pSampleinfo2CLS.lang               = params.python.value

pSSGSEA = proc_factory(
	desc = 'Do single-sample GSEA.',
	config = Diot(annotate = """
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
	"""))
pSSGSEA.input          = "gctfile:file, gmtfile:file"
pSSGSEA.output         = "outdir:file:{{i.gctfile | fn}}-{{i.gmtfile | fn}}-ssGSEA"
pSSGSEA.args.weightexp = 1
pSSGSEA.args.nperm     = 1000
pSSGSEA.args.seed      = -1

pGSEA = proc_factory(
	desc = 'Do GSEA.',
	config = Diot(annotate = """
	@name:
		pGSEA
	@description:
		GSEA
		Refer to http://software.broadinstitute.org/cancer/software/genepattern/file-formats-guide#GCT for GCT file format
		Refer to http://software.broadinstitute.org/cancer/software/genepattern/file-formats-guide#GMT for GMT file format
		Refer to http://software.broadinstitute.org/cancer/software/genepattern/file-formats-guide#CLS for CLS file format
	@input:
		`gctfile:file`: the expression file
		`clsfile:file`: the class file
		`gmtfile:file`: the gmtfile for gene sets
	@output:
		`outdir:file`: the output directory
	@args:
		`weightexp`: Exponential weight employed in calculation of enrichment scores. Default: 0.75
		`nperm`:     Number of permutations. Default: 1000
	"""))
pGSEA.input          = "gctfile:file, clsfile:file, gmtfile:file"
pGSEA.output         = "outdir:dir:{{i.gctfile | fn}}-{{i.gmtfile | fn}}.GSEA"
pGSEA.args.weightexp = 1
pGSEA.args.nperm     = 1000
pGSEA.args.nthread   = 1
pGSEA.args.seed      = -1
pGSEA.lang           = params.Rscript.value

pEnrichr = proc_factory(
	desc   = 'Gene set enrichment analysis using Enrichr APIs',
	config = Diot(annotate = """
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
	input  = 'infile:file',
	output = 'outdir:dir:{{i.infile | stem}}.enrichr',
	lang   = params.python.value,
	errhow = 'retry',
	args   = Diot(
		inopts   = Diot(delimit = '\t', skip = 0, comment = '#'),
		top      = 20,
		cutoff   = 1,
		genecol  = 0,
		include  = None,
		nthread  = 1,
		Rscript  = params.Rscript.value,
		pathview = Diot(), # Diot(fccol = 2)
		libs     = "KEGG_2019_Human",
		devpars  = Diot(res = 300, width = 2000, height = 2000),
		plot     = True
	)
)

pGene2Pathway = proc_factory(
	desc = 'Find pathways that genes are present.',
	config = Diot(annotate = """
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
	"""))
pGene2Pathway.input        = 'infile:file'
pGene2Pathway.output       = 'outfile:file:{{i.infile | fn}}-pw{{i.infile | ext}}'
pGene2Pathway.args.inopts  = Diot(cnames = True)
pGene2Pathway.args.genecol = 0
pGene2Pathway.args.libs    = "KEGG_2019_Human"
pGene2Pathway.lang         = params.python.value

pTargetEnrichr = proc_factory(
	desc = 'Do gene set enrichment analysis for target genes.',
	config = Diot(annotate = """
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
	"""))
pTargetEnrichr.input         = "infile:file"
pTargetEnrichr.output        = "outdir:dir:{{i.infile | fn}}.tenrichr"
pTargetEnrichr.lang          = params.python.value
pTargetEnrichr.args.inopts   = Diot(delimit = '\t', skip = 0, comment = '#', ftype = 'head')
pTargetEnrichr.args.genecol  = "COL2"
pTargetEnrichr.args.dbs      = "KEGG_2016"
pTargetEnrichr.args.norm     = False
pTargetEnrichr.args.rmtags   = True
pTargetEnrichr.args.enrplot  = True
pTargetEnrichr.args.enrn     = 10
pTargetEnrichr.args.netplot  = True
pTargetEnrichr.args.netn     = 5
pTargetEnrichr.args.title    = "Gene enrichment: {db}"
pTargetEnrichr.args.cachedir = params.cachedir.value
pTargetEnrichr.errhow        = 'retry'
