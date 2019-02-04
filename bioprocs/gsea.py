from os import path
from pyppl import Proc, Box
from . import params

"""
@name:
	pGMT2Mat
@description:
	Convert a GMT file to a matrix.
	Rownames of GMT file will be the column names of output matrix.
@input:
	`infile:file`: The input file in GMT format.
@output:
	`outfile:file`: output matrix file
"""
pGMT2Mat        = Proc(desc = 'Convert a GMT file to a matrix.')
pGMT2Mat.input  = "infile:file"
pGMT2Mat.output = "outfile:file:{{i.infile | fn}}.gmat"
pGMT2Mat.lang   = params.python.value
pGMT2Mat.script = "file:scripts/gsea/pGMT2Mat.py"

"""
@name:
	pExprMat2GCT
@description:
	Convert expression matrix to GCT file.
	Refer to http://software.broadinstitute.org/cancer/software/genepattern/file-formats-guide#GCT for file format
@input:
	`expfile:file`: the input expression matrix file. Samples as columns, genes as rows.
@output:
	`outfile:file`: the gct file
"""
pExprMat2GCT        = Proc(desc = 'Convert expression matrix to GCT file.')
pExprMat2GCT.input  = 'expfile:file'
pExprMat2GCT.output = 'outfile:file:{{ i.expfile | fn }}.gct'
pExprMat2GCT.lang   = params.python.value
pExprMat2GCT.script = "file:scripts/gsea/pExprMat2GCT.py"

"""
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
"""
pSampleinfo2CLS                    = Proc(desc = 'Convert sample infomation to cls file.')
pSampleinfo2CLS.input              = 'sifile:file'
pSampleinfo2CLS.output             = 'outfile:file:{{ i.sifile | fn }}.cls'
#pSampleinfo2CLS.envs.txtSampleinfo = txt.sampleinfo.py
pSampleinfo2CLS.lang               = params.python.value
pSampleinfo2CLS.script             = "file:scripts/gsea/pSampleinfo2CLS.py"

"""
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
"""
pSSGSEA = Proc (desc = 'Do single-sample GSEA.')
pSSGSEA.input          = "gctfile:file, gmtfile:file"
pSSGSEA.output         = "outdir:file:{{i.gctfile | fn}}-{{i.gmtfile | fn}}-ssGSEA"
pSSGSEA.args.weightexp = 1
pSSGSEA.args.nperm     = 1000
pSSGSEA.args.seed      = -1
pSSGSEA.lang           = params.Rscript.value
pSSGSEA.script         = "file:scripts/gsea/pSSGSEA.r"

"""
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
"""
pGSEA                = Proc (desc = 'Do GSEA.')
pGSEA.input          = "gctfile:file, clsfile:file, gmtfile:file"
pGSEA.output         = "outdir:dir:{{i.gctfile | fn}}-{{i.gmtfile | fn}}.GSEA"
pGSEA.args.weightexp = 1
pGSEA.args.nperm     = 1000
pGSEA.args.nthread   = 1
pGSEA.args.seed      = -1
pGSEA.lang           = params.Rscript.value
pGSEA.script         = "file:scripts/gsea/pGSEA.r"

"""
@name:
	pEnrichr
@description:
	Use APIs from http://amp.pharm.mssm.edu/Enrichr/help#api&q=1 to analyze a gene list
@input:
	`infile:file`: The gene list, each per line
@output:
	`outdir:dir`:  The output directory, containing the tables and figures.
@args:
	`top`:     Top N pathways used to plot. Default: 10
	`genecol`: The columns index containing the genes. Default: 0
	`inopts`:  The input options.
		- `delimit`: The delimit of input file. Default: `\t`
		- `skip`:    Skip first N lines. Default: `0`
		- `comment`: Line comment mark. Default: `#`
		- Other parameters fit `bioprocs.utils.tsvio.TsvReader`
	`libs`:  The databases to do enrichment against. Default: KEGG_2016
	  - A full list can be found here: http://amp.pharm.mssm.edu/Enrichr/#stats
	  - Multiple dbs separated by comma (,)
	`plot`: Whether to plot the result. Default: True
	`devpars`: Parameters for png. Default: `{'res': 300, 'width': 2000, 'height': 2000}`
"""
pEnrichr               = Proc()
pEnrichr.input         = "infile:file"
pEnrichr.output        = "outdir:dir:{{i.infile | fn}}.enrichr"
pEnrichr.lang          = params.python.value
pEnrichr.args.inopts   = Box(delimit = '\t', skip = 0, comment = '#')
pEnrichr.args.top      = 20
pEnrichr.args.cutoff   = 1
pEnrichr.args.genecol  = ''
pEnrichr.args.nthread  = 1
pEnrichr.args.Rscript  = params.Rscript.value
pEnrichr.args.pathview = Box() # Box(fccol = 2)
pEnrichr.args.libs     = "KEGG_2016"
pEnrichr.args.devpars  = Box(res = 300, width = 2000, height = 2000)
pEnrichr.args.plot     = True
pEnrichr.errhow        = 'retry'
pEnrichr.script        = "file:scripts/gsea/pEnrichr.py"


"""
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
"""
pTargetEnrichr               = Proc(desc = 'Do gene set enrichment analysis for target genes.')
pTargetEnrichr.input         = "infile:file"
pTargetEnrichr.output        = "outdir:dir:{{i.infile | fn}}.tenrichr"
pTargetEnrichr.lang          = params.python.value
pTargetEnrichr.args.inopts   = Box(delimit = '\t', skip = 0, comment = '#', ftype = 'head')
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
#pTargetEnrichr.envs.genenorm = genenorm.py
pTargetEnrichr.errhow        = 'retry'
pTargetEnrichr.script        = "file:scripts/gsea/pTargetEnrichr.py"
