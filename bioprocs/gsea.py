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
pGMT2Mat.output = "outfile:file:{{in.infile | fn}}.gmat"
pGMT2Mat.lang   = params.python.value
pGMT2Mat.script = "file:scripts/gsea/pGMT2Mat.py"

"""
@name:
	pExpmat2Gct
@description:
	Convert expression matrix to GCT file.
	Refer to http://software.broadinstitute.org/cancer/software/genepattern/file-formats-guide#GCT for file format
@input:
	`expfile:file`: the input expression matrix file. Samples as columns, genes as rows.
@output:
	`outfile:file`: the gct file
"""
pExpmat2Gct        = Proc(desc = 'Convert expression matrix to GCT file.')
pExpmat2Gct.input  = 'expfile:file'
pExpmat2Gct.output = 'outfile:file:{{ in.expfile | fn }}.gct'
pExpmat2Gct.lang   = params.python.value
pExpmat2Gct.script = "file:scripts/gsea/pExpmat2Gct.py"

"""
@name:
	pSampleinfo2Cls
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
pSampleinfo2Cls                    = Proc(desc = 'Convert sample infomation to cls file.')
pSampleinfo2Cls.input              = 'sifile:file'
pSampleinfo2Cls.output             = 'outfile:file:{{ in.sifile | fn }}.cls'
#pSampleinfo2Cls.envs.txtSampleinfo = txt.sampleinfo.py
pSampleinfo2Cls.lang               = params.python.value
pSampleinfo2Cls.script             = "file:scripts/gsea/pSampleinfo2Cls.py"

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
pSSGSEA.output         = "outdir:file:{{in.gctfile | fn}}-{{in.gmtfile | fn}}-ssGSEA"
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
	`nperm`:     Number of permutations. Default: 10000
"""
pGSEA                = Proc (desc = 'Do GSEA.')
pGSEA.input          = "gctfile:file, clsfile:file, gmtfile:file"
pGSEA.output         = "outdir:dir:{{in.gctfile | fn}}-{{in.gmtfile | fn}}-GSEA"
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
	`topn`: Top N pathways used to plot. Default: 10
	`col`: The columns index containing the genes. Default: 0
	`delimit`: The delimit of input file. Default: '\\t'
	`dbs`:  The databases to do enrichment against. Default: KEGG_2016
	  - A full list can be found here: http://amp.pharm.mssm.edu/Enrichr/#stats
	  - Multiple dbs separated by comma (,)
	`norm`: Normalize the gene list use [python-mygene](https://pypi.python.org/pypi/mygene/3.0.0)
	`rmtags`: Remove pathway tags in the plot. Default: True
	  - For example: change "Lysine degradation_Homo sapiens_hsa00310" to "Lysine degradation".
	`plot`: Whether to plot the result. Default: True
	`title`: The title for the plot. Default: "Gene enrichment: {db}"
@requires:
	[python-mygene](https://pypi.python.org/pypi/mygene/3.0.0) if `args.norm` is `True`
"""
pEnrichr                = Proc()
pEnrichr.input          = "infile:file"
pEnrichr.output         = "outdir:dir:{{in.infile | fn}}.enrichr"
pEnrichr.lang           = params.python.value
pEnrichr.args.inopts    = Box(delimit = '\t', skip = 0, comment = '#')
pEnrichr.args.top       = 10
pEnrichr.args.genecol   = ''
pEnrichr.args.libs      = "KEGG_2016"
pEnrichr.args.norm      = False
pEnrichr.args.rmtags    = True
pEnrichr.args.plot      = True
pEnrichr.args.title     = "Gene enrichment: {library}"
pEnrichr.args.cachedir  = params.cachedir.value
#pEnrichr.envs.genenorm = genenorm.py
pEnrichr.errhow         = 'retry'
pEnrichr.script         = "file:scripts/gsea/pEnrichr.py"


"""
@name:
	pTargetEnrichr
@description:
	Use APIs from http://amp.pharm.mssm.edu/Enrichr/help#api&q=1 to analyze a gene list
@input:
	`infile:file`: The target genes with regulators
		- Format:
		- Header is not required, but may specified in first line starting with `#`
		- If only 3 columns are there, the 3rd column is anyway the relation!
		- If only 4 columns are there, 3rd is target status, 4th is relation!
		  ```
		  #Regulator	Target	Regulator status	Target status	Relation
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
pTargetEnrichr.output        = "outdir:dir:{{in.infile | fn}}.tenrichr"
pTargetEnrichr.lang          = params.python.value
pTargetEnrichr.args.inopts   = Box(delimit = '\t', skip = 0, comment = '#')
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
