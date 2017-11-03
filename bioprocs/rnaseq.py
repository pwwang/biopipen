from os import path
from glob import glob
from pyppl import Proc, Box
from .utils import plot, txt
from . import params

"""
@name:
	pExpdir2Matrix
@description:
	Convert expression files to expression matrix
	File names will be used as sample names (colnames)
	Each gene and its expression per line.
	Suppose each expression file has the same rownames and in the same order.
@input:
	`expdir:file`:  the directory containing the expression files, could be gzipped
@output:
	`outfile:file`: the expression matrix file
	`outdir:dir`:   the directory containing expr file and plots
@args:
	`pattern` : The pattern to filter files. Default `'*'`
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
"""
pExpdir2Matrix                     = Proc(desc = 'Merge expression files to a matrix.')
pExpdir2Matrix.input               = "expdir:file"
pExpdir2Matrix.output              = [
	"outfile:file:{{in.expdir, args.pattern | fsDirname}}/{{in.expdir, args.pattern | fsDirname}}.expr.txt", 
	"outdir:dir:{{in.expdir, args.pattern | fsDirname}}"
]
pExpdir2Matrix.lang                = params.Rscript.value
pExpdir2Matrix.args.pattern        = '*'
pExpdir2Matrix.args.header         = False
pExpdir2Matrix.args.exrows         = ["^Sample", "^Composite", "^__"]
pExpdir2Matrix.args.boxplot        = False
pExpdir2Matrix.args.heatmap        = False
pExpdir2Matrix.args.heatmapn       = 500
pExpdir2Matrix.args.histplot       = False
pExpdir2Matrix.args.devpars        = Box({'res': 300, 'width': 2000, 'height': 2000})
pExpdir2Matrix.args.boxplotggs     = ['r:ylab("Expression")']
pExpdir2Matrix.args.heatmapggs     = ['r:theme(axis.text.y = element_blank())']
pExpdir2Matrix.args.histplotggs    = ['r:labs(x = "Expression", y = "# Samples")']
pExpdir2Matrix.tplenvs.plotBoxplot = plot.boxplot.r
pExpdir2Matrix.tplenvs.plotHeatmap = plot.heatmap.r
pExpdir2Matrix.tplenvs.plotHist    = plot.hist.r
pExpdir2Matrix.tplenvs.fsDirname   = lambda dir, pat: path.splitext(path.basename(glob(path.join(dir, pat))[0]))[0] + '_etc'
pExpdir2Matrix.script              = "file:scripts/rnaseq/pExpdir2Matrix.r"

"""
@name:
	pBatchEffect
@description:
	Remove batch effect with sva-combat.
@input:
	`expr:file`:  The expression file, generated by pExpdir2Matrix
	`batch:file`: The batch file defines samples and batches.
@output:
	`outfile:file`: the expression matrix file
	`outdir:dir`:   the directory containing expr file and plots
@args:
	`tool`    : The tool used to remove batch effect. Default `'combat'`
	`boxplot` : Whether to plot a boxplot. Default: False
	`heatmap` : Whether to plot a heatmap. Default: False
	`histplot`: Whether to plot a histgram. Default: False
	`devpars` : Parameters for png. Default: `{'res': 300, 'width': 2000, 'height': 2000}`
	`boxplotggs`: The ggplot parameters for boxplot. Default: `['r:ylab("Expression")']`
		- See ggplot2 documentation.
	`heatmapggs`: The ggplot parameters for heatmap. Default: `['r:theme(axis.text.y = element_blank())']`
	`histplotggs`: The ggplot parameters for histgram. Default: `['r:labs(x = "Expression", y = "# Samples")']`	
"""
pBatchEffect                       = Proc(desc = 'Try to remove batch effect of expression data.')
pBatchEffect.input                 = "expr:file, batch:file"
pBatchEffect.output                = "outfile:file:{{in.expr | fn | fn}}/{{in.expr | fn | fn}}.expr.txt, outdir:dir:{{in.expr | fn | fn}}"
pBatchEffect.args.tool             = 'combat'
pBatchEffect.args.boxplot          = False
pBatchEffect.args.heatmap          = False
pBatchEffect.args.heatmapn         = 500
pBatchEffect.args.histplot         = False
pBatchEffect.args.devpars          = Box({'res': 300, 'width': 2000, 'height': 2000})
pBatchEffect.args.boxplotggs       = ['r:ylab("Expression")']
pBatchEffect.args.heatmapggs       = ['r:theme(axis.text.y = element_blank())']
pBatchEffect.args.histplotggs      = ['r:labs(x = "Expression", y = "# Samples")']
pBatchEffect.tplenvs.plotBoxplot   = plot.boxplot.r
pBatchEffect.tplenvs.plotHeatmap   = plot.heatmap.r
pBatchEffect.tplenvs.plotHist      = plot.hist.r
pBatchEffect.tplenvs.txtSampleinfo = txt.sampleinfo.r
pBatchEffect.lang                  = params.Rscript.value
pBatchEffect.script                = "file:scripts/rnaseq/pBatchEffect.r"

"""
@name:
	pRawCounts2
@description:
	Convert raw counts to another unit
@input:
	`expfile:file`: the expression matrix
		- rows are genes, columns are samples
@output:
	`outfile:file`: the converted expression matrix
@args:
	`transpose`: transpose the input matrix? default: False
	`log2`:      whether to take log2? default: False
	`unit`:      convert to which unit? default: cpm (or rpkm, tmm)
	`header`:    whether input file has header? default: True
	`rownames`:  the index of the column as rownames. default: 1
	`glenfile`:  the gene length file, for RPKM
		- no head, row names are genes, have to be exact the same order and length as the rownames of expfile
	`boxplot` : Whether to plot a boxplot. Default: False
	`heatmap` : Whether to plot a heatmap. Default: False
	`histplot`: Whether to plot a histgram. Default: False
	`devpars` : Parameters for png. Default: `{'res': 300, 'width': 2000, 'height': 2000}`
	`boxplotggs`: The ggplot parameters for boxplot. Default: `['r:ylab("Expression")']`
		- See ggplot2 documentation.
	`heatmapggs`: The ggplot parameters for heatmap. Default: `['r:theme(axis.text.y = element_blank())']`
	`histplotggs`: The ggplot parameters for histgram. Default: `['r:labs(x = "Expression", y = "# Samples")']`	
@requires:
	[edgeR](https://bioconductor.org/packages/release/bioc/html/edger.html) if cpm or rpkm is chosen
	[coseq](https://rdrr.io/rforge/coseq/man/transform_RNAseq.html) if tmm is chosen
"""
pRawCounts2                     = Proc(desc = 'Convert raw counts to another unit.')
pRawCounts2.input               = "expfile:file"
pRawCounts2.output              = "outfile:file:{{in.expfile | fn | fn}}/{{in.expfile | fn | fn}}.expr.txt, outdir:dir:{{in.expfile | fn | fn}}"
pRawCounts2.args.unit           = 'cpm'
pRawCounts2.args.header         = True
pRawCounts2.args.log2           = False
pRawCounts2.args.glenfile       = ''
pRawCounts2.args.boxplot        = False
pRawCounts2.args.heatmap        = False
pRawCounts2.args.heatmapn       = 500
pRawCounts2.args.histplot       = False
pRawCounts2.args.devpars        = Box({'res': 300, 'width': 2000, 'height': 2000})
pRawCounts2.args.boxplotggs     = ['r:ylab("Expression")']
pRawCounts2.args.heatmapggs     = ['r:theme(axis.text.y = element_blank())']
pRawCounts2.args.histplotggs    = ['r:labs(x = "Expression", y = "# Samples")']
pRawCounts2.tplenvs.plotBoxplot = plot.boxplot.r
pRawCounts2.tplenvs.plotHeatmap = plot.heatmap.r
pRawCounts2.tplenvs.plotHist    = plot.hist.r
pRawCounts2.lang                = params.Rscript.value
pRawCounts2.script              = "file:scripts/rnaseq/pRawCounts2.r"

p2RawCounts                     = Proc(desc = 'Convert normalized expression back to raw counts.')
p2RawCounts.input               = 'expfile:file'
p2RawCounts.output              = 'outfile:file:{{in.expfile | fn | fn}}/{{in.expfile | fn | fn}}.counts.txt, outdir:dir:{{in.expfile | fn | fn}}'
p2RawCounts.args.unit           = 'fpkm'
p2RawCounts.args.header         = True
p2RawCounts.args.refgene        = params.refgene.value
p2RawCounts.args.boxplot        = False
p2RawCounts.args.heatmap        = False
p2RawCounts.args.heatmapn       = 500
p2RawCounts.args.histplot       = False
p2RawCounts.args.devpars        = Box({'res': 300, 'width': 2000, 'height': 2000})
p2RawCounts.args.boxplotggs     = ['r:ylab("Expression")']
p2RawCounts.args.heatmapggs     = ['r:theme(axis.text.y = element_blank())']
p2RawCounts.args.histplotggs    = ['r:labs(x = "Expression", y = "# Samples")']
p2RawCounts.tplenvs.plotBoxplot = plot.boxplot.r
p2RawCounts.tplenvs.plotHeatmap = plot.heatmap.r
p2RawCounts.tplenvs.plotHist    = plot.hist.r
p2RawCounts.lang                = params.Rscript.value
p2RawCounts.script              = "file:scripts/rnaseq/p2RawCounts.r"


"""
@name:
	pRnaseqDeg
@description:
	Detect DEGs for RNA-seq data
@input:
	`efile:file`: The expression matrix
	`gfile:file`: The group information
		- Like:
		```
		Sample1	Group1
		Sample2	Group1
		Sample3	Group1
		Sample4	group2
		Sample5	group2
		Sample6	group2
		```
@output:
	`outfile:file`: The DEG list
	`outdir:file`:  The output directory containing deg list and plots
@args:
	`tool`      : the tool used to detect DEGs. Default: 'edger' (deseq2)
	`filter`    : filter out low count records. Default: `"1,2"` (At least 2 samples have at least 2 reads)
	`mdsplot`   : whether to plot the MDS plot, default : True
	`volplot`   : whether to plot the volcano plot, default : True
	`maplot`    : whether to plot MA plots within each group, default : False
	`heatmap`   : whether to plot the heatmap using DEGs. Default : False
	`heatmapn`  : How many genes to be used for heatmap. If `heatmapn`, the number will be `heatmapn * # DEGs`. Default: 100
	`heatmapggs`: The ggplots options for heatmap. Default : []
	`maplotggs` : The ggplots options for maplot. Default : []
	`volplotggs`: The ggplots options for volplot. Default : []
	`devpars`   : Parameters for png. Default: `{'res': 300, 'width': 2000, 'height': 2000}`
"""
pRnaseqDeg        = Proc(desc = 'Detect DEGs by RNA-seq data.')
pRnaseqDeg.input  = "efile:file, gfile:file"
pRnaseqDeg.output = [
	"outfile:file:{{in.efile | fn | fn}}-{{in.gfile | fn | fn}}-DEGs/{{in.efile | fn | fn}}-{{in.gfile | fn | fn}}.degs.txt", 
	"outdir:dir:{{in.efile | fn | fn}}-{{in.gfile | fn | fn}}-DEGs"
]
pRnaseqDeg.args.tool             = 'edger' # deseq2
pRnaseqDeg.args.filter           = '1,2'
pRnaseqDeg.args.pval             = 0.05
pRnaseqDeg.args.mdsplot          = True
pRnaseqDeg.args.volplot          = True
pRnaseqDeg.args.maplot           = False
pRnaseqDeg.args.heatmap          = False
pRnaseqDeg.args.heatmapn         = 100
pRnaseqDeg.args.heatmapggs       = ['r:theme(axis.text.y = element_blank())']
pRnaseqDeg.args.maplotggs        = []
pRnaseqDeg.args.volplotggs       = []
pRnaseqDeg.args.devpars          = Box({'res': 300, 'width': 2000, 'height': 2000})
pRnaseqDeg.tplenvs.plotHeatmap   = plot.heatmap.r
pRnaseqDeg.tplenvs.plotMAplot    = plot.maplot.r
pRnaseqDeg.tplenvs.plotVolplot   = plot.volplot.r
pRnaseqDeg.tplenvs.txtSampleinfo = txt.sampleinfo.r
pRnaseqDeg.lang                  = params.Rscript.value
pRnaseqDeg.script                = "file:scripts/rnaseq/pRnaseqDeg.r"

pCoexp             = Proc(desc = "Get co-expression of gene pairs in the expression matrix.")
pCoexp.input       = "infile:file"
pCoexp.output      = "outfile:file:{{in.infile | fn}}.coexp, outpval:file:{{in.infile | fn}}.pval"
pCoexp.args.method = 'pearson'
pCoexp.args.pval   = False
pCoexp.lang        = params.Rscript.value
pCoexp.script      = "file:scripts/rnaseq/pCoexp.r"
