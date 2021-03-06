"""Analysis of expression data from RNA-seq"""
from os import path
from glob import glob
from pyppl import Proc, Diot
from .utils import dirpat2name, fs2name
from . import params, proc_factory

pExprDir2Matrix = proc_factory(
	desc = 'Merge expression files to a matrix.',
	config = Diot(annotate = """
	@name:
		pExprDir2Matrix
	@description:
		Convert expression files to expression matrix
		File names will be used as sample names (colnames)
		Each gene and its expression per line.
		Suppose each expression file has the same rownames and in the same order.
	@input:
		`indir:file`:  the directory containing the expression files, could be gzipped
	@output:
		`outfile:file`: the expression matrix file
		`outdir:dir`:   the directory containing expr file and plots
	@args:
		`pattern` : The pattern to filter files. Default `'*'`
		`namefunc`: Transform filename (no extension) as column name. Default: "function(fn) fn"
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
	"""))
pExprDir2Matrix.input  = "indir:file"
pExprDir2Matrix.output = [
	"outfile:file:{{i.indir, args.pattern | dirpat2name}}.dir/{{i.indir, args.pattern | dirpat2name}}.expr.txt",
	"outdir:dir:{{i.indir, args.pattern | dirpat2name}}.dir"
]
pExprDir2Matrix.lang           = params.Rscript.value
pExprDir2Matrix.args.pattern   = '*'
pExprDir2Matrix.args.fn2sample = 'function(fn) fn'
pExprDir2Matrix.args.exrows    = ["^Sample", "^Composite", "^__"]
pExprDir2Matrix.args.hmrows    = 500
pExprDir2Matrix.args.plot = Diot(
	boxplot = False,
	heatmap = False,
	histogram = False
)
pExprDir2Matrix.args.ggs = Diot(
	boxplot = Diot(ylab = {0: "Expression"}),
	heatmap = Diot(theme = {'axis.text.y': 'r:element_blank()'}),
	histogram = Diot(labs = {'x': 'Expression', 'y': '# Genes'})
)
pExprDir2Matrix.envs.dirpat2name = dirpat2name
# pExprDir2Matrix.envs.rimport = rimport
pExprDir2Matrix.args.devpars = Diot(res = 300, width = 2000, height = 2000)

pExprFiles2Mat = proc_factory(
	desc = 'Merge expression to a matrix from single samples.',
	config = Diot(annotate = """
	@name:
		pExprFiles2Mat
	@description:
		Merge expression to a matrix from single samples.
	@input:
		`infiles:files`: The expression file from single samples, typically with 2 columns: gene and expression value
	@output:
		`outfile:file`: the expression matrix file
	@args:
		`fn2sample`: Transform filename (no extension) as column name. Default: "function(fn) fn"
		`inopts`   : Options to read input files. Default: `Diot(rname = True, cnames = True)`
	"""))
pExprFiles2Mat.input          = 'infiles:files'
pExprFiles2Mat.output         = 'outfile:file:{{i.infiles | fs2name}}.expr.txt'
pExprFiles2Mat.args.inopts    = Diot(cnames = True, rnames = True)
pExprFiles2Mat.args.fn2sample = 'function(fn) unlist(strsplit(fn, ".", fixed=T))[1]'
pExprFiles2Mat.envs.fs2name   = fs2name
# pExprFiles2Mat.envs.rimport   = rimport
pExprFiles2Mat.lang           = params.Rscript.value

pExprStats = proc_factory(
	desc = 'Plot the expression values out.',
	config = Diot(annotate = """
	@name:
		pExprStats
	@description:
		Plot the expression out.
	@input:
		`infile:file`: The expression matrix (rows are genes and columns are samples).
		`gfile:file` : The sample information file. Determines whether to do subgroup stats or not.
			- If not provided, not do for all samples
	@output:
		`outdir:dir`: The directory containing the plots
			- If `args.filter` is given, a filtered expression matrix will be generated in `outdir`.
	@args:
		`inopts`: Options to read `infile`. Default: `Diot(cnames = True, rnames = True)`
		`tsform`: An R function in string to transform the expression matrix (i.e take log).
		`filter`: An R function in string to filter the expression data.
		`plot`  : Which plot to do? Default: `Diot(boxplot = True, histogram = True, qqplot = True)`
		`ggs`   : The ggs for each plot. Default:
			- `boxplot   = Diot(ylab = {0: "Expression"})`,
			- `histogram = Diot(labs = {'x': 'Expression', 'y': '# Genes'})`,
			- `qqplot    = Diot()`
		`params` : The params for each ggplot function.
		`devpars`: Parameters for png. Default: `{'res': 300, 'width': 2000, 'height': 2000}`
	"""))
pExprStats.input       = 'infile:file, gfile:file'
pExprStats.output      = 'outdir:dir:{{i.infile | fn2}}.plots'
pExprStats.args.inopts = Diot(cnames = True, rnames = True)
pExprStats.args.tsform = None
pExprStats.args.filter = None
pExprStats.args.plot   = Diot(boxplot = True, histogram = True, qqplot = True)
pExprStats.args.ggs    = Diot(
	boxplot   = Diot(ylab = {0: "Expression"}),
	histogram = Diot(labs = {'x': 'Expression', 'y': '# Genes'}),
	qqplot    = Diot()
)
pExprStats.args.params   = Diot(
	boxplot   = Diot(),
	histogram = Diot(),
	qqplot    = Diot()
)
pExprStats.args.devpars = Diot(res = 300, width = 2000, height = 2000)
# pExprStats.envs.rimport = rimport
pExprStats.lang         = params.Rscript.value

pBatchEffect = proc_factory(
	desc = 'Try to remove batch effect of expression data.',
	config = Diot(annotate = """
	@name:
		pBatchEffect
	@description:
		Remove batch effect with sva-combat.
	@input:
		`expr:file`:  The expression file, generated by pExprDir2Matrix
		`batch:file`: The batch file defines samples and batches.
	@output:
		`outfile:file`: the expression matrix file
		`outdir:dir`:   the directory containing expr file and plots
	@args:
		`tool`    : The tool used to remove batch effect. Default `'combat'`
		`hmrows`  : How many rows to be used to plot heatmap
		`plot`: Whether to plot
			- `boxplot`   : Whether to plot a boxplot. Default: False
			- `heatmap`   : Whether to plot a heatmap. Default: False
			- `histogram` : Whether to plot a histgram. Default: False
		`devpars`    : Parameters for png. Default: `{'res': 300, 'width': 2000, 'height': 2000}`
		`ggs`: The ggplot parameters
			- `boxplot`  : The ggplot parameters for boxplot. Default: `Diot(ylab = {0: "Log2 Intensity"})`
			- `heatmap`  : The ggplot parameters for heatmap. Default: `Diot(theme = {'axis.text.y': 'r:element_blank()'})`
			- `histogram`: The ggplot parameters for histgram. Default: `Diot(labs = {'x': "Log2 Intensity", "y": "Density"})`
	"""))
pBatchEffect.input        = "expr:file, batch:file"
pBatchEffect.output       = "outfile:file:{{i.expr | fn2}}/{{i.expr | fn2}}.expr.txt, outdir:dir:{{i.expr | fn2}}"
pBatchEffect.args.tool    = 'combat'
pBatchEffect.args.plot    = Diot(boxplot = False, heatmap = False, histogram = False)
pBatchEffect.args.hmrows  = 500
pBatchEffect.args.devpars = Diot(res = 300, width = 2000, height = 2000)
pBatchEffect.args.ggs     = Diot(
	boxplot   = Diot(ylab = {0: "Expression"}),
	heatmap   = Diot(theme = {'axis.text.y': 'r:element_blank()'}),
	histogram = Diot(labs  = {'x': "Expression", "y": "Frequency"})
)
# pBatchEffect.envs.rimport = rimport
pBatchEffect.lang   = params.Rscript.value

pUnitConversion = proc_factory(
	desc   = 'Convert RNAseq data in different units back and forth',
	config = Diot(annotate = """
	@description:
		Convert expression between units
		See here: https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/ and
		https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/#fpkm
		Available converstions:
		- `count -> cpm, fpkm/rpkm, fpkm-uq/rpkm-rq, tpm, tmm`
		- `fpkm/rpkm -> count, tpm, cpm`
		- `tpm -> count, fpkm/rpkm, cpm`
		- `cpm -> count, fpkm/rpkm, tpm`
		NOTE that during some conversions, `sum(counts/effLen)` is approximated to `sum(counts)/sum(effLen) * length(effLen))`
	@input:
		`infile`: the expression matrix
			- rows are genes, columns are samples
	@output:
		`outfile`: the converted expression matrix
	@args:
		`inunit` : The unit of input expression values.
		`outunit`: The unit of output expression values.
		`meanfl` : A file containing the mean fragment length for each sample by rows, without header.
			- Or a fixed universal estimated number (1 used by TCGA).
		`nreads`:    The estimatied total number of reads for each sample.
			- Or you can pass a file with the number for each sample by rows, without header.
			- In converting `fpkm/rpkm -> count`: it should be total reads of that sample
			- In converting `cpm -> count`: it should be total reads of that sample
			- In converting `tpm -> count`: it should be total reads of that sample
			- In converting `tpm -> cpm`: it should be total reads of that sample
			- In converting `tpm -> fpkm/rpkm`: it should be `sum(fpkm)` of that sample
		`inform`:    Transform the input to the unit specified. (sometimes the values are log transformed)
			- For example, if the `inunit` is `tpm`, but it's actually `log2(expr+1)` transformed,
			- Then this should be `function(expr) 2^(expr - 1)`.
		`outform`: Transform to be done on the output expression values.
		`refexon`: the exome gff file, for RPKM/FPKM
			- `gene_id` is required for gene names
	@requires:
		[edgeR](https://bioconductor.org/packages/release/bioc/html/edger.html) if cpm or rpkm is chosen
		[coseq](https://rdrr.io/rforge/coseq/man/transform_RNAseq.html) if tmm is chosen
	"""),
	lang   = params.Rscript.value,
	input  = 'infile:file',
	output = "outfile:file:{{i.infile | fn2}}.{{args.outunit}}{{args.outform | lambda x: '_t' if x else ''}}.txt",
	args   = Diot(
		inunit  = 'count',
		outunit = 'tpm',
		inopts  = Diot(rnames = True, cnames = True),
		meanfl  = 1,
		nreads  = 50000000,
		refexon = params.refexon.value,
		inform  = None,
		outform = None
	)
)

pRNASeqDEG = proc_factory(
	desc   = 'Detect DEGs from RNA-seq data.',
	config = Diot(annotate = """
	@input:
		efile:file: The expression matrix
			- Columns other than samples in gfile will be used as annotations
			- See `args.mapping`
		gfile:file: The group information
			- For example:
				```
				Sample	[Patient	]Group
				sample1	[patient1	]group1
				sample2	[patient1	]group1
				sample3	[patient2	]group1
				sample4	[patient2	]group2
				sample5	[patient3	]group2
				sample6	[patient3	]group2
				```
			- If it is not paired comparison, you should not include Patient.
	@output:
		outfile:file: The DEG list
		outdir:file:  The output directory containing deg list and plots
	@args:
		tool  : The tool used to detect DEGs. Default: 'deseq2' (edger is also available).
		inopts: Options to read `infile`. Default: `Diot(cnames = True, rnames = True)`
		cutoff: The cutoff used to filter the results. Default: `0.05`
			- `0.05` implies `{"by": "p", "value": "0.05", "sign": "<"}`
		ggs   : The ggs for each plot. Default:
			- For heatmap: should be the `draw` argument from `plot.heatmap2` in `plot.r`
			- Not available for `mdsplot`.
			- Others are empty `Diot()`s
			- To disable a plot: `ggs.heatmap = FALSE`
		params: Parameters for each plot. Default:
			- `volplot`: `Diot(pcut = 0.05, logfccut = 2)`
			- `maplot` : `Diot(pcut = 0.05)`
			- `heatmap`: `Diot(ngenes = None, <other arguments for plot.heatmap2's params>)`, all genes in heatmap or a number for up/down genes in heatmap
		devpars: Parameters for png. Default: `{'res': 300, 'width': 2000, 'height': 2000}`
		mapping: Probe to gene mapping file. If not provided, assume genes are used as rownames. This could be:
			- A column name or a number (index without samples, starting from 1) in `i.efile` to specify which column is to use a gene names
			- A file path with probe-to-gene mapping, could also include other annotations
				- A suffix(`file:1` or `file:Gene`) is available to specify the gene column. An integer indicates no header for the file, while a columname indicates `header = TRUE`
				- Without suffix, rownames will not be replaced and header is TRUE anyway.
				- Columns are ignored if this is provided.
	"""),
	lang   = params.Rscript.value,
	input  = "efile:file, gfile:file",
	output = """outfile:file:{{i.efile | stem | stem}}-{{i.gfile | stem
		}}.DEGs/{{i.efile | stem | stem}}-{{i.gfile | stem
		}}.degs.xls, outdir:dir:{{i.efile | stem | stem}}-{{i.gfile | stem}}.DEGs""",
	args = Diot(
		tool    = 'deseq2',
		inopts  = Diot(cnames = True, rnames = True, dup = 'drop'),
		mapping = "",
		cutoff  = 0.05,
		ggs = Diot(
			mdsplot = Diot(),
			volplot = Diot(),
			maplot  = Diot(),
			heatmap = Diot(),
			qqplot  = Diot(labs = {'x': 'Expected', 'y': 'Observed -log10(PValue)'})
		),
		params  = Diot(
			volplot = Diot(logfccut = 2),
			maplot = Diot(),
			heatmap = Diot(ngenes = None, show_row_names = False)),
		devpars = Diot(res = 300, width = 2000, height = 2000)
	)
)

pCoexp = proc_factory(
	desc = "Get co-expression of gene pairs in the expression matrix.",
	config = Diot(annotate = """
	@name:
		pCoexp
	@description:
		Get co-expression of gene pairs in the expression matrix.
	"""))
pCoexp.input       = "infile:file"
pCoexp.output      = "outfile:file:{{i.infile | fn}}.coexp, outpval:file:{{i.infile | fn}}.pval"
pCoexp.args.method = 'pearson'
pCoexp.args.pval   = False
pCoexp.lang        = params.Rscript.value

pExprSimulate = proc_factory(
	desc = "Simulate expression values",
	config = Diot(annotate = """
	@name:
		pExprSimulate
	"""))
pExprSimulate.input         = 'seed:var'
pExprSimulate.output        = 'outfile:file:exprsim.{{i.seed if isinstance(i.seed, int) else "noseed"}}.txt'
pExprSimulate.args.nsamples = 100
pExprSimulate.args.ngenes   = 1000
pExprSimulate.args.slabel   = 'Sample'
pExprSimulate.args.glabel   = 'Gene'
pExprSimulate.args.params   = Diot()
pExprSimulate.lang          = params.Rscript.value
