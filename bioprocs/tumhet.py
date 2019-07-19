"""A set of processes for Tumor heterogeneity analysis"""
from modkit import Modkit
from pyppl import Box, Proc
from . import rimport, params, delefactory, procfactory
from .utils import fs2name
Modkit().delegate(delefactory())

@procfactory
def _pSciClone():
	"""
	@input:
		`vfvcfs:files`: The VCF files of mutations of each sample
			- Single-sample/paired-sample vcf files.
			- If it is paired, you have to specify which sample (`args.vfsamcol`) is the target.
		`cnvcfs:files`: The VCF files of copy number variations of each sample
	@output:
		`outdir:dir`: The output directory.
	@args:
		`params`  : Other parameters for original `sciClone` function. Default: `Box()`
		`exfile`  : The regions to be excluded. In BED3 format
		`vfsamcol`: The index of the target sample in mutation VCF file, 1-based. Default: `1`
		`cnsamcol`: The index of the target sample in copy number VCF file, 1-based. Default: `1`
		`varcount`: An R function string to define how to get the variant allele count.
			- If this function returns `NULL`, record will be skipped.
			- It can use the sample calls (`fmt`) and also the record info (`info`)
			- Both `function(fmt) ...` and `function(fmt, info) ...` can be used.
			- Don't include `info` if not necessary. This saves time.
			- This function can return the variant count directly, or
			- an R `list` like: `list(count = <var count>, depth = <depth>)`.
			- By default, the `depth` will be read from `fmt$DP`
		`cncount` : An R function string to define how to get the copy number.
			- Similar as `varcount`
			- Returns copy number directly, or
			- an R `list` like: `list(cn = <copy number>, end = <end>, probes = <probes>)`
			- `end` defines where the copy number variation stops
			- `probes` defines how many probes cover this copy number variantion.
	"""
	return Box(
		desc   = "Clonality analysis using SciClone.",
		input  = "vfvcfs:files, cnvcfs:files",
		output = "outdir:dir:{{i.vfvcfs | fs2name}}.sciclone",
		envs   = Box(fs2name  = fs2name),
		lang   = params.Rscript.value,
		args   = Box(
			params   = Box(),
			exfile   = "",
			vfsamcol = 1,     # the first sample is the target sample in variant vcf
			cnsamcol = 1,     # the first sample is the target sample in copy number vcf
			# how to get the var count
			varcount = 'function(fmt) as.integer(unlist(strsplit(fmt$AD, ","))[2])',
			cncount  = 'function(fmt) fmt$CN' # how to get the copy number
		)
	)

@procfactory
def _pTMBurden():
	"""
	@input:
		infile: The input MAF file
	@output:
		outfile: The tumor mutation burden file
	@args:
		type: The type of mutation burden.
			- `nonsyn`: Counting nonsynonymous mutations
	"""
	return Box(
		desc   = 'Calculation of tumor mutation burden.',
		input  = 'infile:file',
		output = 'outfile:file:{{i.infile | stem}}.tmb.txt',
		args   = Box(type = 'nonsyn'),
		lang   = params.python.value
	)

@procfactory
def _pPyClone():
	"""
	@input:
		`vfvcfs:files`: The VCF files of mutations of each sample
			- Single-sample/paired-sample vcf files.
			- If it is paired, you have to specify which sample (`args.vfsamcol`) is the target.
		`cnvcfs:files`: The VCF files of copy number variations of each sample
	@output:
		`outdir:dir`: The output directory.
	@args:
		params (dict): Other parameters for original `PyClone run_analysis_pipeline` function. Default: `Box()`
		vfsamcol (int): The index of the target sample in mutation VCF file, 1-based. Default: `1`
		cnsamcol (int): The index of the target sample in copy number VCF file, 1-based. Default: `1`
		varcount (str): A python lambda string to define how to get #ref, #alt and allele frequency ( `ref`, `alt`, `af`).
			- If this function returns `None`, record will be skipped.
			- It can use the sample calls (`fmt`) and also the record info (`info`), meaning both `function(fmt) ...` and `function(fmt, info) ...` can be used.
			- This function should return a dict with keys `ref`, `alt`, `dp` and `af`.
			- By default, `ref` and `alt` will get from `FORMAT[AD]`; `af` from `FORMAT[AF]`, if not provided, then `alt/sum(FORMAT[AD])`
		cncount (str): An python lambda string to define how to get the copy number.
			- Similar as `varcount`, returns a `dict` of copy number of ref and alt alleles (key: `major` and `minor`), as well as end of the copy number region (`end`).
			- By default, `major` and `minor` will get from `FORMAT[CN]` if it is a `tuple/list`. If it is an integer, then it's `major` and `minor` will be `0`. Otherwise, `major` and `minor` will be `2` and `0`, respectively.
		pyclone (str): The path of `PyClone`.
	"""
	return Box(
		desc   = "Run PyClone for clonality analysis.",
		input  = "vfvcfs:files, cnvcfs:files",
		output = "outdir:dir:{{i.vfvcfs | fs2name}}.pyclone",
		lang   = params.python.value,
		envs   = Box(fs2name  = fs2name),
		args   = Box(
			params   = Box(),
			vfsamcol = 1, # 1-based
			cnsamcol = 1,
			varcount = None,
			cncount  = None,
			pyclone  = params.pyclone.value,
		)
	)

@procfactory
def _pQuantumClone():
	"""
	@description:
		Clonality analysis using QuantumClone:
		https://academic.oup.com/bioinformatics/article/34/11/1808/4802225
	@input:
		`vfvcfs:files`: The input vcf files with mutations
	@output:
		`outdir:dir`: The output directory
	@args:
		`params`  : other parameters for `QuantumClone`'s `One_step_clustering`
		`vfsamcol`: The index of the target sample in mutation VCF file, 1-based. Default: `1`
		`varcount`: An R function string to define how to get the variant allele count. Default: `function(fmt) as.integer(unlist(strsplit(fmt$AD, ","))[2])`
			- If this function returns `NULL`, record will be skipped.
			- It can use the sample calls (`fmt`) and also the record info (`info`)
			- Both `function(fmt) ...` and `function(fmt, info) ...` can be used.
			- Don't include `info` if not necessary. This saves time.
			- This function can return the variant count directly, or
			- an R `list` like: `list(count = <var count>, depth = <depth>)`.
			- By default, the `depth` will be read from `fmt$DP`
		`nthread` : # threads to use. Default: `1`
	"""
	pQuantumClone               = Proc(desc = "Clonality analysis using QuantumClone")
	pQuantumClone.input         = 'vfvcfs:files'
	pQuantumClone.output        = "outdir:dir:{{i.vfvcfs | fs2name}}.qclone"
	pQuantumClone.envs.fs2name  = fs2name
	pQuantumClone.args.params   = Box()
	pQuantumClone.args.vfsamcol = 1 # 1-based
	pQuantumClone.args.varcount = 'function(fmt) as.integer(unlist(strsplit(fmt$AD, ","))[2])'
	pQuantumClone.args.nthread  = 1
	pQuantumClone.lang          = params.Rscript.value
	return pQuantumClone

@procfactory
def _pTheta():
	"""
	@description:
		Run THetA2 for tumor purity calculation
		Set lower MIN_FRAC if interval is not enough and NO_CLUSTERING if it raises
		"No valid Copy Number Profiles exist", but have to pay attention to the results.
		(see: https://groups.google.com/forum/#!topic/theta-users/igrEUol3sZo)
	@args:
		`affysnps`: The affymetrix Array snps, or other candidate snp list, in BED6-like format
			- The first 6 columns should be in BED6 format
			- The 7th column is reference allele, and 8th column is mutation allele.
	@install:
		`conda install -c bioconda theta2`
		`conda install -c bioconda bam-readcount`
	"""
	pTheta                    = Proc(desc = "Tumor purity calculation using THetA2")
	pTheta.input              = 'itvfile:file, tumbam:file, normbam:file'
	pTheta.output             = 'outdir:dir:{{i.itvfile | fn2}}.theta'
	pTheta.args.params        = Box(BAF = True, FORCE = True, n = 2)
	pTheta.args.bam_readcount = params.bam_readcount.value
	pTheta.args.ref           = params.ref.value
	pTheta.args.theta         = params.theta2.value
	pTheta.args.nthread       = 1
	pTheta.args.affysnps      = params.affysnps.value
	pTheta.lang               = params.python.value
	return pTheta

@procfactory
def _pSuperFreq():
	pSuperFreq              = Proc(desc = "Subclonal analysis with superFreq")
	pSuperFreq.input        = "indir:dir, gfile:file"
	pSuperFreq.output       = "outdir:dir:{{i.indir | fn2}}-{{i.gfile | fn2}}.superfreq"
	pSuperFreq.args.nthread = 1
	pSuperFreq.args.baits   = '' # target regions
	pSuperFreq.args.ref     = params.ref.value
	pSuperFreq.args.resdir  = params.superfreq_res.value
	pSuperFreq.args.genome  = params.genome.value
	pSuperFreq.args.params  = Box(
		systematicVariance = .02, maxCov = 150, BQoffset = 33,
		mode = 'exome', splitRun = True
	)
	pSuperFreq.lang         = params.Rscript.value
	return pSuperFreq

@procfactory
def _pClonEvol():
	"""
	@input:
		mutfile: The mutation file.
		saminfo: The sample information file.
	@output:
		outdir: The output directory.
	@args:
		inopts: The input options to read the mutation file.
		params: The parameters for individual `ClonEvol` functions.
	"""
	pClonEvol             = Proc(desc = "Inferring and visualizing clonal evolution in multi-sample cancer sequencing.")
	pClonEvol.input       = 'mutfile:file, saminfo:file'
	pClonEvol.output      = 'outdir:dir:{{i.mutfile | stem}}.clonevol'
	pClonEvol.lang        = params.Rscript.value
	pClonEvol.args.inopts = Box(rnames = False, cnames = True)
	pClonEvol.args.params = Box({
		'plot.variant.clusters': Box(
			# see https://rdrr.io/github/hdng/clonevol/man/plot.variant.clusters.html
		),
		'plot.cluster.flow': Box(
			# see https://rdrr.io/github/hdng/clonevol/man/plot.cluster.flow.html
		),
		'infer.clonal.models': Box({
			# see https://rdrr.io/github/hdng/clonevol/man/infer.clonal.models.html
			"founding.cluster": 1,
			"cluster.center": "mean",
			"sum.p.cutoff": 0.05,
			"alpha": 0.05
		}),
		'transfer.events.to.consensus.trees': Box({
			# see https://rdrr.io/github/hdng/clonevol/man/transfer.events.to.consensus.trees.html
			"event.col.name": "gene"
		}),
		'convert.consensus.tree.clone.to.branch': Box({
			# see https://rdrr.io/github/hdng/clonevol/man/convert.consensus.tree.clone.to.branch.html
			"branch.scale": "sqrt"
		}),
		'plot.clonal.models': Box({
			# see https://rdrr.io/github/hdng/clonevol/man/plot.clonal.models.html
			# "clone.shape"                     : 'bell',
			# "bell.event"                      : True,
			# "bell.event.label.color"          : 'blue',
			# "bell.event.label.angle"          : 60,
			# "clone.time.step.scale"           : 1,
			# "bell.curve.step"                 : 2,
			# "merged.tree.plot"                : True,
			# "tree.node.label.split.character" : None,
			# "tree.node.shape"                 : 'circle',
			# "tree.node.size"                  : 30,
			# "tree.node.text.size"             : 0.5,
			# "merged.tree.node.size.scale"     : 1.25,
			# "merged.tree.node.text.size.scale": 2.5,
			# "merged.tree.cell.frac.ci"        : False,
			# "mtcab.event.sep.char"            : ',',
			"mtcab.branch.text.size"          : .8,
			# "mtcab.branch.width"              : 0.75,
			# "mtcab.node.size"                 : 3,
			# "mtcab.node.label.size"           : 1,
			"mtcab.node.text.size"            : .8,
			# "cell.plot"                       : True,
			# "num.cells"                       : 100,
			# "cell.border.size"                : 0.25,
			# "cell.border.color"               : 'black',
			# "clone.grouping"                  : 'horizontal',
			# "show.score"                      : False,
			# "cell.frac.ci"                    : True,
			# "disable.cell.frac"               : False
		})
	})
	pClonEvol.args.devpars = Box(width = 2000, height = 2000, res = 300)
	return pClonEvol

@procfactory
def _pPyClone2ClonEvol():
	pPyClone2ClonEvol               = Proc(desc = "Convert PyClone results to ClonEvol input format.")
	pPyClone2ClonEvol.input         = 'indir:dir'
	pPyClone2ClonEvol.output        = 'outfile:file:{{i.indir | fn}}.clonevol.txt'
	pPyClone2ClonEvol.args.refgene  = params.refgene.value
	pPyClone2ClonEvol.args.drivers  = []
	pPyClone2ClonEvol.args.bedtools = params.bedtools.value
	pPyClone2ClonEvol.lang          = params.python.value
	return pPyClone2ClonEvol

