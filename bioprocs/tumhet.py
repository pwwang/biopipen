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
		muts: Somatic mutations. Could be one of:
			- Single(paired)-sample VCF files, separared by comma
				- Must have `FORMAT/AF` and `FORMAT/AD`
			- Multi-sample VCF file
			- TCGA MAF file
				- Must have `t_alt_count` and `t_ref_count`
			- SciClone TSV file
				- Chr, Pos, RefCount, VarCount, VAF*100
		cnvs: Copy numbers. Could be one of:
			- Single-sample VCF files, separared by comma
			- Multi-sample VCF files, must have `INFO/END` and `FORMAT/CN`
			- SciClone TSV file
				- Chr, Start, End, CopyNumber
	@output:
		outdir: The output directory.
	@args:
		params (Box)  : Other parameters for original `sciClone` function. Default: `Box()`
		exfile (file)  : The regions to be excluded. In BED3 format
		mutctrl (int|NoneType): Index of the control sample in paired-sample mutation VCF files or multi-sample VCF file, 0-based.
			- `None` indicates no control sample.
			- For paired-sample VCF files, if `None` specified, second sample(1) will be used as control.
		cnctrl (int|NoneType): Index of the control sample in copy number VCF files or multi-sample VCF file, 0-based.
	"""
	return Box(
		desc   = "Clonality analysis using SciClone.",
		input  = "muts:var, cnvs:var",
		output = "outdir:dir:{{i.muts | bn | .split('.')[0]}}.sciclone",
		lang   = params.Rscript.value,
		args   = Box(
			params  = Box(),
			exfile  = "",
			mutctrl = None,
			cnctrl  = None
		)
	)

@procfactory
def _pPyClone():
	"""
	@input:
		muts: Somatic mutations. Could be one of:
			- Single(paired)-sample VCF files, separated by comma (i.e. a.vcf,b.vcf)
				- Must have `FORMAT/AF` and `FORMAT/AD`
			- Multi-sample VCF file (i.e. a.vcf)
			- TCGA MAF file (i.e. a.maf)
				- Must have `t_alt_count` and `t_ref_count` columns
			- PyClone TSV file (i.e. a.tsv or a.txt)
				- See https://github.com/aroth85/pyclone/blob/master/examples/mixing/tsv/SRR385941.tsv without `*_cn` columns
				- `mutation_id` should be either `<sample>:<gt>:<chr>:<pos>` or `<chr>:<pos>`
		cnvs: Copy numbers. Could be one of:
			- Single-sample VCF files, separated by comma
			- Multi-sample VCF file
				- Must have `INFO/END` and `FORMAT/CN`
			- PyClone TSV file
				- See https://github.com/aroth85/pyclone/blob/master/examples/mixing/tsv/SRR385941.tsv without `*_counts` and `variant_freq` columns
	@output:
		`outdir:dir`: The output directory.
	@args:
		pyclone (str): Path to `PyClone`.
		bcftools(str): Path to `bcftools`, used to get information from VCF file.
		bedtools(str): Path to `bedtools`, used to intersect mutations with CNVs.
		params  (Box): Other parameters for original `PyClone run_analysis_pipeline` function.
		nthread (int): Number of threads to use by openblas.
		mutctrl (int|NoneType): Index of the control sample in paired-sample mutation VCF files or multi-sample VCF file, 0-based.
			- `None` indicates no control sample.
			- For paired-sample VCF files, if `None` specified, second sample(1) will be used as control.
		cnctrl (int|NoneType): Index of the control sample in copy number VCF files or multi-sample VCF file, 0-based.
	"""
	return Box(
		desc   = "Clonality analysis using PyClone",
		input  = "muts:var, cnvs:var",
		output = "outdir:dir:{{i.muts | bn | .split('.')[0]}}.pyclone",
		lang   = params.python.value,
		args   = Box(
			pyclone  = params.pyclone.value,
			bcftools = params.bcftools.value,
			bedtools = params.bedtools.value,
			refgene  = params.refgene.value,
			nthread  = 1,
			params   = Box(),
			mutctrl  = None,
			cnctrl   = None,
		)
	)

@procfactory
def _pAllFIT():
	"""
	@description:
		All-FIT - Allele-Frequency-based Imputation of Tumor Purity infers specimen purity from tumor-only samples sequenced with deep sequencing. It is developed in Khiabanian Lab by Jui Wan Loh and Hossein Khiabanian.
		See: https://github.com/KhiabanianLab/All-FIT
	@input:
		infile: The input file, could be one of:
			- Single(paired)-sample VCF file with `FORMAT/AF` and `FORMAT/AD`
			- All-FIT format (see: https://github.com/KhiabanianLab/All-FIT/blob/master/test/input/sampleFile1.xls)
		cnfile: The copy number variation file, could be one of:
			- Omitted (then Ploidy in `infile` must be provided)
			- Single(paired)-sample VCF file with `INFO/END` and `FORMAT/CN`
			- Bed3 file with 4th column as the copy number.
	@output:
		outfile: The output file with purity
		outdir: The output directory with output file and figures.
	@args:
		allfit (str): Path to All-FIT.py
		bcftools (str): Path to `bcftools`, used to get information from VCF file.
		bedtools (str): Path to `bedtools`, used to intersect mutations with CNVs.
		params (Box): Other parameters for All-FIT.py
		mutctrl (int|NoneType): Index of the control sample in paired-sample mutation VCF file, 0-based.
			- `None` indicates no control sample.
			- For paired-sample VCF files, if `None` specified, second sample(1) will be used as control.
		cnctrl (int|NoneType): Index of the control sample in copy number VCF file, 0-based.
		nthread (int): Number of threads to use for openblas.
	"""
	return Box(
		desc = "Allele-Frequency-based Imputation of Tumor Purity Inference",
		input = "infile:file, cnfile:var",
		output = [
			"outfile:file:{{i.infile | stem | @append: '.allfit'}}/{{i.infile | stem | @append: '.purity.txt'}}",
			"outdir:dir:{{  i.infile | stem | @append: '.allfit'}}"
		],
		lang = params.python.value,
		args = Box(
			allfit   = params.allfit.value,
			bcftools = params.bcftools.value,
			bedtools = params.bedtools.value,
			params   = Box(t = 'somatic'),
			mutctrl  = None,
			cnctrl   = None,
			nthread  = 1
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
		desc    = 'Calculation of tumor mutation burden.',
		input   = 'infile:file',
		output  = 'outfile:file:{{i.infile | stem}}.tmb.txt',
		args    = Box(type = 'nonsyn'),
		lang    = params.python.value,
		mutctrl = None,
		cnctrl  = None,
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

		THetA2 needs an interval file with the coverages of both tumor and normal and SNP files
		with counts of mutation allele and reference allele for both tumor and normal as well.
	@input:
		infile: A bed3 file, or a THetA2 interval file.
			- See: https://github.com/raphael-group/THetA/blob/master/example/Example.intervals
		tumbam: bam file for tumor or SNP file for tumor.
			- See: https://raw.githubusercontent.com/raphael-group/THetA/master/example/TUMOR_SNP.formatted.txt
		normbam: bam file for normal or SNP file for normal.
			- See: https://raw.githubusercontent.com/raphael-group/THetA/master/example/NORMAL_SNP.formatted.txt
	@output:
		outfile: The output file with purity.
		outdir: The output directory with output files and figures.
	@args:
		theta (str): Path to THetA2
		bam_readcount (str): Path to bam_readcount.
		bedtools (str): Path to bedtools, used to extract coverage of regions.
		samtools (str): Path to samtools, used to index bam file if possible.
		ref (file): The reference genome file.
		params (Box): Other parameters for THetA2.
		nthread (int): The number of threads to use.
		affysnps (str): The affymetrix Array snps, or other candidate snp list, in BED6-like format
			- The first 6 columns should be in BED6 format
			- The 7th column is reference allele, and 8th column is mutation allele.
	"""
	return Box(
		desc   = "Tumor purity calculation using THetA2",
		lang   = params.python.value,
		input  = 'infile:file, tumbam:file, normbam:file',
		output = [
			'outfile:file:{{i.infile | stem}}.theta/{{i.infile | stem}}.purity.txt',
			'outdir:dir:{{i.infile | stem}}.theta'
		],
		args = Box(
			theta         = params.theta2.value,
			bedtools      = params.bedtools.value,
			bam_readcount = params.bam_readcount.value,
			samtools      = params.samtools.value,
			params        = Box(),
			ref           = params.ref.value,
			nthread       = 1,
			affysnps      = params.affysnps.value,
		)
	)

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
		mutfile: The mutation file or output directory from PyClone.
			- https://github.com/hdng/clonevol/issues/4#issuecomment-280997440
			- VAF, CCF should x100
		samfile: The sample information file.
		drivers: The driver genes.
			- One per line or
			- Mutsig file with top 10 genes
	@output:
		outdir: The output directory.
	@args:
		inopts: The input options to read the mutation file.
		params: The parameters for individual `ClonEvol` functions.
	"""
	return Box(
		desc   = "Inferring and visualizing clonal evolution in multi-sample cancer sequencing",
		input  = 'mutfile:file, samfile:file, drivers:var',
		output = 'outdir:dir:{{i.mutfile | stem}}.clonevol',
		lang   = params.Rscript.value,
		args   = Box(
			# only for clonevol input format
			inopts   = Box(rnames = False,  cnames = True),
			drivers  = [],
			refgene  = params.refgene.value,
			bedtools = params.bedtools.value,
			devpars  = Box(width = 2000, height = 2000, res = 300),
			params   = Box({
				'plot.variant.clusters': Box(
					# see https://rdrr.io/github/hdng/clonevol/man/plot.variant.clusters.html
				),
				'plot.cluster.flow': Box({
					# see https://rdrr.io/github/hdng/clonevol/man/plot.cluster.flow.html
				}),
				'infer.clonal.models': Box({
					# see https://rdrr.io/github/hdng/clonevol/man/infer.clonal.models.html
					"founding.cluster": 1,
					"cluster.center": "mean",
					"sum.p.cutoff": 0.05,
					"alpha": 0.05,
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
		)
	)

@procfactory
def _pSchism():
	"""
	@description:
		Infer subclonal hierarchy and the tumor evolution from somatic mutations using SCHISM.
		See: https://github.com/KarchinLab/SCHISM
	@input:
		infile: The input mutation file
			- In format of:
				```
				sampleID	mutationID	referenceReads	variantReads	copyNumber	clusterID
				S1	0	368	132	2	0
				S2	1	381	119	2	0
				S3	2	367	133	3	1
				...
				```
			- Or `pPyClone` output directory
			- Or a MAF file
		purity: The purity of the samples, two columns:
			- Sample[tab]Purity, without header.
			- If omitted, will assume 100% purity for all samples.
	@output:
		outdir: The output directory.
	@args:
		schism   (str) : Path to runSchism
		params   (Box) : Other parameters for runSchism in yaml config.
		fishplot (bool): Whether generate a fishplot or not, requires fishplot install with R.
			- Not implemented yet
		Rscript  (str) : Path to Rscript to run fishplot
	"""
	return Box(
		desc   = 'Infer subclonal hierarchy and the tumor evolution from somatic mutations using SCHISM',
		lang   = params.python.value,
		input  = 'infile:file, purity:file',
		output = 'outdir:dir:{{i.infile | stem}}.schism',
		args   = Box(
			schism   = params.schism.value,
			dot      = params.dot.value,
			fishplot = True,
			Rscript  = params.Rscript.value,
			devpars  = Box(res = 100),
			params   = Box(
				cellularity_estimator = Box(),
				hypothesis_test       = Box(),
				genetic_algorithm     = Box()
			)
		)
	)

@procfactory
def _pLichee():
	"""
	@description:
		Fast and scalable inference of multi-sample cancer lineages using LICHeE.
		See: https://github.com/pwwang/lichee
	@input:
		infile: The input file, could be one of:
			- LICHeE input file. Set `args.params.cp = True` if input data represents cell prevalence instead of VAF.
			- A pPyClone output directory.
			- A MAF file.
	@output:
		outdir: The output directory
	@args:
		dot (str): Path to dot to generate figures
		params (Box): Other parameters for `lichee`.
			- Set a larger `e` if you have noisy data. Default is 0.1.
	"""
	return Box(
		desc   = 'Fast and scalable inference of multi-sample cancer lineages using LICHeE',
		input  = 'infile:file',
		output = 'outdir:dir:{{i.infile | stem}}.lichee',
		lang   = params.python.value,
		args   = Box(
			dot      = params.dot.value,
			lichee   = params.lichee.value,
			fishplot = True,
			Rscript  = params.Rscript.value,
			devpars  = Box(res = 100),
			params   = Box(maxVAFAbsent = 0.005, minVAFPresent = 0.005)
		)
	)

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



