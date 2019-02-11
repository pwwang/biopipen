# A set of processes for Tumor heterogeneity analysis
from pyppl import Box, Proc
from . import rimport, params
from .utils import fs2name

"""
@name:
	pSciClone
@description:
	Run sciClone for subclonal analysis.
@input:
	`vfvcfs:files`: The VCF files of mutations of each sample
	`cnvcfs:files`: The VCF files of copy number variations of each sample
@output:
	`outdir:dir`: The output directory.
@args:
	`params`  : Other parameters for original `sciClone` function. Default: `Box()`
	`exfile`  : The regions to be excluded. In BED3 format
	`vfsamcol`: The index of the target sample in mutation VCF file, 1-based. Default: `1`
	`cnsamcol`: The index of the target sample in copy number VCF file, 1-based. Default: `1`
	`varcount`: An R function string to define how to get the variant allele count. Default: `function(fmt) as.integer(unlist(strsplit(fmt$AD, ","))[2])`
		- If this function returns `NULL`, record will be skipped.
		- It can use the sample calls (`fmt`) and also the record info (`info`)
		- Both `function(fmt) ...` and `function(fmt, info) ...` can be used.
		- Don't include `info` if not necessary. This saves time.
		- This function can return the variant count directly, or 
		- an R `list` like: `list(count = <var count>, depth = <depth>)`.
		- By default, the `depth` will be read from `fmt$DP`
	`cncount` : An R function string to define how to get the copy number. Default: `function(fmt) fmt$CN`
		- Similar as `varcount`
		- Returns copy number directly, or
		- an R `list` like: `list(cn = <copy number>, end = <end>, probes = <probes>)`
		- `end` defines where the copy number variation stops
		- `probes` defines how many probes cover this copy number variantion.
"""
pSciClone               = Proc(desc = "Run sciClone")
pSciClone.input         = "vfvcfs:files, cnvcfs:files"
pSciClone.output        = "outdir:dir:{{i.vfvcfs | fs2name}}.sciclone"
pSciClone.envs.fs2name  = fs2name
pSciClone.envs.rimport  = rimport
pSciClone.args.params   = Box()
pSciClone.args.exfile   = ""
pSciClone.args.vfsamcol = 1 # the first sample is the target sample in variant vcf
pSciClone.args.cnsamcol = 1 # the first sample is the target sample in copy number vcf
pSciClone.args.varcount = 'function(fmt) as.integer(unlist(strsplit(fmt$AD, ","))[2])' # how to get the var count
pSciClone.args.cncount  = 'function(fmt) fmt$CN' # how to get the copy number 
pSciClone.lang          = params.Rscript.value
pSciClone.script        = "file:scripts/tumhet/pSciClone.r"

"""
@name:
	pPyClone
@description:
	Run PyClone for subclonal analysis
@input:
	`vfvcfs:files`: The VCF files of mutations of each sample
	`cnvcfs:files`: The VCF files of copy number variations of each sample
@output:
	`outdir:dir`: The output directory.
@args:
	`params`  : Other parameters for original `PyClone run_analysis_pipeline` function. Default: `Box()`
	`vfsamcol`: The index of the target sample in mutation VCF file, 1-based. Default: `1`
	`cnsamcol`: The index of the target sample in copy number VCF file, 1-based. Default: `1`
	`varcount`: A python lambda string to define how to get the variant allele count. Default: `lambda fmt: fmt.get("AD") and fmt.get("AD")[1]`
		- If this function returns `None`, record will be skipped.
		- It can use the sample calls (`fmt`) and also the record info (`info`)
		- Both `function(fmt) ...` and `function(fmt, info) ...` can be used.
		- This function can return the variant count directly, or 
		- a `dict` like: `dict(count = <var count>, depth = <depth>)`.
		- By default, the `depth` will be read from `fmt.DP`
	`cncount` : An python lambda string to define how to get the copy number. Default: `lambda fmt: fmt.get("CN")`
		- Similar as `varcount`
		- Returns copy number directly, or
		- a `dict` like: `dict(cn = <copy number>, end = <end>)`
		- `end` defines where the copy number variation stops
"""
pPyClone               = Proc(desc = "Run pyclone.")
pPyClone.input         = "vfvcfs:files, cnvcfs:files"
pPyClone.output        = "outdir:dir:{{i.vfvcfs | fs2name}}.pyclone"
pPyClone.envs.fs2name  = fs2name
pPyClone.args.params   = Box()
pPyClone.args.vfsamcol = 1 # 1-based
pPyClone.args.cnsamcol = 1
pPyClone.args.varcount = 'lambda fmt: fmt.get("AD") and fmt.get("AD")[1]'
pPyClone.args.cncount  = 'lambda fmt: fmt.get("CN")'
pPyClone.args.pyclone  = params.pyclone.value
pPyClone.lang          = params.python.value
pPyClone.script        = "file:scripts/tumhet/pPyClone.py"

"""
@name:
	pQuantumClone
@description:
	Run QuantumClone: https://academic.oup.com/bioinformatics/article/34/11/1808/4802225
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
pQuantumClone               = Proc(desc = "Run QuantumClone")
pQuantumClone.input         = 'vfvcfs:files'
pQuantumClone.output        = "outdir:dir:{{i.vfvcfs | fs2name}}.qclone"
pQuantumClone.envs.fs2name  = fs2name
pQuantumClone.envs.rimport  = rimport
pQuantumClone.args.params   = Box()
pQuantumClone.args.vfsamcol = 1 # 1-based
pQuantumClone.args.varcount = 'function(fmt) as.integer(unlist(strsplit(fmt$AD, ","))[2])'
pQuantumClone.args.nthread  = 1
pQuantumClone.lang          = params.Rscript.value
pQuantumClone.script        = "file:scripts/tumhet/pQuantumClone.r"

"""
@name:
	pTheta
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
pTheta                    = Proc(desc = "Run THetA2 for tumor purity calculation")
pTheta.input              = 'itvfile:file, tumbam:file, normbam:file'
pTheta.output             = 'outdir:dir:{{i.itvfile | fn2}}.theta'
pTheta.args.params        = Box(BAF = True, FORCE = True, n = 2)
pTheta.args.bam_readcount = params.bam_readcount.value
pTheta.args.ref           = params.ref.value
pTheta.args.theta         = params.theta2.value
pTheta.args.nthread       = 1
pTheta.args.affysnps      = params.affysnps.value
pTheta.lang               = params.python.value
pTheta.script             = "file:scripts/tumhet/pTheta.py"
