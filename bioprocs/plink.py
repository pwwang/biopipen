"""
Procs for plink 1.9
"""
from pyppl import Proc, Box
from . import params, rimport, bashimport
from .vcf import pVcf2Plink
from .tcgamaf import pGTMat2Plink
from . import delefactory, procfactory
from modkit import Modkit
Modkit().delegate(delefactory())

@procfactory
def _pPlink2Vcf():
	"""
	@name:
		pPlink2Vcf
	"""
	pPlinkFromVcf   = pVcf2Plink.copy()
	return pPlink2Vcf

@procfactory
def _pPlinkFromGTMat():
	"""
	@name:
		pPlinkFromGTMat
	"""
	pPlinkFromGTMat = pGTMat2Plink.copy()
	return pPlinkFromGTMat

@procfactory
def _pPlinkStats():
	"""
	@name:
		pPlinkStats
	"""
	pPlinkStats               = Proc(desc = 'Do basic statistics with plink 1.9')
	pPlinkStats.input         = 'indir:dir'
	pPlinkStats.output        = 'outdir:dir:{{i.indir | fn}}.plinkStats'
	pPlinkStats.args.plink    = params.plink.value
	pPlinkStats.args.nthread  = False
	pPlinkStats.args.params   = Box({
		'hardy'    : True,
		'het'      : True,
		'freq'     : True,
		'missing'  : True,
		'check-sex': True,
	})
	pPlinkStats.args.cutoff   = Box({
		'hardy.hwe'     : 1e-5,
		'hardy.mingt'   : None,
		'het'           : 3,
		'freq'          : 0.01,
		'missing.sample': .95,
		'missing.snp'   : .95,
	})
	pPlinkStats.args.plot     = Box({
		'hardy.hwe'     : True,
		'hardy.mingt'   : True,
		'het'           : True,
		'freq'          : True,
		'missing.sample': True,
		'missing.snp'   : True,
	})
	pPlinkStats.args.devpars  = Box(res=300, width=2000, height=2000)
	pPlinkStats.envs.rimport  = rimport
	pPlinkStats.lang          = params.Rscript.value
	pPlinkStats.script        = "file:scripts/plink/pPlinkStats.r"
	return pPlinkStats

@procfactory
def _pPlinkSampleFilter():
	"""
	@name:
		pPlinkSampleFilter
	"""
	pPlinkSampleFilter              = Proc(desc = 'Do sample filtering or extraction using `--keep[-fam]` or `--remove[-fam]`')
	pPlinkSampleFilter.input        = 'indir:dir, samfile:file'
	pPlinkSampleFilter.output       = 'outdir:dir:{{i.indir | bn}}'
	pPlinkSampleFilter.args.plink   = params.plink.value
	pPlinkSampleFilter.args.keep    = True
	pPlinkSampleFilter.args.samid   = 'iid' # both or fid
	pPlinkSampleFilter.args.fam     = False
	pPlinkSampleFilter.args.params  = Box()
	pPlinkSampleFilter.args.nthread = False
	pPlinkSampleFilter.lang         = params.python.value
	pPlinkSampleFilter.script       = "file:scripts/plink/pPlinkSampleFilter.py"
	return pPlinkSampleFilter

@procfactory
def _pPlinkMiss():
	"""
	@name:
		pPlinkMiss
	@description:
		Find samples and snps with missing calls, calculate the call rates and plot them.
	@input:
		`indir:dir`: The input directory containing .bed/.bim/.fam files
	@output:
		`outdir:dir`: The output directory. Default: `{{i.indir | fn}}.miss`
			- `.imiss`: The miss calls for samples
			- `.lmiss`: The miss calls for snps
			- `.samplecr.fail`: The samples fail sample call rate cutoff (`args.samplecr`)
			- `.snpcr.fail`: The SNPs fail snp call rate cutoff (`args.snpcr`)
	@args:
		`plink`: The path to plink.
		`samplecr`: The sample call rate cutoff. Default: `.95`
		`snpcr`: The SNP call rate cutoff. Default: `.95`
		`plot`: Whether plot the distribution of the call rates? Default: `True`
		`devpars`: The device parameters for the plot. Default: `Box(res=300, width=2000, height=2000)`
	"""
	pPlinkMiss               = Proc(desc = 'Find samples and snps with missing calls')
	pPlinkMiss.input         = 'indir:dir'
	pPlinkMiss.output        = 'outdir:dir:{{i.indir | fn}}.miss'
	pPlinkMiss.args.plink    = params.plink.value
	pPlinkMiss.args.samplecr = .95
	pPlinkMiss.args.snpcr    = .95
	pPlinkMiss.args.plot     = True
	pPlinkMiss.args.devpars  = Box(res=300, width=2000, height=2000)
	pPlinkMiss.envs.rimport  = rimport
	pPlinkMiss.lang          = params.Rscript.value
	pPlinkMiss.script        = "file:scripts/plink/pPlinkMiss.r"
	return pPlinkMiss

@procfactory
def _pPlinkFreq():
	"""
	@name:
		pPlinkFreq
	"""
	pPlinkFreq              = Proc(desc = 'Filter snps with minor allele frequency.')
	pPlinkFreq.input        = 'indir:dir'
	pPlinkFreq.output       = 'outdir:dir:{{i.indir | fn}}.freq'
	pPlinkFreq.args.plink   = params.plink.value
	pPlinkFreq.args.cutoff  = 0.01
	pPlinkFreq.args.plot    = True
	pPlinkFreq.args.devpars = Box(res=300, width=2000, height=2000)
	pPlinkFreq.lang         = params.python.value
	pPlinkFreq.script       = "file:scripts/plink/pPlinkFreq.py"
	return pPlinkFreq

@procfactory
def _pPlinkSexcheck():
	"""
	@name:
		pPlinkSexcheck
	@description:
		Check inconsistency between sex denoted and from genotypes.
	@input:
		`indir:dir`: The input directory containing .bed/.bim/.fam files
	@output:
		`outdir:dir`: The output directory. Default: `{{i.indir | fn}}.sexcheck`
			- `.sexcheck`: The orginal sex check report from `plink`
			- `.sex.fail`: The samples that fail sex check.
	@args:
		`plink`: The path to plink.
	"""
	pPlinkSexcheck               = Proc(desc = 'Check inconsistency between sex denoted and from genotypes.')
	pPlinkSexcheck.input         = 'indir:dir'
	pPlinkSexcheck.output        = 'outdir:dir:{{i.indir | fn}}.sexcheck'
	pPlinkSexcheck.args.plink    = params.plink.value
	pPlinkSexcheck.envs.rimport  = rimport
	pPlinkSexcheck.lang          = params.Rscript.value
	pPlinkSexcheck.script        = "file:scripts/plink/pPlinkSexcheck.r"
	return pPlinkSexcheck

@procfactory
def _pPlinkHet():
	"""
	@name:
		pPlinkHet
	@description:
		Calculate the heterozygosity of each sample.
	@input:
		`indir:dir`: The input directory containing .bed/.bim/.fam files
	@output:
		`outdir:dir`: The output directory. Default: `{{i.indir | fn}}.het`
			- `.het`: The heterozygosity file generated by `plink`.
			- `.het.fail`: The samples fail sample heterozygosity cutoff (`args.cutoff`)
	@args:
		`plink`: The path to plink.
		`cutoff`: The sample heterozygosity cutoff. Default: `3` (mean-3SD ~ mean+3SD)
		`plot`: Whether plot the distribution of the heterozygosity? Default: `True`
		`devpars`: The device parameters for the plot. Default: `Box(res=300, width=2000, height=2000)`
	"""
	pPlinkHet              = Proc(desc = 'Calculate the heterozygosity of each sample')
	pPlinkHet.input        = 'indir:dir'
	pPlinkHet.output       = 'outdir:dir:{{i.indir | fn}}.het'
	pPlinkHet.args.plink   = params.plink.value
	pPlinkHet.args.cutoff  = 3
	pPlinkHet.args.plot    = True
	pPlinkHet.args.devpars = Box(res=300, width=2000, height=2000)
	pPlinkHet.envs.rimport = rimport
	pPlinkHet.lang         = params.Rscript.value
	pPlinkHet.script       = "file:scripts/plink/pPlinkHet.r"
	return pPlinkHet

@procfactory
def _pPlinkHWE():
	"""
	@name:
		pPlinkHWE
	@description:
		Hardy-Weinberg Equilibrium report and filtering.
	@input:
		`indir:dir`: The input directory containing .bed/.bim/.fam files
	@output:
		`outdir:dir`: The output directory. Default: `{{i.indir | fn}}.hwe`
			- `.hwe`: The HWE report by `plink`
			- `.hardy.fail`: The SNPs fail HWE test
	@args:
		`plink`: The path to plink.
		`cutoff`: The HWE p-value cutoff. Default: `1e-5`
		`plot`: Whether plot the distribution of the HWE p-values? Default: `True`
		`devpars`: The device parameters for the plot. Default: `Box(res=300, width=2000, height=2000)`
	"""
	pPlinkHWE = Proc(desc = "Hardy-Weinberg Equilibrium report and filtering.")
	pPlinkHWE.input        = 'indir:dir'
	pPlinkHWE.output       = 'outdir:dir:{{i.indir | fn}}.hwe'
	pPlinkHWE.args.plink   = params.plink.value
	pPlinkHWE.args.cutoff  = 1e-5
	pPlinkHWE.args.mingt   = .05
	pPlinkHWE.args.plot    = True
	pPlinkHWE.args.devpars = Box(res=300, width=2000, height=2000)
	pPlinkHWE.envs.rimport = rimport
	pPlinkHWE.lang         = params.Rscript.value
	pPlinkHWE.script       = "file:scripts/plink/pPlinkHWE.r"
	return pPlinkHWE

@procfactory
def _pPlinkIBD():
	"""
	@name:
		pPlinkIBD
	@description:
		Estimate the degree of recent shared ancestry individual pairs,
		the identity by descent (IBD)
	@input:
		`indir:dir`: The input directory containing .bed/.bim/.fam files
	@output:
		`outdir:dir`: The output directory. Default: `{{i.indir | fn}}.ibd`
			- `.genome`: The original genome report from `plink`
			- `.ibd.png`: The heatmap of PI_HAT
	@args:
		`plink`: The path to plink.
		`indep`: To give a concrete example: the command above that specifies 50 5 0.2 would a) consider a window of 50 SNPs, b) calculate LD between each pair of SNPs in the window, b) remove one of a pair of SNPs if the LD is greater than 0.5, c) shift the window 5 SNPs forward and repeat the procedure.
		`pihat`: The PI_HAT cutoff. Default: 0.1875 (see: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5007749/)
		`plot` : Whether plot the PI_HAT heatmap? Default: `True`
		`devpars`: The device parameters for the plot. Default: `Box(res=300, width=2200, height=1600)`
		`samid`: Sample ids on the heatmap. Default: `iid`
			- Could also be `fid` or `fid<sep>iid`, or an R function: `function(fid, iid)`
		`anno` : The annotation file for the samples. Names must match the ones that are transformed by `args.samid`. Default: `''`
	"""
	pPlinkIBD              = Proc(desc = "Estimate the identity by descent (IBD)")
	pPlinkIBD.input        = 'indir:dir'
	pPlinkIBD.output       = 'outdir:dir:{{i.indir | fn}}.ibd'
	pPlinkIBD.args.plink   = params.plink.value
	pPlinkIBD.args.highld  = params.highld.value
	pPlinkIBD.args.samid   = 'iid' # fid or a function (fid, iid)
	pPlinkIBD.args.indep   = [50, 5, .2]
	pPlinkIBD.args.pihat   = 0.1875 # ref: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5007749/
	pPlinkIBD.args.plot    = True
	pPlinkIBD.args.anno    = ''
	pPlinkIBD.args.seed    = None
	pPlinkIBD.args.devpars = Box(res=300, width=2000, height=2000)
	pPlinkIBD.envs.rimport = rimport
	pPlinkIBD.lang         = params.Rscript.value
	pPlinkIBD.script       = "file:scripts/plink/pPlinkIBD.r"
	return pPlinkIBD

@procfactory
def _pPlinkRemove():
	"""
	@name:
		pPlinkRemove
	@description:
		Remove failed samples and/or SNPs
		The samples/SNPs to be removed should be generated by one of:
		`pPlinkHet`, `pPlinkHWE`, `pPlinkIBD` or `pPlinkMiss`
	@input:
		`indir:dir`: The input directory containing .bed/.bim/.fam files
		`pdir:dir` : The output directory from one of the processes listed in description
			- It could also be the `.fail` file generated by those processes
	@output:
		`outdir:dir`: The output directory containing the `.bed/.bim/.fam` after filtering.
	@args:
		`plink`: The path to plink.
	"""
	pPlinkRemove              = Proc(desc = "Remove failed samples and SNPs")
	pPlinkRemove.input        = 'indir:dir, pdir:dir'
	pPlinkRemove.output       = 'outdir:dir:{{i.indir | fn}}'
	pPlinkRemove.args.plink   = params.plink.value
	pPlinkRemove.envs.rimport = rimport
	pPlinkRemove.lang         = params.Rscript.value
	pPlinkRemove.script       = "file:scripts/plink/pPlinkRemove.r"
	return pPlinkRemove

@procfactory
def _pPlink2Vcf():
	"""
	@name:
		pPlink2Vcf
	@description:
		Convert plink binaries into VCF file.
	@input:
		`indir:dir`: The input directory containing .bed/.bim/.fam files
	@output:
		`outfile:file`: The output vcf file.
	@args:
		`plink`: The path to plink.
		`gz`   : Whether bgzip the output vcf file. Default: `False`
		`samid`: What to use as sample ID. Default: `both`
			- `both`: use `<FID>_<IID>` as sample id
			- `fid` : use `<FID>` as sample id
			- `iid` : use `<IID>` as sample id
	"""
	pPlink2Vcf             = Proc(desc = "Convert plink binary files to VCF file.")
	pPlink2Vcf.input       = 'indir:dir'
	pPlink2Vcf.output      = 'outfile:file:{{i.indir | bn}}.vcf{% if args.gz %}.gz{% endif %}'
	pPlink2Vcf.args.plink  = params.plink.value
	pPlink2Vcf.args.gz     = False
	pPlink2Vcf.args.samid  = 'both' # fid, iid
	pPlink2Vcf.args.chroms = {"23": "X", "24": "Y", "25": "XY", "26": "M"}
	pPlink2Vcf.lang        = params.python.value
	pPlink2Vcf.script      = "file:scripts/plink/pPlink2Vcf.py"
	return pPlink2Vcf

@procfactory
def _pPlink2GTMat():
	"""
	@name:
		pPlink2GTMat
	@description:
		Convert plink binaries into genotype matrix.
	@input:
		`indir:dir`: The input directory containing .bed/.bim/.fam files
	@output:
		`outfile:file`: The output genotype matrix file.
	@args:
		`plink`: The path to plink.
		`samid`: What to use as sample ID. Default: `both`
			- `both`: use `<FID>_<IID>` as sample id
			- `fid` : use `<FID>` as sample id
			- `iid` : use `<IID>` as sample id
	"""
	pPlink2GTMat             = Proc(desc = "Convert plink binary files to genotype matrix")
	pPlink2GTMat.input       = 'indir:dir'
	pPlink2GTMat.output      = 'outfile:file:{{i.indir | bn}}.gtmat.txt'
	pPlink2GTMat.args.plink  = params.plink.value
	pPlink2GTMat.args.samid  = 'both' # fid, iid
	pPlink2GTMat.args.addchr = True
	pPlink2GTMat.args.snpid  = '{chr}_{pos}_{rs}_{ref}_{alt}' # or raw
	pPlink2GTMat.args.chroms = {"23": "X", "24": "Y", "25": "XY", "26": "M"}
	pPlink2GTMat.args.nors   = "NOVEL"
	pPlink2GTMat.lang        = params.python.value
	pPlink2GTMat.script      = "file:scripts/plink/pPlink2GTMat.py"
	return pPlink2GTMat

@procfactory
def _pPlinkPCA():
	"""
	@name:
		pPlinkPCA
	@description:
		Do PCA on genotype data with PLINK
	@input:
		`indir`: The input directory with .bed/.bim/.fam files
	@output:
		`outfile:file`: The output file of selected PCs, Default: `{{i.indir | fn}}.plinkPCA/{{i.indir | fn}}.pcs.txt`
		`outdir:dir`: The output directory with output file and plots. Default: `{{i.indir | fn}}.plinkPCA`
	@args:
		`plink`: The path to `plink`, Default: `<params.plink>`
		`samid`: Which IDs to report in results, Default: `both`
			- `both`: Both family ID and individual ID connected with `_`
			- `iid`:  Individual ID
			- `fid`:  Family ID
		`nthread`: # threads to use, Default: `False`
			- `False`: Don't put `--threads` in plink command
		`indep`: `indep` used to prune LD SNPs. Default: `[50, 5, .2]`
		`highld`: High LD regions. Default: `<params.highld>`
		`params`: Other parameters for `plink --pca`. Default: `Box(mind = .95)`
		`select`: Select first PCs in the output file. Default: `0.2`
			- `select < 1`: select PCs with contribution greater than `select`
			- `select >=1`: select first `select` PCs
		`plots` : Output plots. Default: `Box(scree = Box(ncp = 20))`
		`devpars`: The parameters for ploting device. Default: `Box(height = 2000, width = 2000, res = 300)`
	"""
	pPlinkPCA          = Proc(desc = "Perform PCA on genotype data and covariates.")
	pPlinkPCA.input    = 'indir:dir'
	pPlinkPCA.output   = [
		'outfile:file:{{i.indir | fn}}.plinkPCA/{{i.indir | fn}}.pcs.txt',
		'outdir:dir:{{i.indir | fn}}.plinkPCA'
	]
	pPlinkPCA.args.plink   = params.plink.value
	pPlinkPCA.args.samid   = 'both' # fid, iid
	pPlinkPCA.args.nthread = False
	pPlinkPCA.args.indep   = [50, 5, .2] # used to prune LD SNPs
	pPlinkPCA.args.highld  = params.highld.value
	pPlinkPCA.args.params  = Box(mind = .95)
	pPlinkPCA.args.select  = .2
	pPlinkPCA.args.plots   = Box(
		scree = Box(ncp = 20),
		# rownames of anno should be consistent with `args.samid`
		pairs = Box(anno = '', ncp = 4, params = Box(upper = Box(continuous = 'density')), ggs = Box(theme = {"axis.text.x": "r:ggplot2::element_text(angle = 60, hjust = 1)" })),
		# more to add
	)
	pPlinkPCA.args.devpars = Box(height = 2000, width = 2000, res = 300)
	pPlinkPCA.envs.rimport = rimport
	pPlinkPCA.lang         = params.Rscript.value
	pPlinkPCA.script       = "file:scripts/plink/pPlinkPCA.r"
	return pPlinkPCA

@procfactory
def _pPlinkSimulate():
	"""
	@name:
		pPlinkSimulate
	"""
	pPlinkSimulate              = Proc(desc = "Simulate a set of SNPs")
	pPlinkSimulate.input        = 'seed'
	pPlinkSimulate.output       = 'outdir:dir:simsnps.{{i.seed if isinstance(i.seed, int) else "noseed"}}.plink'
	pPlinkSimulate.args.plink   = params.plink.value
	pPlinkSimulate.args.ncases  = 1000
	pPlinkSimulate.args.nctrls  = 1000
	pPlinkSimulate.args.nsnps   = 100
	pPlinkSimulate.args.label   = 'SNP'
	pPlinkSimulate.args.dprev   = .01
	pPlinkSimulate.args.minfreq = 0
	pPlinkSimulate.args.maxfreq = 1
	pPlinkSimulate.args.hetodds = 1
	pPlinkSimulate.args.homodds = 1
	pPlinkSimulate.args.params  = Box()
	pPlinkSimulate.lang         = params.python.value
	pPlinkSimulate.script       = "file:scripts/plink/pPlinkSimulate.py"
	return pPlinkSimulate

