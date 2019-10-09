"""TCGA data analysis"""
from pyppl import Proc, Box
from . import params, rimport
from .utils import fs2name
from . import delefactory, procfactory
from modkit import Modkit
Modkit().delegate(delefactory())

@procfactory
def _pTCGADownload():
	"""
	@name:
		pTCGADownload
	@description:
		Download TCGA use `gdc-client` and a manifest file
	@input:
		`manifile:file`: the manifest file
	@output:
		`outdir:file`:   the directory containing downloaded file
	@args:
		`params`    : other params for `gdc-client download`, default: `{'no-file-md5sum': True}`
		`gdc_client`: the executable file of `gdc-client`,    default: "gdc-client"
		`nthread`   : Number of threads to use. Default     : `1`
		`token`     : The token file if needed.
	"""
	pTCGADownload                 = Proc (desc = 'Download data with gdc-client and a menifest file.')
	pTCGADownload.input           = "manifile:file"
	pTCGADownload.output          = "outdir:dir:{{i.manifile | fn}}"
	pTCGADownload.echo            = Box(jobs = 0, type = 'stderr')
	pTCGADownload.args.params     = Box({'no-file-md5sum': True})
	pTCGADownload.args.nthread    = 1
	pTCGADownload.args.token      = None
	pTCGADownload.args.gdc_client = params.gdc_client.value
	pTCGADownload.lang            = params.python.value
	pTCGADownload.script          = "file:scripts/tcga/pTCGADownload.py"
	return pTCGADownload

@procfactory
def _pSample2SubmitterID():
	"""
	@name:
		pSample2SubmitterID
	@description:
		convert TCGA sample names with submitter id with metadata and sample containing folder
	@input:
		`indir:file`:    the directory containing the samples
		`mdfile:file`: the metadata file
	@output:
		`outdir:file`: the directory containing submitter-id named files
	@args:
		`method`: How the deal with the files. Default: `symlink`
			- We can also do `copy`
		`nthread`: Number threads to use. Default: `1`
	"""
	pSample2SubmitterID              = Proc(desc = 'convert TCGA sample names with submitter id with metadata and sample containing folder')
	pSample2SubmitterID.input        = "indir:file, mdfile:file"
	pSample2SubmitterID.output       = "outdir:dir:{{i.indir | fn}}"
	pSample2SubmitterID.args.method  = 'symlink' # or copy
	pSample2SubmitterID.args.nthread = 1
	pSample2SubmitterID.lang         = params.python.value
	pSample2SubmitterID.script       = "file:scripts/tcga/pSample2SubmitterID.py"
	return pSample2SubmitterID

@procfactory
def _pGtFiles2Mat():
	"""
	@name:
		pGtFiles2Mat
	@description:
		Convert TCGA genotype files to a matrix.
	@input:
		`infiles:files`: The input genotypes files
	@output:
		`outfile:file`: The output matrix file
	@args:
		`rsmap`    : The rsid probe mapping file. If not provided, will use the probe id for matrix rownames.
		`fn2sample`: How to convert filename(without extension) to sample name. Default: `None`
		`confcut`  : The confidence cutoff. Genotype will be NA for lower confidence snps. Default: `0.05`
	"""
	pGtFiles2Mat                = Proc(desc = 'Convert genotype files to a matrix.')
	pGtFiles2Mat.input          = 'infiles:files'
	pGtFiles2Mat.output         = 'outfile:file:{{i.infiles | fs2name}}.gt.txt'
	pGtFiles2Mat.args.rsmap     = params.rsmap_gwsnp6.value
	pGtFiles2Mat.args.confcut   = 0.05
	pGtFiles2Mat.args.fn2sample = "lambda fn: fn.split('.')[0]"
	pGtFiles2Mat.envs.fs2name   = fs2name
	pGtFiles2Mat.lang           = params.python.value
	pGtFiles2Mat.script         = "file:scripts/tcga/pGtFiles2Mat.py"
	return pGtFiles2Mat

@procfactory
def _pClinic2Survival():
	"""
	@name:
		pClinic2Survival
	@description:
		Convert TCGA clinic data to survival data
		The clinic data should be downloaded as "BCR Biotab" format
	@input:
		`infile:file`: The clinic data file downloaded from TCGA
	@output:
		`outfile:file`: The output file
		`covfile:file`: The covariate file
	@args:
		`cols`: The column names:
			- `time_lastfollow`: The column names of last follow up. Default: `['days_to_last_followup']`
			- `time_death`: The column names of time to death. Default: `['days_to_death']`
			- `status`: The columns of vital status. Default: `['vital_status']`
			- `age`: The columns of days to birth. Default: `['days_to_birth']`
		`covs`: The covariates to output. Default:
			- `gender`, `race`, `ethnicity`, `age`
		`mat`: An expression or genotype matrix with samples as column names, used to get sample names for patient instead of short ones. Default: `None`
	"""
	pClinic2Survival           = Proc(desc = 'Convert TCGA clinic data to survival data.')
	pClinic2Survival.input     = 'infile:file'
	pClinic2Survival.output    = 'outfile:file:{{i.infile | stem}}.survdata.txt, covfile:file:{{i.infile | stem}}.survcov.txt'
	pClinic2Survival.args.cols = Box(
		time_lastfollow = ['days_to_last_followup'],
		time_death      = ['days_to_death'],
		status          = ['vital_status'],
		age             = ['days_to_birth'],
		patient         = ['bcr_patient_barcode']
	)
	pClinic2Survival.args.covs = ['gender', 'race', 'ethnicity', 'age']
	pClinic2Survival.args.mat  = None
	pClinic2Survival.lang      = params.python.value
	pClinic2Survival.script    = "file:scripts/tcga/pClinic2Survival.py"
	return pClinic2Survival

@procfactory
def _pClinic2PlinkMeta():
	"""
	@name:
		pClinic2PlinkMeta
	@description:
		Convert TCGA clinic data to plink meta file.
		Used with a `GTMat` file to convert to plink binary file (see vcfnext.pGTMat2Plink).
	@input:
		`infile:file`: The input TCGA patiant clinic file.
	@output:
		`outfile:file`: The meta file.
	@args:
		`suffix`: The suffix added to the barcode. Default: `''`
			- By default the barcode is like `TCGA-55-8087`, but in analysis, we want to be specific to samples
			- For example, for solid tumor, it'll be `TCGA-55-8087-01`, then suffix `-01` will be added
			- You can add multiple suffices by `-01, -10` or `['-01', '-10']`
	"""
	pClinic2PlinkMeta             = Proc(desc = 'Convert TCGA clinic data to plink meta file.')
	pClinic2PlinkMeta.input       = 'infile:file'
	pClinic2PlinkMeta.output      = 'outfile:file:{{i.infile | fn}}.meta.txt'
	pClinic2PlinkMeta.args.suffix = ''
	pClinic2PlinkMeta.lang        = params.python.value
	pClinic2PlinkMeta.script      = "file:scripts/tcga/pClinic2PlinkMeta.py"
	return pClinic2PlinkMeta

@procfactory
def _pClinic2Cov():
	"""
	@name:
		pClinic2Cov
	"""
	pClinic2Cov              = Proc(desc = 'Extract information from clinic data as covariance matrix')
	pClinic2Cov.input        = 'infile:file'
	pClinic2Cov.output       = 'outfile:file:{{i.infile | fn}}.cov.txt'
	pClinic2Cov.args.covs    = ['gender', 'race', 'pathologic_stage', 'age_at_initial_pathologic_diagnosis']
	pClinic2Cov.args.asnum   = True
	pClinic2Cov.args.suffix  = ''
	pClinic2Cov.envs.rimport = rimport
	pClinic2Cov.lang         = params.Rscript.value
	pClinic2Cov.script       = "file:scripts/tcga/pClinic2Cov.r"
	return pClinic2Cov

