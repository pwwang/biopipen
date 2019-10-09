"""Processes to process files with SNPs"""
from pyppl import Proc, Box
from . import params
from . import delefactory, procfactory
from modkit import Modkit
Modkit().delegate(delefactory())

@procfactory
def _pRs2Bed():
	"""
	@name:
		pRs2Bed
	@description:
		Find coordinates for SNPs in BED format.
	@input:
		`snpfile:file`: the snp file, each snp per line
	@output:
		`outfile:file`: the result file, could be a 3 or 6-col bed file.
	@args:
		`tool`    : The tool used to retrieve the coordinates. Default: `cruzdb`
			- Or retrieve from a dbsnp VCF file (dbsnp or local)
		`dbsnp`   : The dbsnp vcf file
		`notfound`: What to do if the snp is not found. Default: skip
		`inopts`  : The input options for input file. Default: `Box(delimit = '\t', skip = 0, comment = '#')`
		`snpcol`  : The column where the snp is. Could be index (0-based) or the column name with `args.inopts.cnames = True`
		`ncol`    : How many columns to output. Possible values: 6, 8, 9. Default: `8`
			- `6`: The BED6 format
			- `8`: BED6 + `ref` + `alt`
			- `9`: BED6 + `ref` + `alt` + `allele frequences from 1000 genome project`
			- `args.ncol = 9` only avaliable when `args.tool = "dbsnp"`
		`vcftools`: The path to vcftools to extract snp from `args.dbsnp`
		`genome`  : The genome to locate the right database from `cruzdb`
		`dbsnpver`: The version of dbsnp from `cruzdb`
		`sortby`  : Sort the output file by coordinates (coord) or name. Default: `coord`
			- `False`: don't sort.
	"""
	pRs2Bed               = Proc(desc = 'Find coordinates for SNPs in BED format.')
	pRs2Bed.input         = "snpfile:file"
	pRs2Bed.output        = "outfile:file:{{i.snpfile | fn}}.bed"
	pRs2Bed.args.tool     = 'cruzdb' # or local/dbsnp
	pRs2Bed.args.notfound = 'skip' # error
	pRs2Bed.args.inopts   = Box(delimit = '\t', skip = 0, comment = '#')
	pRs2Bed.args.ncol     = 6
	pRs2Bed.args.snpcol   = ''
	pRs2Bed.args.dbsnp    = params.dbsnp_all.value
	pRs2Bed.args.vcftools = params.vcftools.value
	pRs2Bed.args.sortby   = 'coord' # name/False
	pRs2Bed.args.genome   = params.genome.value
	pRs2Bed.args.dbsnpver = params.dbsnpver.value
	pRs2Bed.lang          = params.python.value
	pRs2Bed.script        = "file:scripts/snp/pRs2Bed.py"
	return pRs2Bed

@procfactory
def _pCoord2SnpBedx():
	"""
	@name:
		pCoord2SnpBedx
	"""
	pCoord2SnpBedx               = Proc(desc = 'Find snps with coordinates.')
	pCoord2SnpBedx.input         = 'infile:file'
	pCoord2SnpBedx.output        = 'outfile:file:{{i.infile | fn}}-snps.bed'
	pCoord2SnpBedx.args.genome   = params.genome.value
	pCoord2SnpBedx.args.dbsnpver = params.dbsnpver.value
	pCoord2SnpBedx.errhow        = 'retry'
	pCoord2SnpBedx.args.inmeta   = ['CHR', 'START', 'END']
	pCoord2SnpBedx.args.inopts   = Box(delimit = '\t', skip = 0, comment = '#')
	pCoord2SnpBedx.args.xcols    = ['refUCSC', 'alleles', 'alleleFreqs', 'alleleFreqCount']
	#pCoord2SnpBedx.envs.read     = read
	#pCoord2SnpBedx.envs.write    = write
	pCoord2SnpBedx.lang          = params.python.value
	pCoord2SnpBedx.script        = "file:scripts/snp/pCoord2SnpBedx.py"
	return pCoord2SnpBedx

@procfactory
def _pSnp2Avinput():
	"""
	@name:
		pSnp2Avinput
	@description:
		Convert SNP list to avinput to ANNOVAR.
	@input:
		`snpfile:file`: the snp file, each snp per line
	@output:
		`outfile:file`: the result avinput file
	@args:
		`genome`: default: hg19
		`snpver`: default: snp147
	@requires:
		[`python-cruzdb`](https://github.com/brentp/cruzdb)
	"""
	pSnp2Avinput               = Proc (desc = 'Convert SNP list to avinput to ANNOVAR.')
	pSnp2Avinput.input         = "snpfile:file"
	pSnp2Avinput.output        = "outfile:file:{{i.snpfile | fn}}.avinput"
	pSnp2Avinput.lang          = params.python.value
	pSnp2Avinput.args.genome   = params.genome.value
	pSnp2Avinput.args.dbsnpver = params.dbsnpver.value
	pSnp2Avinput.errhow        = 'retry'
	pSnp2Avinput.args.notfound = 'skip' # error
	pSnp2Avinput.args.header   = False
	pSnp2Avinput.args.skip     = 0
	pSnp2Avinput.args.comment  = '#'
	pSnp2Avinput.args.delimit  = '\t'
	pSnp2Avinput.args.col      = 0
	pSnp2Avinput.script        = "file:scripts/snp/pSnp2Avinput.py"
	return pSnp2Avinput

