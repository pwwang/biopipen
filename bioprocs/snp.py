from pyppl import Proc, Box
from . import params
from .utils import read, write

"""
@name:
	pSnp2Bedx
@description:
	Find coordinates for SNPs in BEDX format.
@input:
	`snpfile:file`: the snp file, each snp per line
@output:
	`outfile:file`: the result file, columns are:
		- chrom, start(0-based), end, name, score, strand, ref, allele
@args:
	`genome`: default: hg19
	`snpver`: default: snp147
	`notfound`: What to do if the snp is not found. Default: skip
	`inmeta`: The metadata for input file to determine which column is rsID
	`xcols`: The extra columns to extract and output to extra columns in output file.
	`indem`: The input delimit. Default: '\\t'
	`incom`: The input comment. Default: '#'
	`skip`: The lines to skip for input file. Default: 0
@requires:
	[`python-cruzdb`](https://github.com/brentp/cruzdb)
"""
pSnp2Bedx               = Proc(desc = 'Find coordinates for SNPs in BED format.')
pSnp2Bedx.input         = "snpfile:file"
pSnp2Bedx.output        = "outfile:file:{{in.snpfile | fn}}.bed"
pSnp2Bedx.args.genome   = params.genome.value
pSnp2Bedx.args.dbsnpver = params.dbsnpver.value
pSnp2Bedx.errhow        = 'retry'
pSnp2Bedx.args.notfound = 'skip' # error
pSnp2Bedx.args.inopts   = Box(delimit = '\t', skip = 0, comment = '#')
pSnp2Bedx.args.inmeta   = ['SNP']
pSnp2Bedx.args.xcols    = ['refUCSC', 'alleles', 'alleleFreqs', 'alleleFreqCount']
pSnp2Bedx.envs.read     = read
pSnp2Bedx.envs.write    = write
pSnp2Bedx.lang          = params.python.value
pSnp2Bedx.script        = "file:scripts/snp/pSnp2Bedx.py"

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
pSnp2Avinput.output        = "outfile:file:{{in.snpfile | fn}}.avinput"
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

