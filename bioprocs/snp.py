from pyppl import Proc
from . import params

"""
@name:
	pSnp2Bed
@description:
	Find coordinates for SNPs in BED format.
@input:
	`snpfile:file`: the snp file, each snp per line
@output:
	`outfile:file`: the result file, columns are:
		- chrom, start(0-based), end, name, score, strand, ref, allele
@args:
	`genome`: default: hg19
	`snpver`: default: snp147
@requires:
	[`python-cruzdb`](https://github.com/brentp/cruzdb)
"""
pSnp2Bed               = Proc(desc = 'Find coordinates for SNPs in BED format.')
pSnp2Bed.input         = "snpfile:file"
pSnp2Bed.output        = "outfile:file:{{in.snpfile | fn}}.bed"
pSnp2Bed.args.genome   = params.genome.value
pSnp2Bed.args.dbsnpver = params.dbsnpver.value
pSnp2Bed.errhow        = 'retry'
pSnp2Bed.args.notfound = 'skip' # error
pSnp2Bed.args.header   = False
pSnp2Bed.args.skip     = 0
pSnp2Bed.args.comment  = '#'
pSnp2Bed.args.delimit  = '\t'
pSnp2Bed.args.col      = 0
pSnp2Bed.lang          = params.python.value
pSnp2Bed.script        = "file:scripts/snp/pSnp2Bed.py"

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

