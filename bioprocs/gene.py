from pyppl import Proc
from . import params
from .seq import pPromoters
from .utils import genenorm

pGenePromoters = pPromoters.copy()

"""
@name:
	pGeneNameNorm
@description:
	Normalize gene names using MyGeneinfo.
@input:
	`infile:file`: The input file
@output:
	`outfile:file`: The output file
@args:
	`notfound`: What if a symbol is not found. Default: ignore
		- skip  : skip the record(don't write it to output file)
		- ignore: use the original name;
		- error : report erro
	`col`: the column index containing the gene names
	`from`: the original format. Default: 'symbol, alias'
	`to`: the output gene name format. Default: 'symbol'
	`genome`: the genome. Default: 'hg19'
"""
pGeneNameNorm               = Proc(desc = 'Normalize gene names using MyGeneinfo.')
pGeneNameNorm.input         = 'infile:file'
pGeneNameNorm.output        = 'outfile:file:{{in.infile | bn}}'
pGeneNameNorm.errhow        = 'retry'
pGeneNameNorm.args.notfound = 'ignore'
pGeneNameNorm.args.header   = False
pGeneNameNorm.args.skip     = 0
pGeneNameNorm.args.comment  = '#'
pGeneNameNorm.args.delimit  = '\t'
pGeneNameNorm.args.col      = 0
pGeneNameNorm.args.frm      = 'symbol, alias'
pGeneNameNorm.args.to       = 'symbol'
pGeneNameNorm.args.tmpdir   = params.tmpdir.value
pGeneNameNorm.args.genome   = params.genome.value
pGeneNameNorm.envs.genenorm = genenorm.py
pGeneNameNorm.lang          = params.python.value
pGeneNameNorm.script        = "file:scripts/gene/pGeneNameNorm.py"

"""
@name:
	pGeneTss
@description:
	Get gene TSS in BEd format.
@input:
	`infile:file`: The input file containing genes
@output:
	`outfile:file`: The output BED file
@args:
	`notfound`: What if the gene is not found. Default: skip.
		- error: report error
	`header`: Whether the input file contains header. Default: False
	`skip`: Skip N lines of input file. Default: 0 
		- This has highest priority of header and comment
	`comment`: The comment line start sign. Default: #
	`delimit`: The delimit of input file if it has multiple column. Default: `\\t`
	`col`: The column index contains the genes. Default: 0
	`frm`: The format of the genes. Default: `symbol, alias`
	`tmpdir`: The tmpdir used to store mygene cache files. 
	`genome`: In which genome to fetch the coordinates. Default: hg19
"""
pGeneTss               = Proc(desc = 'Get gene TSS in BED format')
pGeneTss.input         = 'infile:file'
pGeneTss.output        = 'outfile:file:{{in.infile | fn}}-tss.bed'
pGeneTss.errhow        = 'retry'
pGeneTss.args.notfound = 'skip' # error
pGeneTss.args.header   = False
pGeneTss.args.skip     = 0
pGeneTss.args.comment  = '#'
pGeneTss.args.delimit  = '\t'
pGeneTss.args.col      = 0
pGeneTss.args.frm      = 'symbol, alias'
pGeneTss.args.tmpdir   = params.tmpdir.value
pGeneTss.args.genome   = params.genome.value
pGeneTss.envs.genenorm = genenorm.py
pGeneTss.lang          = params.python.value
pGeneTss.script        = "file:scripts/gene/pGeneTss.py"
