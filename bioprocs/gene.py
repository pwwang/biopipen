from pyppl import Proc, Box
from . import params
from .seq import pPromoters
#from .utils import genenorm, write

"""
@name:
	pGenePromoters
@description:
	Alias of `seq.pPromoters`.
"""
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
		- error : report error
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
pGeneNameNorm.args.inopts   = Box(skip = 0, comment = '#', delimit = '\t')
pGeneNameNorm.args.outopts  = Box(delimit = '\t', headDelimit = '\t', headPrefix = '', headTransform = None, head = True, query = False)
pGeneNameNorm.args.genecol  = ''
pGeneNameNorm.args.frm      = 'symbol, alias'
pGeneNameNorm.args.to       = 'symbol'
pGeneNameNorm.args.cachedir = params.cachedir.value
pGeneNameNorm.args.genome   = params.genome.value
pGeneNameNorm.args.cachedir = params.cachedir.value
#pGeneNameNorm.envs.genenorm = genenorm.py
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
pGeneTss                = Proc(desc = 'Get gene TSS in BED format')
pGeneTss.input          = 'infile:file'
pGeneTss.output         = 'outfile:file:{{in.infile | fn}}-tss.bedx'
pGeneTss.errhow         = 'retry'
pGeneTss.args.notfound  = 'skip' # error
pGeneTss.args.genecol   = ''
pGeneTss.args.inopts    = Box(skip = 0, comment = '#', delimit = '\t')
pGeneTss.args.outopts   = Box(delimit = '\t', headDelimit = '\t', headPrefix = '', headTransform = None, head = False, query = False, ftype = 'bed')
pGeneTss.args.frm       = 'symbol, alias'
pGeneTss.args.cachedir  = params.cachedir.value
pGeneTss.args.genome    = params.genome.value
#pGeneTss.envs.genenorm  = genenorm.py
#pGeneTss.envs.writeBedx = write.bedx.py
pGeneTss.lang           = params.python.value
pGeneTss.script         = "file:scripts/gene/pGeneTss.py"

"""
@name:
	pGeneBody
@description:
	Get gene body region in BED format
@input:
	`infile:file`: The input file containing genes
@output:
	`outfile:file`: The gene body region
@args:
	`notfound`: What if a gene is not found when transfer the gene names to gene symbols
		- error: report error
		- skip (default): skip it
	`inmeta`:   The metadata for input file, mainly to indicate where the GENE column is.
	`inopts`:   Input options for reading input file.
		- skip: number of lines to skip. Default: 0
		- comment: the starting string for comment lines. Default: #
		- delimit: The delimit for the input file. Default: '\\t'
	frm: The gene name format in the input file. Default: 'symbol, alias'
	tmpdir: The tmpdir to cache the gene name conversion.
	genome: The genome used to do the conversion.
"""
pGeneBody               = Proc(desc = 'Get gene body in BED format')
pGeneBody.input         = 'infile:file'
pGeneBody.output        = 'outfile:file:{{in.infile | fn}}-body.bedx'
pGeneBody.errhow        = 'retry'
pGeneBody.args.notfound = 'skip' # error
pGeneBody.args.genecol  = ''
pGeneBody.args.inopts   = Box(skip = 0, comment = '#', delimit = '\t')
pGeneBody.args.outopts   = Box(delimit = '\t', headDelimit = '\t', headPrefix = '', headTransform = None, head = False, query = False, ftype = 'bed')
pGeneBody.args.frm      = 'symbol, alias'
pGeneBody.args.cachedir  = params.cachedir.value
pGeneBody.args.genome   = params.genome.value
#pGeneBody.envs.genenorm = genenorm.py
pGeneBody.lang          = params.python.value
pGeneBody.script        = "file:scripts/gene/pGeneBody.py"
