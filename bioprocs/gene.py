"""Gene related processes"""
from pyppl import Proc
from diot import Diot
from . import params, proc_factory

pGenePromoters = pPromoters.copy()

pGeneNameNorm = proc_factory(
	desc = 'Normalize gene names using MyGeneinfo.',
	config = Diot(annotate = """
	@name:
		pGeneNameNorm
	@description:
		Normalize gene names using MyGeneinfo.
	@input:
		`infile:file`: The input file
	@output:
		`outfile:file`: The output file
	@args:
		`inopts`: options for reading input file.
		`outopts`: options for writing output file.
			- `query` : Output the original query column? Default: `False`
			- `cnames`: Output headers? Default: `True`
		`notfound`: What if a symbol is not found. Default: ignore
			- skip  : skip the record(don't write it to output file)
			- ignore: use the original name;
			- error : report error
		`genecol` : the column index containing the gene names
		`frm`     : the original format. Default: 'symbol, alias'
		`to`      : the output gene name format. Default: 'symbol'
		`genome`  : the genome. Default: 'hg19'
		`cachedir`: The cache directory
	"""))
pGeneNameNorm.input         = 'infile:file'
pGeneNameNorm.output        = 'outfile:file:{{i.infile | bn}}'
pGeneNameNorm.errhow        = 'retry'
pGeneNameNorm.args.notfound = 'ignore'
pGeneNameNorm.args.inopts   = Diot(skip = 0, comment = '#', delimit = '\t')
pGeneNameNorm.args.outopts  = Diot(delimit = '\t', cnames = True, query = False)
pGeneNameNorm.args.genecol  = ''
pGeneNameNorm.args.frm      = 'symbol, alias'
pGeneNameNorm.args.to       = 'symbol'
pGeneNameNorm.args.genome   = params.genome.value
pGeneNameNorm.args.cachedir = params.cachedir.value
pGeneNameNorm.lang          = params.python.value

pIPI = proc_factory(
	desc = 'Convert gene symbol to IPI protein accession and vice versa.',
	config = Diot(annotate = """
	@name:
		pIPI
	@description:
		Convert gene symbol to IPI protein accession and vice versa.
		One gene symbol could map to multiple IPIs, which will be separated by pipe (|)
	@input:
		`infile:file` : The input file
	@output:
		`outfile:file`: The output file
	@args:
		`notfound`: What if a record is not found: Default: `ignore`
			- `skip`  : skip the record(don't write it to output file)
			- `ignore`: use the original name;
			- `error` : report error
		`genecol`: The column index containing the gene/protein record
		`ipidb`: The IPI xref database (see http://ftp.ebi.ac.uk/pub/databases/IPI/last_release/current/).
		`fromipi`: Whether the input is IPI or genes
		`inopts`: The options for input file
		`outopts`: The options for output file
	"""))
pIPI.input         = 'infile:file'
pIPI.output        = 'outfile:file:{{i.infile | bn}}'
pIPI.errhow        = 'retry'
pIPI.args.notfound = 'ignore'
pIPI.args.inopts   = Diot(skip = 0, comment = '#', delimit = '\t')
pIPI.args.outopts  = Diot(delimit = '\t', headDelimit = '\t', headPrefix = '', headTransform = None, head = True, query = False)
pIPI.args.genecol  = None
pIPI.args.fromipi  = True
pIPI.args.ipidb    = params.ipidb.value
pIPI.lang          = params.python.value

pGeneTss = proc_factory(
	desc = 'Get gene TSS in BED format',
	config = Diot(annotate = """
	@name:
		pGeneTss
	@description:
		Get gene TSS in BED format.
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
	"""))
pGeneTss.input          = 'infile:file'
pGeneTss.output         = 'outfile:file:{{i.infile | fn}}-tss.bedx'
pGeneTss.errhow         = 'retry'
pGeneTss.args.notfound  = 'skip' # error
pGeneTss.args.genecol   = ''
pGeneTss.args.inopts    = Diot(skip = 0, comment = '#', delimit = '\t')
pGeneTss.args.outopts   = Diot(delimit = '\t', headDelimit = '\t', headPrefix = '', headTransform = None, head = False, query = False, ftype = 'bed')
pGeneTss.args.frm       = 'symbol, alias'
pGeneTss.args.cachedir  = params.cachedir.value
pGeneTss.args.genome    = params.genome.value
#pGeneTss.envs.genenorm  = genenorm.py
#pGeneTss.envs.writeBedx = write.bedx.py
pGeneTss.lang           = params.python.value

pGeneBody = proc_factory(
	desc = 'Get gene body in BED format',
	config = Diot(annotate = """
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
	"""))
pGeneBody.input         = 'infile:file'
pGeneBody.output        = 'outfile:file:{{i.infile | fn}}-body.bed'
pGeneBody.args.inopts   = Diot(cnames = False)
pGeneBody.args.notfound = 'skip' # error
pGeneBody.args.genecol  = ''
pGeneBody.args.refgene  = params.refgene.value
pGeneBody.lang          = params.python.value
