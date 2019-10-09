"""Some misc processes"""
from pyppl import Proc, Box
from bioprocs import params, rimport
from . import delefactory, procfactory
from modkit import Modkit
Modkit().delegate(delefactory())

@procfactory
def _pGEP70():
	"""
	@name:
		pGEP70
	@description:
		Calculate GEP70 scores for multiple mylenoma 70-gene-signatures and do survival plot.
		Or add one gene to it to see the survival plot.
		see: https://www.ncbi.nlm.nih.gov/pubmed/17105813
	@input:
		`exprfile:file`: The input file with expression matrix
			- Columns are samples, rows are genes
			- make sure the expression values are log2-scale normalized
		`survfile:file`: The survival data file (header is required).
			- col1: rownames
			- col2: the survival time
			- col3: the status. 0/1 for alive/dead or 1/2 for alive dead
		`gene`: An extra gene to be added to the plot.
			- If not provided, only do the GEP70 plot.
	@output:
		`outdir:file`: The output directory, containing:
			- The survival data file.
			- The GEP70 plot and results
			- The GEP70 + gene plot and results if `i.gene` is provided.
	@args:
		`gep70`: The GEP70 genes.
			- Column 1: genes
			- Column 2: how is the regulated (up or down)
		`inunit`  : The time unit in input file. Default: days
		`outunit` : The output unit for plots. Default: days
		`params`  : The params for `ggsurvplot`. Default: `Box({'risk.table': True, 'conf.int': True, 'font.legend': 13, 'pval': 'Log-rank p = {pval}'})`
			- You may do `ylim.min` to set the min ylim. Or you can set it as 'auto'. Default: 0.
		`ggs`     : Extra ggplot2 elements for main plot. `ggs.table` is for the risk table.
		`devpars` : The device parameters for png. Default: `{res:300, height:2000, width:2000}`
	"""
	pGEP70              = Proc(desc = 'Do GEP70 plot for multiple mylenoma 70-gene-signatures')
	pGEP70.input        = 'exprfile:file, survfile:file, gene'
	pGEP70.output       = 'outdir:dir:{{i.survfile | fn2}}.{{args.name}}{{i.gene}}'
	pGEP70.args.name    = 'GEP70'
	pGEP70.args.gep70   = params.gep70.value
	pGEP70.args.inunit  = 'days' # months, weeks, years
	pGEP70.args.outunit = 'days'
	pGEP70.args.params  = Box({'font.legend': 13, 'pval': 'Log-rank p = {pval}', 'risk.table': True})
	pGEP70.args.devpars = Box(res = 300, height = 2000, width = 2000)
	pGEP70.args.ggs     = Box(table = Box())
	pGEP70.envs.rimport = rimport
	pGEP70.lang         = 'Rscript'
	pGEP70.script       = "file:scripts/misc/pGEP70.r"
	return pGEP70

@procfactory
def _pNCBI():
	"""
	@name:
		pNCBI
	@description:
		The NCBI E-Utils
	@input:
		`term`: The term or the id argument for esearch or efetch
	@output:
		`outfile:file`: The output file
	@args:
		`prog`   : The program to use, esearch (Default) or efetch
		`apikey` : The api key for E-utils
			- Without API key, we can only query 3 time in a second
			- With it, we can do 10.
		`sleep`  : Sleep sometime after job done. Default: `0.15`
			- Because of the limit of # queries/sec, we need to sleep sometime after the job is done
			- At the same time, we also have to limit # jobs to run at the same time. typically: `pNCBI.forks = 10`
		`db`     : The database to query. Default: `pubmed`. Available databases:
			- annotinfo, assembly, biocollections, bioproject, biosample, biosystems, blastdbinfo, books,
			- cdd, clinvar, clone, dbvar, gap, gapplus, gds, gencoll, gene, genome, geoprofiles, grasp, gtr,
			- homologene, ipg, medgen, mesh, ncbisearch, nlmcatalog, nuccore, nucest, nucgss, nucleotide,
			- omim, orgtrack, pcassay, pccompound, pcsubstance, pmc, popset, probe, protein, proteinclusters,
			- pubmed, pubmedhealth, seqannot, snp, sparcle, sra, structure, taxonomy, unigene
		`joiner` : The delimit to use if the field is a list
		`record` : A function to transform the record.
	@requires:
		[python-eutils](https://github.com/biocommons/eutils)
	"""
	pNCBI             = Proc(desc = 'The NCBI E-Utils')
	pNCBI.input       = 'term'
	pNCBI.output      = 'outfile:file:{{i.term | lambda x: __import__("re").sub(r"[^\\w_]", "_", x)[:255]}}.{{args.prog}}.txt'
	pNCBI.errhow      = 'retry'
	pNCBI.args.prog   = 'esearch'
	pNCBI.args.apikey = params.ncbikey.value
	pNCBI.args.sleep  = .15
	# annotinfo, assembly, biocollections, bioproject, biosample, biosystems, blastdbinfo, books,
	# cdd, clinvar, clone, dbvar, gap, gapplus, gds, gencoll, gene, genome, geoprofiles, grasp, gtr,
	# homologene, ipg, medgen, mesh, ncbisearch, nlmcatalog, nuccore, nucest, nucgss, nucleotide,
	# omim, orgtrack, pcassay, pccompound, pcsubstance, pmc, popset, probe, protein, proteinclusters,
	# pubmed, pubmedhealth, seqannot, snp, sparcle, sra, structure, taxonomy, unigene
	pNCBI.args.db     = 'pubmed'
	pNCBI.args.joiner = '|'
	pNCBI.args.record = None
	pNCBI.lang        = params.python.value
	pNCBI.script      = 'file:scripts/misc/pNCBI.py'
	return pNCBI

