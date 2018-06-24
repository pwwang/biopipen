from pyppl import Proc, Box
from bioprocs import params, rimport

"""
@name:
	pGEP70
@description:
	Calculate GEP70 scores for multiple mylenoma 70-gene-signatures
@input:
	`infile:file`: The input file with expression matrix
		- Columns are samples, rows are genes
@output:
	`outfile:file`: The output files with gep70 scores for each sample.
		- Samples become rows, just one column is in the file.
@args:
	`inopts`: The input options.
		- `cnames`: Whether the input file has column names. Default: `True`
	`gep70`: The GEP70 genes. 
		- Column 1: up-regulated genes (51)
		- Column 2: down-regulated genes (19)
"""
pGEP70              = Proc(desc = 'Calculate GEP70 scores for multiple mylenoma 70-gene-signatures')
pGEP70.input        = 'infile:file' # make sure the expression values are log2-scale normalized
pGEP70.output       = 'outfile:file:{{in.infile | fn2}}.gep70.txt'
pGEP70.args.gep70   = params.gep70.value
pGEP70.args.inopts  = Box(cnames = True)
pGEP70.args.lang    = 'Rscript'
pGEP70.envs.rimport = rimport
pGEP70.script       = "file:scripts/misc/pGEP70.r"

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
pNCBI.output      = 'outfile:file:{{in.term | lambda x: __import__("re").sub(r"[^\\w_]", "", x)[:20]}}.{{args.prog}}.txt'
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
