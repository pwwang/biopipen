"""A set of resources to download via API or URL links"""
from pyppl import Proc, Box
from . import params
#from .utils import txt, download, runcmd
from . import delefactory, procfactory
from modkit import Modkit
Modkit().delegate(delefactory())

@procfactory
def _pTxt():
	"""
	@name:
		pTxt
	@description:
		Download CSV format files.
	@input:
		`in`: The name of the resource
	@output:
		`outfile:file`: The output file
	@args:
		`cols`:      Select the columns to keep. Default: '' (all cols)
		`rowfilter`: Filter rows. For example, to filter out rows not start with 'Chr':
			- `"lambda x: not x[0].startswith('Chr')"`
			- Note that rowfilter applied before cols filter.
		`urls`:      Available resources and their urls.
		`gz`:        Whether to gzip the output file.
	@requires:
		[`curl`](https://en.wikipedia.org/wiki/CURL)
	"""
	pTxt                      = Proc(desc = 'Download CSV format files.')
	pTxt.input                = "in"
	pTxt.output               = "outfile:file:{{i.in}}.txt{{args.gz | lambda x: '.gz' if x else ''}}"
	pTxt.args.gz              = False
	pTxt.args.delimit         = "\t"
	pTxt.args.skip            = 0
	pTxt.args.cols            = ''
	pTxt.args.header          = True
	pTxt.args.rowfilter       = ''
	pTxt.args.transform       = ''
	pTxt.args.username        = ''
	pTxt.args.password        = ''
	pTxt.args.curl            = params.curl.value
	# pTxt.tplenvs.txtFilter    = txt.filter.py
	# pTxt.tplenvs.txtTransform = txt.transform.py
	# pTxt.tplenvs.downloadCurl = download.curl.py
	# pTxt.tplenvs.runcmd       = runcmd.py
	pTxt.args.urls            = Box({
		'drugbank-target-all': 'https://www.drugbank.ca/releases/5-0-7/downloads/target-all-uniprot-links',
		'drugbank-target-approved': 'https://www.drugbank.ca/releases/5-0-7/downloads/target-approved-uniprot-links',
		'ccle-sample-info': 'https://data.broadinstitute.org/ccle_legacy_data/cell_line_annotations/CCLE_sample_info_file_2012-10-18.txt',
		'ccle-rseq-rpkm': 'https://data.broadinstitute.org/ccle/CCLE_RNAseq_081117.rpkm.gct',
		'ccle-rseq-reads': 'https://data.broadinstitute.org/ccle/CCLE_RNAseq_081117.reads.gct',
		'KEGG_2016_gmt': 'http://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=KEGG_2016',
		'GO_Molecular_Function_2017_gmt': 'http://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=GO_Molecular_Function_2017',
		'GO_Cellular_Component_2017_gmt': 'http://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=GO_Cellular_Component_2017',
		'GO_Biological_Process_2017_gmt': 'http://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=GO_Biological_Process_2017',
		'TargetScan_microRNA': 'http://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=TargetScan_microRNA',
		'TRANSFAC_and_JASPAR_PWMs': 'http://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=TRANSFAC_and_JASPAR_PWMs',
	})
	pTxt.lang   = params.python.value
	pTxt.script = "file:scripts/resource/pTxt.py"
	return pTxt

@procfactory
def _pGtf():
	"""
	@name:
		pGtf
	"""
	pGtf                      = Proc(desc = 'Download GTF files.')
	pGtf.input                = "in"
	pGtf.output               = "outfile:file:{{in}}.gtf{{args.gz | lambda x: '.gz' if x else ''}}"
	pGtf.args.gz              = False
	pGtf.args.curl            = 'curl'
	pGtf.args.username        = ''
	pGtf.args.password        = ''
	pGtf.args.genepredtogtf   = 'genePredToGtf'
	# pGtf.tplenvs.runcmd       = runcmd.py
	# pGtf.tplenvs.downloadCurl = download.curl.py
	pGtf.args.urls            = {
		'hg19-refgene': 'http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz',
		'hg19-knowngene': 'http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/knownGene.txt.gz',
		'hg38-refgene': 'http://hgdownload.cse.ucsc.edu/goldenPath/hg38database/refGene.txt.gz',
		'hg38-knowngene': 'http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/knownGene.txt.gz',
	}
	pGtf.lang   = params.python.value
	pGtf.script = """
	import os, shutil
	from subprocess import check_output
	{{downloadCurl}}
	{{runcmd}}
	url = {{args.urls | json}}["{{in}}"]
	tmpdir = "{{job.outdir}}/tmp"
	if not os.path.exists(tmpdir):
		os.makedirs(tmpdir)

	downfile = os.path.join(tmpdir, 'downloaded')
	downloadCurl(url, downfile, {{args.username | quote}}, {{args.password | quote}}, {{args.curl | quote}})
	outfile = "{{outfile}}"[:-3] if {{args.gz | bool}} else "{{outfile}}"
	output  = check_output(['file', downfile])
	if 'gzip' in output:
		ugfile = downfile + '.ungz'
		with open(ugfile, 'w') as f:
			f.write(check_output(['gunzip', downfile, '-c']))
		downfile = ugfile
	elif 'Zip' in output:
		zipdir = os.path.join(tmpdir, '_unzipped')
		import zipfile, glob
		zipref = zipfile.ZipFile(downfile, 'r')
		zipref.extractall(zipdir)
		zipref.close()
		downfile = glob.glob(os.path.join(zipdir, '*'))[0]
	cutfile = downfile + '.cutf2'
	with open(downfile) as fin, open(cutfile, 'w') as fout:
		for line in fin:
			if not line.strip(): continue
			fout.write("\\t".join(line.split("\\t")[1:]))
	runcmd ('{{args.genepredtogtf}} file "%s" "%s" -source="{{in}}"' % (cutfile, outfile))
	if {{args.gz | bool}}:
		runcmd ('gz "%s"' % outfile)
	shutil.rmtree(tmpdir)
	"""
	return pGtf

