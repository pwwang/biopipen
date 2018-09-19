from pyppl import Proc, Box
from . import params

"""
@name:
	pDownload
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
pDownload                 = Proc (desc = 'Download data with gdc-client and a menifest file.')
pDownload.input           = "manifile:file"
pDownload.output          = "outdir:dir:{{i.manifile | fn}}"
pDownload.args.params     = Box({'no-file-md5sum': True})
pDownload.args.nthread    = 1
pDownload.args.token      = None
pDownload.args.gdc_client = params.gdc_client.value
pDownload.lang            = params.python.value
pDownload.script          = "file:scripts/tcga/pDownload.py"

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
	`rsmap`  : The rsid probe mapping file. If not provided, will use the probe id for matrix rownames.
	`fn2sam` : How to convert filename(without extension) to sample name. Default: `None`
"""
pGtFiles2Mat              = Proc(desc = 'Convert genotype files to a matrix.')
pGtFiles2Mat.input        = 'infiles:files'
pGtFiles2Mat.output       = 'outfile:file'
pGtFiles2Mat.args.rsmap   = params.rsmap_gwsnp6.value
pGtFiles2Mat.args.fn2sam  = None
pGtFiles2Mat.lang         = params.python.value
pGtFiles2Mat.script       = "file:scripts/tcga/pGtFiles2Mat.py"

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

"""
@name:
	pConvertExpFiles2Matrix
@description:
	convert TCGA expression files to expression matrix, and convert sample name to submitter id
@input:
	`dir:file`:    the directory containing the samples
	`mdfile:file`: the metadata file
@output:
	`outfile:file`:the output matrix
@requires:
	[python-mygene](https://pypi.python.org/pypi/mygene/3.0.0)
"""
pConvertExpFiles2Matrix = Proc()
pConvertExpFiles2Matrix.input     = "dir:file, mdfile:file"
pConvertExpFiles2Matrix.output    = "outfile:file:expmat-{{#}}.txt"
pConvertExpFiles2Matrix.lang = "python"
pConvertExpFiles2Matrix.script    = """
import json, glob, os
import gzip
import mygene
mg = mygene.MyGeneInfo()

sam_meta = None
sample_ids = {}
with open ("{{mdfile}}") as f:
	sam_meta = json.load (f)

fout = open ("{{outfile}}", 'w')
ext = ''
for sam in sam_meta:
	if not ext:
		ext = os.path.splitext (sam['file_name'])[1]
	sample_ids[sam['file_name']] = sam['associated_entities'][0]['entity_submitter_id'][:15]

gz = ext == '.gz'
exps  = {}
bns   = ['sample']
genes = {}
for samfile in glob.glob (os.path.join("{{dir}}", "*", "*" + ext)):
	bn = os.path.basename (samfile)
	if not sample_ids.has_key(bn): continue
	bns.append (sample_ids[bn])
	f = gzip.open (samfile) if gz else open (samfile)
	if not genes:
		origenes = [line.split("\t")[0].split('.')[0] for line in f if line.strip()]
		f.seek(0)
		origenes = mg.getgenes(origenes, fields='symbol', species='human')
		genes2 = {}
		for ret in origenes:
			if not ret.has_key('symbol'): continue
			genes2[ret['symbol']] = ret['query']
		genes = {v:k for k,v in genes2.iteritems()}
				
	for line in f:
		line = line.strip()
		if not line or "\\t" not in line: continue
		(gene, exp) = line.split("\\t")
		exp = "%.3f" % float(exp)
		gene = gene.split('.')[0]
		if not genes.has_key(gene):continue
		gene = genes[gene]
		if exps.has_key(gene):
			exps[gene].append(exp)
		else:
			exps[gene] = [exp]
	f.close()

fout.write ("\\t".join(bns) + "\\n")
for gene, expvec in exps.iteritems():
	fout.write ("%s\\t%s\\n" % (gene, "\\t".join(expvec)))
fout.close()
"""

"""
@name:
	pConvertMutFiles2Matrix
@description:
	convert TCGA mutation files (vcf.gz) to mut matrix, and convert sample name to submitter id
@input:
	`dir:file`:    the directory containing the samples
	`mdfile:file`: the metadata file
@output:
	`outfile:file`:the output matrix
"""
pConvertMutFiles2Matrix = Proc()
pConvertMutFiles2Matrix.input     = "dir:file, mdfile:file"
pConvertMutFiles2Matrix.output    = "outfile:file:mutmat-{{#}}.txt"
pConvertMutFiles2Matrix.lang = "python"
pConvertMutFiles2Matrix.args      = {'minCount': 2}
pConvertMutFiles2Matrix.script    = """
import json, glob, os
import gzip
import mygene

sam_meta = None
sample_ids = {}
with open ("{{mdfile}}") as f:
	sam_meta = json.load (f)

fout = open ("{{outfile}}", 'w')
ext = ''
for sam in sam_meta:
	if not ext:
		ext = os.path.splitext (sam['file_name'])[1]
	assoc_entities = sam['associated_entities']
	sid1           = assoc_entities[0]['entity_submitter_id'][:15]
	sid2           = assoc_entities[1]['entity_submitter_id'][:15]
	sid            = sid1 if sid1[-2] == '0' else sid2
	sample_ids[sam['file_name']] = sid

gz = ext == '.gz'
muts  = {}
bns   = ['sample']

for samfile in glob.glob (os.path.join("{{dir}}", "*", "*" + ext)):
	bn = os.path.basename (samfile)
	if not sample_ids.has_key(bn): continue
	sid = sample_ids[bn]
	bns.append (sid)
	f = gzip.open (samfile) if gz else open (samfile)
				
	for line in f:
		line = line.strip()
		if line.startswith ('#') or not line: continue
		
		items = line.split("\\t")
		mut   = "%s:%s" % (items[0], items[1])

		if not muts.has_key(mut):
			muts[mut] = {sid: 1}
		else:
			muts[mut][sid] = 1

	f.close()

fout.write ("\\t".join(bns) + "\\n")
for mut, sids in muts.iteritems():
	if len(sids) < {{args.minCount}}: continue
	fout.write (mut)
	for s in bns:
		if sids.has_key(s):
			fout.write ("\\t1")
		else:
			fout.write ("\\t0")
	fout.write ("\\n")
fout.close()
"""

