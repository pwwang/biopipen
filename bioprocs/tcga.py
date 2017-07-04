from pyppl import proc

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
	`params`:        other params for `gdc-client`, default: "--no-related-files --no-file-md5sum -n 20"
	`bin-gdc`:       the executable file of `gdc-client`, default: "gdc-client"
"""
pDownload = proc ()
pDownload.input     = "manifile:file"
pDownload.output    = "outdir:dir:{{manifile | fn}}"
pDownload.args      = {"params": " --no-file-md5sum -n 20 ", "bin-gdc": "gdc-client"}
pDownload.script    = """
{{args.bin-gdc}} download -m "{{manifile}}" -d "{{outdir}}" {{args.params}}
"""

"""
@name:
	pSample2SubmitterID
@description:
	convert TCGA sample names with submitter id with metadata and sample containing folder
@input:
	`dir:file`:    the directory containing the samples
	`mdfile:file`: the metadata file
@output:
	`outdir:file`: the directory containing submitter-id named files
"""
pSample2SubmitterID = proc ()
pSample2SubmitterID.input     = "dir:file, mdfile:file"
pSample2SubmitterID.output    = "outdir:dir:out-{{#}}"
pSample2SubmitterID.defaultSh = "python"
pSample2SubmitterID.script    = """
import json, glob, os

sam_meta = None
sample_ids = {}
with open ("{{mdfile}}") as f:
	sam_meta = json.load (f)

ext = ''
for sam in sam_meta:
	if not ext:
		ext = os.path.splitext (sam['file_name'])[1]
	sample_ids[sam['file_name']] = sam['associated_entities'][0]['entity_submitter_id'][:15]

for samfile in glob.glob (os.path.join(os.path.abspath("{{dir}}"), "*" + ext)):
	bn = os.path.basename (samfile)
	if not sample_ids.has_key (bn): continue
	newfile = os.path.join ("{{outdir}}", sample_ids[bn] + ext)
	if os.path.exists (newfile):
		os.remove(newfile)
	os.symlink (samfile, newfile)
"""

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
pConvertExpFiles2Matrix = proc()
pConvertExpFiles2Matrix.input     = "dir:file, mdfile:file"
pConvertExpFiles2Matrix.output    = "outfile:file:expmat-{{#}}.txt"
pConvertExpFiles2Matrix.defaultSh = "python"
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
pConvertMutFiles2Matrix = proc()
pConvertMutFiles2Matrix.input     = "dir:file, mdfile:file"
pConvertMutFiles2Matrix.output    = "outfile:file:mutmat-{{#}}.txt"
pConvertMutFiles2Matrix.defaultSh = "python"
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

