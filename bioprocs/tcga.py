from pyppl import proc

"""
@name:
	pSample2SubmitterID
@description:
	convert TCGA sample names with submitter id with metadata and sample containing folder
@input:
	`dir:file`:    the directory containing the samples
	`mdfile:file`: the metadata file
@output:
	`outdir`:      the directory containing submitter-id named files
"""
pSample2SubmitterID = proc ()
pSample2SubmitterID.input     = "dir:file, mdfile:file"
pSample2SubmitterID.output    = "outdir:file:out"
pSample2SubmitterID.defaultSh = "python"
pSample2SubmitterID.script    = """
import json, glob, os
if not os.path.exists ("{{outdir}}"):
	os.makedirs("{{outdir}}")
sam_meta = None
sample_ids = {}
with open ("{{mdfile}}") as f:
	sam_meta = json.load (f)

ext = ''
for sam in sam_meta:
	if not ext:
		ext = os.path.splitext (sam['file_name'])[1]
	sample_ids[sam['file_name']] = sam['associated_entities'][0]['entity_submitter_id'][:15]

for samfile in glob.glob (os.path.join("{{dir}}", "*", "*" + ext)):
	bn = os.path.basename (samfile)
	newfile = os.path.join ("{{outdir}}", sample_ids[bn] + ext)
	if os.path.exists (newfile):
		os.remove(newfile)
	os.symlink (samfile, newfile)
"""