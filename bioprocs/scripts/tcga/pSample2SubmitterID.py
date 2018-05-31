import json, glob, os
from bioprocs.utils.parallel import Parallel
from threading import Lock

# load varialbes
indir   = {{in.indir | quote}}
mdfile  = {{in.mdfile | quote}}
outdir  = {{out.outdir | quote}}
nthread = {{args.nthread | repr}}

sam_meta = None
sample_ids = {}
with open (mdfile) as f:
	sam_meta = json.load (f)

ext = ''
for sam in sam_meta:
	if not ext:
		ext = os.path.splitext (sam['file_name'])[1]
	sample_ids[sam['file_name']] = sam['associated_entities'][0]['entity_submitter_id'][:15]

samfiles = glob.glob (os.path.join(os.path.abspath(indir), "*" + ext))
# or direct dir from TCGA download
samfiles += glob.glob (os.path.join(os.path.abspath(indir), "*", "*" + ext))

lock = Lock()
def single(samfile):
	bn = os.path.basename (samfile)
	if not bn in sample_ids: return
	newfile = os.path.join (outdir, sample_ids[bn] + ext)
	with lock:
		if os.path.exists (newfile):
			os.remove(newfile)
		os.symlink (samfile, newfile)

p = Parallel(nthread = nthread, backend = 'threading', raiseExc = True)
p.run(single, [(samfile,) for samfile in samfiles])