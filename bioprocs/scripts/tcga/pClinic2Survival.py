from pyppl import Box
from bioprocs.utils.tsvio import TsvReader, TsvWriter, TsvRecord

infile  = {{i.infile | quote}}
outfile = {{o.outfile | quote}}
covfile = {{o.covfile | quote}}
cols    = {{args.cols}}
covs    = {{args.covs}}
mat     = {{args.mat | repr}}

reader = TsvReader(infile, ftype = 'head', skip = 1)
next(reader) # skip CDE_ID line

writer    = TsvWriter(outfile)
covwriter = TsvWriter(covfile)
writer.meta.add('patient', 'time', 'status')
writer.writeHead()
covwriter.meta.add('patient', *covs)
covwriter.writeHead()

realcols = {}
for key, val in cols.items():
	for v in val:
		if v in reader.meta:
			realcols[key] = v
			break

# get the sample names
samples = []
if mat:
	with open(mat) as f:
		for line in f:
			if line.startswith('#'): continue
			samples = line.strip().split('\t')
			break

def getSampleName(patient):
	ret = []
	for sample in samples:
		if sample.startswith(patient):
			ret.append(sample)
	return ret or [patient]

for r in reader:
	patient  = r[realcols.get('patient', 'patient')]
	patients = getSampleName(patient)
	status   = r[realcols.get('status', 'status')]

	if status == 'Alive':
		status = 0
	elif status == 'Dead':
		status = 1
	else:
		continue
	if status:
		time = r[realcols.get('time_death', 'time_death')]
	else:
		time = r[realcols.get('time_lastfollow', 'time_lastfollow')]
	if time == '[Not Available]' or time == '[Completed]' or time == '0':
		continue
	
	for pat in patients:
		rout         = TsvRecord()
		rout.patient = pat
		rout.time    = time
		rout.status  = status
		writer.write(rout)

		rcov         = TsvRecord()
		rcov.patient = pat
		for cov in covs:
			rcov[cov] = r[realcols.get(cov, cov)]
			if cov == 'age':
				if rcov[cov].startswith('-'):
					rcov[cov] = '{:.2f}'.format(abs(int(rcov[cov]))/365.25)
		covwriter.write(rcov)