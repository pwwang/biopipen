from os import path
from diot import Diot
from bioprocs.utils import shell2 as shell
from bioprocs.utils.tsvio2 import TsvReader, TsvWriter

fq1      = {{ i.fqfile1 | quote}}
fq2      = {{ i.fqfile2 | quote}}
outfile  = {{ o.outfile | quote}}
outdir   = {{ o.outdir | quote}}
params   = {{ args.params | repr}}
idxfile  = {{ args.idxfile | quote}}
kallisto = {{ args.kallisto | quote}}
nthread  = {{ args.nthread | repr}}

shell.load_config(kallisto = kallisto)

params.i = idxfile
params.o = outdir
params.t = nthread
params._ = [fq1, fq2]

shell.kallisto.quant(**params).fg

imfile        = path.join(outdir, 'abundance.tsv')
reader        = TsvReader(imfile)
writer        = TsvWriter(outfile)
writer.cnames = ['target_id', 'est_counts']
writer.writeHead()

for r in reader:
	r.target_id = r.target_id.split('::')[0]
	try:
		r.est_counts = int(round(float(r.est_counts)))
	except TypeError:
		r.est_counts = 0
	writer.write(r)
writer.close()
