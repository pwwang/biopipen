from pyppl import Box
from bioprocs.utils import runcmd, cmdargs
from bioprocs.utils.tsvio import TsvReader, TsvWriter

params = Box()
params.fi = {{args.ref | quote}}
params.update({{args.params}})

infile  = {{i.infile | quote}}
infile2 = {{i._infile | quote}} + '.names'
reader  = TsvReader(infile, ftype = 'bed')
writer  = TsvWriter(infile2)
writer.meta.update(reader.meta)
for r in reader:
	if not params.name:
		r.NAME = ''
	r.NAME += '::{}:{}-{}'.format(r.CHR, r.START, r.END)
	writer.write(r)
writer.close()
params.name = True
params.bed  = infile2

cmd = "{{args.bedtools}} getfasta %s > {{o.outfile | squote}}" % cmdargs(params, dash='-', equal=' ')
print cmd
runcmd(cmd)
