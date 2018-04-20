from pyppl import Box
from collections import OrderedDict
from bioprocs.utils import runcmd, cmdargs
from bioprocs.utils.tsvio import TsvReader, TsvWriter

{% if args.sorted %}
cmd = "ln -s {{in.infile | squote}} {{out.outfile | squote}}"

{% else %}
params = Box()
params['T'] = {{args.tmpdir | quote}}
params['S'] = {{args.mem | quote}}
params['u'] = {{args.unique}}
params['t'] = {{args.inopts.delimit | quote}}
params.update({{args.params}})
kopts = {k:v for k,v in params.items() if k.startswith('k')}
for i, k in enumerate(sorted(kopts.keys())):
	del params[k]
	params['k%s' % (' '*i)] = kopts[k]

infile   = {{in.infile | quote}}
outfile  = {{out.outfile | quote}}
tmpfile  = outfile + '.tmp'
skip     = {{args.inopts | lambda x: x.skip if 'skip' in x else 0}}
delimit  = {{args.inopts | lambda x: x.delimit if 'delimit' in x else '\t' | quote}}
comment  = {{args.inopts | lambda x: x.comment if 'comment' in x else '#' | quote}}

if not skip and not comment:
	tmpfile = infile
else:
	readerSkip = TsvReader(infile, delimit = delimit, comment = '', skip = 0)
	readerSkip.autoMeta()
	writerSkip = TsvWriter(outfile, delimit = delimit)
	writerSkip.meta.update(readerSkip.meta)
	for i, r in enumerate(readerSkip):
		if i < skip:
			writerSkip.write(r)
	writerSkip.close()
	readerSkip.close()
	readerTmp = TsvReader(infile, delimit = delimit, comment = comment, skip = skip)
	readerTmp.autoMeta()
	writerTmp = TsvWriter(tmpfile, delimit = delimit)
	writerTmp.meta.update(readerTmp.meta)
	for r in readerTmp:
		writerTmp.write(r)
	writerTmp.close()
	
{% if args.case %}
case = "LANG=C"
{% else %}
case = "LANG=en_US.UTF-8"
{% endif %}

cmd = '%s sort %s "%s" >> {{out.outfile | quote}}' % (case, cmdargs(params), tmpfile)
{% endif %}

runcmd(cmd)