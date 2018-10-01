from pyppl import Box
from collections import OrderedDict
from bioprocs.utils import runcmd, cmdargs
from bioprocs.utils.tsvio import TsvReader, TsvWriter

{% if args.sorted %}
cmd = "ln -s {{i.infile | squote}} {{o.outfile | squote}}"

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

infile   = {{i.infile | quote}}
outfile  = {{o.outfile | quote}}
tmpfile  = outfile + '.tmp'
skip     = {{args.inopts | lambda x: x.get('skip', 0)}}
delimit  = {{args.inopts | lambda x: x.get('delimit', '\t') | quote}}
comment  = {{args.inopts | lambda x: x.get('comment', '#') | quote}}

if not skip and not comment:
	tmpfile = infile
else:
	with open(infile) as readerSkip, open(outfile, 'w') as writerSkip:
		for i, line in enumerate(readerSkip):
			if i >= skip: break
			writerSkip.write(line)

	readerTmp = TsvReader(infile, delimit = delimit, comment = comment, skip = skip, ftype = 'nometa', head = False)
	#readerTmp.autoMeta()
	writerTmp = TsvWriter(tmpfile, delimit = delimit, ftype = 'nometa')
	#writerTmp.meta.update(readerTmp.meta)
	for r in readerTmp:
		writerTmp.write(r)
	writerTmp.close()
	
{% if args.case %}
case = "LANG=C"
{% else %}
case = "LANG=en_US.UTF-8"
{% endif %}

cmd = '%s sort %s "%s" >> {{o.outfile | quote}}' % (case, cmdargs(params), tmpfile)
{% endif %}

runcmd(cmd)