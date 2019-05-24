from os import path
from pyppl import Box
from bioprocs.utils import regionOverlap, log2pyppl
from bioprocs.utils.tsvio2 import TsvJoin

infiles = {{i.infiles}}
inopts  = {{args.inopts | repr}}
outopts = {{args.outopts | repr}}
debug   = {{args.debug}}
outfile = {{o.outfile | quote}}
sefile  = {{job.errfile | quote}}
helper  = {{args.helper | repr}}
if not isinstance(helper, list):
	helper = [helper]

compare = TsvJoin.compare
helper  = [line for line in helper if line]
exec('\n'.join(helper), globals())

do      = {{args.do}}
match   = {{args.match}}

tj = TsvJoin(*infiles, **inopts)
tj.debug = debug
outopts['headCallback'] = {{args.outopts.get('headCallback', True)}}
tj.join(do, outfile, match, outopts)

with open(sefile, 'r') as f:
	f.flush()
	i = 0
	for line in f:
		if line.startswith('pyppl.log'):
			continue
		log2pyppl(line.rstrip(), 'warning')
		i += 1
		if i >= 20:
			break
	if i > 0:
		log2pyppl('See more in {}'.format(sefile), 'warning')
