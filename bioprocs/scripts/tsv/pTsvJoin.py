from os import path
from pyppl import Box
from bioprocs.utils import regionOverlap
from bioprocs.utils.tsvio2 import TsvJoin

infiles = {{i.infiles}}
inopts  = {{args.inopts}}
outopts = {{args.outopts}}
outfile = {{o.outfile | quote}}

{{args.helper}}
do      = {{args.do}}
match   = {{args.match}}

tj = TsvJoin(*infiles, **inopts)
tj.join(do, outfile, match, outopts)
