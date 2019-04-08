from pyppl import Box
from bioprocs.utils.shell2 import sh, _update_args

infile  = {{i.infile | quote}}
outfile = {{o.outfile | quote}}
ref     = {{args.ref | quote}}
picard  = {{args.picard | quote}}
params  = {{args.params | repr}}

_update_args(picard = picard)

sh.picard.CollectOxoGMetrics(I = infile, R = ref, O = outfile, **params)