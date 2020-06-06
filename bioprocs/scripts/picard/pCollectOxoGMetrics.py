from diot import Diot
from bioprocs.utils import shell2 as shell

infile  = {{i.infile | quote}}
outfile = {{o.outfile | quote}}
ref     = {{args.ref | quote}}
picard  = {{args.picard | quote}}
params  = {{args.params | repr}}

shell.load_config(picard = picard)

shell.picard.CollectOxoGMetrics(I = infile, R = ref, O = outfile, **params).fg
