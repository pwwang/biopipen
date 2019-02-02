from pyppl import Box
from bioprocs.utils import shell

cnvkit   = {{args.cnvkit | quote}}
infile   = {{i.tgfile | quote}}
atgfile  = {{i.atgfile | quote}}
ref      = {{args.ref | quote}}
outfile  = {{o.outfile | quote}}
params   = {{args.params}}

shell.TOOLS['cnvkit'] = cnvkit
ckshell = shell.Shell(subcmd = True, equal = ' ').cnvkit

params.o = outfile
params.f = ref
params.t = infile
params.a = atgfile
ckshell.reference(**params).run()
