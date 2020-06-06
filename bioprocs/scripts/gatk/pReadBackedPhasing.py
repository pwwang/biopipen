from diot import Diot
from bioprocs.utils import shell2 as shell, mem2

infile   = {{i.infile | quote}}
bamfile  = {{i.bamfile | quote}}
outfile  = {{o.outfile | quote}}
gatk3    = {{args.gatk | quote}}
ref      = {{args.ref | quote}}
tmpdir   = {{args.tmpdir | quote}}
mem      = {{args.mem | quote}}
params   = {{args.params | repr}}
interval = {{args.interval | repr}}

shell.load_config(gatk3 = gatk3)

params.T = 'ReadBackedPhasing'
params.R = ref
params.I = bamfile
params.o = outfile
params.variant = infile
params.L = interval or infile

shell.gatk3(f'-Djava.io.tmpdir={tmpdir}',
            *mem2(mem, 'java').split(), **params).fg
