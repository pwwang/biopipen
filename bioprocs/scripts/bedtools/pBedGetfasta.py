from pyppl import Box
from bioprocs.utils import shell

infile   = {{ i.infile | quote}}
outfile  = {{ o.outfile | quote}}
bedtools = {{ args.bedtools | quote}}
params   = {{ args.params | repr }}
ref      = {{ args.ref | quote}}

shell.TOOLS.bedtools = bedtools

params.fi      = ref
params.bed     = infile
params._stdout = outfile
shell.Shell(subcmd = True, dash = '-', equal = ' ').bedtools.getfasta(**params).run()
