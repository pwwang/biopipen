from pyppl import Box
from bioprocs.utils import runcmd, cmdargs

infile   = {{ i.infile | quote}}
outfile  = {{ o.outfile | quote}}
bedtools = {{ args.bedtools | quote}}
params   = {{ args.params | repr }}
ref      = {{ args.ref | quote}}

params.fi = ref
params.bed = infile
cmd = '{bedtools!r} getfasta {params} > {outfile!r}'.format(
	bedtools = bedtools,
	params   = cmdargs(params, dash = '-', equal = ' '),
	outfile  = outfile
)
runcmd(cmd)
