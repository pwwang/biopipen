from pyppl import Box
from bioprocs.utils.helpers import runcmd, cmdargs

params = Box()
params['fi']  = {{args.ref | quote}}
params['bed'] = {{in.infile | quote}}
params.update({{args.params}})

cmd = '{{args.bedtools}} getfasta %s > {{out.outfile | quote}}' % cmdargs(params, dash='-', equal=' ')
runcmd(cmd)
