from pyppl import Box
from bioprocs.utils.helpers import runcmd, cmdargs

params = Box()
params['i'] = {{in.infile | quote}}

params.update({{args.params}})

cmd  = "{{args.bedtools}} merge %s > {{out.outfile | squote}}" % cmdargs(params, dash = '-', equal = ' ')
runcmd(cmd)