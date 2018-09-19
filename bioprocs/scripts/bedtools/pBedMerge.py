from pyppl import Box
from bioprocs.utils import runcmd, cmdargs

params = Box()
params['i'] = {{i.infile | quote}}

params.update({{args.params}})

cmd  = "{{args.bedtools}} merge %s > {{o.outfile | squote}}" % cmdargs(params, dash = '-', equal = ' ')
runcmd(cmd)
