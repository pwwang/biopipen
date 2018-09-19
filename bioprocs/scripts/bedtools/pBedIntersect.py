from pyppl import Box
from bioprocs.utils import runcmd, cmdargs

params = Box()
params['a']           = {{i.afile | quote}}
params['b']           = {{i.bfile | quote}}
params['wao']         = True
params['nonamecheck'] = True
       
params.update({{args.params}})

cmd = '{{args.bedtools}} intersect %s > {{o.outfile | quote}}' % cmdargs(params, dash='-', equal=' ')
runcmd(cmd)
