from pyppl import Box
from bioprocs.utils.helpers import runcmd, cmdargs

params = Box()
params['a']           = {{in.afile | quote}}
params['b']           = {{in.bfile | quote}}
params['wao']         = True
params['nonamecheck'] = True
       
params.update({{args.params}})

cmd = '{{args.bedtools}} intersect %s > {{out.outfile | quote}}' % cmdargs(params, dash='-', equal=' ')
runcmd(cmd)
