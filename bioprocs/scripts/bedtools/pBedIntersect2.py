from pyppl import Box
from bioprocs.utils.helpers import runcmd, cmdargs

params = Box()
params['a']           = {{in.afile | quote}}
params['wao']         = True
params['nonamecheck'] = True
       
params.update({{args.params}})

cmd = '{{args.bedtools}} intersect -b {{in.bfiles | lambda x: ' '.join(x)}} %s > {{out.outfile | quote}}' % cmdargs(params, dash='-', equal=' ')
runcmd(cmd)
