from os import path
from bioprocs.utils import runcmd, cmdargs

params = {}

if path.isfile({{in.region | quote}}):
	params['R'] = {{in.region | quote}}
	params.update({{args.params}})
	cmd = '{{args.tabix}} %s "{{in.infile}}" > "{{out.outfile}}"' % cmdargs(params, equal=' ')
else:
	params.update({{args.params}})
	cmd = '{{args.tabix}} %s "{{in.infile}}" {{in.region}} > "{{out.outfile}}"' % cmdargs(params, equal=' ')
runcmd(cmd)