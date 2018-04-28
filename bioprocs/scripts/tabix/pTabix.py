from os import path
from pyppl import Box
from bioprocs.utils import runcmd, cmdargs

params = {}
region = {{in.region | quote}}
outfile = {{out.outfile | quote}}

if path.isfile(region):
	if path.getsize(region) > 0:
		params['R'] = {{in.region | quote}}
		params.update({{args.params}})
		cmd = '{{args.tabix}} %s "{{in.infile}}" > "{{out.outfile}}"' % cmdargs(params, equal=' ')
		runcmd(cmd)
	else:
		open(outfile, 'w').close()
else:
	if region:
		params.update({{args.params}})
		cmd = '{{args.tabix}} %s "{{in.infile}}" {{in.region}} > "{{out.outfile}}"' % cmdargs(params, equal=' ')
		runcmd(cmd)
	else:
		open(outfile, 'w').close()
