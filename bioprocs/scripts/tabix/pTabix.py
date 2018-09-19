from os import path
from pyppl import Box
from bioprocs.utils import runcmd, cmdargs

params = {}
region = {{i.region | quote}}
outfile = {{o.outfile | quote}}

if path.isfile(region):
	if path.getsize(region) > 0:
		params['R'] = {{i.region | quote}}
		params.update({{args.params}})
		cmd = '{{args.tabix}} %s "{{i.infile}}" > "{{o.outfile}}"' % cmdargs(params, equal=' ')
		runcmd(cmd)
	else:
		open(outfile, 'w').close()
else:
	if region:
		params.update({{args.params}})
		cmd = '{{args.tabix}} %s "{{i.infile}}" {{i.region}} > "{{o.outfile}}"' % cmdargs(params, equal=' ')
		runcmd(cmd)
	else:
		open(outfile, 'w').close()
