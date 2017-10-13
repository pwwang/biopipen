from os import path
{{runcmd}}
{{params2CmdArgs}}

params = {{args.params}}

if path.isfile({{in.region | quote}}):
	params['R'] = {{in.region | quote}}
	cmd = '{{args.tabix}} %s "{{in.infile}}" > "{{out.outfile}}"' % params2CmdArgs(params, equal=' ')
else:
	cmd = '{{args.tabix}} %s "{{in.infile}}" {{in.region}} > "{{out.outfile}}"' % params2CmdArgs(params, equal=' ')
runcmd(cmd)