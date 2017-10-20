{{runcmd}}
{{params2CmdArgs}}

params = {{args.params}}

params['fi']  = {{args.ref | quote}}
params['bed'] = {{in.infile | quote}}

cmd = '{{args.bedtools}} getfasta %s > {{out.outfile | quote}}' % params2CmdArgs(params, dash='-', equal=' ')
runcmd(cmd)
