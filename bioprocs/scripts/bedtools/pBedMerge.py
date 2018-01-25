{{runcmd}}
{{params2CmdArgs}}

params = {}
params['i'] = {{in.infile | quote}}

params.update({{args.params}})

cmd  = '{{args.bedtools}} merge %s' % params2CmdArgs(params, dash = '-', equal = ' ')
runcmd(cmd)