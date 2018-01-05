{{runcmd}}
{{params2CmdArgs}}

params = {}
params['a']           = {{in.afile | quote}}
params['b']           = {{in.bfile | quote}}
params['wao']         = True
params['nonamecheck'] = True
       
params.update({{args.params}})

cmd = '{{args.bedtools}} intersect %s > {{out.outfile | quote}}' % params2CmdArgs(params, dash='-', equal=' ')
runcmd(cmd)
