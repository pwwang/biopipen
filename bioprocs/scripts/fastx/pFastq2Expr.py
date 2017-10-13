from collections import OrderedDict

{{params2CmdArgs}}
{{runcmd}}
{{txtFilter}}
{{txtTransform}}

params = OrderedDict({{args.params}})
params['i'] = {{args.idxfile | quote}}
params['o'] = {{out.outdir | quote}}
params['t'] = {{args.nthread}}

cmd = '{{args.kallisto}} quant %s "{{in.fqfile1}}" "{{in.fqfile2}}"' % (params2CmdArgs(params))
runcmd (cmd)

retfile = "{{out.outdir}}/abundance.tsv"
txtFilter(retfile, "{{out.outfile}}.dec", cols=[0,3], header=False, skip=1)
txtTransform("{{out.outfile}}.dec", "{{out.outfile}}", transformer = lambda parts: [parts[0], str(int(float(parts[1])))], header=False)