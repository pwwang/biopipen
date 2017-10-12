from collections import OrderedDict

{{params2CmdArgs}}
{{runcmd}}
{{txtFilter}}

params = OrderedDict({{args.params}})
params['i'] = {{args.idxfile | quote}}
params['o'] = {{out.outdir | quote}}
params['t'] = {{args.nthread}}

cmd = '{{args.kallisto}} quant %s "{{in.fqfile1}}" "{{in.fqfile2}}"' % (params2CmdArgs(params))
runcmd (cmd)

retfile = "{{out.outdir}}/abundance.tsv"
txtFilter(retfile, "{{out.outfile}}", cols=[0,3])