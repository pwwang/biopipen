import shlex
from os import path, rename
{{runcmd}}
{{params2CmdArgs}}

params = {{args.params}}
{% if i.infile.endswith('avinput') %}
params['annovar-file'] = {{i.infile | quote}}
{% else %}
params['vcf-file'] = {{i.infile | quote}}
{% endif %}
params['regulatory-causing-predict'] = 'all'
params['cell']                       = {{args.cell | quote}}
params['db-score']                   = 'dbncfp'
params['out']                        = {{o.outfile | prefix | prefix | quote}}

jarfile = [x for x in shlex.split({{args.cepip | quote}}) if x.endswith('.jar')][0]
jardir  = path.dirname(jarfile)

cmd = 'cd "%s"; {{args.cepip}} %s' % (jardir, params2CmdArgs(params, dash='--', equal=' '))
runcmd(cmd)

rename(params['out'] + '.flt.txt', {{o.outfile | quote}})