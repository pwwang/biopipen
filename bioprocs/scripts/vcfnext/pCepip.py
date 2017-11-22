import shlex
from os import path, rename
{{runcmd}}
{{params2CmdArgs}}

params = {{args.params}}
{% if in.infile | lambda x: x.endswith('avinput') %}
params['annovar-file'] = {{in.infile | quote}}
{% else %}
params['vcf-file'] = {{in.infile | quote}}
{% endif %}
params['regulatory-causing-predict'] = 'all'
params['cell']                       = {{args.cell | quote}}
params['db-score']                   = 'dbncfp'
params['out']                        = {{out.outfile | prefix | prefix | quote}}

jarfile = [x for x in shlex.split({{args.cepip | quote}}) if x.endswith('.jar')][0]
jardir  = path.dirname(jarfile)

cmd = 'cd "%s"; {{args.cepip}} %s' % (jardir, params2CmdArgs(params, dash='--', equal=' '))
runcmd(cmd)

rename(params['out'] + '.flt.txt', {{out.outfile | quote}})