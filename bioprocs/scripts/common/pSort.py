{{runcmd}}
{{params2CmdArgs}}

from collections import OrderedDict
params = {{args.params}}
params = OrderedDict(sorted(params.items()))

{% if args.case %}
case = "LANG=C"
{% else %}
case = "LANG=en_US.UTF-8"
{% endif %}

{% if args.skip %}
cmd = '(head -n {{args.skip}} "{{in.infile}}" && tail -n +{{args.skip | lambda x: x+1}} "{{in.infile}}" | %s sort %s) > "{{out.outfile}}"' % (case, params2CmdArgs(params, dash='-', equal=' '))
{% else %}
cmd = '%s sort %s "{{in.infile}}" > "{{out.outfile}}"' % (case, params2CmdArgs(params, dash='-', equal=' '))
{% endif %}
runcmd(cmd)