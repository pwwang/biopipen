{{runcmd}}
{{params2CmdArgs}}

from collections import OrderedDict
params = {{args.params}}
params = OrderedDict(sorted(params.items()))

{% if args.skip %}
cmd = '(head -n {{args.skip}} "{{in.infile}}" && tail -n +{{args.skip | lambda x: x+1}} "{{in.infile}}" | sort %s) > "{{out.outfile}}"' % params2CmdArgs(params, dash='-', equal=' ')
{% else %}
cmd = 'sort %s "{{in.infile}}" > "{{out.outfile}}"' % params2CmdArgs(params, dash='-', equal=' ')
{% endif %}
runcmd(cmd)