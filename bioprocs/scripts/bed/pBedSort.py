{{runcmd}}
{{params2CmdArgs}}

from collections import OrderedDict
params = OrderedDict(sorted({{args.params}}))

####### sort
{% if args.tool | lambda x: x == 'sort' %}
params['T']  = {{args.tmpdir | quote}}
params['S']  = {{args.mem | quote}}
params['u']  = {{args.unique}}
params['k']  = '1,1'
params[' k'] = '2,2n'
cmd = 'grep "^#" {{in.infile | quote}} > {{out.outfile | quote}}; grep -v "^#" {{in.infile | quote}} | sort %s >> {{out.outfile | quote}}' % params2CmdArgs(params)

####### bedops
{% elif args.tool | lambda x: x == 'bedops' %}
params['max-mem'] = {{args.mem | quote}}
params['tmpdir']  = {{args.tmpdir | quote}}
cmd = 'grep "^#" {{in.infile | quote}} > {{out.outfile | quote}}; {{args.bedops}} %s {{in.infile | quote}} {% if args.unique %}|uniq{% endif %} >> {{out.outfile | quote}}' % (params2CmdArgs(params, equal = ' '))

####### bedtools
{% elif args.tool | lambda x: x == 'bedtools' %}
params['i'] = {{in.infile | quote}}
cmd = 'grep "^#" {{in.infile | quote}} > {{out.outfile | quote}}; {{args.bedtools}} sort %s {% if args.unique %}|uniq{% endif %} >> {{out.outfile | quote}}' % params2CmdArgs(params, dash = '-', equal = ' ')
{% endif %}

runcmd(cmd)