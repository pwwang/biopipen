from collections import OrderedDict
from pyppl import Box
from bioprocs.utils.helpers import runcmd, cmdargs
params = Box()

####### sort
{% if args.tool | lambda x: x == 'sort' %}
params['T']  = {{args.tmpdir | quote}}
params['S']  = {{args.mem | quote}}
params['u']  = {{args.unique}}
params['k']  = '1,1'
params[' k'] = '2,2n'
params.update({{args.params}})
cmd = 'grep "^#" {{in.infile | quote}} > {{out.outfile | quote}}; grep -v "^#" {{in.infile | quote}} | sort %s >> {{out.outfile | quote}}' % cmdargs(params)

####### bedops
{% elif args.tool | lambda x: x == 'bedops' %}
params['max-mem'] = {{args.mem | quote}}
params['tmpdir']  = {{args.tmpdir | quote}}
params.update({{args.params}})
cmd = 'grep "^#" {{in.infile | quote}} > {{out.outfile | quote}}; {{args.bedops}} %s {{in.infile | quote}} {% if args.unique %}|uniq{% endif %} >> {{out.outfile | quote}}' % (cmdargs(params, equal = ' '))

####### bedtools
{% elif args.tool | lambda x: x == 'bedtools' %}
params['i'] = {{in.infile | quote}}
params.update({{args.params}})
cmd = 'grep "^#" {{in.infile | quote}} > {{out.outfile | quote}}; {{args.bedtools}} sort %s {% if args.unique %}|uniq{% endif %} >> {{out.outfile | quote}}' % cmdargs(params, dash = '-', equal = ' ')
{% endif %}

runcmd(cmd)