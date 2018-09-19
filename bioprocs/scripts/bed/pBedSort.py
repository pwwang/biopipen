from collections import OrderedDict
from pyppl import Box
from bioprocs.utils import runcmd, cmdargs
params = Box()

by = {{args.by | quote}}

####### sort
{% case args.tool %}
{% when 'sort' %}
params['T']  = {{args.tmpdir | quote}}
params['S']  = {{args.mem | quote}}
params['u']  = {{args.unique}}
if by == 'coord':
	params['k']  = '1,1'
	params[' k'] = '2,2n'
else:
	params['k'] = '4'
params.update({{args.params}})
cmd = 'grep "^#" {{i.infile | quote}} > {{o.outfile | quote}}; grep -v "^#" {{i.infile | quote}} | sort %s >> {{o.outfile | quote}}' % cmdargs(params)

####### bedops
{% when 'bedops' %}
params['max-mem'] = {{args.mem | quote}}
params['tmpdir']  = {{args.tmpdir | quote}}
params.update({{args.params}})
cmd = 'grep "^#" {{i.infile | quote}} > {{o.outfile | quote}}; {{args.bedops}} %s {{i.infile | quote}} {% if args.unique %}|uniq{% endif %} >> {{o.outfile | quote}}' % (cmdargs(params, equal = ' '))

####### bedtools
{% when 'bedtools' %}
params['i'] = {{i.infile | quote}}
params.update({{args.params}})
cmd = 'grep "^#" {{i.infile | quote}} > {{o.outfile | quote}}; {{args.bedtools}} sort %s {% if args.unique %}|uniq{% endif %} >> {{o.outfile | quote}}' % cmdargs(params, dash = '-', equal = ' ')
{% endcase %}

runcmd(cmd)
