{{runcmd}}
{{params2CmdArgs}}

from collections import OrderedDict

params = {}
params['T'] = {{args.tmpdir | quote}}
params['S'] = {{args.mem | quote}}
params['u'] = {{args.unique}}
params['t'] = {{args.delimit | quote}}
params.update({{args.params}})
params = OrderedDict(sorted(params.items()))

{% if args.case %}
case = "LANG=C"
{% else %}
case = "LANG=en_US.UTF-8"
{% endif %}

infile = {{in.infile | quote}}
{% if args.noeline %}
noelinefile = '{{out.outfile}}.noeline'
with open(infile) as f, open(noelinefile, 'w') as fout:
	for line in f:
		if not line.strip(): continue
		fout.write(line)
infile = noelinefile
{% endif %}

{% if args.skip %}
inskipfile = '{{out.outfile}}.skip'
with open(infile) as f, \
	open(inskipfile, 'w') as fin, \
	open({{out.outfile | quote}}, 'w') as fout:
	for _ in range({{args.skip}}):
		fout.write(f.readline())
	for line in f:
		fin.write(line)
infile = inskipfile
{% endif %}

cmd = '%s sort %s "%s" >> {{out.outfile | quote}}' % (case, params2CmdArgs(params), infile)
runcmd(cmd)