import re
from os import makedirs, path
from shutil import rmtree
from sys import stderr
from pyppl import Box
from bioprocs.utils import cmdargs, runcmd, mem2

# detemine default read group
rg = {{ args.rg }}
rg = {key.upper():val for key, val in rg.items()}
if not rg['ID']:
	g = re.search (r'[^a-zA-Z0-9]+(L\\d+)[^a-zA-Z0-9]+', "{{out.outfile | fn}}")
	rg['ID'] = g.group(1) if g else "{{out.outfile | fn}}.L{{job.index}}"
if not rg['SM']:
	rg['SM'] = "{{out.outfile | fn}}"

tmpdir    = path.join ("{{args.tmpdir}}", "{{proc.id}}.{{in.infile | fn}}.{{job.index}}")
if not path.exists (tmpdir):
	makedirs (tmpdir)

params = {{args.params}}
try:
	############## picard
	{% if args.tool | lambda x: x == 'picard' %}
	mem = mem2({{ args.mem | quote }})
	params['-Djava.io.tmpdir'] = tmpdir
	params['TMP_DIR'] = tmpdir
	params['I'] = {{in.infile | quote}}
	params['O'] = {{out.outfile | quote}}
	for k,v in rg.items():
		params['RG' + k] = v

	runcmd ('{{args.picard}} AddOrReplaceReadGroups %s %s' % (mem, cmdargs(params, dash='', equal='=')))

	############## bamutil
	{% elif args.tool | lambda x: x == 'bamutil' %}
	params['RG'] = "@RG\\tID:%s\\t%s" % (rg['ID'], "\\t".join([k + ":" + v for k,v in rg.items() if k!='ID']))
	params['in'] = {{in.infile | quote}}
	params['out'] = {{out.outfile | quote}}

	runcmd ('{{args.bamutil}} polishBam %s' % cmdargs(params, equal = ' '))

	{% endif %}
except Exception as ex:
	stderr.write ("Job failed: %s" % str(ex))
	raise
finally:
	rmtree (tmpdir)
