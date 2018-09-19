from os import makedirs, path
from shutil import rmtree
from pyppl import Box
from bioprocs.utils import cmdargs, runcmd, mem2

tmpdir    = path.join ("{{args.tmpdir}}", "{{proc.id}}.{{i.infile | fn}}.{{job.index}}")
if not path.exists (tmpdir): makedirs (tmpdir)

mem    = mem2({{ args.mem | quote }}, 'java')
ref    = {{args.ref | quote}}
params = {{args.params}}
try:
	params['TMP_DIR'] = tmpdir
	params['I']       = {{i.infile | quote}}
	params['O']       = {{o.outfile | quote}}
	params['R']       = ref
	runcmd ('{{args.picard}} ReorderSam %s -Djava.io.tmpdir="%s" %s' % (mem, tmpdir, cmdargs(params, dash = '', equal = '=')))
except Exception as ex:
	stderr.write ("Job failed: %s" % str(ex))
	raise
finally:
	rmtree (tmpdir)
