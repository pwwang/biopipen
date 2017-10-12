from os import makedirs, path
from shutil import rmtree
	
tmpdir    = path.join ("{{args.tmpdir}}", "{{proc.id}}.{{in.infile | fn}}.{{job.index}}")
if not path.exists (tmpdir): makedirs (tmpdir)
	
{{ runcmd }}
{{ mem2 }}
{{ params2CmdArgs }}

mem    = mem2({{ args.mem | quote }}, 'java')
ref    = {{args.ref | quote}}
params = {{args.params}}
try:
	params['TMP_DIR'] = tmpdir
	params['I']       = {{in.infile | quote}}
	params['O']       = {{out.outfile | quote}}
	params['R']       = ref
	runcmd ('{{args.picard}} ReorderSam %s -Djava.io.tmpdir="%s" %s' % (mem, tmpdir, params2CmdArgs(params, dash = '', equal = '=')))
except Exception as ex:		
	stderr.write ("Job failed: %s" % str(ex))
	raise
finally:
	rmtree (tmpdir)