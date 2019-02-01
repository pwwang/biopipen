from pyppl import Box
from bioprocs.utils import shell
from bioprocs.utils.shell import Shell

infile   = {{ i.infile | quote}}
outfile  = {{ o.outfile | quote }}
outidx   = {{ o.outidx | quote }}
tool     = {{ args.tool | quote}}
samtools = {{ args.samtools | quote}}
params   = {{ args.params | repr}}
nthread  = {{ args.nthread | repr}}

shell.TOOLS['samtools'] = samtools

shell.ln_s(infile, outfile)

def run_samtools():
	shellst = Shell(subcmd = True).samtools
	params['@'] = nthread
	shellst.index(_ = [infile, outidx]).run()

tools = dict(
	samtools = run_samtools
)

try:
	tools[tool]()
except KeyError:
	raise KeyError('Tool {!r} is not supported.'.format(tool))
except:
	raise
finally:
	pass
