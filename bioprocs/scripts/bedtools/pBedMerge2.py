from os import remove
from pyppl import Box
from bioprocs.utils import shell

infiles = {{i.infile | repr}}
params  = {{args.params | repr}}

shell.TOOLS.bedtools = {{args.bedtools | quote}}

mergedtmp       = '{{job.outdir}}/mergedtmp.bed'
mergedtmpsorted = '{{job.outdir}}/mergedtmp.sorted.bed'

shell.touch(mergedtmp)

for infile in infiles:
	shell.grep(v = '^#', _ = infile, __stdout = mergedtmp)

shell.sort(k = ['1,1', '2,2n'], _ = mergedtmp, _stdout = mergedtmpsorted)

params.i = mergedtmpsorted
params._stdout = outfile
shell.Shell(subcmd = True, dash = '-', equal = ' ').bedtools.merge(**params).run()

shell.rm_rf(mergedtmp)
shell.rm_rf(mergedtmpsorted)
