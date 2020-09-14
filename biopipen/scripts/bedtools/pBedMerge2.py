from os import remove
from diot import Diot
from bioprocs.utils import shell2 as shell

infiles = {{i.infile | repr}}
params  = {{args.params | repr}}

shell.load_config(bedtools = {{args.bedtools | quote}})

mergedtmp       = '{{job.outdir}}/mergedtmp.bed'
mergedtmpsorted = '{{job.outdir}}/mergedtmp.sorted.bed'

shell.touch(mergedtmp)

for infile in infiles:
	shell.grep(v = '^#', _ = infile).r >> mergedtmp

shell.sort(k = ['1,1', '2,2n'], _ = mergedtmp).r > mergedtmpsorted

params.i = mergedtmpsorted
# params._stdout = outfile
shell.bedtools.merge(**params).r > outfile

shell.rm_rf(mergedtmp)
shell.rm_rf(mergedtmpsorted)
