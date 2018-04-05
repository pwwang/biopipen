from os import remove
from pyppl import Box
from bioprocs.utils.helpers import runcmd, cmdargs

mergedtmp       = '{{job.outdir}}/mergedtmp.bed'
mergedtmpsorted = '{{job.outdir}}/mergedtmp.sorted.bed'
cmd             = 'echo "" > "%s"' % mergedtmp
for infile in {{in.infiles}}:
	cmd += '; grep -v "^#" "%s" >> "%s"' % (infile, mergedtmp)

params = Box()
params['i'] = mergedtmpsorted

params.update({{args.params}})
cmd += '; sort -k1,1 -k2,2n "%s" > "%s"; {{args.bedtools}} merge %s > {{out.outfile | quote}}' % (mergedtmp, mergedtmpsorted, cmdargs(params, dash = '-', equal = ' '))
runcmd(cmd)

remove (mergedtmp)
remove (mergedtmpsorted)