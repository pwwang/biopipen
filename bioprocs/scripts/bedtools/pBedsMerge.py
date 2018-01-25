{{runcmd}}
{{params2CmdArgs}}
from os import remove

mergedtmp       = '{{job.outdir}}/mergedtmp.bed'
mergedtmpsorted = '{{job.outdir}}/mergedtmp.sorted.bed'
cmd             = 'echo "" > "%s"' % mergedtmp
for infile in {{in.infiles}}:
	cmd += '; grep -v "^#" "%s" >> "%s"' % (infile, mergedtmp)

params = {}
params['i'] = mergedtmpsorted

params.update({{args.params}})
cmd += '; sort -k1,1 -k2,2n "%s" > "%s"; {{args.bedtools}} merge %s > {{out.outfile | quote}}' % (mergedtmp, mergedtmpsorted, params2CmdArgs(params, dash = '-', equal = ' '))
runcmd(cmd)

remove (mergedtmp)
remove (mergedtmpsorted)