
from pyppl import Box
from bioprocs.utils import shell2 as shell

infile  = {{i.infile | quote}}
outdir  = {{o.outdir | quote}}
maf2vcf = {{args.maf2vcf | quote}}
ref     = {{args.ref | quote}}
params  = {{args.params | repr}}

shell.load_config(maf2vcf = maf2vcf)

params['input-maf']  = infile
params['output-dir'] = outdir
params['ref-fasta']  = ref

shell.fg.maf2vcf(**params)
