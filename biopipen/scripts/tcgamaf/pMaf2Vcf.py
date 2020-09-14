
from diot import Diot
from bioprocs.utils import shell2 as shell
from tempfile import gettempdir
from filelock import FileLock


infile  = {{i.infile | quote}}
outdir  = {{o.outdir | quote}}
maf2vcf = {{args.maf2vcf | quote}}
ref     = {{args.ref | repr}}
params  = {{args.params | repr}}

params['input-maf']  = infile
params['output-dir'] = outdir
params['ref-fasta']  = ref

# due to https://github.com/mskcc/vcf2maf/issues/234, now the safest way to do it is to patch
# maf2vcf.pl and change that 5000 to 1.
maf2vcf_path = shell.which(maf2vcf)

maf2vcf_real = gettempdir() + '/maf2vcf_patched_by_pMaf2vcf.pl'
# thread-safe
with FileLock(maf2vcf_real):
	with open(maf2vcf_path) as fin, open(maf2vcf_real, 'w') as fout:
		fout.write(fin.read().replace('splice( @regions, 0, 5000 )', 'splice( @regions, 0, 1 )'))
	shell.chmod('+x', maf2vcf_real)

shell.load_config(maf2vcf = maf2vcf_real)

shell.maf2vcf(**params).fg
