from os import path
from sys import stderr
from pyppl import Box
from bioprocs.utils import parallel, shell2 as shell
from bioprocs.utils.reference import vcfIndex

infile   = {{i.infile | quote}}
prefix   = {{i.infile | fn2 | quote}}
outdir   = {{o.outdir | quote}}
samples  = {{i.samples | quote}}
tool     = {{args.tool | quote}}
bcftools = {{args.bcftools | quote}}
gatk     = {{args.gatk | quote}}
tabix    = {{args.tabix | quote}}
ref      = {{args.ref | quote}}
params   = {{args.params | repr}}
nthread  = {{args.nthread | repr}}

shell.load_config(bcftools = bcftools, gatk = gatk)
vcfIndex(infile, tabix = tabix)

allsamples = shell.bcftools.query(l = infile).splitlines()
allsamples = [s.strip() for s in allsamples if s.strip()]
if samples:
	with open(samples) as f:
		samples = f.readlines()
	samples = list(set(allsamples) & set(samples))
else:
	samples = allsamples

def run_bcftools_one(sample):
	shell.fg.bcftools.view(_ = infile, s = sample, o = path.join(outdir, '{}-{}.vcf'.format(prefix, sample)), **params)

def run_bcftools():
	parallel.Parallel(nthread).run(run_bcftools_one, [(sample,) for sample in samples])

def run_awk_one(sample, index, awkfile):
	shell.awk(
		v = ["sample={!r}".format(sample), "index={}".format(index + 10)],
		_stderr = stderr,
		f = awkfile,
		_ = infile,
		_out = path.join(outdir, '{}-{}.vcf'.format(prefix, sample)))

def run_awk():
	# write the awk script
	awkfile = path.join(outdir, 'vcfsample.awk')
	awkfh   = open(awkfile, 'w')
	awkfh.write("""
BEGIN {
	OFS="\\t"
}
$0 ~ "^##" {
	print
}
$0 ~ "^#CHROM" {
	print "#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\tFORMAT\\t"sample
}
$0 !~ "^#" {
	print $1,$2,$3,$4,$5,$6,$7,$8,$9,$index
}
""")
	awkfh.close()
	parallel.Parallel(nthread).run(run_awk_one, [(sample, i, awkfile) for i, sample in enumerate(samples)])

def run_gatk_one(sample):
	shell.fg.gatk(
		R = ref,
		V = infile,
		o = path.join(outdir, '{}-{}.vcf'.format(prefix, sample)),
		sample_name = sample,
		T = 'SelectVariants',
		excludeNonVariants = True,
		**params
	)

def run_gatk():
	parallel.Parallel(nthread).run(run_gatk_one, [(sample, ) for sample in samples])

tools = dict(bcftools = run_bcftools, awk = run_awk, gatk = run_gatk)

try:
	tools[tool]()
except KeyError:
	raise ValueError('Tool {!r} not supported yet.'.format(tool))
