from os import path
from pyppl import Box
from bioprocs.utils import shell, parallel

infile    = {{i.infile | quote}}
outfile   = {{o.outfile | quote}}
tool      = {{args.tool | quote}}
vcf2maf   = {{args.vcf2maf | quote}}
vep       = {{args.vep | quote}}
vepDb     = {{args.vepDb | quote}}
filtervcf = {{args.filtervcf | quote}}
ref       = {{args.ref | quote}}
bcftools  = {{args.bcftools | quote}}
tumoridx  = {{args.tumor | repr}}
nthread   = {{args.nthread | repr}}
params    = {{args.params | repr}}

shell.TOOLS.vcf2maf  = vcf2maf
shell.TOOLS.bcftools = bcftools

veppath  = path.dirname(shell.which(vep))
vcf2maf  = shell.Shell(equal = ' ').vcf2maf
bcftools = shell.Shell(subcmd = True, equal = ' ').bcftools

def run_vcf2maf_one(vcf, maf, tumor, normal = None, forks = nthread):
	params['input-vcf']  = vcf
	params['output-maf'] = maf
	params['vep-data']   = vepDb
	params['vep-forks']  = forks
	params['filter-vcf'] = filtervcf
	params['ref-fasta']  = ref
	params['vep-path']   = veppath
	params['tumor-id']   = tumor
	params['normal-id']  = normal
	vcf2maf(**params).run()

def extract_sample_from_vcf(vcf, sample, outvcf):
	bcftools.view(s = sample, o = outvcf + '.tmp', _ = vcf).run()
	# exclude sites without genotypes
	# because this is split from a merged vcf file
	bcftools.view(_ = outvcf + '.tmp', g = '^miss', o = outvcf).run()
	shell.rmrf(outvcf + '.tmp')

def run_vcf2maf():
	vcfsams  = bcftools.query(l = infile).stdout.splitlines()
	vcfsams  = [s for s in vcfsams if s.strip()]
	
	if tumoridx < 2: # single sample
		tumor  = vcfsams[tumoridx]
		normal = vcfsams[1-tumoridx]
		run_vcf2maf_one(infile, outfile, tumor, normal)
	else:
		# split vcf file
		splitdir = path.join(path.dirname(outfile), "splits")
		shell.mkdir(p = splitdir)
		mafdir = path.join(path.dirname(outfile), "mafs")
		shell.mkdir(p = mafdir)

		para = parallel.Parallel(nthread)
		para.run(extract_sample_from_vcf, [
			(infile, s, path.join(splitdir, "split{}.vcf".format(i+1))) 
			for i, s in enumerate(vcfsams)
		])
		restThreads = int(float(nthread)/float(len(vcfsams)) - 1.0)
		restThreads = max(restThreads, 1)
		para.run(run_vcf2maf_one, [
			(path.join(splitdir, "split" + str(i+1) + ".vcf"), path.join(mafdir, "split{}.maf".format(i+1)), s, None, restThreads) 
			for i, s in enumerate(vcfsams)
		])
		del para

		# merge mafs
		shell.head(n = 1, _ = path.join(mafdir, "split1.maf"), _stdout = outfile)
		for i, s in enumerate(vcfsams):
			shell.tail(n = '+2', _ = path.join(mafdir, "split{}.maf".format(i+1)), _stdout_ = outfile)

tools = dict(
	vcf2maf = run_vcf2maf
)

try:
	tools[tool]()
except KeyError:
	raise ValueError('Tool {!r} not supported.'.format(tool))