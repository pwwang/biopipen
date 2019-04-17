from os import path, chdir
from pyppl import Box
from bioprocs.utils import parallel, logger, shell2 as shell

infile       = {{i.infile | quote}}
outfile      = {{o.outfile | quote}}
tool         = {{args.tool | quote}}
vcf2maf      = {{args.vcf2maf | quote}}
vep          = {{args.vep | quote}}
vepDb        = {{args.vepDb | quote}}
filtervcf    = {{args.filtervcf | quote}}
ref          = {{args.ref | quote}}
bcftools     = {{args.bcftools | quote}}
tumoridx     = {{args.tumor | repr}}
nthread      = {{args.nthread | repr}}
params       = {{args.params | repr}}
oncotator    = {{args.oncotator | quote}}
oncotator_db = {{args.oncotator_db | quote}}

shell.load_config(vcf2maf = vcf2maf, bcftools = bcftools, oncotator = oncotator)

veppath   = path.dirname(shell.which(vep).strip())
vcf2maf   = shell.fg.vcf2maf
bcftools  = shell.bcftools

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
	vcf2maf(**params)

def extract_sample_from_vcf(vcf, sample, outvcf):
	bcftools.view(s = sample, o = outvcf + '.tmp', _ = vcf)
	# exclude sites without genotypes
	# because this is split from a merged vcf file
	bcftools.view(_ = outvcf + '.tmp', g = '^miss', o = outvcf)
	shell.rm_rf(outvcf + '.tmp')

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
		shell.out.head(n = 1, _ = path.join(mafdir, "split1.maf")) > outfile
		for i, s in enumerate(vcfsams):
			shell.out.tail(n = '+2', _ = path.join(mafdir, "split{}.maf".format(i+1))) >> outfile

def run_oncotator():
	params._ = [infile, outfile, 'hg19']
	params.v = True
	params.input_format = 'VCF'
	params.output_format = 'TCGAMAF'
	params.log_name = outfile + '.oncotator.log'
	# db_dir
	params['db-dir'] = oncotator_db
	# c
	from glob import glob
	txfiles = glob(path.join(oncotator_db, 'tx_exact_uniprot_matches*'))
	if txfiles:
		params.c = txfiles[0]

	# configs
	outdir = path.dirname(outfile)
	chdir(outdir)
	confs = glob(path.join(oncotator_db, '*.config'))
	for conf in confs:
		shell.ln_s(conf, path.join(outdir, path.basename(conf)))

	# don't use **params, otherwise input_format will be turned into input-format
	logger.info('See %s for oncotator logs.', outfile + '.oncotator.log')
	shell.fg.oncotator(params)

tools = dict(
	vcf2maf   = run_vcf2maf,
	oncotator = run_oncotator
)

try:
	tools[tool]()
except KeyError:
	raise ValueError('Tool {!r} not supported.'.format(tool))