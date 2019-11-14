from os import path, chdir
from glob import glob
from pyppl import Box
from bioprocs.utils import parallel, logger, shell2 as shell
from bioprocs.utils.tsvio2 import TsvWriter, TsvReader
from bioprocs.utils.reference import vcfIndex

infile       = {{i.infile | quote}}
outfile      = {{o.outfile | quote}}
tool         = {{args.tool | quote}}
vcf2maf      = {{args.vcf2maf | quote}}
vep          = {{args.vep | quote}}
tabix        = {{args.tabix | quote}}
vepDb        = {{args.vepDb | quote}}
filtervcf    = {{args.filtervcf | quote}}
ref          = {{args.ref | quote}}
bcftools     = {{args.bcftools | quote}}
genome       = {{args.genome | quote}}
tumoridx     = {{args.tumor | repr}}
withchr      = {{args.withchr | repr}}
nthread      = {{args.nthread | repr}}
params       = {{args.params | repr}}
oncotator    = {{args.oncotator | quote}}
oncotator_db = {{args.oncotator_db | quote}}

shell.load_config(vcf2maf = vcf2maf, bcftools = bcftools, oncotator = oncotator)

def extract_sample_from_vcf(vcf, sample, outvcf, nomiss = True):
	shell.fg.bcftools.view(s = sample, o = outvcf + '.tmp', _ = vcf)
	# exclude sites without genotypes
	# because this is split from a merged vcf file
	if nomiss:
		shell.fg.bcftools.view(_ = outvcf + '.tmp', g = '^miss', o = outvcf)
	else:
		shell.fg.bcftools.view(_ = outvcf + '.tmp', o = outvcf)
	shell.rm_rf(outvcf + '.tmp')

def run_vcf2maf_one(vcf, maf, tumor, normal = None, forks = nthread):
	params_one = params.copy()
	params_one['input-vcf']  = vcf
	params_one['output-maf'] = maf
	params_one['vep-data']   = vepDb
	params_one['vep-forks']  = forks
	params_one['filter-vcf'] = filtervcf
	params_one['ref-fasta']  = ref
	params_one['vep-path']   = path.dirname(shell.which(vep).strip())
	params_one['tumor-id']   = tumor
	params_one['normal-id']  = normal
	shell.fg.vcf2maf(**params_one)

def run_oncotator_one(vcf, maf, tumor, normal = None, forks = nthread):
	vcf = vcfIndex(vcf, tabix = tabix)
	vcf0 = vcf
	if normal: # split tumor and normal
		vcf_tumor  = vcf + '.' + tumor
		extract_sample_from_vcf(vcf, tumor, vcf_tumor, nomiss = False)
		vcf = vcf_tumor

	openblas_threads = max(nthread - forks, 1)
	params_one               = params.copy()
	params_one._             = [vcf, maf, 'hg19']
	params_one.v             = True
	params_one.input_format  = 'VCF'
	params_one.output_format = 'TCGAMAF'
	params_one.log_name      = maf + '.oncotator.log'
	# db_dir
	params_one['db-dir'] = oncotator_db
	# c
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
	logger.info('See %s for oncotator logs.', params_one.log_name)
	shell.fg.oncotator(params_one, _env =  dict(
		OPENBLAS_NUM_THREADS = str(openblas_threads),
		OMP_NUM_THREADS      = str(openblas_threads),
		NUMEXPR_NUM_THREADS  = str(openblas_threads),
		MKL_NUM_THREADS      = str(openblas_threads)
	))

	# oncotator cannot put Matched_Norm_Sample_Barcode and allele information
	# add them if normal is specified
	if normal:
		# I suppose there isn't many somatic mutations, so I read the normal information first
		# into memory
		from pysam import VariantFile
		normal_infos = {}
		vfile        = VariantFile(vcf0)
		tum_index    = list(vfile.header.samples).index(tumor)
		norm_index   = list(vfile.header.samples).index(normal)
		for record in vfile:
			tum_alt_index1  = record.samples[tum_index]['GT'][0]
			tum_alt_index2  = record.samples[tum_index]['GT'][-1]
			tum_alt1        = record.alleles[tum_alt_index1]
			tum_alt2        = record.alleles[tum_alt_index2]
			chrom = record.chrom[3:] if record.chrom.startswith('chr') else record.chrom

			key = '{3}_{0.pos}_{0.ref}_{1}_{2}'.format(
				record, tum_alt1, tum_alt2, chrom)
			try:
				norm_alt_index1 = record.samples[norm_index]['GT'][0]
				norm_alt_index2 = record.samples[norm_index]['GT'][-1]
				normal_infos[key] = (record.alleles[norm_alt_index1], record.alleles[norm_alt_index2])
			except (TypeError, IndexError):
				pass

	reader = TsvReader(maf, cnames = True)
	writer = TsvWriter(maf + '.tmp')
	writer.cnames = reader.cnames

	# correct Start_position to Start_Position, End_position to End_Position
	#               ^                               ^
	writer.writeHead(lambda cnames: [
		'Start_Position' if cname == 'Start_position' else
		'End_Position' if cname == 'End_position' else cname
		for cname in cnames])
	for r in reader:
		key = '{r.Chromosome}_{r.Start_position}_{r.Reference_Allele}_{r.Tumor_Seq_Allele1}_{r.Tumor_Seq_Allele2}'.format(r = r)
		r.NCBI_Build = genome if r.NCBI_Build == '__UNKNOWN__' else r.NCBI_Build
		r.Tumor_Sample_Barcode = tumor
		r.Matched_Norm_Sample_Barcode = normal or 'NORMAL'
		if normal:
			try:
				r.Match_Norm_Seq_Allele1 = normal_infos[key][0]
				r.Match_Norm_Seq_Allele2 = normal_infos[key][1]
			except KeyError:
				r.Match_Norm_Seq_Allele1 = r.Reference_Allele
				r.Match_Norm_Seq_Allele2 = r.Reference_Allele
		if r.Variant_Type == 'INS' and r.Match_Norm_Seq_Allele1 == '-' and r.Tumor_Seq_Allele2:
			r.Match_Norm_Seq_Allele1 = r.Tumor_Seq_Allele2[0]
		if r.Variant_Type == 'INS' and r.Match_Norm_Seq_Allele2 == '-' and r.Tumor_Seq_Allele2:
			r.Match_Norm_Seq_Allele2 = r.Tumor_Seq_Allele2[0]
		if withchr and r.Chromosome[:3] != 'chr':
			r.Chromosome = 'chr' + r.Chromosome
		writer.write(r)

	shell.mv(maf + '.tmp', maf)

ones = dict(
	vcf2maf   = run_vcf2maf_one,
	oncotator = run_oncotator_one
)

def run(tool):
	global tumoridx

	one = ones[tool]

	vcfsams  = shell.bcftools.query(l = infile).splitlines()
	vcfsams  = [s for s in vcfsams if s.strip()]

	if len(vcfsams) == 1:
		one(infile, outfile, vcfsams[0])
	elif len(vcfsams) >= 2 and tumoridx in (0, 1, 'auto'):
		if len(vcfsams) == 2 and tumoridx == 'auto':
			from difflib import SequenceMatcher
			similarity = lambda a, b: SequenceMatcher(None, a, b).ratio()
			sim0       = similarity(path.basename(infile), vcfsams[0])
			sim1       = similarity(path.basename(infile), vcfsams[1])
			tumoridx   = int(sim1 > sim0)
			logger.info('Using sample %s(%s) as tumor', tumoridx + 1, vcfsams[tumoridx])

		tumor  = vcfsams[tumoridx]
		normal = vcfsams[1-tumoridx]
		one(infile, outfile, tumor, normal)
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
		restThreads = int(float(nthread)/float(len(vcfsams)))
		restThreads = max(restThreads, 1)
		para.run(one, [
			(path.join(splitdir, "split" + str(i+1) + ".vcf"),
			 path.join(mafdir, "split{}.maf".format(i+1)),
			 s, None, restThreads)
			for i, s in enumerate(vcfsams)
		])
		del para

		# merge mafs
		shell.out.head(n = 1, _ = path.join(mafdir, "split1.maf")) > outfile
		for i, s in enumerate(vcfsams):
			shell.out.tail(n = '+2', _ = path.join(mafdir, "split{}.maf".format(i+1))) >> outfile

try:
	run(tool)
except KeyError:
	raise ValueError('Tool {!r} not supported.'.format(tool))
