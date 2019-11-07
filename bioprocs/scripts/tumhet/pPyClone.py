import cmdy
from pysam import VariantFile
from pyppl import Box
from os import path, environ
from bioprocs.utils import logger, shell2 as shell
from bioprocs.utils.tsvio2 import TsvReader, TsvWriter, TsvRecord

inmuts   = {{i.muts | quote}}
incnvs   = {{i.cnvs | quote}}
outdir   = {{o.outdir | quote}}
pyclone  = {{args.pyclone | quote}}
bcftools = {{args.bcftools | quote}}
bedtools = {{args.bedtools | quote}}
params   = {{args.params | repr}}
nthread  = {{args.nthread | int }}
refgene  = {{args.refgene | quote }}
mutctrl  = {{args.mutctrl | repr}}
cnctrl   = {{args.cnctrl | repr}}

env = environ.copy()
env.update(dict(
	OPENBLAS_NUM_THREADS = str(nthread),
	OMP_NUM_THREADS      = str(nthread),
	NUMEXPR_NUM_THREADS  = str(nthread),
	MKL_NUM_THREADS      = str(nthread)
) if nthread else {})

shell.load_config(pyclone = dict(
	_exe = pyclone,
	_cwd = outdir,
	_env = env
), bcftools = bcftools, bedtools = bedtools)

def defaultCNs(chrom, gt):
	# normal_cn, minor_cn, major_cn
	ret = [2, 0, 2]
	if gt == 'AA':
		ret[0] = 1 if chrom == 'chrY' else 2 # sex information needed for chrX
		ret[1] = 0
		ret[2] = ret[0]
	elif gt == 'AB':
		ret[0] = 2
		ret[1] = 1
		ret[2] = 1
	else: # BB
		ret[0] = 1 if chrom == 'chrY' else 2
		ret[1] = ret[0]
		ret[2] = 0
	return ret

def singleSampleMutVCF2BED(vcffiles):
	"""
	Convert single(paired)-sample VCFs to a BED file like:
	#CHROM  START   END	    NAME      GENOTYPE REFCOUNTS VARCOUNTS CASE    VARFREQ
	chr1    70820   70820   chr:70820 BB       1083      1996      NA12156 0.648
	"""
	ret = {}
	for vcffile in vcffiles:
		samples = shell.bcftools.query(l = vcffile).strip().splitlines()
		assert len(samples) <= 2
		samidx  = 0 if len(samples) == 1 or mutctrl in (1, -1, None) else 1
		sample  = samples[samidx]
		ret[sample] = path.join(outdir, sample + '.muts.bed')
		writer = TsvWriter(ret[sample])
		writer.cnames = ['CHROM', 'START', 'END', 'NAME', 'GENOTYPE', 'REFCOUNTS', 'VARCOUNTS', 'CASE', 'VARFREQ']
		for line in shell.bcftools.query(
			_ = vcffile,
			s = sample,
			f = '%CHROM\t%POS\t%POS\t%CHROM:%POS\t%GT\t%FORMAT/AD{0}\t%FORMAT/AD{1}\t%SAMPLE\t%FORMAT/AF\n'):

			items = line.split('\t')
			if items[4] in ('0/0', '0|0'):
				items[4] = 'AA'
			elif items[4] in ('0/1', '0|1'):
				items[4] = 'AB'
			elif items[4] in ('1/1', '1|1'):
				items[4] = 'BB'
			writer.write(items)

		writer.close()
	return ret

def multiSampleMutVCF2BED(vcffile):
	ret = {}
	writers = {}
	samples = shell.bcftools.query(l = vcffile).strip().splitlines()
	samples = [sample for i, sample in samples
		if mutctrl is None or i != (mutctrl if mutctrl >= 0 else len(samples) + mutctrl)]
	for sample in samples:
		ret[sample] = path.join(outdir, sample + '.muts.bed')
		writer = TsvWriter(sample)
		writer.cnames = ['CHROM', 'START', 'END', 'NAME', 'GENOTYPE', 'REFCOUNTS', 'VARCOUNTS', 'CASE', 'VARFREQ']
		writers[sample] = writer
	for line in shell.bcftools.query(
		_ = vcffile,
		s = ','.join(samples),
		f = '[%CHROM\t%POS\t%POS\t%CHROM:%POS\t%GT\t%FORMAT/AD{0}\t%FORMAT/AD{1}\t%SAMPLE\t%FORMAT/AF\n]'):
		items = line.split('\t')

		writer = writers[items[7]]
		if items[4] in ('0/0', '0|0'):
			items[4] = 'AA'
		elif items[4] in ('0/1', '0|1'):
			items[4] = 'AB'
		elif items[4] in ('1/1', '1|1'):
			items[4] = 'BB'
		writer.write(items)

	for writer in writers.values():
		writer.close()
	return ret

def MAF2BED(maffile):

	reader = TsvReader(maffile)
	if 't_alt_count' not in reader.cnames:
		raise ValueError('t_alt_count not found in MAF file.')
	if 't_ref_count' not in reader.cnames:
		raise ValueError('t_ref_count not found in MAF file.')

	ret = {}
	writers = {}
	for r in reader:
		r.CHROM = r.Chromosome
		r.START = r.Start_Position
		r.END   = r.End_Position
		r.NAME  = '%s:%s' % (r.CHROM, r.START)
		r.GENOTYPE = 'AA' if r.Tumor_Seq_Allele1 == r.Tumor_Seq_Allele2 \
						and r.Tumor_Seq_Allele1 == r.Reference_Allele else \
					 'BB' if r.Tumor_Seq_Allele1 == r.Tumor_Seq_Allele2 \
						and r.Tumor_Seq_Allele1 != r.Reference_Allele else 'AB'
		r.REFCOUNTS = r.t_ref_count
		r.VARCOUNTS = r.t_alt_count
		r.CASE      = r.Tumor_Sample_Barcode
		try:
			varcount  = float(r.VARCOUNTS)
			refcount  = float(r.REFCOUNTS)
			depth     = float(r.get('t_depth', 0))
			if depth == 0:
				depth = varcount + refcount
			r.VARFREQ = varcount / depth
		except (ValueError, TypeError, ZeroDivisionError):
			logger.warning('Variant %s drop due to unknown t_ref_count(%s) or t_alt_count(%s)' % (
				r.NAME,
				r.t_ref_count,
				r.t_alt_count
			))
			continue

		if r.CASE not in ret:
			ret[r.CASE] = path.join(outdir, r.CASE + '.muts.bed')
			writers[r.CASE] = TsvWriter(ret[r.CASE])
			writers[r.CASE].cnames = ['CHROM', 'START', 'END', 'NAME', 'GENOTYPE', 'REFCOUNTS', 'VARCOUNTS', 'CASE', 'VARFREQ']
		writers[r.CASE].write(r)

	for writer in writers.values():
		writer.close()
	reader.close()
	return ret

def PyCloneMutTSV2BED(pcfile):
	ret = {}
	writers = {}
	reader = TsvReader(pcfile)
	for r in reader:
		if ':' not in r.mutation_id:
			raise "`mutation_id` should end with `<chr>:<pos>`"
		chr, pos = r.mutation_id.split(':')[:-2]

		r.CHROM     = chr
		r.START     = pos
		r.END       = pos
		r.NAME      = '%s:%s' % (chr, pos)
		r.GENOTYPE  = r.genotype
		r.REFCOUNTS = r.ref_counts
		r.VARCOUNTS = r.var_counts
		r.CASE      = r.variant_case
		r.VARFREQ   = r.variant_freq

		if r.CASE not in ret:
			ret[r.CASE] = path.join(outdir, r.CASE + '.muts.bed')
			writers[r.CASE] = TsvWriter(ret[r.CASE])
			writers[r.CASE].cnames = ['CHROM', 'START', 'END', 'NAME', 'GENOTYPE', 'REFCOUNTS', 'VARCOUNTS', 'CASE', 'VARFREQ']
		writers[r.CASE].write(r)
	for writer in writers.values():
		writer.close()
	reader.close()
	return ret

def singleSampleCnVCF2BED(vcffiles):
	"""
	Convert single(paired)-sample VCFs to a BED file like:
	#CHROM  START   END	    GENOTYPE NORMAL_CN MINOR_CN MAJOR_CN CASE
	chr1    70820   70820   BB       2         0        2        NA12156
	"""

	ret = {}
	for vcffile in vcffiles:
		samples = shell.bcftools.query(l = vcffile).strip().splitlines()
		assert len(samples) <= 2
		samidx  = 0 if len(samples) == 1 or cnctrl in (1, -1, None) else 1
		sample  = samples[samidx]
		ret[sample] = path.join(outdir, sample + '.cn.bed')
		writer = TsvWriter(ret[sample])
		writer.cnames = ['CHROM', 'START', 'END', 'GENOTYPE', 'NORMAL_CN', 'MINOR_CN', 'MAJOR_CN', 'CASE']
		for line in shell.bcftools.query(
			_ = vcffile,
			s = sample,
			f = '%CHROM\t%POS\t%POS\t%GT\t2\t0\t2\t%SAMPLE\n'):

			items = line.split('\t')
			if items[3] in ('0/0', '0|0'):
				items[3] = 'AA'
			elif items[3] in ('0/1', '0|1'):
				items[3] = 'AB'
			elif items[3] in ('1/1', '1|1'):
				items[3] = 'BB'

			items[4], items[5], items[6] = defaultCNs(items[0], items[3])

			writer.write(items)
		writer.close()
	return ret

def multiSampleCnVCF2BED(vcffile):
	ret = {}
	writers = {}
	samples = shell.bcftools.query(l = vcffile).strip().splitlines()
	samples = [sample for i, sample in samples
		if cnctrl is None or i != (cnctrl if cnctrl >= 0 else len(samples) + cnctrl)]
	for sample in samples:
		ret[sample] = path.join(outdir, sample + '.cn.bed')
		writer = TsvWriter(sample)
		writer.cnames = ['CHROM', 'START', 'END', 'GENOTYPE', 'NORMAL_CN', 'MINOR_CN', 'MAJOR_CN', 'CASE']
		writers[sample] = writer
	for line in shell.bcftools.query(
		_ = vcffile,
		s = ','.join(samples),
		f = '[%CHROM\t%POS\t%POS\t%GT\t2\t0\t2\t%SAMPLE\n]'):

		items = line.split('\t')
		if items[3] in ('0/0', '0|0'):
			items[3] = 'AA'
		elif items[3] in ('0/1', '0|1'):
			items[3] = 'AB'
		elif items[3] in ('1/1', '1|1'):
			items[3] = 'BB'

		items[4], items[5], items[6] = defaultCNs(items[0], items[3])
		writer = writers[items[7]]
		writer.write(items)
	for writer in writers.values():
		writer.close()
	return ret

def PyCloneCnTSV2BED(pcfile):
	reader = TsvReader(pcfile)
	ret = {}
	writers = {}
	for r in reader:
		if ':' not in r.mutation_id:
			raise "`mutation_id` should end with `<chr>:<pos>`"
		chr, pos = r.mutation_id.split(':')[:-2]

		r.CHROM     = chr
		r.START     = pos
		r.END       = pos
		r.GENOTYPE  = r.genotype
		r.CASE      = r.variant_case
		r.NORMAL_CN = r.normal_cn
		r.MINOR_CN  = r.minor_cn
		r.MAJOR_CN  = r.major_cn
		if r.CASE not in ret:
			ret[r.CASE] = path.join(outdir, r.CASE + '.cn.bed')
			writers[r.CASE] = TsvWriter(ret[r.CASE])
			writers[r.CASE].cnames = ['CHROM', 'START', 'END', 'GENOTYPE', 'NORMAL_CN', 'MINOR_CN', 'MAJOR_CN', 'CASE']
		writers[r.CASE].write(r)
	for writer in writers.values():
		writer.close()
	reader.close()
	return ret

def GetPyCloneTsv(mutfile, outfile, cnfile = None):
	# mutbed:
	# #CHROM  START   END	    NAME      GENOTYPE REFCOUNTS VARCOUNTS CASE    VARFREQ
	# cnbed:
	# #CHROM  START   END	    GENOTYPE NORMAL_CN MINOR_CN MAJOR_CN CASE
	# outfile:
	# mutation_id	ref_counts	var_counts	normal_cn	minor_cn	major_cn	variant_case	variant_freq	genotype
	writer = TsvWriter(outfile)
	writer.cnames = ['mutation_id', 'ref_counts', 'var_counts', 'normal_cn', 'minor_cn', 'major_cn', 'variant_case', 'variant_freq', 'genotype']
	writer.writeHead()
	if not cnfile:
		reader = TsvReader(mutbed)
		reader.cnames = ['CHROM', 'START', 'END', 'NAME', 'GENOTYPE', 'REFCOUNTS', 'VARCOUNTS', 'CASE', 'VARFREQ']
		for r in reader:
			r.mutation_id = r.NAME
			r.ref_counts = r.REFCOUNTS
			r.var_counts = r.VARCOUNTS
			r.normal_cn, r.minor_cn, r.major_cn = defaultCNs(r.CHROM, r.GENOTYPE)
			r.variant_case = r.CASE
			r.variant_freq = r.VARFREQ
			r.genotype = r.GENOTYPE
			writer.write(r)
		writer.close()
		return

	for line in shell.bedtools.intersect(
		a = mutfile,
		b = cnfile,
		loj = True):

		# CHROM  START END NAME GENOTYPE REFCOUNTS VARCOUNTS CASE VARFREQ
		# 0      1     2   3    4        5         6         7    8
		# CHROM START END GENOTYPE NORMAL_CN MINOR_CN MAJOR_CN CASE
		# 9     10    11  12       13        14       15       16
		parts = line.split('\t')

		rec = TsvRecord()
		rec.mutation_id = parts[3]
		rec.ref_counts = parts[5]
		rec.var_counts = parts[6]
		rec.variant_case = parts[7]
		rec.variant_freq = parts[8]
		rec.genotype = parts[4]
		if parts[9] == '.':
			rec.normal_cn, rec.minor_cn, rec.major_cn = defaultCNs(parts[0], parts[4])
		else:
			rec.normal_cn, rec.minor_cn, rec.major_cn = parts[13:16]
		writer.write(rec)
	writer.close()

params  = Box(in_files = [])
mutbeds = {}
cnbeds  = {}
if path.isfile(inmuts):
	if inmuts.endswith('.vcf') or inmuts.endswith('.vcf.gz'):
		mutbeds = multiSampleMutVCF2BED(inmuts)
	elif inmuts.endswith('.maf') or inmuts.endswith('.maf.gz'):
		mutbeds = MAF2BED(inmuts)
	else:
		mutbeds = PyCloneMutTSV2BED(inmuts)
else:
	inmuts = inmuts.split(',')
	mutbeds = singleSampleMutVCF2BED(inmuts)

if incnvs and path.isfile(incnvs):
	if incnvs.endswith('.vcf') or incnvs.endswith('.vcf.gz'):
		cnbeds = multiSampleCnVCF2BED(incnvs)
	else:
		cnbeds = PyCloneCnTSV2BED(incnvs)
elif incnvs:
	incnvs = incnvs.split(',')
	cnbeds = singleSampleCnVCF2BED(incnvs)

for sample, mutbed in mutbeds.items():
	pcfile = path.join(outdir, sample + '.tsv')
	params.in_files.append(pcfile)
	GetPyCloneTsv(mutbed, pcfile, cnbeds.get(sample))

#PyClone run_analysis_pipeline --in_files A.tsv B.tsv C.tsv --working_dir pyclone_analysis
# create matplotlibrc file
with open(path.join(outdir, 'matplotlibrc'), 'w') as f:
	f.write('backend: Agg\n')
params.working_dir      = outdir
params.prior            = 'total_copy_number'
params.plot_file_format = 'svg'
params._env             = env
params._raise = False

c = shell.fg.pyclone.run_analysis_pipeline(**params)
if c.rc != 0:
	# Let it go if 'No mutations found in common across samples'
	exit(0)

# annotate tables/loci.tsv with genes
"""
mutation_id	sample_id	cluster_id	cellular_prevalence	cellular_prevalence_std	variant_allele_frequency
chr10:102776870	AK_nevus_1	0	0.656	0.127	0.417
chr10:102776870	AK_nevus_2	0	0.432	0.126	0.333
"""
locifile = path.join(outdir, 'tables', 'loci.tsv')
locifilebed = path.join(outdir, 'tables', 'loci.bed')
reader = TsvReader(locifile)
writer = TsvWriter(locifilebed)
writer.cnames = ['CHR', 'START', 'END'] + reader.cnames
for r in reader:
	r.CHR, r.START = r.mutation_id.split(':')
	r.END = r.START
	writer.write(r)
reader.close()
writer.close()

writer = TsvWriter(locifile)
writer.cnames = reader.cnames + ['gene']
writer.writeHead()
for line in shell.bedtools.intersect(a = locifilebed, b = refgene, loj = True, _iter = True):
	parts = line.split('\t')
	r = parts[3:9] + [parts[17].split('; ')[0][9:-1]]
	writer.write(r)
writer.close()


