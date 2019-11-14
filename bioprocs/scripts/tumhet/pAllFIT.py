"""
usage: All-FIT.py [-h] -i INPUTDIRFILE -d OUTPUTDIR -o OUTPUTNAME
				  [-s STANDARDDEVIATION] [-t {somatic,all}]

This is a program that estimates specimen purity from tumor-only sample
sequenced with deep sequencing, called All-FIT (Allele Frequency based
Imputation of Tumor Purity). We do not provide model for a variant that is
being called on chrX (male) or chrY. Users may consider the option of removing
all variants on chrX (male) and on chrY, if there is any in the input sample.
Input sample should be a tab-delimited file with 4 columns and a header of ID
Allele_Freq Depth Ploidy

ID	Allele_Freq	Depth	Ploidy
Var1	39.0	378	2
Var2	19.0	996	2
Var3	64.0	849	1
Var4	77.0	487	2
Var5	49.0	639	2
Var6	52.0	891	2
Var7	52.0	841	2

optional arguments:
  -h, --help            show this help message and exit
  -i INPUTDIRFILE, --inputDirFile INPUTDIRFILE
						Input file name with path
  -d OUTPUTDIR, --outputDir OUTPUTDIR
						Output directory
  -o OUTPUTNAME, --outputName OUTPUTNAME
						Output file name prefix
  -s STANDARDDEVIATION, --standardDeviation STANDARDDEVIATION
						How many standard deviation or confidence interval of
						estimated purity
  -t {somatic,all}, --typeOfInput {somatic,all}
						Types of variants input whether germline variants are
						removed(somatic) or not(all)
"""
from os import path
from pyppl import Box
from bioprocs.utils import logger, shell2 as shell
from bioprocs.utils.tsvio2 import TsvReader, TsvWriter, TsvRecord

infile   = {{i.infile | quote}}
cnfile   = {{i.cnfile | quote}}
outstem  = {{o.outfile | stem | quote}}
outdir   = {{o.outdir | quote}}
allfit   = {{args.allfit | quote}}
bcftools = {{args.bcftools | quote}}
bedtools = {{args.bedtools | quote}}
params   = {{args.params | repr}}
mutctrl  = {{args.mutctrl | repr}}
cnctrl   = {{args.cnctrl | repr}}
nthread  = {{args.nthread | repr}}

shell.load_config(allfit = dict(
	_exe = allfit,
	_env = dict(
		OPENBLAS_NUM_THREADS = str(nthread),
		OMP_NUM_THREADS      = str(nthread),
		NUMEXPR_NUM_THREADS  = str(nthread),
		MKL_NUM_THREADS      = str(nthread)
	)),
	bcftools = bcftools, bedtools = bedtools)

def MutVCF2BED(vcffile):
	"""
	Convert single(paired)-sample VCFs to a BED file like:
	#CHROM  START   END	    NAME      AF  Depth GT
	chr1    70820   70820   chr:70820 .39 1083  1
	"""
	ret = None
	samples = shell.bcftools.query(l = vcffile).strip().splitlines()
	assert len(samples) <= 2
	samidx  = 0 if len(samples) == 1 or mutctrl in (1, -1, None) else 1
	sample  = samples[samidx]
	ret = path.join(outdir, sample + '.muts.bed')
	writer = TsvWriter(ret)
	writer.cnames = ['CHROM', 'START', 'END', 'NAME', 'AF', 'Depth', 'CN']
	for line in shell.bcftools.query(
		_iter = True,
		_ = vcffile,
		s = sample,
		f = '%CHROM\t%POS\t%POS\t%CHROM:%POS\t[%AF]\t[%AD{1}]\t[%GT]\n'):
		items = line.strip().split('\t')
		try:
			gt = int(items[-1][0]) + int(items[-1][-1])
		except (TypeError, ValueError):
			logger.warning('Variant %s dropped failure fetching genotype' % (items[3]))
			continue
		items[-1] = gt
		writer.write(items)
	writer.close()
	return ret

def CnVCF2BED(vcffile):
	"""
	Convert single(paired)-sample VCFs to a BED file like:
	#CHROM  START   END	    CN
	chr1    70820   70820   5
	"""

	ret = {}
	samples = shell.bcftools.query(l = vcffile).strip().splitlines()
	assert len(samples) <= 2
	samidx  = 0 if len(samples) == 1 or cnctrl in (1, -1, None) else 1
	sample  = samples[samidx]
	ret = path.join(outdir, sample + '.cn.bed')
	writer = TsvWriter(ret)
	writer.cnames = ['CHROM', 'START', 'END', 'CN']
	for line in shell.bcftools.query(
		_iter = True,
		_ = vcffile,
		s = sample,
		f = '%CHROM\t%POS\t[%END]\t[%CN]\n'):
		items = line.strip().split('\t')
		if items[3] == '0':
			logger.warning('Record skipped with CN=0: %s', line)
			continue
		writer.write(items)
	writer.close()
	return ret

if infile.endswith('.vcf') or infile.endswith('.vcf.gz'):
	mutbed = MutVCF2BED(infile)
	if cnfile.endswith('.vcf') or cnfile.endswith('.vcf.gz'):
		cnbed = CnVCF2BED(cnfile)
	elif cnfile.endswith('.bed'):
		cnbed = cnfile
	else:
		raise ValueError('Unsupported CNV file.')

	infile = path.join(outdir, path.splitext(mutbed)[0] + '.allfit.txt')
	writer = TsvWriter(infile)
	writer.cnames = ['ID', 'Allele_Freq', 'Depth', 'Ploidy']
	writer.writeHead()
	for line in shell.bedtools.intersect(
		_iter = True,
		a = mutbed,
		b = cnbed,
		loj = True):
		# CHROM  START END NAME AF   Depth       GT
		# 0      1     2   3    4        5        6
		# CHROM START END CN
		# 7     8     9   10
		parts = line.strip().split('\t')
		writer.write([
			parts[3],
			float(parts[4]) * 100,
			parts[5],
			2 if parts[10] == '.' else parts[10]])
	writer.close()

elif cnfile:
	logger.warning('All-FIT format detected, cnfile will be ignored!')

params.i = infile
params.d = outdir
params.o = outstem

shell.fg.allfit(**params)

