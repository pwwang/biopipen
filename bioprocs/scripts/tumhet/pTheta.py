
# region THetA example
# Generates Sample_T.BEST.results:
# /path/to/theta2/bin/RunTHetA Sample_T.interval_count \
#     --TUMOR_FILE Sample_T.tumor.snp_formatted.txt \
#     --NORMAL_FILE Sample_T.normal.snp_formatted.txt \
#     --BAF --NUM_PROCESSES `nproc` --FORCE
# endregion

# region THetA usage
# usage: RunTHetA.py [-h] [--TUMOR_FILE TUMOR_FILE] [--NORMAL_FILE NORMAL_FILE]
#                    [-n N] [-k K] [-t TAU] [-d DIR] [-p PRE] [-m MAX_NORMAL]
#                    [--NUM_PROCESSES NUM_PROCESSES]
#                    [--NUM_INTERVALS NUM_INTERVALS] [--BOUND_HEURISTIC BH]
#                    [--NORMAL_BOUND_HEURISTIC NBH] [--HEURISTIC_LB LB]
#                    [--HEURISTIC_UB UB] [--BOUNDS_ONLY] [--NO_MULTI_EVENT]
#                    [--RESULTS filename] [--FORCE] [--GET_VALUES]
#                    [--NO_INTERVAL_SELECTION] [--READ_DEPTH_FILE FILENAME]
#                    [--GRAPH_FORMAT GRAPH_FORMAT] [--BAF]
#                    [--RATIO_DEV RATIO_DEV] [--MIN_FRAC MIN_FRAC]
#                    [--NO_CLUSTERING]
#                    QUERY_FILE

# positional arguments:
#   QUERY_FILE            Interval file

# optional arguments:
#   -h, --help            show this help message and exit
#   --TUMOR_FILE TUMOR_FILE
#                         File containing allelic counts for tumor sample SNPs.
#   --NORMAL_FILE NORMAL_FILE
#                         File containing allelic counts for normal samlpe SNPs.
#   -n N, --N N           Number of subpopulations
#   -k K, --MAX_K K       The maximum value allowed for entries in C
#   -t TAU, --TAU TAU     Expected number of copies in normal genome
#   -d DIR, --DIR DIR     Directory where result file is written to
#   -p PRE, --OUTPUT_PREFIX PRE
#                         Prefix for output files created. By default, it will
#                         be the beginning of the input file name (i.e.if input
#                         filename were example.input, output filed would be
#                         example.output and example.withbounds)
#   -m MAX_NORMAL, --MAX_NORMAL MAX_NORMAL
#                         The maximum fraction to consider for normal. Only
#                         enforced for n=2
#   --NUM_PROCESSES NUM_PROCESSES
#                         The number of processes to be used
#   --NUM_INTERVALS NUM_INTERVALS
#                         The maximum number of intervals used by automatic
#                         interval selection.
#   --BOUND_HEURISTIC BH
#   --NORMAL_BOUND_HEURISTIC NBH
#   --HEURISTIC_LB LB
#   --HEURISTIC_UB UB
#   --BOUNDS_ONLY
#   --NO_MULTI_EVENT
#   --RESULTS filename
#   --FORCE
#   --GET_VALUES
#   --NO_INTERVAL_SELECTION
#   --READ_DEPTH_FILE FILENAME
#   --GRAPH_FORMAT GRAPH_FORMAT
#                         Options are .pdf, .jpg, .png, .eps
#   --BAF                 Option to run the BAF model.
#   --RATIO_DEV RATIO_DEV
#                         The deviation away from 1.0 that a ratio must be to be
#                         considered a potential copy number aberration.
#   --MIN_FRAC MIN_FRAC   The minimum fraction of the genome that must have a
#                         potential copy number aberration to be a valid sample
#                         for THetA analysis.
#   --NO_CLUSTERING       Option to run THetA without clustering.
# endregion

# region bam_readcount
# Usage: bam-readcount [OPTIONS] <bam_file> [region]
# Generate metrics for bam_file at single nucleotide positions.
# Example: bam-readcount -f ref.fa some.bam

# Available options:
#   -h [ --help ]                         produce this message
#   -v [ --version ]                      output the version number
#   -q [ --min-mapping-quality ] arg (=0) minimum mapping quality of reads used
#                                         for counting.
#   -b [ --min-base-quality ] arg (=0)    minimum base quality at a position to
#                                         use the read for counting.
#   -d [ --max-count ] arg (=10000000)    max depth to avoid excessive memory
#                                         usage.
#   -l [ --site-list ] arg                file containing a list of regions to
#                                         report readcounts within.
#   -f [ --reference-fasta ] arg          reference sequence in the fasta format.
#   -D [ --print-individual-mapq ] arg    report the mapping qualities as a comma
#                                         separated list.
#   -p [ --per-library ]                  report results by library.
#   -w [ --max-warnings ] arg             maximum number of warnings of each type
#                                         to emit. -1 gives an unlimited number.
#   -i [ --insertion-centric ]            generate indel centric readcounts.
#                                         Reads containing insertions will not be
#                                         included in per-base counts
# endregion

# region snp formted file
# (1) Chrm
# (2) pos
# (3) A -  # reads with A at this position
# (4) C - # reads with C at this position
# (5) G - # reads with G at this position
# (6) T - # reads with T at this position
# (7) Total - the # reads at this position
# (8) refCount - # reads with the reference allele (indicated in SNP file)
# (9) mutCount - # reads with the mutant allele (indicated in SNP file)
# endregion

from pyppl import Box
from os import path, makedirs
from bioprocs.utils import logger, shell2 as shell
from bioprocs.utils.tsvio2 import TsvReader, TsvWriter, TsvRecord
from bioprocs.utils.parallel import Parallel, distribute
from bioprocs.utils.reference import bamIndex
from bioprocs.utils.shell2 import ln_s, wc_l

{% python from os import path %}
itvfile  = {{ i.infile | quote}}
prefix   = {{ o.outdir | path.join: fn(i.infile) | quote }}
tumbam   = {{ i.tumbam | quote}}
tumsnp   = {{ o.outdir | path.join: fn(i.tumbam) | @append: '.bamrc' | quote}}
normbam  = {{ i.normbam | quote}}
normsnp  = {{ o.outdir | path.join: fn(i.normbam) | @append: '.bamrc' | quote}}
outfile  = {{ o.outfile | quote}}
outdir   = {{ o.outdir | quote}}
params   = {{ args.params | repr}}
ref      = {{ args.ref | quote}}
bedtools = {{ args.bedtools | quote }}
bamrc    = {{ args.bam_readcount | quote }}
nthread  = {{ args.nthread | repr}}
affysnps = {{ args.affysnps | quote}}
theta    = {{ args.theta | quote}}
samtools = {{ args.samtools | quote}}

shell.load_config(
	theta2        = Box(_exe = theta, _raw = True),
	samtools      = samtools,
	bedtools      = bedtools,
	bam_readcount = bamrc,
)
if tumbam.endswith('.bam'):
	logger.info('- Try to index tumor bam file...')
	bamIndex(tumbam,  samtools = samtools, nthread = nthread)
if normbam.endswith('.bam'):
	logger.info('- Try to index normal bam file...')
	bamIndex(normbam, samtools = samtools, nthread = nthread)

if not path.isfile(affysnps) and (tumbam.endswith('.bam') or normbam.endswith('.bam')):
	raise ValueError('Formatted SNPs for tumor or normal not provided, so snp list is required.')

if itvfile.endswith('.bed') and not (tumbam.endswith('.bam') and normbam.endswith('.bam')):
	raise ValueError('Bam files of tumor and normal required as interval file is a BED.')

def getCoverage(bamfiles, region, outfile):
	"""
	Get coverage of regions
	"""
	logger.info('- Get read counts from %s ...' % bamfiles)
	if not region.endswith('.bed'): # read counts already provided
		logger.info('  Read counts already provided in region file, skip.')
		ln_s(bamfile, outfile)
		return

	writer = TsvWriter(outfile)
	writer.cnames = ['ID', 'Chrm', 'Start', 'End', 'numTumor', 'numNormal']
	writer.writeHead(lambda ns: '#' + '\t'.join(ns))
	i = 0
	for line in shell.bedtools.multicov(
		bams  = bamfiles,
		bed   = region,
		_iter = True):
		# chr1	1	12444	104	87
		i += 1
		parts = line.strip().split('\t')
		parts.insert(0, i)
		writer.write(parts[:4] + parts[-2:])
	writer.close()

def getAlleleCount(bamfile, snpfile, outfile):
	logger.info('- Get allele counts from %s ...' % bamfile)
	if not bamfile.endswith('.bam'): # it's THetA2 formatted SNP file already
		logger.info('  THetA2 formatted SNP file already, skip.')
		ln_s(bamfile, outfile)
		return
	brcparams       = Box()
	brcparams.f     = ref
	brcparams.w     = 0
	brcparams.l     = snpfile
	brcparams._iter = True

	brcparams[''] = bamfile

	writer = TsvWriter(outfile)
	writer.cnames = ['Chrm', 'pos', 'A', 'C', 'G', 'T', 'Total', 'refCount', 'mutCount']
	writer.writeHead(lambda ns: '#' + '\t'.join(ns))
	for line in shell.bam_readcount(**brcparams):
		#chr1	564773	C	14	=:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00	A:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00	C:14:...	G:0:...	T:0:...	N:0:...
		r = line.split('\t')
		counts = dict(
			A = r[5].split(':', 2)[1],
			C = r[6].split(':', 2)[1],
			G = r[7].split(':', 2)[1],
			T = r[8].split(':', 2)[1]
		)
		rec    = TsvRecord()
		rec.Chrm  = r[0]
		rec.pos   = r[1]
		rec.Total = r[3]
		rec.A = counts['A']
		rec.C = counts['C']
		rec.G = counts['G']
		rec.T = counts['T']
		# if reference allele is unknown, assuming all are mut alleles
		refallele = r[2].upper()
		rec.refCount = counts.get(refallele, 0)
		counts = [int(count) for allele, count in counts.items() if allele != refallele and count.isdigit()]
		rec.mutCount = max(counts) if counts else 0
		writer.write(rec)
	writer.close()

if nthread == 1:
	getAlleleCount(tumbam, affysnps, tumsnp)
	getAlleleCount(normbam, affysnps, normsnp)
else:
	# try to split the affysnps into N files and distribute the jobs to nthreads
	# get number of lines of affysnps file
	total  = int(wc_l(affysnps).split()[0])
	dists  = distribute(total, nthread)
	reader = TsvReader(affysnps, cnames = False)
	# dir to save the split file and result file
	thdir  = path.join(outdir, 'bamrc.nthreads')
	if not path.exists(thdir):
		makedirs(thdir)

	asbname = path.basename(affysnps).split('.')[0]
	for i, dist in enumerate(dists):
		writer = TsvWriter(path.join(thdir, '{bname}.thread{i}.snp'.format(
			bname = asbname, i = i
		)))
		for _ in range(dist):
			writer.write(next(reader))
		writer.close()

	para   = Parallel(nthread, raiseExc = True)
	para.run(getAlleleCount, [
		(tumbam, path.join(
			thdir, '{bname}.thread{i}.snp'.format(bname = asbname, i = i)
		), path.join(
			thdir, '{tumbn}.thread{i}.bamrc'.format(tumbn = path.basename(tumbam), i = i)
		)) for i in range(nthread)
	])
	# merge to tumsnp
	writer = TsvWriter(tumsnp)
	writer.cnames = ['Chrm', 'pos', 'A', 'C', 'G', 'T', 'Total', 'refCount', 'mutCount']
	writer.writeHead(lambda cn: "#" + "\t".join(cn))
	for i in range(nthread):
		subrc = path.join(
			thdir, '{tumbn}.thread{i}.bamrc'.format(tumbn = path.basename(tumbam), i = i)
		)
		reader = TsvReader(subrc, cnames = False)
		for r in reader:
			writer.write(r.values())
		reader.close()
	writer.close()

	# normal
	para.run(getAlleleCount, [
		(normbam, path.join(
			thdir, '{bname}.thread{i}.snp'.format(bname = asbname, i = i)
		), path.join(
			thdir, '{normbn}.thread{i}.bamrc'.format(normbn = path.basename(normbam), i = i)
		)) for i in range(nthread)
	])
	# merge to normsnp
	writer = TsvWriter(normsnp)
	writer.cnames = ['Chrm', 'pos', 'A', 'C', 'G', 'T', 'Total', 'refCount', 'mutCount']
	writer.writeHead(lambda cn: "#" + "\t".join(cn))
	for i in range(nthread):
		subrc = path.join(
			thdir, '{normbn}.thread{i}.bamrc'.format(normbn = path.basename(normbam), i = i)
		)
		reader = TsvReader(subrc, cnames = False)
		for r in reader:
			writer.write(r.values())
		reader.close()
	writer.close()

# check type of itvfile
if itvfile.endswith('.bed'):
	tmpfile = path.join(outdir, path.basename(itvfile).split('.')[0] + '.interval.txt')
	getCoverage([tumbam, normbam], itvfile, tmpfile)
	itvfile = tmpfile


params.NUM_PROCESSES = nthread
params.TUMOR_FILE    = tumsnp
params.NORMAL_FILE   = normsnp
params.GRAPH_FORMAT  = '.png'
params.OUTPUT_PREFIX = prefix
#params.RATIO_DEV     = params.get('RATIO_DEV', .5)
#params.MIN_FRAC      = params.get('MIN_FRAC', .1)
params['']           = itvfile
shell.fg.theta2(**params)

ln_s(prefix + '.BEST.results', outfile)

