from os import path, makedirs
from shutil import rmtree
from sys import stderr
from pyppl import Box
from bioprocs.utils import cmdargs, runcmd, mem2
from bioprocs.utils.reference import bamIndex

{% python from os import path %}

ref       = {{args.ref | quote}}
jobindex  = {{job.index | repr}}
tmpdir    = {{args.tmpdir | quote}}
procid    = {{proc.id | quote}}
tumor     = {{i.tumor | quote}}
tprefix   = {{i.tumor | fn2 | quote}}
normal    = {{i.normal | quote}}
nprefix   = {{i.normal | fn2 | quote}}
tool      = {{args.tool | quote}}
gz        = {{args.gz | bool}}
outfile   = {{o.outfile | quote}}
params    = {{args.params | repr}}
samtools  = {{args.samtools | quote}}
mem       = {{args.mem | quote}}
nthread   = {{args.nthread | repr}}
gatk      = {{args.gatk | quote}}
ssniper   = {{args.somaticsniper | quote}}
ssniffer  = {{args.snvsniffer | quote}}
strelka   = {{args.strelka | quote}}
virmid    = {{ args.virmid | quote}}
vardict   = {{ args.vardict | quote}}
joboutdir = {{job.outdir | quote}}
bamIndex(tumor,  samtools = None, nthread = nthread)
bamIndex(normal, samtools = None, nthread = nthread)

if tprefix == nprefix:
	tprefix = tprefix + '_TUMOR'
	nprefix = nprefix + '_NORMAL'

if gz:	
	# remove ending .gz
	outfile = outfile[:-3]

if not path.exists (ref):
	stderr.write ("Reference file not specified")
	exit (1)

tmpdir = path.join (tmpdir, "{procid}.{tprefix}.{nprefix}.{jobindex}".format(
	jobindex = jobindex,
	procid   = procid,
	tprefix  = tprefix,
	nprefix  = nprefix
))
if not path.exists (tmpdir):
	makedirs (tmpdir)

def run_gatk():
	# generate interval list file
	intvfile = {{job.outdir | path.join: "interval.list" | quote}}
	cmd = '{samtools} idxstats {tumor!r} | head -1 | cut -f1 > {intvfile!r}'.format(
		samtools = samtools,
		tumor    = tumor,
		intvfile = intvfile
	)
	runcmd(cmd)

	mem = mem2(mem, 'java')

	params['I:tumor']  = tumor
	params['I:normal'] = normal
	
	params.R   = ref
	params.o   = outfile
	params.nct = nthread
	params.L   = intvfile

	cmd = '{gatk} -T MuTect2 {mem} -Djava.io.tmpdir={tmpdir!r} {args}'.format(
		gatk   = gatk,
		mem    = mem,
		tmpdir = tmpdir,
		args   = cmdargs(params, dash = '-', equal = ' ')
	)
	runcmd(cmd)
	if gz: runcmd(['gzip', outfile])

def run_somaticsniper():
	params.f   = ref
	params.F   = 'vcf'
	params[''] = [tumor, normal, outfile]
	cmd = '{ssniper} {args}'.format(ssniper = ssniper, args = cmdargs(params))
	runcmd(cmd)
	if gz: runcmd(['gzip', outfile])

def run_snvsniffer():
	# generate a header file
	theader = {{job.outdir | path.join: bn(i.tumor)  | @append: '.header' | quote}}
	nheader = {{job.outdir | path.join: bn(i.normal) | @append: '.header' | quote}}
	cmd = '{samtools} view -H {infile!r} > {hfile!r}'
	runcmd(cmd.format(samtools = samtools, infile = tumor, hfile = theader))
	runcmd(cmd.format(samtools = samtools, infile = normal, hfile = nheader))

	params.g = ref
	params.o = outfile

	params[''] = [theader, nheader, tumor, normal]
	cmd = '{ssniffer} somatic {args}'.format(ssniffer, cmdargs(params))
	runcmd(cmd)
	if gz: runcmd(['gzip', outfile])

def _mergeAndAddGT(snvvcf, indvcf, outfile):
	from pysam import VariantFile
	snv = VariantFile(snvvcf)
	ind = VariantFile(indvcf)
	
	snv.header.info.add('TYPE', 1, 'String', 'Type of somatic mutation')
	ind.header.info.add('TYPE', 1, 'String', 'Type of somatic mutation')
	snv.header.info.add('QSI', 1, 'Integer', 'Quality score for any somatic variant, ie. for the ALT haplotype to be present at a significantly different frequency in the tumor and normal')
	snv.header.info.add('TQSI', 1, 'Integer', 'Data tier used to compute QSI')
	snv.header.info.add('QSI_NT', 1, 'Integer', 'Quality score reflecting the joint probability of a somatic variant and NT')
	snv.header.info.add('TQSI_NT', 1, 'Integer', 'Data tier used to compute QSI_NT')
	snv.header.info.add('IC', 1, 'Integer', 'Number of times RU repeats in the indel allele')
	snv.header.info.add('IHP', 1, 'Integer', 'Largest reference interrupted homopolymer length intersecting with the indel')
	snv.header.info.add('OVERLAP', 0, 'Flag', 'Somatic indel possibly overlaps a second indel.')
	snv.header.info.add('RC', 1, 'Integer', 'Number of times RU repeats in the reference allele')
	snv.header.info.add('RU', 1, 'String', 'Smallest repeating sequence unit in inserted or deleted sequence')
	snv.header.formats.add('GT', 1, 'String', 'Possible genotype')
	ind.header.formats.add('GT', 1, 'String', 'Possible genotype')
	snv.header.formats.add('BCN50', 1, 'Float', 'Fraction of filtered reads within 50 bases of the indel.')
	snv.header.formats.add('DP2', 1, 'Integer', 'Read depth for tier2')
	snv.header.formats.add('DP50', 1, 'Float', 'Average tier1 read depth within 50 bases')
	snv.header.formats.add('FDP50', 1, 'Float', 'Average tier1 number of basecalls filtered from original read depth within 50 bases')
	snv.header.formats.add('SUBDP50', 1, 'Float', 'Average number of reads below tier1 mapping quality threshold aligned across sites within 50 bases')
	snv.header.formats.add('TAR', 2, 'Integer', 'Reads strongly supporting alternate allele for tiers 1,2')
	snv.header.formats.add('TIR', 2, 'Integer', 'Reads strongly supporting indel allele for tiers 1,2')
	snv.header.formats.add('TOR', 2, 'Integer', 'Other reads (weak support or insufficient indel breakpoint overlap) for tiers 1,2')

	contigs = list(snv.header.contigs.keys())
	out = open(outfile, 'w')
	#Can't change sample names with VariantFile
	#out = VariantFile(outfile, 'w', header = snv.header)
	headers = str(snv.header).splitlines()
	cnames  = headers[-1].split("\t")

	cnames [-2] = nprefix
	cnames [-1] = tprefix
	headers[-1] = "\t".join(cnames)
	out.write("\n".join(headers) + "\n")
	
	r1 = r2 = None
	indel_gts = {
		'ref': (0, 0),
		'het': (0, 1),
		'hom': (1, 1)
	}
	while True:
		if not r1:
			try:
				r1 = next(snv)
				r1.info['TYPE'] = 'SNV'
				alleles = (r1.ref, ) + r1.alts
				gts = r1.info['SGT'].split('->')
				try:
					r1.samples['NORMAL']['GT'] = tuple(sorted(alleles.index(gt) for gt in list(gts[0])))
					r1.samples['TUMOR']['GT'] = tuple(sorted(alleles.index(gt) for gt in list(gts[1])))
				except ValueError:
					r1 = None
					continue
			except StopIteration:
				r1 = None
		if not r2:
			try:
				r2 = next(ind)
				r2.info['TYPE'] = 'INDEL'
				gts = r2.info['SGT'].split('->')
				r2.samples['NORMAL']['GT'] = indel_gts[gts[0]]
				r2.samples['TUMOR']['GT']  = indel_gts[gts[1]]
			except StopIteration:
				r2 = None
		
		if r1 and r2:
			if (contigs.index(r1.chrom), r1.pos) < (contigs.index(r2.chrom), r2.pos):
				out.write(str(r1))
				r1 = None
			else:
				out.write(str(r2))
				r2 = None
		elif r1:
			out.write(str(r1))
			r1 = None
		elif r2:
			out.write(str(r2))
			r2 = None
		else:
			break
	out.close()

def run_strelka():
	cparams                = {{args.configParams | repr}}
	cparams.normalBam      = normal
	cparams.tumorBam       = tumor
	cparams.referenceFasta = ref
	cparams.runDir         = joboutdir
	runcmd('{strelka} {args}'.format(strelka = strelka, args = cmdargs(cparams)))

	params.m = 'local'
	params.g = mem2(mem, 'G')[:-1]
	params.j = nthread
	runcmd('{joboutdir}/runWorkflow.py {args}'.format(joboutdir = joboutdir, args = cmdargs(params)))

	snvvcf = path.join(joboutdir, 'results', 'variants', 'somatic.snvs.vcf.gz')
	indvcf = path.join(joboutdir, 'results', 'variants', 'somatic.indels.vcf.gz')
	_mergeAndAddGT(snvvcf, indvcf, outfile)
	if gz: runcmd(['gzip', outfile])

def run_virmid():
	params.R = ref
	params.D = tumor
	params.N = normal
	params.w = joboutdir
	cmd = '{virmid} {mem} -Djava.io.tmpdir={tmpdir!r} {args}'.format(
		virmid = virmid,
		mem    = mem2(mem, 'java'),
		tmpdir = tmpdir,
		args   = cmdargs(params)
	)
	runcmd(['mv', path.join(joboutdir, '*.virmid.som.passed.vcf'), outfile])
	if gz: runcmd(['gzip', outfile])

def run_vardict():
	params.v = True
	params.G = ref
	params.b = '{}|{}'.format(tumor, normal)
	cmd = '{vardict} {args} > {outfile!r}'.format(
		vardict = vardict,
		args    = cmdargs(params),
		outfile = outfile
	)
	runcmd(cmd)
	if gz: runcmd(['gzip', outfile])

tools = {
	'gatk'         : run_gatk,
	'somaticsniper': run_somaticsniper,
	'snvsniffer'   : run_snvsniffer,
	'strelka'      : run_strelka,
	'virmid'       : run_virmid,
	'vardict'      : run_vardict
}

try:
	tools[tool]()
except Exception:
	raise
finally:
	rmtree (tmpdir)
