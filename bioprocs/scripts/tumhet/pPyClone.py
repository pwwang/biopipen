from pysam import VariantFile
from pyppl import Box
from os import path
from bioprocs.utils import funcargs, logger, shell2 as shell

vfvcfs   = {{ i.vfvcfs | repr}}
cnvcfs   = {{ i.cnvcfs | repr}}
outdir   = {{ o.outdir | quote}}
params   = {{ args.params | repr}}
pyclone  = {{ args.pyclone | quote}}
vfsamcol = {{ args.vfsamcol | repr}} - 1
cnsamcol = {{ args.cnsamcol | repr}} - 1
varcount = {{ args.varcount}}
cncount  = {{ args.cncount}}

shell.load_config(pyclone = dict(
	_exe = pyclone,
	_cwd = outdir
))

varcount_default = lambda fmt, info = None: dict(
	ref = fmt['AD'][0],
	alt = fmt['AD'][1],
	af  = fmt.get('AF', float(fmt['AD'][1]) / float(fmt.get('DP', sum(fmt['AD']))))
)
if not callable(varcount):
	varcount = varcount_default
else:
	def varcount2(fmt, info = None):
		def_varc = varcount_default(fmt, info)
		varc = varcount(fmt) if len(funcargs(varcount)) == 1 else varcount(fmt, info)
		def_varc.update(varc)
		return def_varc
	varcount = varcount2

cncount_default = lambda fmt, info: dict(
	major = fmt.get('CN') if isinstance(fmt.get('CN'), int) else \
		fmt.get('CN')[0] if isinstance(fmt.get('CN'), (tuple,list)) else None,
	minor = fmt.get('CN')[1] if isinstance(fmt.get('CN'), (tuple,list)) else 0,
	end   = (info or {}).get('END')
)
if not callable(cncount):
	cncount = cncount_default
else:
	def cncount2(fmt, info = None):
		info = info or {}
		def_cn = cncount_default(fmt, info)
		cn = cncount(fmt) if len(funcargs(cncount)) == 1 else cncount(fmt, info)
		def_cn.update(cn)
		return def_cn
	cncount = cncount2

def vcf2vaf(vcf):
	for record in vcf:
		sample = record.header.samples[vfsamcol]
		fmt  = Box(record.samples[vfsamcol].items())
		info = Box(record.info.items())
		try:
			vc = varcount(fmt, info)
		except Exception as ex:
			logger.warning('[{sample}] Failed to get count for {chr}:{pos}: {ex}'.format(
				sample = sample, chr = record.chrom, pos = record.pos, ex = ex))
			continue
		gtsum = 0 if fmt['GT'][0] is None else sum(fmt['GT'])
		gt = 'BB' if gtsum == 2 else 'AB' if gtsum == 1 else 'AA'
		yield dict(
			mutation_id = "{chr}:{pos}".format(
				chr = record.chrom, pos = record.pos),
			chrom       = record.chrom,
			pos         = record.pos,
			ref_counts  = vc['ref'],
			var_counts  = vc['alt'],
			# need sex information to determine for X
			normal_cn    = 1 if 'Y' in record.chrom else 2,
			minor_cn     = 0,
			major_cn     = 1 if 'Y' in record.chrom else 2,
			variant_case = sample,
			variant_freq = vc['af'],
			genotype     = gt
		)

def vcf2cn(vcf):
	for record in vcf:
		fmt  = Box(record.samples[cnsamcol].items())
		info = Box(record.info.items())
		cn   = cncount(fmt, info)
		if cn['major'] is None:
			continue
		yield dict(
			chrom = record.chrom,
			start = record.pos,
			end   = cn.get('end') and int(cn.get('end')) or int(info.get('END', record.pos + 1)),
			major = int(cn['major']),
			minor = int(cn['minor'])
		)

def combineVafAndCn(vafs, outfile, cns, contigs):
	# vafs and cns are generators
	keys = [
		'mutation_id', 'ref_counts', 'var_counts', 'normal_cn',
		'minor_cn', 'major_cn', 'variant_case', 'variant_freq', 'genotype']
	with open(outfile, 'w') as out:
		out.write('\t'.join(keys) + '\n')
		if not cns:
			for vaf in vafs:
				out.write('\t'.join(str(rec[key]) for key in keys) + '\n')
			return
		for cn in cns:
			flag = ''
			while True:
				try:
					vaf  = next(vafs)
					if vaf['chrom'] == cn['chrom']:
						if vaf['pos'] < cn['start']:
							out.write('\t'.join(str(vaf[key]) for key in keys) + '\n')
							continue
						elif vaf['pos'] > cn['end']:
							flag = 'continue'
							break
						else:
							vaf['major_cn'] = cn['major']
							vaf['minor_cn'] = cn['minor']
							out.write('\t'.join(str(vaf[key]) for key in keys) + '\n')
							continue
					elif contigs.index(vaf['chrom']) < contigs.index(cn['chrom']):
						out.write('\t'.join(str(vaf[key]) for key in keys) + '\n')
						continue
					else:
						flag = 'continue'
						break
				except StopIteration:
					flag = 'break'
					break
			if flag == 'continue':
				continue
			if flag == 'break':
				break

if cnvcfs and len(cnvcfs) != len(vfvcfs):
	raise ValueError('Unequal length of variant and copy number vcf files.')

params = Box(in_files = [])

for i, vfvcf in enumerate(vfvcfs):
	vcf     = VariantFile(vfvcf)
	contigs = list(vcf.header.contigs.keys())
	sample  = vcf.header.samples[vfsamcol]
	vafs    = vcf2vaf(vcf)
	cns     = None
	if cnvcfs:
		cnvcf = VariantFile(cnvcfs[i])
		#if contigs != list(cnvcf.header.contigs.keys()):
		#	raise ValueError('Inconsistency contigs of variant and copy number vcf files.\n- First:\n%s\nSecond:\n%s\n' % (contigs, list(cnvcf.header.contigs.keys())))
		cns = vcf2cn(cnvcf)
	mutfile = path.join(outdir, sample + '.tsv')
	combineVafAndCn(vafs, mutfile, cns, contigs)
	params.in_files.append(mutfile)

#PyClone run_analysis_pipeline --in_files A.tsv B.tsv C.tsv --working_dir pyclone_analysis
# create matplotlibrc file
with open(path.join(outdir, 'matplotlibrc'), 'w') as f:
	f.write('backend: Agg\n')
params.working_dir = outdir
params.prior = 'total_copy_number'
shell.fg.pyclone.run_analysis_pipeline(**params)

assert path.exists(path.join(outdir, 'tables'))
assert path.exists(path.join(outdir, 'plots'))