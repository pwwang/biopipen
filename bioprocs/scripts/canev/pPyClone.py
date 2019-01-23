from pysam import VariantFile
from pyppl import Box
from os import path
from bioprocs.utils import funcargs, cmdargs, runcmd

vfvcfs   = {{ i.vfvcfs | repr}}
cnvcfs   = {{ i.cnvcfs | repr}}
outdir   = {{ o.outdir | quote}}
params   = {{ args.params | repr}}
pyclone  = {{ args.pyclone | quote}}
vfsamcol = {{ args.vfsamcol | repr}} - 1
cnsamcol = {{ args.cnsamcol | repr}} - 1
varcount = {{ args.varcount}}
cncount  = {{ args.cncount}}
if not callable(varcount):
	varcount = lambda fmt: fmt.get('AD')
if not callable(cncount):
	cncount = lambda fmt: fmt.get('CN')

def vcf2vaf(vcf):
	for record in vcf:
		fmt  = Box(record.samples[vfsamcol].items())
		info = Box(record.info.items())
		vc   = varcount(fmt) if len(funcargs(varcount)) == 1 else varcount(fmt, info)
		if not isinstance(vc, dict):
			vc = {'count': vc}
		if vc['count'] is None:
			continue
		yield dict(
			mutation_id = "{chr}:{pos}".format(chr = record.chrom, pos = record.pos),
			chrom       = record.chrom,
			pos         = record.pos,
			ref_counts  = vc.get('depth', fmt.DP),
			var_counts  = vc['count'],
			# need sex information to determine for X
			normal_cn   = 1 if 'Y' in record.chrom else 2,
			minor_cn    = 0,
			major_cn    = 2
		)

def vcf2cn(vcf):
	for record in vcf:
		fmt  = Box(record.samples[cnsamcol].items())
		info = Box(record.info.items())
		cn   = cncount(fmt) if len(funcargs(cncount)) else cncount(fmt, info)
		if not isinstance(cn, dict):
			cn = {'cn': cn}
		if cn['cn'] is None:
			continue
		yield dict(
			chrom      = record.chrom,
			start      = record.pos,
			end        = int(cn.get('end', info.get('END', record.pos + 1))),
			cn         = cn['cn']
		)

def combineVafAndCn(vafs, outfile, cns = None, contigs = None):
	# vafs and cns are generators
	# write header
	out = open(outfile, 'w')
	out.write('mutation_id\tref_counts\tvar_counts\tnormal_cn\tminor_cn\tmajor_cn\n')
	if not cns:
		for vaf in vafs:
			out.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(
				vaf['mutation_id'],
				vaf['ref_counts'],
				vaf['var_counts'],
				vaf['normal_cn'],
				vaf['minor_cn'],
				vaf['major_cn']
			))
		out.close()
		return
	for cn in cns:
		flag = ''
		while True:
			try:
				vaf  = next(vafs)
				if vaf['chrom'] == cn['chrom']:
					if vaf['pos'] < cn['start']:
						out.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(
							vaf['mutation_id'],
							vaf['ref_counts'],
							vaf['var_counts'],
							vaf['normal_cn'],
							vaf['minor_cn'],
							vaf['major_cn']
						))
						continue
					elif vaf['pos'] > cn['end']:
						flag = 'continue'
						break
					else:
						vaf['major_cn'] = cn['cn']
						out.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(
							vaf['mutation_id'],
							vaf['ref_counts'],
							vaf['var_counts'],
							vaf['normal_cn'],
							vaf['minor_cn'],
							vaf['major_cn']
						))
						continue
				elif contigs.index(vaf['chrom']) < contigs.index(cn['chrom']):
					out.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(
						vaf['mutation_id'],
						vaf['ref_counts'],
						vaf['var_counts'],
						vaf['normal_cn'],
						vaf['minor_cn'],
						vaf['major_cn']
					))
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
	out.close()

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
params.prior = True
cmd = 'cd {outdir}; {pyclone} run_analysis_pipeline {args}'.format(outdir = outdir, pyclone = pyclone, args = cmdargs(params, equal = ' '))
runcmd(cmd)

