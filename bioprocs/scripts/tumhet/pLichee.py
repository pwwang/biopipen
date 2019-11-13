from pyppl import Box
from os import path
from bioprocs.utils.tsvio2 import TsvReader, TsvWriter
from bioprocs.utils import shell2 as shell, logger

infile   = {{i.infile | quote}}
outdir   = {{o.outdir | quote}}
dot      = {{args.dot | quote}}
lichee   = {{args.lichee |quote}}
params   = {{args.params | repr}}
fishplot = {{args.fishplot |repr}}
Rscript  = {{args.Rscript |quote}}
devpars  = {{args.devpars |repr}}

shell.load_config(lichee = {'_exe': lichee, '_prefix': '-'}, Rscript = Rscript, dot = dot)

# expected input:
# #chr    position    description    Normal    S1    S2    S3    S4
# 17      123456      A/T DUSP19     0.0       0.1   0.2   0.25  0.15
# 11      341567      C/G MUC16      0.0       0.4   0.09  0.38  0.24
# 9       787834      A/C OR2A14     0.0       0.35  0.14  0.17  0.48

# let's see if the input file is pyclone directory
if path.isdir(infile) and infile.endswith('.pyclone'):
	# tables/loci.tsv
	#mutation_id	sample_id	cluster_id	cellular_prevalence	cellular_prevalence_std	variant_allele_frequency	gene
	#chr10:17875816	NA12878	2	0.828865827591	0.0355432431918	0.904883227176	Gene1
	#chr10:17875816	NA12156	2	0.741264908504	0.11050033704	0.560311284047	Gene2

	indata = {} # {mutation: {gene: gene, sample1: {cluster_id: cluster_id, cp = cp}, sample2: ...}}
	reader = TsvReader(path.join(infile, 'tables/loci.tsv'))
	samples = set()
	for i, r in enumerate(reader):
		indata.setdefault(r.mutation_id, Box())
		if not indata[r.mutation_id].get('_gene'):
			indata[r.mutation_id]._gene = r.gene
		indata[r.mutation_id]._cluster     = r.cluster_id
		indata[r.mutation_id][r.sample_id] = r.cellular_prevalence
		samples.add(r.sample_id)

	profile = lambda vafs: ''.join(str(int(bool(vaf))) for vaf in vafs)

	samples = list(samples)
	writer = TsvWriter(path.join(outdir, 'input.txt'))
	writer.cnames = ['chr', 'position', 'description', 'Normal'] + samples
	writer.writeHead(lambda cnames: '#' + '\t'.join(cnames))
	clusters = {}
	i = 1
	for mutation, info in indata.items():
		vafs = [0.0] + [info[sample] for sample in samples]
		prof = profile(vafs)
		clusters.setdefault(info.cluster, [])
		clusters[info.cluster].append(Box(
			profile = prof,
			vafs    = vafs,
			index   = str(i)
		))
		i += 1
		writer.write(mutation.split(':') + [info._gene] + vafs)
	writer.close()

	# 01111     0.0   0.1   0.23  0.23  0.13    1,2
	# 01111     0.0   0.4   0.4   0.45  0.45    3
	# 01001     0.0   0.4   0.0   0.0   0.24    4
	writer = TsvWriter(path.join(outdir, 'clusters.txt'))
	for cluster_id, info in clusters.items():
		writer.write([info[0].profile] + [
			sum(float(inf.vafs[i]) for inf in info)/float(len(info)) for i in range(len(info[0].vafs))
		] + [','.join(inf.index for inf in info)])
	writer.close()
	params.cp = True
	infile = path.join(outdir, 'input.txt')
	params.clustersFile = path.join(outdir, 'clusters.txt')
elif infile.endswith('.maf') or infile.endswith('.maf.gz'):
	reader = TsvReader(infile)
	if 't_alt_count' not in reader.cnames:
		raise ValueError('t_alt_count not found in MAF file.')
	if 't_ref_count' not in reader.cnames:
		raise ValueError('t_ref_count not found in MAF file.')

	# get all mutations
	mutations = {}
	samples = set()
	for r in reader:
		mut = r.Chromosome + ':' + r.Start_Position
		mutations.setdefault(mut, Box())
		try:
			mutations[mut][r.Tumor_Sample_Barcode] = float(r.t_alt_count) / (float(r.t_alt_count) + float(r.t_ref_count))
		except ZeroDivisionError:
			logger.warning(
				f'Failed to get allele frequency for variant: {mut} with alt/ref counts: {r.t_alt_count}, {r.t_ref_count}, setting VAF to 0')
			mutations[mut][r.Tumor_Sample_Barcode] = 0.0
		mutations[mut]._gene = r.Hugo_Symbol

		samples.add(r.Tumor_Sample_Barcode)
	samples = list(samples)

	writer = TsvWriter(path.join(outdir, 'input.txt'))
	writer.cnames = ['chr', 'position', 'description', 'Normal'] + samples
	writer.writeHead(lambda cnames: '#' + '\t'.join(cnames))
	for mut, info in mutations.items():
		writer.write(mut.split(':') + [info._gene, 0.0] + [
			info.get(sample, 0.0) for sample in samples
		])
	writer.close()
	infile = path.join(outdir, 'input.txt')
else:
	shell.ln_s(infile, path.join(outdir, path.basename(infile)))

prefix = path.join(outdir, path.basename(infile).split('.')[0])
outfile = prefix + '.trees.txt'

params.i = infile
params.n = 0
params.o = outfile
params.color = True
params.v = True
params.dot = True
params.dotFile = prefix + '.dot'

shell.fg.lichee('-build', **params)

shell.fg.dot(params.dotFile, f'-Gdpi={devpars.res}', T = 'png', o = prefix + '.png')

# save clusters and snv infos in different file for reporting
snvfile = prefix + '.snvinfo.txt'
snvinfos = {}
nodes_started = False
info_started = False
with open(outfile, 'r') as fout:
	for line in fout:
		line = line.strip()
		if line == 'Nodes:':
			nodes_started = True
			continue
		if line == 'SNV info:':
			info_started = True
			continue
		if not line:
			nodes_started = False
			info_started = False
			continue
		if nodes_started:
			parts = line.split('\t')
			cluster = int(parts[0]) - 1
			snvs = parts[3:]
			for snv in snvs:
				snvinfos.setdefault(snv, {})
				snvinfos[snv]['cluster'] = cluster
			continue
		if info_started:
			snv, info = line.split(': ', 1)
			snvinfos[snv]['info'] = info
with open(snvfile, 'w') as fsnv:
	fsnv.write('SNV\tClusterID\tDescription\n')
	for snv, info in snvinfos.items():
		fsnv.write(f'{snv}\t{info["cluster"]}\t{info["info"]}\n')

