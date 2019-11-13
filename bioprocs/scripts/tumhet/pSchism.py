import random
import yaml
from os import path
from glob import glob
from pyppl import Box
from bioprocs.utils import shell2 as shell
from bioprocs.utils.tsvio2 import TsvReader, TsvWriter

infile   = {{i.infile | quote}}
purity   = {{i.purity | quote}}
outdir   = {{o.outdir | quote}}
schism   = {{args.schism | quote}}
fishplot = {{args.fishplot | repr}}
Rscript  = {{args.Rscript | quote}}
dot      = {{args.dot | quote}}
params   = {{args.params |repr}}
devpars  = {{args.devpars |repr}}

shell.load_config(schism = schism, Rscript = Rscript, dot = dot)

# default configurations
config                        = Box()
config.working_dir            = outdir
config.output_prefix          = {{i.infile | stem | quote}}
config.cellularity_estimation = "schism"
config.cellularity_estimator  = Box(coverage_threshold = 50.0, absent_mode = 1)
config.tumor_sample_purity    = Box()
config.hypothesis_test        = Box(test_level = "mutations", significance_level = 0.05, store_pvalues = True)
config.genetic_algorithm      = Box(
	instance_count         = 10,
	generation_count       = 50,
	generation_size        = 100,
	random_object_fraction = 0.2,
	mutation_probability   = .9,
	crossover_probability  = .25,
	fitness_coefficient    = 5.0,
	verbose                = False
)

# purity
if purity:
	purites = {}
	reader = TsvReader(purity, cnames = False)
	for r in reader:
		purites[r[0]] = r[1]
	reader.close()
	config.tumor_sample_purity.update(purites)

# expected mutation_to_cluster_assignment file:
# mutationID	clusterID
# 0	0
# 1	0
# 2	0

# expected mutation_raw_input file:
# sampleID	mutationID	referenceReads	variantReads	copyNumber
# S1	0	368	132	2
# S1	1	381	119	2
# S1	2	367	133	2

# let's see if the input file is pyclone directory
if path.isdir(infile) and infile.endswith('.pyclone'):
	# tables/loci.tsv
	#mutation_id	sample_id	cluster_id	cellular_prevalence	cellular_prevalence_std	variant_allele_frequency	gene
	#chr10:17875816	NA12878	2	0.828865827591	0.0355432431918	0.904883227176	Gene1
	#chr10:17875816	NA12156	2	0.741264908504	0.11050033704	0.560311284047	Gene2

	# <sample>.muts.bed
	##CHR   START   END     NAME            GT  REF VAR CASE    VAF
	#chr2	3482633	3482633	chr2:3482633	AB	103	2	sample1	0.019
	#chr2	9661972	9661972	chr2:9661972	AB	60	6	sample1	0.091

	# get cluster file first:
	writer = TsvWriter(path.join(outdir, 'mutation_to_cluster_assignment.tsv'))
	writer.cnames = ['mutationID', 'clusterID']
	writer.writeHead()
	reader = TsvReader(path.join(infile, 'tables/loci.tsv'))
	mutclusters = {}
	for r in reader:
		mutclusters[r.mutation_id] = r.cluster_id
	for mutid, clustid in mutclusters.items():
		writer.write([mutid, clustid])
	reader.close()
	writer.close()

	# get the mutation_raw_input file
	writer = TsvWriter(path.join(outdir, 'mutation_raw_input.tsv'))
	writer.cnames = ['sampleID', 'mutationID', 'referenceReads', 'variantReads', 'copyNumber']
	writer.writeHead()
	for mutbed in glob(path.join(infile, '*.muts.bed')):
		sample = path.basename(mutbed)[:-9]
		# update default purity
		if sample not in config.tumor_sample_purity:
			config.tumor_sample_purity[sample] = 1
		reader = TsvReader(mutbed)
		mclusters = mutclusters.copy()
		for r in reader:
			if r[3] not in mclusters:
				continue
			del mclusters[r[3]] # remove duplicates
			writer.write([sample, r[3], r[5], r[6], 2])
	writer.close()
	reader.close()

	config.mutation_to_cluster_assignment = 'mutation_to_cluster_assignment.tsv'
	config.mutation_raw_input = 'mutation_raw_input.tsv'
# let's see if the input file is a mutation file, then split it into 2 files,
# mutation_to_cluster_assignment and mutation_raw_input
else:
	writer = TsvWriter(path.join(outdir, 'mutation_to_cluster_assignment.tsv'))
	writer.cnames = ['mutationID', 'clusterID']
	writer.writeHead()
	reader = TsvReader(infile)
	mutclusters = {}
	for r in reader:
		mutclusters[r.mutationID] = r.clusterID
	for mutid, clustid in mutclusters.items():
		writer.write([mutid, clustid])
	writer.close()
	reader.rewind()

	writer = TsvWriter(path.join(outdir, 'mutation_raw_input.tsv'))
	writer.cnames = ['sampleID', 'mutationID', 'referenceReads', 'variantReads', 'copyNumber']
	writer.writeHead()
	for r in reader:
		if r[0] not in config.tumor_sample_purity:
			config.tumor_sample_purity[r[0]] = 1
		writer.write(r[:-1])
	writer.close()
	reader.close()
	config.mutation_to_cluster_assignment = 'mutation_to_cluster_assignment.tsv'
	config.mutation_raw_input = 'mutation_raw_input.tsv'

configfile = path.join(outdir, 'config.yaml')
with open(configfile, 'w') as fconf:
	yaml.dump(config.to_dict(), fconf)

shell.fg.schism.analyze(c = configfile)

# schism always returns 0
# check if consensus plot generated
if not glob(path.join(outdir, '*.GA.consensusTree.pdf')):
	raise RuntimeError

ctreefile = glob(path.join(outdir, '*.GA.consensusTree'))[0]
# parent	child	frequency	label
# 1	0	0.5	1/2
# 2	0	0.5	1/2
# 3	1	1.0	2/2
# 3	2	1.0	2/2

# make it to dot

dotfile = ctreefile[:-17] + '.dot'
reader = TsvReader(ctreefile)
nodes = set()
edges = []
childs = set()
for r in reader:
	nodes.add(r.parent)
	nodes.add(r.child)
	edges.append(Box(parent = r.parent, child = r.child, freq = float(r.frequency), label = r.label))
	childs.add(r.child)

dotstr = ['digraph G {']
dotstr.append('forcelabels=true')
dotstr.append('999 [ shape=plaintext label="GL" ]')
for node in nodes:
	fillcolor = '#' + ''.join(str(hex(random.randint(150, 255)))[2:] for _ in range(3))
	dotstr.append(f'{node} [ shape=circle style=filled label="{node}" fillcolor="{fillcolor}"]'.format(node))

random.seed(8525)
for edge in edges:
	if edge.freq < 1:
		dotstr.append(f'{edge.parent} -> {edge.child} [ label="{edge.label}" style=dashed ]')
	else:
		dotstr.append(f'{edge.parent} -> {edge.child} [ label="{edge.label}" style=solid ]')
# add GL to roots
for child in (nodes - childs):
	dotstr.append(f'999 -> {child} [ label="shared" style=solid ]')
dotstr.append('}')
with open(dotfile, 'w') as f:
	f.write('\n'.join(dotstr))

shell.fg.dot(dotfile, f'-Gdpi={devpars.res}', T = 'png', o = ctreefile[:-17] + '.tree.png')
