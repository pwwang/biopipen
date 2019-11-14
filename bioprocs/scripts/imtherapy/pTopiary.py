"""
./topiary \
  --vcf somatic.vcf \
  --mhc-predictor netmhcpan \
  --mhc-alleles HLA-A*02:01,HLA-B*07:02 \
  --ic50-cutoff 500 \
  --percentile-cutoff 2.0 \
  --mhc-epitope-lengths 8-11 \
  --rna-gene-fpkm-tracking-file genes.fpkm_tracking \
  --rna-min-gene-expression 4.0 \
  --rna-transcript-fpkm-tracking-file isoforms.fpkm_tracking \
  --rna-min-transcript-expression 1.5 \
  --output-csv epitopes.csv
"""
import re
from os import environ
from pathlib import Path
from pyppl import Box
from bioprocs.utils import shell2 as shell
from bioprocs.utils.tsvio2 import TsvReader, TsvWriter

infile        = {{i.infile | quote}}
mhcallele     = {{i.mhcallele | quote}}
gexpr         = {{i.gexpr | quote}}
texpr         = {{i.texpr | quote}}
outfile       = Path({{o.outfile | quote}})
neatfile      = outfile.with_suffix('.neat.txt')
topiary       = {{args.topiary | quote}}
netmhc        = {{args.netmhc | quote}}
netmhcpan     = {{args.netmhcpan | quote}}
netmhciipan   = {{args.netmhciipan | quote}}
netmhccons    = {{args.netmhccons | quote}}
smm           = {{args.smm | quote}}
smm_pmbec     = {{args.smm_pmbec | quote}}
mhc_predictor = {{args.mhc_predictor | quote}}
genome        = {{args.genome | quote}}
params        = {{args.params | repr}}
wildtype      = {{args.wildtype | repr}}
tmpdir        = Path({{args.tmpdir | quote}}) / '.'.join([
	{{proc.id | quote}}, {{proc.tag | quote}}, {{proc.suffix | quote}}, {{job.index | quote}}])
tmpdir.mkdir(exist_ok = True, parents = True)

# check if we have downloaded annotation data for the genome
gmaps = {'hg19': 'GRCh37', 'hg38': 'GRCh38'}
datadir = Path.home().joinpath('.cache', 'pyensembl')
if not datadir.joinpath(genome).is_dir() and not datadir.joinpath(gmaps.get(genome, genome)).is_dir():
	raise RuntimeError("You don't have annotation data for genome {}{} installed. "
		"Either you run 'pyensembl install' first or "
		"specify 'params.download_reference_genome_data = True'. "
		"If you have it installed somewhere else, make a symbolic link to {}".format(genome, ('/' + gmaps[genome]) if genome in gmaps else '', datadir))
if not datadir.joinpath(genome).is_dir() and datadir.joinpath(gmaps.get(genome, genome)).is_dir():
	genome = gmaps[genome]
datadir = datadir / genome

shell.load_config(topiary = topiary)
if infile.endswith('.vcf') or infile.endswith('.vcf.gz'):
	params.vcf = infile
else:
	params.maf = infile

if Path(mhcallele).is_file():
	mhcallelefile = outfile.parent.joinpath('mhcalleles.txt')
	with open(mhcallele, 'r') as fin, open(mhcallelefile, 'w') as fout:
		for line in fin:
			fout.write(line.replace('*', ''))
	params['mhc-alleles-file'] = mhcallelefile
else:
	params['mhc-alleles'] = mhcallele.replace('*', '')

if gexpr:
	params['rna-gene-fpkm-tracking-file'] = gexpr
if texpr:
	params['rna-transcript-fpkm-tracking-file'] = texpr

params.genome = genome
params['output-csv'] = outfile
params['mhc-predictor'] = mhc_predictor

# make sure those mhc-predictors are in PATH
PATHs = set()
for mhcpred in (netmhc, netmhcpan, netmhciipan, netmhccons, smm, smm_pmbec):
	if '/' in mhcpred:
		PATHs.add(str(Path(mhcpred).parent))

params._env = Box(PATH = environ['PATH'] + ':' + ':'.join(PATHs))

shell.fg.topiary(**params)

if not wildtype:
	reader = TsvReader(outfile, comment = '###')
	writer = TsvWriter(neatfile)
	writer.cnames = ['HLA_allele', 'mt_peptide', 'mt_affinity', 'gene']
	writer.writeHead()
	for r in reader:
		writer.write([r.allele, r.peptide, r.affinity, r.gene])
else:
	if mhc_predictor not in ('netmhc', 'netmhcpan', 'netmhciipan', 'netmhccons', 'smm', 'smm_pmbec'):
		raise ValueError("Calculations for wildtype peptides are only available with local predictors.")

	tmpfile = outfile.with_suffix('.topiary.txt')
	shell.mv(outfile, tmpfile)

	wildpeps = set()
	mutpeps = {}
	tpreader = TsvReader(tmpfile, comment = '###')
	for r in tpreader:
		if r.effect_type != 'Substitution':
			# I don't know how to get the wildtype peptides if it is not a substitution
			continue
		# parse effect: p.N84S
		m = re.match(r'^p\.([A-Z])\d+([A-Z])$', r.effect)
		if not m:
			continue
		wildpep = r.peptide
		index = int(r.mutation_start_in_peptide)
		if wildpep[index] != m.group(2):
			continue
		wildpep = wildpep[:index] + m.group(1) + wildpep[(index+1):]
		mutpeps[r.peptide + '\t' + r.allele] = wildpep
		wildpeps.add(wildpep)

	def run_netmhc():
		pass

	def run_netmhcpan():
		shell.load_config(netmhcpan = netmhcpan)
		mhcallele2 = params['mhc-alleles'] if 'mhc-alleles' in params else ','.join(
			allele for allele in Path(params['mhc-alleles-file']).read_text().splitlines() if allele
		)

		wildfile = outfile.parent / 'wildtype.peptides.txt'
		wildfile.write_text('\n'.join(wildpeps))
		xlsfile = outfile.parent / 'wildtype.binding.txt'
		nparams = Box(
			a = mhcallele2, v = True, BA = True, inptype = 1, f = wildfile, _prefix = '-',
			xls = True, xlsfile = xlsfile)
		shell.fg.netmhcpan(**nparams)

		if not xlsfile.is_file():
			raise RuntimeError("Failed to run netmhcpan, output file not generated.")
		# read the output
		"""
					HLA-A24:02					HLA-A29:02
		Pos	Peptide	ID	core	icore	1-log50k	nM	Rank	core	icore	1-log50k	nM	Rank
		0	LYLPALWFH	PEPLIST	LYLPALWFH	LYLPALWFH	0.1353	11560.5488	6.1138	LYLPALWFH	LYLPALWFH	0.4137	568.6087	1.1231
		0	RRQRRQRRW	PEPLIST	RRQRRQRRW	RRQRRQRRW	0.0788	21311.8301	12.3392	RRQRRQRRW	RRQRRQRRW	0.0308	35829.9805	47.6206
		"""
		with xlsfile.open('r') as f:
			alleles = [allele.replace('HLA-A', 'HLA-A*').replace('HLA-B', 'HLA-B*').replace('HLA-C', 'HLA-C*')
				for allele in f.readline().strip().split('\t') if allele]
		reader = TsvReader(xlsfile, comment = '\t\t\t')
		wildbindings = {}
		for r in reader:
			peptide = r[1]
			for i, hla in enumerate(alleles):
				wildbindings[peptide + '\t' + hla] = float(r[7 + i*5])

		writer = TsvWriter(outfile)
		writer.cnames = tpreader.cnames + ['wildpeptide', 'wildaffinity', 'deltaaffinity']
		writer.writeHead()
		nwriter = TsvWriter(neatfile)
		nwriter.cnames = ['HLA_allele', 'mt_peptide', 'mt_affinity', 'wt_peptide', 'wt_affinity', 'delta_affinity', 'gene']
		nwriter.writeHead()
		tpreader.rewind()
		for r in tpreader:
			r.wildpeptide = mutpeps.get(r.peptide + '\t' + r.allele, '-')
			r.wildaffinity = wildbindings.get(r.wildpeptide + '\t' + r.allele, '-')
			if r.wildaffinity != '-':
				r.deltaaffinity = float(r.affinity) - r.wildaffinity
			else:
				r.deltaaffinity = '-'
			nwriter.write([	r.allele, r.peptide, r.affinity, r.wildpeptide,
							r.wildaffinity, r.deltaaffinity, r.gene])
			writer.write(r)

	def run_netmhciipan():
		pass

	def run_netmhccons():
		pass

	def run_smm():
		pass

	def run_smm_pmbec():
		pass

	runner = {
		'netmhc'      : run_netmhc,
		'netmhcpan'   : run_netmhcpan,
		'netmhciipan' : run_netmhciipan,
		'netmhccons'  : run_netmhccons,
		'smm'         : run_smm,
		'smm-pmbec'   : run_smm_pmbec,
	}
	runner.get(mhc_predictor)()
