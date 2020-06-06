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
from cyvcf2 import VCF
from gff import Gff
from diot import Diot
from cmdy import CmdyReturnCodeError
from bioprocs.utils import shell2 as shell, logger
from bioprocs.utils.tsvio2 import TsvReader, TsvWriter, TsvRecord
{% from os import path%}
infile        = {{i.infile | quote}}
afile         = {{i.afile | ?path.isfile | =readlines | !alwaysList | repr}}
outfile       = Path({{o.outfile | quote}})
outdir        = Path({{o.outdir | quote}})
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
refall        = {{args.refall | quote}}
tmpdir        = Path({{args.tmpdir | quote}}) / '.'.join([
	{{proc.id | quote}}, {{proc.tag | quote}}, {{proc.suffix | quote}}, {{job.index | quote}}])
tmpdir.mkdir(exist_ok = True, parents = True)

# check if we have downloaded annotation data for the genome
# topiary will use it to annotate the data
gmaps = {'hg19': 'GRCh37', 'hg38': 'GRCh38'}
datadir = Path.home().joinpath('.cache', 'pyensembl')
if not datadir.joinpath(genome).is_dir() and not datadir.joinpath(gmaps.get(genome, genome)).is_dir():
	raise RuntimeError("You don't have annotation data for genome {}{} installed. "
		"Either you run 'pyensembl install' first or "
		"specify 'params.download_reference_genome_data = True'. "
		"If you have it installed somewhere else, make a symbolic link to {}".format(genome, ('/' + gmaps[genome]) if genome in gmaps else '', datadir))
# if not datadir.joinpath(genome).is_dir() and datadir.joinpath(gmaps.get(genome, genome)).is_dir():
# 	genome = gmaps[genome]

# extract expression from VCF file
vcf      = VCF(infile)
gxfile   = txfile = False
features = set()
if vcf.contains('GX'):
	if not vcf.contains('CSQ'):
		raise ValueError('VCF file has to be annotated with by VEP')
	# tracking_id class_code  nearest_ref_id  gene_id gene_short_name tss_id  locus   length  coverage    FPKM    FPKM_conf_lo    FPKM_conf_hi    FPKM_status
	# ENSG00000240361 -   -   ENSG00000240361 OR4G11P -   chr1:62947-63887    -   -   0   0   0   OK
	# ENSG00000268020 -   -   ENSG00000268020 AL627309.1  -   chr1:53048-54936    -   -   0   0   0   OK
	gxfile = outfile.with_suffix('.gx_nopos')
	writer = TsvWriter(gxfile)
	writer.cnames = ['tracking_id', 'class_code', 'nearest_ref_id', 'gene_id', 'gene_short_name', 'tss_id', 'locus', 'length', 'coverage', 'FPKM', 'FPKM_conf_lo', 'FPKM_conf_hi', 'FPKM_status']
	writer.writeHead()
	for variant in vcf:
		# try..except
		try:
			gx = variant.format('GX')[0]
		except (KeyError, TypeError):
			continue
		csqs = variant.INFO['CSQ'].split(',')
		gxs = gx.split(',')
		for gx in gxs:
			gene, expr = gx.split('|', 1)
			csq = [csq for csq in csqs if f'|{gene}|' in csq][0].split('|')
			r                 = TsvRecord()
			r.tracking_id     = csq[4]
			r.class_code      = '-'
			r.nearest_ref_id  = '-'
			r.gene_id         = csq[4]
			r.gene_short_name = csq[3]
			r.tss_id          = '-'
			r.locus           = '<pos>'
			r.length          = '-'
			r.coverage        = '-'
			r.FPKM            = expr
			r.FPKM_conf_lo    = 0
			r.FPKM_conf_hi    = 1000
			r.FPKM_status     = 'OK'
			writer.write(r)
			features.add(r.tracking_id)
	writer.close()

if vcf.contains('TX'):
	if not vcf.contains('CSQ'):
		raise ValueError('VCF file has to be annotated with by VEP')
	# tracking_id class_code  nearest_ref_id  gene_id gene_short_name tss_id  locus   length  coverage    FPKM    FPKM_conf_lo    FPKM_conf_hi    FPKM_status
	# ENSG00000240361 -   -   ENSG00000240361 OR4G11P -   chr1:62947-63887    -   -   0   0   0   OK
	# ENSG00000268020 -   -   ENSG00000268020 AL627309.1  -   chr1:53048-54936    -   -   0   0   0   OK
	txfile = outfile.with_suffix('.tx_nopos')
	writer = TsvWriter(txfile)
	writer.cnames = ['tracking_id', 'class_code', 'nearest_ref_id', 'gene_id', 'gene_short_name', 'tss_id', 'locus', 'length', 'coverage', 'FPKM', 'FPKM_conf_lo', 'FPKM_conf_hi', 'FPKM_status']
	writer.writeHead()
	for variant in vcf:
		# try..except
		try:
			tx = variant.format('TX')[0]
		except (KeyError, TypeError):
			continue
		csqs = variant.INFO['CSQ'].split('|')
		txs = tx.split(',')
		for tx in txs:
			transcript, expr = tx.split('|', 1)
			csq = [csq for csq in csqs if f'|{transcript}|' in csq][0].split('|')
			r                 = TsvRecord()
			r.tracking_id     = csq[6]
			r.class_code      = '-'
			r.nearest_ref_id  = '-'
			r.gene_id         = csq[4]
			r.gene_short_name = csq[3]
			r.tss_id          = '-'
			r.locus           = '<pos>'
			r.length          = '-'
			r.coverage        = '-'
			r.FPKM            = expr
			r.FPKM_conf_lo    = 0
			r.FPKM_conf_hi    = 1000
			r.FPKM_status     = 'OK'
			writer.write(r)
			features.add(r.tracking_id)
	writer.close()

if gxfile or txfile:
	allpos = {}
	for gff in Gff(refall):
		if gff['type'] == 'gene':
			feature = gff['attributes']['gene_id']
		elif gff['type'] == 'transcript':
			feature = gff['attributes']['transcript_id']
		else:
			continue
		if feature not in features:
			continue
		allpos[feature] ='{}:{}-{}'.format(gff['seqid'], gff['start'], gff['end'])
if gxfile:
	gxfile2 = outfile.with_suffix('.gx')
	with open(gxfile) as fin, open(gxfile2, 'w') as fout:
		for line in fin:
			if '<pos>' not in line:
				fout.write(line)
			else:
				feature_id = line.split('\t', 1)[0]
				if feature_id not in allpos:
					logger.warning('Cannot find position information for: %s, skipping', feature_id)
				else:
					fout.write(line.replace('<pos>', allpos[feature_id]))
	gxfile = gxfile2
if txfile:
	txfile2 = outfile.with_suffix('.tx')
	with open(txfile) as fin, open(txfile2, 'w') as fout:
		for line in fin:
			if '<pos>' not in line:
				fout.write(line)
			else:
				feature_id = line.split('\t', 1)[0]
				if feature_id not in allpos:
					logger.warning('Cannot find position information for: %s, skipping', feature_id)
				else:
					fout.write(line.replace('<pos>', allpos[feature_id]))
	txfile = txfile2

params['rna-gene-fpkm-tracking-file'] = gxfile
params['rna-transcript-fpkm-tracking-file'] = txfile

shell.load_config(topiary = topiary)
if infile.endswith('.vcf') or infile.endswith('.vcf.gz'):
	params.vcf = infile
else:
	params.maf = infile

alleles = [allele.replace('*', '') for allele in afile]
params['mhc-alleles'] = ','.join(alleles)
params.genome = genome
params['output-csv'] = outfile.with_suffix('.nowt')
params['mhc-predictor'] = mhc_predictor

# make sure those mhc-predictors are in PATH
PATHs = set()
for mhcpred in (netmhc, netmhcpan, netmhciipan, netmhccons, smm, smm_pmbec):
	try:
		PATHs.add(str(Path(shell.which(mhcpred).str()).parent))
	except CmdyReturnCodeError:
		continue

params._env = Diot(PATH = environ['PATH'] + ':' + ':'.join(PATHs))

shell.topiary(**params).fg

# add wildtype binding
# #,variant,peptide_offset,peptide,allele,affinity,percentile_rank,prediction_method_name,peptide_length,gene,gene_id,transcript_id,transcript_name,effect,effect_type,contains_mutant_residues,mutation_start_in_peptide,mutation_end_in_peptide,gene_expression
# 0,chr6 g.31237146C>A,353,AACSNSAHG,HLA-A*02:01,35651.3,65.0,netMHC,9,HLA-C,ENSG00000204525,ENST00000383329,HLA-C-002,p.Q361H,Substitution,True,7,8,0.0
# 1,chr6 g.33037619G>T,40,AAFVQTHRT,HLA-A*02:01,22758.73,32.0,netMHC,9,HLA-DPA1,ENSG00000231389,ENST00000419277,HLA-DPA1-001,p.P49T,Substitution,True,8,9,0.0

if mhc_predictor in ('netmhc', 'netmhcpan', 'netmhciipan', 'netmhccons', 'smm', 'smm_pmbec'):

	wildpeps = set()
	mutpeps = {}
	tpreader = TsvReader(outfile.with_suffix('.nowt'), comment = '###', delimit = ',')
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
		shell.load_config(netmhc = netmhc)
		mhcallele2 = params['mhc-alleles'].replace(':', '').replace('*', '')

		wildfile = outfile.parent / 'wildtype.peptides.txt'
		wildfile.write_text('\n'.join(wildpeps))

		nparams = Diot(
			a = mhcallele2, v = True, inptype = 1, f = wildfile, _prefix = '-',
			_iter = True, _debug = True)
		res = shell.netmhc(**nparams).iter()
		pos_hit = False
		wildbindings = {allele: {} for allele in mhcallele2.split(',')}
		for line in res:
			if 'PEPLIST' not in line or line.startswith('Protein'):
				continue
			parts = line.split()
			wildbindings[parts[1]][parts[2]] = parts[12]

		writer = TsvWriter(outfile)
		writer.cnames = ['HLA_allele', 'Peptide', 'Affinity', 'Gene', 'ENSG', 'ENST', 'Ref_peptide', 'Ref_affinity', 'Mutation', 'AAChange']
		writer.writeHead()
		writerall = TsvWriter(outfile.with_suffix('.all.txt'))
		writerall.cnames = writer.cnames
		writerall.writeHead()
		tpreader.rewind()
		for r in tpreader:
			out              = TsvRecord()
			out.HLA_allele   = r.allele
			out.Peptide      = r.peptide
			out.Affinity     = r.affinity
			out.Gene         = r.gene
			out.ENSG         = r.gene_id
			out.ENST         = r.transcript_id
			wtpep            = mutpeps.get(r.peptide + '\t' + r.allele, '-')
			out.Ref_peptide  = wtpep
			out.Ref_affinity = wildbindings[r.allele.replace(':', '').replace('*', '')].get(wtpep, '>500')
			out.Mutation     = r.variant
			out.AAChange     = r.effect
			writerall.write(out)
			if float(out.Affinity) < 500 and ('>' in out.Ref_affinity or float(out.Ref_affinity) >= 2000):
				writer.write(out)

	def run_netmhcpan():
		shell.load_config(netmhcpan = netmhcpan)
		mhcallele2 = params['mhc-alleles'] if 'mhc-alleles' in params else ','.join(
			allele for allele in Path(params['mhc-alleles-file']).read_text().splitlines() if allele
		)

		wildfile = outfile.parent / 'wildtype.peptides.txt'
		wildfile.write_text('\n'.join(wildpeps))
		xlsfile = outfile.parent / 'wildtype.binding.txt'
		nparams = Diot(
			a = mhcallele2, v = True, BA = True, inptype = 1, f = wildfile, _prefix = '-',
			xls = True, xlsfile = xlsfile)
		shell.netmhcpan(**nparams).fg

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
