from os import path, makedirs
from pyppl import Box
from bioprocs.utils import shell, mem2

infile          = {{i.infile | quote}}
outfile         = {{o.outfile | quote}}
outdir          = {{o.outdir | quote}}
tmpdir          = {{args.tmpdir | quote}}
tool            = {{args.tool | quote}}
snpeff          = {{args.snpeff | quote}}
vep             = {{args.vep | quote}}
gz              = {{args.gz | repr}}
vcfanno         = {{args.vcfanno | quote}}
annovar         = {{args.annovar | quote}}
annovar_convert = {{args.annovar_convert | quote}}
genome          = {{args.genome | quote}}
nthread         = {{args.nthread | repr}}
dbs             = {{args.dbs | repr}}
dbs             = dbs[tool]
snpeffStats     = {{args.snpeffStats | repr}}
params          = {{args.params | repr}}
mem             = {{args.mem | repr}}
tmpdir          = path.join(tmpdir, "{{proc.id}}.{{i.infile | fn}}.{{job.index}}")
if not path.exists(tmpdir):
	makedirs(tmpdir)
if gz:
	outfile = outfile[:-3]

shell.TOOLS.vcfanno         = vcfanno
shell.TOOLS.annovar         = annovar
shell.TOOLS.annovar_convert = annovar_convert
shell.TOOLS.vep             = vep
shell.TOOLS.snpeff          = snpeff

def run_snpeff():
	if not path.exists(dbs):
		raise ValueError('Database does not exist: {}'.format(dbs))
	snpeff = shell.Shell(subcmd = True, dash = '-', equal = ' ').snpeff
	params.dataDir = dbs
	if snpeffStats:
		params.csvStats = path.join(outdir, path.basename(infile) + '.stats.csv')
	else:
		params.noStats = True
	params.i = 'vcf'
	params.o = 'vcf'
	params["Djava.io.tmpdir=%s" % tmpdir] = True
	params._ = [genome, infile]
	params._stdout = outfile
	snpeff.ann(**params).run()
	if gz: shell.bgzip(outfile)

def run_vep():
	if not path.exists(dbs):
		raise ValueError('Database does not exist: {}'.format(dbs))
	vep             = shell.Shell(equal = ' ').vep
	params.i        = infile
	params.o        = outfile
	params.format   = 'vcf'
	params.vcf      = True
	params.cache    = True
	params.dir      = dbs
	params.assembly = genome
	vep(**params).run()
	if gz: shell.bgzip(outfile)

def run_annovar():
	if not path.exists(dbs):
		raise ValueError('Database does not exist: {}'.format(dbs))
	import vcf
	avinput              = path.join(outdir, infile + '.avinput')
	avparams             = Box()
	avparams.includeinfo = True
	avparams.allsample   = True
	avparams.withfreq    = True
	avparams.format      = 'vcf4'
	avparams.outfile     = avinput
	avparams._           = infile
	shell.Shell(equal = ' ').annovar_convert(**avparams).run()

	params.buildver = genome
	params._        = [avinput, dbs]
	params.outfile  = path.splitext(outfile)[0]
	shell.Shell(equal = ' ').annovar(**params).run()
	# convert .variant_function to vcf
	# like snpEff ann, add ANN field to vcf
	def info2dict(info):
		ret   = {}
		infos = info.split(';')
		for i in infos:
			if '=' in i:
				ii = i.split('=')
				if ',' in ii[1]:
					ret[ii[0]] = ii[1].split(',')
				else:
					ret[ii[0]] = ii[1]
			else:
				ret[i] = True
		return ret
			
	def rs2record(rs):
		r       = rs[0]
		CHROM   = r[2]
		POS     = int(r[3])
		ID      = r[12]
		REF     = r[5]
		ALT     = vcf.model._Substitution(r[6])
		QUAL    = r[15]
		FILTER  = r[16]
		INFO    = info2dict(r[17])
		FORMAT  = r[18]
		samples = r[19:]
		
		anns = []
		alts = []
		for i, lr in enumerate(rs):
			alts.append(vcf.model._Substitution(lr[6]))
			ann = [lr[6]]
			ann.append(lr[0])  #Annotation 
			ann.append('')     #impact
			gene = ''
			dist = '0'
			for genes in lr[1].split(','):
				if genes.startswith('NONE'): continue
				genes = genes.split('(dist=')
				gene  = genes[0]
				if len(genes)>1:
					dist  = genes[1][:-1]
				break
			ann.append(gene)
			ann.append('')    # geneid
			ann.append('')    # Feature type
			ann.append('')    # Feature ID
			ann.append('')    # Transcript biotype
			ann.append('1/1') # Rank / total
			ann.append('')    # HGVS.c
			ann.append('')    # HGVS.p
			ann.append('')    # cDNA_position
			ann.append('')    # CDS_position
			ann.append('')    # Protein_position
			ann.append(dist)  # Distance to feature
			ann.append('')  # Errors, Warnings or Information messages
			anns.append ('|'.join(ann))
		INFO['ANN'] = anns
		
		record  = vcf.model._Record(CHROM, POS, ID, REF, alts, QUAL, FILTER, INFO, FORMAT, snames)
		record.samples = reader._parse_samples (samples, FORMAT, record)
		return record
			
	reader = vcf.Reader(filename="{{i.infile}}")
	reader.infos["ANN"] = vcf.parser._Info("ANN", 1, "String", "Annotation by ANNOVAR", "", "")
	snames = {v:k for k,v in enumerate(reader.samples)}
	writer = vcf.Writer(open(outfile, 'w'), reader)
	f2conv = "{{o.outfile | prefix}}.variant_function"
	lastvid= ''
	lastr  = []
	with open (f2conv) as f:
		for line in f:
			line = line.strip("\r\n")
			if not line: continue
			
			parts = line.split("\t")
			varid = parts[2] + '|' + parts[3] + '|' + parts[12] + '|' + parts[5]
			
			if lastvid != varid and lastvid:
				record  = rs2record(lastr)
				writer.write_record(record)
				lastvid = varid
				lastr   = [parts]
			else:
				lastvid = varid
				lastr.append (parts)
				
	record  = rs2record(lastr)
	writer.write_record(record)		
	writer.close()

	if gz: shell.bgzip(outfile)

def run_vcfanno():
	# compose toml file
	toml = path.join(outdir, 'config.toml')
	with open(toml, 'w') as f:
		for config in dbs:
			f.write("\n[[annotation]]\n")
			for key, val in config.items():
				f.write('{key}={val!r}\n'.format(key = key, val = val))
	params.p       = nthread
	params._       = [toml, infile]
	params._stdout = outfile
	shell.Shell(dash = '-', equal = ' ').vcfanno(**params).run()
	if gz: shell.bgzip(outfile)

tools = dict(
	snpeff  = run_snpeff,
	vep     = run_vep,
	annovar = run_annovar,
	vcfanno = run_vcfanno
)


try:
	tools[tool]()
except KeyError:
	raise ValueError('Tool {!r} not supported yet.'.format(tool))
except:
	raise
finally:
	shell.rmrf(tmpdir)
