from os import path, makedirs
from shutil import rmtree
from sys import stderr, exit
import vcf

{{ runcmd }}
{{ mem2 }}
{{ params2CmdArgs }}
tmpdir    = path.join ("{{args.tmpdir}}", "{{proc.id}}.{{in.infile | fn}}.{{job.index}}")
if not path.exists (tmpdir): makedirs (tmpdir)

gz      = {{args.gz | lambda x: bool(x)}}
outfile = {{out.outfile | [:-3] | quote}} if gz else {{out.outfile | quote}}
tool    = {{args.tool | quote}}
dbpath  = {{args.dbpath}}[tool]
params  = {{args.params}}
try:
	###### snpeff
	{% if args.tool | lambda x: x == 'snpeff' %}
	mem = mem2({{args.mem | quote}}, 'java')
	params['dataDir'] = dbpath
	{% if args.snpeffStats %}
	params['csvStats'] = "{{out.outdir}}/{{in.infile | fn}}.stats.csv"
	{% else %}
	params['noStats'] = True
	{% endif %}
	params['i'] = 'vcf'
	params['o'] = 'vcf'
		
	cmd = '{{args.snpeff}} %s -Djava.io.tmpdir="%s" %s {{args.genome}} "{{in.infile}}" > "%s"' % (mem, tmpdir, params2CmdArgs(params, dash='-', equal=' '), outfile)
	runcmd (cmd)
	
	###### vep
	{% elif args.tool | lambda x: x == 'vep' %}
	params['i']        = {{in.infile | quote}}
	params['o']        = outfile
	params['format']   = 'vcf'
	params['vcf']      = True
	params['cache']    = True
	params['dir']      = dbpath
	params['assembly'] = {{args.genome | quote}}
	cmd = '{{args.vep}} %s' % (params2CmdArgs(params, equal=' '))
	runcmd (cmd)

	###### annovar
	{% elif args.tool | lambda x: x == 'annovar' %}
	# convert vcf to avinput
	avinput  = "{{out.outdir}}/{{in.infile | fn}}.avinput"
	avparams = {}

	avparams['includeinfo'] = True
	avparams['allsample']   = True
	avparams['withfreq']    = True
	avparams['format']      = 'vcf4'
	avparams['outfile']     = avinput

	cmd = '{{args.annovar_convert}} %s "{{in.infile}}"' % params2CmdArgs(avparams, equal=' ')
	runcmd (cmd)

	params['outfile']  = {{out.outfile | prefix | quote}}
	params['buildver'] = {{args.genome | quote}}

	cmd = '{{args.annovar}} %s "%s" "%s"' % (params2CmdArgs(params, equal=' '), avinput, dbpath) 
	runcmd (cmd)
		
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
			
	reader = vcf.Reader(filename="{{in.infile}}")
	reader.infos["ANN"] = vcf.parser._Info("ANN", 1, "String", "Annotation by ANNOVAR", "", "")
	snames = {v:k for k,v in enumerate(reader.samples)}
	writer = vcf.Writer(open(outfile, 'w'), reader)
	f2conv = "{{out.outfile | prefix}}.variant_function"
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
	{% endif %}

	if gz: runcmd ('gzip "%s"' % outfile)
except Exception as ex:
	stderr.write ("Job failed: %s" % str(ex))
	raise
finally:
	rmtree(tmpdir)