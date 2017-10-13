from os import path, makedirs, remove
from shutil import rmtree, copyfile, move
from sys import stderr
import vcf

{{ runcmd }}
{{ mem2 }}
{{ buildrefIndex }}
{{ params2CmdArgs }}

filters   = {{args.filters}}
selectors = {{args.selectors}}
params    = {{args.params}}
tmpdir    = {{args.tmpdir | quote}}
tmpdir    = path.join (tmpdir, "{{proc.id}}.{{in.infile | fn}}.{{job.index}}")
if not path.exists (tmpdir):
	makedirs (tmpdir)

if "type" in selectors:
	if not isinstance(selectors['type'], list):
		selectors['type'] = [selectors['type']]
	for i, t in enumerate(selectors['type']):
		selectors['type'][i] = t.lower()

if 'genotype' in selectors:
	for i, gt in selectors['genotype'].items():
		if not isinstance(gt, list):
			selectors['genotype'][i] = [gt]
		
		for j, t in enumerate(selectors['genotype'][i]):
			selectors['genotype'][i][j] = t.replace('|', '/')

if 'filter' in selectors:
	if not isinstance(selectors['filter'], list):
		selectors['filter'] = [selectors['filter']]
	
infile = "{{in.infile}}"
if selectors:
	mdfile = "{{job.outdir}}/{{in.infile | fn}}.selected.vcf"
	reader = vcf.Reader(filename=infile)
	writer = vcf.Writer(open(mdfile, 'w'), reader)
	for record in reader:
		reflen = len (record.REF)
		if "type" in selectors and sorted(selectors["type"]) != sorted(['snp', 'indel']):
			t = "snp"
			if reflen != 1:
				t = "indel"
			else:
				for a in record.ALT:
					if a is None or reflen != len(a):
						t = "indel"
						break
			if t not in selectors["type"]:
				continue
		if "genotype" in selectors:
			foundBadRecord = False
			for i, gt in selectors["genotype"].items():
				if str(i).isdigit():
					i = int(i)
					rgt = record.samples[i]['GT']
				else:
					rgt = record.genotype(i)['GT']
				a = rgt[0]
				b = rgt[-1]
				rgt = '%s/%s' % (a,b)	
				if not rgt in gt:
					foundBadRecord = True
					break
			if foundBadRecord:
				continue
		if "qual" in selectors:
			if record.QUAL < selectors['qual']:
				continue
		if "filter" in selectors:
			if record.FILTER and not (set(record.FILTER) & set(selectors['filter'])):
				continue
		writer.write_record(record)
	infile = mdfile
	writer.close()
	
try:	
	gz      = {{args.gz | lambda x: bool(x)}}
	tool    = {{args.tool | quote}}
	keep    = {{args.keep | lambda x: bool(x)}}
	outfile = {{out.outfile | [:-3] | quote}} if gz else {{out.outfile | quote}}
	if not filters:
		if not selectors:
			copyfile(infile, outfile)
		else:
			move(infile, outfile)
	else:
		if not isinstance(filters, dict):
			filters = {tool: filters}
		{% if args.tool | lambda x: x == 'gatk' %}
		ref = buildrefIndex (
			{{job.index}}, 
			"{{args.ref}}", 
			["{{args.ref}}.fai", "{{ args.ref | prefix }}.dict"], 
			'{{args.samtools}} faidx "{{args.ref}}"; {{args.picard}} CreateSequenceDictionary R="{{args.ref}}" O="{{args.ref | prefix}}.dict"',
			"{{ proc.workdir }}",
			'{{args.samtools}} faidx "{{job.outdir}}/{{args.ref | bn}}"; {{args.picard}} CreateSequenceDictionary R="{{job.outdir}}/{{args.ref | bn}}" O="{{job.outdir}}/{{args.ref | fn}}.dict"')
		
		# get filters
		i = 0
		for key, val in filters.items():
			argkey = 'filterExpression' + ' ' * i
			argval = 'filterName' + ' ' * i
			params[argkey] = key
			params[argval] = val
			i += 1

		mem       = mem2({{args.mem | quote}}, 'java')
		mdoutfile = outfile + 'WithFilter.vcf'

		params['num_threads'] = {{args.nthread}}
		params['R']  = {{args.ref | quote}}
		params['V']  = infile
		params['o']  = mdoutfile

		cmd = '{{args.gatk}} -T VariantFiltration %s -Djava.io.tmpdir="%s" %s' % (mem, tmpdir, params2CmdArgs(params, equal=' '))
		runcmd (cmd)

		if not keep:
			svparams = {}
			svparams['excludeFiltered'] = True
			svparams['num_threads']     = {{args.nthread}}
			svparams['R'] = {{args.ref | quote}}
			svparams['V'] = mdoutfile
			svparams['o'] = outfile

			cmd = '{{args.gatk}} -T SelectVariants %s -Djava.io.tmpdir="%s" %s' % (mem, tmpdir, params2CmdArgs(params, equal=' '))
			runcmd (cmd)
			remove (mdoutfile)
		else:
			move(mdoutfile, outfile)

		{% elif args.tool | lambda x: x == 'snpsift' %}
		# only one filter allowed
		key, val = filters.items()[0]
		mem = mem2({{ args.mem | quote }}, 'java')
		mdoutfile = outfile + 'WithFilter.vcf'
		params['n'] = True
		params['f'] = infile
		params['i'] = key

		cmd = '{{args.snpsift}} filter %s -Djava.io.tmpdir="%s" %s "%s" > "%s"' % (mem, tmpdir, params2CmdArgs(params, equal=' '), val, mdoutfile)
		runcmd (cmd)

		if not keep:
			cmd = '{{args.bcftools}} view -O v -f .,PASS "%s" > "%s"' % (mdoutfile, outfile)
			runcmd (cmd)
			remove (mdoutfile)
		else:
			move(mdoutfile, outfile)

		{% elif args.tool | lambda x: x == 'bcftools' %}
		params['O'] = 'v'
		params['o'] = outfile
		params['e'] = filters.values()[0]

		cmd = '{{args.bcftools}} view %s "%s"' % (params2CmdArgs(params), infile)
		runcmd (cmd)

		{% elif args.tool | lambda x: x == 'vcflib' %}
		i = 0
		for _, val in filters.items():
			argkey = 'f' + ' ' * i
			params[argkey] = "! ( %s )" % val
			i += 1

		cmd = '{{args.vcflib_vcffilter}} %s "%s" > "%s"' % (params2CmdArgs(params), infile, outfile)
		runcmd (cmd)
		{% endif %}
		
		if filters:
			pass
			#remove (infile)
	if gz: runcmd ('gzip "%s"' % outfile)
except Exception as ex:
	stderr.write ("Job failed: %s" % str(ex))
	raise
finally:
	rmtree(tmpdir)