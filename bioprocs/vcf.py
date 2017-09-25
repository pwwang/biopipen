"""
A set of processes to generate/process vcf files
"""

from pyppl import Proc, Box
from .utils import mem, runcmd, buildArgIndex, buildArgsFastaFai, buildArgsFastaDict, checkArgsRef, cbindFill, plot

"""
@name:
	pVcfFilter
@description:
	Filter records in vcf file.
@input:
	`infile:file`: The input file
@output:
	`outfile:file`: The output file
@args:
	`tool`: Which tool to use for filtering. Default: 'vcflib'
	`vcflib_vcffilter`: The path of vcffilter from vcflib. Default: 'vcffilter'
	`gatk`            : The path of gatk. Default: 'gatk'
	`snpsift`         : The path of snpsift. Default: 'SnpSift'
	`bcftools`        : The path of bcftools. Default: 'bcftools'
	`samtools`        : The path of samtools. Default: 'samtools' (used by gatk to generate reference index)
	`picard`          : The path of picard. Default: 'picard' (used by picard to generate reference dict) 
	`params`          : Other params of `tool`. Default: ""
	`mem`             : The memory to be used. Default: "4G" (only for snpsift and gatk)
	`gz`              : Whether to gzip the output file. Default: False
	`keep`            : Whether to keep the filtered records. Default: True. (only for gatk, snpsift at filter step)
	`ref`             : The path of reference. Default: "" (for gatk)
	`tmpdir`          : The path of tmpdir. Default: <system tmpdir> (only used by gatk and snpsift)
	`nthread`         : The path of Default: 1	
	`selectors`:   Select records by:
	  - type (snp, indel), sample genotypes (0, 1, 2), min genotype quality, filter (PASS, .)
	  - for example:
	    ```
		{"type": "snp", "genotype": {0: '0/0'}, "qual": 30}
		to select snps and whose genotype is '0/0' in 1st sample with quality >= 30
		{"genotype": {0: ['1/1', '0|1']}, "filter": ["PASS"]}
		to select records with PASS and genotype in 1st sample is '1/1' or '0/1'
		```
	`filters`:     Filters depend on the tool you use on INFO filelds
	  - format: `{"name1": "expression1", ...}`
	  - If a string is specified, will convert to `{<tool name>: <expression>}`
	  - Remember it filters OUT the records when ANY of the expression is true
@requires:
	[`pyvcf`](https://github.com/jamescasbon/PyVCF)
	[`gatk`](https://software.broadinstitute.org/gatk)
	[`bcftools`](http://www.htslib.org/doc/bcftools-1.2.html)
	[`snpsift`](http://snpeff.sourceforge.net/SnpSift.version_4_0.html)
	[`samtools`](https://github.com/samtools/samtools) if `args.ref` is not indexed, or bamutil is used for bam index file generation.
	[`picard`](https://broadinstitute.github.io/picard/command-line-overview.html) if `args.ref is not dicted.`
"""
pVcfFilter                       = Proc(desc = 'Filter records in vcf file.')
pVcfFilter.input                 = "infile:file"
pVcfFilter.output                = "outfile:file:{{in.infile | fn}}.vcf{{args.gz | lambda x: '.gz' if x else ''}}"
pVcfFilter.args.tool             = 'vcflib' # vcflib,    gatk, snpsift
pVcfFilter.args.vcflib_vcffilter = 'vcffilter'
pVcfFilter.args.gatk             = 'gatk'
pVcfFilter.args.snpsift          = 'SnpSift'
pVcfFilter.args.bcftools         = 'bcftools'
pVcfFilter.args.samtools         = 'samtools'
pVcfFilter.args.picard           = 'picard'
pVcfFilter.args.selectors        = Box()
pVcfFilter.args.filters          = Box()
pVcfFilter.args.params           = ""
pVcfFilter.args.mem              = "4G"
pVcfFilter.args.gz               = False
pVcfFilter.args.keep             = True # only for gatk, snpsift at filter step
pVcfFilter.args.ref              = "" # gatk
pVcfFilter.args.tmpdir           = __import__('tempfile').gettempdir()
pVcfFilter.args.nthread          = 1
pVcfFilter.tplenvs.memtoJava     = mem.toJava.python
pVcfFilter.tplenvs.runcmd        = runcmd.python
pVcfFilter.tplenvs.buildArgIndex = buildArgIndex.python
pVcfFilter.lang                  = "python"
pVcfFilter.script                = """
from os import path, makedirs, remove
from shutil import rmtree, copyfile, move
from sys import stderr, exit
import vcf

{{ runcmd }}
{{ memtoJava }}
{{ buildArgIndex }}
filters = {{args.filters | json}}
selectors = {{args.selectors | json}}
tmpdir    = {{ args.tmpdir | quote }}
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
			ref = buildArgIndex (
				{{job.index}}, 
				"{{args.ref}}", 
				["{{args.ref}}.fai", "{{ args.ref | prefix }}.dict"], 
				'{{args.samtools}} faidx "{{args.ref}}"; {{args.picard}} CreateSequenceDictionary R="{{args.ref}}" O="{{args.ref | prefix}}.dict"',
				"{{ proc.workdir }}",
				'{{args.samtools}} faidx "{{job.outdir}}/{{args.ref | bn}}"; {{args.picard}} CreateSequenceDictionary R="{{job.outdir}}/{{args.ref | bn}}" O="{{job.outdir}}/{{args.ref | fn}}.dict"')
			# get filters
			gatkfilters = ['--filterExpression "%s" --filterName "%s"' % (val, key) for key, val in filters.items()]
			gatkfilters = ' '.join(gatkfilters)
			mem = memtoJava({{ args.mem | quote }})
			mdoutfile = outfile + 'WithFilter.vcf'
			cmd = '{{args.gatk}} -T VariantFiltration %s -Djava.io.tmpdir="%s" -nt {{args.nthread}} -R "{{args.ref}}" -V "%s" -o "%s" %s {{args.params}}' % (mem, tmpdir, infile, mdoutfile, gatkfilters)
			runcmd (cmd)
			if not keep:
				cmd = '{{args.gatk}} -T SelectVariants --excludeFiltered %s -Djava.io.tmpdir="%s" -nt {{args.nthread}} -R "{{args.ref}}" -V "%s" -o "%s"' % (mem, tmpdir, mdoutfile, outfile)
				runcmd (cmd)
				remove (mdoutfile)
			else:
				move(mdoutfile, outfile)
		{% elif args.tool | lambda x: x == 'snpsift' %}
			# only one filter allowed
			(key, val) = filters.items()[0]
			mem = memtoJava({{ args.mem | quote }})
			mdoutfile = outfile + 'WithFilter.vcf'
			cmd = '{{args.snpsift}} filter %s -Djava.io.tmpdir="%s" -n -f "%s" -i "%s" "%s" > "%s"' % (mem, tmpdir, infile, key, val, mdoutfile)
			runcmd (cmd)
			if not keep:
				cmd = '{{args.bcftools}} view -O v -f .,PASS "%s" > "%s"' % (mdoutfile, outfile)
				runcmd (cmd)
				remove (mdoutfile)
			else:
				move(mdoutfile, outfile)
		{% elif args.tool | lambda x: x == 'bcftools' %}
			cmd = '{{args.bcftools}} view -O v -o "%s" -e "%s" "%s"' % (outfile, filters.values()[0], infile)
			runcmd (cmd)
		{% elif args.tool | lambda x: x == 'vcflib' %}
			vcffilters = ['-f "! ( %s )"' % val for key, val in filters.items()]
			cmd = '{{args.vcflib_vcffilter}} %s "%s" > "%s"' % (' '.join(vcffilters), infile, outfile)
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
"""

"""
@name:
	pVcfAnno
@description:
	Annotate the variants in vcf file.
	You have to prepare the databases for each tool.
@input:
	`infile:file`: The input vcf file
@output:
	`outfile:file`: The output file (output file of annovar will also be converted to vcf)
	`outdir`: The output directory, used to fetch some stat/summary files
@args:
	`tool`:            The tool used to do annotation. Default: snpeff
	`snpeff`:          The path of snpeff. Default: snpEff
	`vep`:             The path to vep. Default: vep
	`gz`:              Whether to gzip the result file. Default: False
	`annovar`:         The path of annovar. Default: annotate_variation.pl
	`annovar_convert`: The path of convert2annovar.pl, used to convert vcf to annovar input file. Default: convert2annovar.pl
	`genome`:          The genome for annotation. Default: hg19
	`tmpdir`:          The tmpdir, mainly used by snpeff. Default: <system tmpdir>
	`dbpath`:          The path of database for each tool. Required by 'annovar' and 'vep'
	`params`:          Other params for tool. Default: ''
	`snpeffStats`:     Whether to generate stats file when use snpeff. Default: False
	`mem`:             The memory used by snpeff. Default: '4G'
@requires:
	[`annovar`](http://doc-openbio.readthedocs.io/projects/annovar/en/latest/)
	[`snpeff`](http://snpeff.sourceforge.net/SnpEff_manual.html#intro)
	[`vep`](http://www.ensembl.org/info/docs/tools/vep/script/vep_tutorial.html)
"""
pVcfAnno                      = Proc(desc = 'Annotate the variants in vcf file.')
pVcfAnno.input                = "infile:file"
pVcfAnno.output               = "outfile:file:{{in.infile | fn}}.{{args.tool}}.vcf{{args.gz | lambda x: '.gz' if x else ''}}, outdir:{{job.outdir}}"
pVcfAnno.args.tool            = 'snpeff'
pVcfAnno.args.snpeff          = 'snpEff'
pVcfAnno.args.vep             = 'vep'
pVcfAnno.args.gz              = False
pVcfAnno.args.annovar         = 'annotate_variation.pl'
pVcfAnno.args.annovar_convert = 'convert2annovar.pl'
pVcfAnno.args.genome          = 'hg19'
pVcfAnno.args.tmpdir          = __import__('tempfile').gettempdir()
pVcfAnno.args.dbpath          = ''
pVcfAnno.args.snpeffStats     = False
pVcfAnno.args.params          = ''
pVcfAnno.args.mem             = '4G'
pVcfAnno.tplenvs.runcmd       = runcmd.python
pVcfAnno.tplenvs.memtoJava    = mem.toJava.python
pVcfAnno.beforeCmd            = """
# check dbpath
if [[ "{{args.tool}}" != "snpeff" ]] && [[ ! -e "{{args.dbpath}}" ]]; then
	echo "You have to specify args.dbpath." 1>&2 
	echo "  - For vep: /path/to/cache" 1>&2
	echo "  - For snpEff: /path/to/datadir" 1>&2
	echo "  - For annovar: /path/to/db" 1>&2
	exit 1
fi
"""
pVcfAnno.lang                 = 'python'
pVcfAnno.script               = """
from os import path, makedirs
from shutil import rmtree
from sys import stderr, exit
import vcf

{{ runcmd }}
{{ memtoJava }}
tmpdir    = path.join ("{{args.tmpdir}}", "{{proc.id}}.{{in.infile | fn}}.{{job.index}}")
if not path.exists (tmpdir): makedirs (tmpdir)

gz      = {{args.gz | lambda x: bool(x)}}
outfile = {{out.outfile | [:-3] | quote}} if gz else {{out.outfile | quote}}
tool    = {{args.tool | quote}}
try:
	if tool == 'snpeff':
		datadir = ''
		if path.exists ("{{args.dbpath}}"):
			datadir = '-dataDir "{{args.dbpath}}"'
		mem = memtoJava({{ args.mem | quote }})
		stats = '-csvStats "{{job.outdir}}/{{in.infile | fn}}.stats.csv"'
		if not {{args.snpeffStats | lambda x: bool(x)}}:
			stats = '-noStats'
			
		cmd = '{{args.snpeff}} %s -Djava.io.tmpdir="%s" "%s" {{args.params}} -i vcf -o vcf %s {{args.genome}} "{{in.infile}}" > "%s"' % (mem, tmpdir, stats, datadir, outfile)
		runcmd (cmd)
		
	elif tool == 'vep':
		cmd = '{{args.vep}} -i "{{in.infile}}" -o "%s" --format vcf --vcf --cache --dir "{{args.dbpath}}" --assembly {{args.genome}} {{args.params}}' % (outfile)
		runcmd (cmd)
	elif tool == 'annovar':
		# convert vcf to avinput
		avinput = "{{job.outdir}}/{{in.infile | fn}}.avinput"
		cmd = '{{args.annovar_convert}} --includeinfo --allsample --withfreq --format vcf4 "{{in.infile}}" --outfile "{{job.outdir}}/{{in.infile | fn}}.avinput"'
		runcmd (cmd)
		cmd = '{{args.annovar}} {{args.params}} --outfile "{{out.outfile | prefix}}" --buildver {{args.genome}} "%s" "{{args.dbpath}}"' % (avinput) 
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
				line = line.strip("\\r\\n")
				if not line: continue
				
				parts = line.split("\\t")
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
		
	if gz: runcmd ('gzip "%s"' % outfile)
except Exception as ex:
	stderr.write ("Job failed: %s" % str(ex))
	raise
finally:
	rmtree(tmpdir)
"""

"""
@name:
	pCallRate
@description:
	Calculate sample/snp call rate from single sample vcfs
@input:
	`indir:file`:     The dir containing the vcfs
@output:
	`outsample:file`: The report of call rate for each sample
	`figsample:file`: The bar chat of sample call rates
	`outsnp:file`:    The report of call rate for each snp
	`figsnp:file`:    The bar chat of snp call rates
"""
pCallRate        = Proc()
pCallRate.input  = "indir:file"
pCallRate.output = [
	"outsample:file:{{in.indir | fn}}.sampleCallRate.txt",
	"figsample:file:{{in.indir | fn}}.sampleCallRate.png",
	"outsnp:file:{{in.indir | fn}}.snpCallRate.txt",
	"figsnp:file:{{in.indir | fn}}.snpCallRate.png"
]
pCallRate.tplenvs.runcmd    = runcmd.r
pCallRate.tplenvs.cbindFill = cbindFill.r
pCallRate.tplenvs.plotHist  = plot.hist.r
pCallRate.lang            = "Rscript"
pCallRate.script          = """
{{runcmd}}
{{cbindFill}}
{{plotHist}}

setwd("{{in.indir}}")
files = list.files(pattern = "*.vcf")

data  = NULL

for (file in files) {
	sample  = tools::file_path_sans_ext(tools::file_path_sans_ext(file))
	tmp     = runcmd (paste('awk', '\\'BEGIN{OFS="\\\\t"} $0 !~ /^#/ {print $1":"$2,1}\\'', file), intern=T)
	con     =  textConnection(tmp)
	tmp     = read.table (con, header=F, row.names=1, check.names=F)
	close(con)
	colnames(tmp) = c(sample)
	data    = cbindFill(data, tmp)
}

samplecr = colSums(data) / nrow(data)
snpcr    = rowSums(data) / ncol(data)

# sample/snp call rate
write.table (samplecr, "{{out.outsample}}", quote=F, sep="\\t", col.names=F)
plotHist (samplecr, "{{out.figsample}}", main="Sample call rate", ylab = "# samples")
write.table (snpcr, "{{out.outsnp}}", quote=F, sep="\\t", col.names=F)
plotHist (snpcr, "{{out.figsnp}}", main="SNP call rate", ylab = "# snps")
"""
