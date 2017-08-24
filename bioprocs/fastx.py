"""
A set of processes to generate/process fastq/fasta files
"""

from pyppl import proc
from .utils import mem, runcmd, buildArgIndex, checkArgsRef

"""
@name:
	pFastqPESim
@description:
	Simulate reads
@input:
	`in`: Index of the job/simulation, typically use range(10) for 10-time simulations
@output:
	`fq1:file`: The first pair read file
	`fq2:file`: The second pair read file
@args:
	`tool`:  The tool used for simulation. Default: wgsim (dwgsim)
	`len1`:  The length of first pair read. Default: 100
	`len2`:  The length of second pair read. Default: 100
	`num`:   The number of read PAIRs. Default: 1000000
	`seed`:  The seed for randomization. Default: None
	`gz`:    Whether generate gzipped read file. Default: True
	`wgsim`: The path of wgsim. Default: wgsim
	`dwgsim`:The path of wgsim. Default: dwgsim
	`ref`:   The reference genome. Required
	`params`:Other params for `tool`. Default: ""
@requires:
	[`wgsim`](https://github.com/lh3/wgsim)
"""
pFastqPESim              = proc (desc = 'Simulate reads')
pFastqPESim.input        = "in"
pFastqPESim.output       = "fq1:file:read{{in}}_1.fq{{args.gz | lambda x: '.gz' if x else ''}}, fq2:file:read{{in}}_2.fq{{args.gz | lambda x: '.gz' if x else ''}}"
pFastqPESim.args.tool    = 'wgsim'
pFastqPESim.args.wgsim   = 'wgsim'
pFastqPESim.args.dwgsim  = 'dwgsim'
pFastqPESim.args.len1    = 100
pFastqPESim.args.len2    = 100
pFastqPESim.args.num     = 1000000
pFastqPESim.args.gz      = True
pFastqPESim.args.seed    = None
pFastqPESim.args.params  = ""
pFastqPESim.args.ref     = ""
pFastqPESim.lang         = "python"
pFastqPESim.args._runcmd = runcmd.python
pFastqPESim.script       = """
from os import path
from shutil import move
from sys import stderr, exit

tool = {{args.tool | quote}}
ref  = {{args.ref | quote}}
if not ref or not path.exists (ref):
	stderr.write ("Reference file is not spcified or not exists.")
	exit (1)
	
{{args._runcmd}}

fq1 = {{fq1 | quote}}
fq2 = {{fq2 | quote}}
try:
	if tool == 'wgsim':
		if {{args.gz}}:
			fq1 = "{{fq1 | [:-3]}}"
			fq2 = "{{fq2 | [:-3]}}"
		cmd = '{{args.wgsim}} {{args.params}} -N {{args.num}} -1 {{args.len1}} -2 {{args.len2}} -S {{args.seed | lambda x: -1 if x is None else x}} "%s" "%s" "%s"' % (ref, fq1, fq2)
		runcmd (cmd)
		if {{args.gz}}:
			runcmd ('gzip "%s"' % fq1)
			runcmd ('gzip "%s"' % fq2)
	elif tool == 'dwgsim':
		prefix = {{fq1 | [:-5] | quote}}
		if {{args.gz}}:
			fq1 = "{{fq1 | [:-3]}}"
			fq2 = "{{fq2 | [:-3]}}"
			prefix = "{{fq1 | [:-8]}}"
		cmd = '{{args.dwgsim}} {{args.params}} -N {{args.num}} -1 {{args.len1}} -2 {{args.len2}} {{args.seed | lambda x: '' if x is None else "-z " + str(x)}} "%s" "%s"' % (ref, prefix)
		runcmd (cmd)
		if path.exists (prefix + '.bwa.read1.fastq.gz'):
			move (prefix + '.bwa.read1.fastq.gz', prefix + '_1.fq.gz')
			move (prefix + '.bwa.read2.fastq.gz', prefix + '_2.fq.gz')
			if not {{args.gz}}:
				runcmd ('gunzip "%s"' % (prefix + '_1.fq.gz'))
				runcmd ('gunzip "%s"' % (prefix + '_2.fq.gz'))
		else:
			move (prefix + '.bwa.read1.fastq', prefix + '_1.fq')
			move (prefix + '.bwa.read2.fastq', prefix + '_2.fq')
			if {{args.gz}}:
				runcmd ('gzip "%s"' % (prefix + '_1.fq'))
				runcmd ('gzip "%s"' % (prefix + '_2.fq'))
except Exception as ex:
	stderr.write ("Job failed: %s" % str(ex))
	raise
"""

"""
@name:
	pFastQC
@description:
	QC report for fastq file
@input:
	`fq:file`:    The fastq file (also fine with gzipped)
@output:
	`outdir:dir`: The output direcotry
@args:
	`tool`:    The tool used for simulation. Default: fastqc 
	`fastqc`:  The path of fastqc. Default: fastqc
	`nthread`: Number of threads to use. Default: 1
	`params`:Other params for `tool`. Default: ""
@requires:
	[`fastqc`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
"""
pFastQC              = proc (desc = 'QC report for fastq file')
pFastQC.input        = "fq:file"
pFastQC.output       = "outdir:dir:{{fq | fn | fn}}" # if it has .gz 
pFastQC.args.tool    = 'fastqc'
pFastQC.args.fastqc  = 'fastqc'
pFastQC.args.nthread = 1
pFastQC.args.params  = ""
pFastQC.lang         = "python"
pFastQC.args._runcmd = runcmd.python
pFastQC.script       = """
from sys import stderr, exit
	
{{args._runcmd}}

fq   = {{fq | quote}}
tool = {{args.tool | quote}}
try:
	if tool == 'fastqc':
		cmd = '{{args.fastqc}} "{{fq}}" -o "{{outdir}}" -p'
		runcmd(cmd)
	else:
		raise Exception('Tool %s not supported.' % (tool))
except Exception as ex:
	stderr.write ("Job failed: %s" % str(ex))
	raise
"""

"""
@name:
	pFastMC
@description:
	Multi-QC based on pFastQC
@input:
	`qcdir:file`:  The direcotry containing QC files
@output:
	`outdir:dir`: The output direcotry
@args:
	`tool`:    The tool used for simulation. Default: multiqc 
	`multiqc`: The path of fastqc. Default: multiqc
	`params`:  Other params for `tool`. Default: ""
@requires:
	[`multiqc`](http://multiqc.info/)
"""
pFastMC              = proc (desc = 'Multi-QC based on pFastQC')
pFastMC.input        = "qcdir:file"
pFastMC.output       = "outdir:dir:{{qcdir | fn}}_multiqc_{{#}}"
pFastMC.args.tool    = 'multiqc'
pFastMC.args.multiqc = 'multiqc'
pFastMC.args.params  = ""
pFastMC.lang         = "python"
pFastMC.args._runcmd = runcmd.python
pFastMC.script       = """
from sys import stderr, exit

	
{{args._runcmd}}

qcdir = {{qcdir | quote}}
tool  = {{args.tool | quote}}
try:
	if tool == 'multiqc':
		cmd = '{{args.multiqc}} "{{qcdir}}" -o "{{outdir}}" -p'
		runcmd(cmd)
	else:
		raise Exception('Tool %s not supported.' % (tool))
except Exception as ex:
	stderr.write ("Job failed: %s" % str(ex))
	raise
"""

"""
@name:
	pFastqPETrim
@description:
	Trim pair-end FASTQ reads
@input:
	`fq1:file`:  The input fastq file
	`fq2:file`:  The input fastq file
@output:
	`outfq1:file`: The trimmed fastq file
	`outfq2:file`: The trimmed fastq file
@args:
	`tool`        : The tools used for trimming. Default: trimmomatic (cutadapt|skewer)
	`cutadapt`    : The path of seqtk. Default: cutadapt
	`skewer`      : The path of fastx toolkit trimmer. Default: skewer
	`trimmomatic` : The path of trimmomatic. Default: trimmomatic
	`params`      : Other params for `tool`. Default: ""
	`nthread`     : Number of threads to be used. Default: 1
	- Not for cutadapt
	`gz`          : Whether gzip output files. Default: True
	`mem`         : The memory to be used. Default: 4G
	- Only for trimmomatic
	`minlen`      : Discard trimmed reads that are shorter than `minlen`. Default: 18
	- For trimmomatic, the number will be `minlen`*2 for MINLEN, as it filters before trimming
	`minq`        : Minimal mean qulity for 4-base window or leading/tailing reads. Default: 3
	`cut5`        : Remove the 5'end reads if they are below qulity. Default: 3
	`cut3`        : Remove the 3'end reads if they are below qulity. Default: 3
	- Not for skewer
	`adapter1`    : The adapter for sequence. Default: AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
	`adapter2`    : The adapter for pair-end sequence. Default: AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA
@requires:
	[`cutadapt`](http://cutadapt.readthedocs.io/en/stable/guide.html)
	[`skewer`](https://github.com/relipmoc/skewer)
	[`trimmomatic`](https://github.com/timflutre/trimmomatic)
"""
pFastqPETrim                  = proc (desc = 'Trim pair-end reads in fastq file.')
pFastqPETrim.input            = "fq1:file, fq2:file"
pFastqPETrim.output           = [
	"outfq1:file:{{fq1 | fn | lambda x: x.rpartition('.')[0] if x.endswith('fastq') or x.endswith('fq') else x}}.fq{{args.gz | lambda x: '.gz' if x else ''}}",
	"outfq2:file:{{fq2 | fn | lambda x: x.rpartition('.')[0] if x.endswith('fastq') or x.endswith('fq') else x}}.fq{{args.gz | lambda x: '.gz' if x else ''}}"
]
pFastqPETrim.lang             = "python"
pFastqPETrim.args.tool        = 'skewer'
pFastqPETrim.args.cutadapt    = 'cutadapt'
pFastqPETrim.args.skewer      = 'skewer'
pFastqPETrim.args.trimmomatic = 'trimmomatic'
pFastqPETrim.args.params      = ''
pFastqPETrim.args.nthread     = 1
pFastqPETrim.args.gz          = True
pFastqPETrim.args.mem         = '4G'
pFastqPETrim.args.minlen      = 18
pFastqPETrim.args.minq        = 3
pFastqPETrim.args.cut5        = 3
pFastqPETrim.args.cut3        = 3
pFastqPETrim.args.adapter1    = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'
pFastqPETrim.args.adapter2    = 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA'
pFastqPETrim.args._memtoJava  = mem.toJava.python
pFastqPETrim.args._runcmd     = runcmd.python
pFastqPETrim.script           = """
from sys import stderr, exit
from shutil import move

{{args._runcmd}}
{{args._memtoJava}}

tool = {{args.tool | quote}}
try:
	if tool == 'trimmomatic':
		def seqrev (seq):
			d = {
				'A':'T', 'T':'A', 'G':'C', 'C':'G', 
				'a':'t', 't':'a', 'g':'c', 'c':'g'
			}
			return ''.join([d[s] for s in seq])
		
		mem    = memtoJava ({{args.mem | quote}})
		minlen = str({{args.minlen}} * 2)
		adfile = "{{job.outdir}}/adapters.fa"
		with open (adfile, "w") as ad:
			ad.write (">PE1\\n")
			ad.write (seqrev({{args.adapter1 | quote}}) + "\\n")
			ad.write (">PE1_rc\\n")
			ad.write ({{args.adapter1 | quote}} + "\\n")
			ad.write (">PE2\\n")
			ad.write (seqrev({{args.adapter2 | quote}}) + "\\n")
			ad.write (">PE2_rc\\n")
			ad.write ({{args.adapter2 | quote}} + "\\n")
		cmd    = '{{args.trimmomatic}} %s PE -threads {{args.nthread}} "{{fq1}}" "{{fq2}}" "{{outfq1}}" /dev/null "{{outfq2}}" /dev/null ILLUMINACLIP:%s:2:30:10 LEADING:{{args.cut5}} TRAILING:{{args.cut3}} SLIDINGWINDOW:4:{{args.minq}} MINLEN:%s {{args.params}} ' % (mem, adfile, minlen)
		runcmd (cmd)
	elif tool == 'cutadapt':
		cmd = '{{args.cutadapt}} -a {{args.adapter1}} -A {{args.adapter2}} -u {{args.cut5}} -u -{{args.cut3}} -U {{args.cut5}} -U -{{args.cut3}} -m {{args.minlen}} -q {{args.minq}},{{args.minq}} -o "{{outfq1}}" -p "{{outfq2}}" "{{fq1}}" "{{fq2}}"'
		runcmd (cmd)
	elif tool == 'skewer':
		cmd = '{{args.skewer}} -m pe -t {{args.nthread}} -x {{args.adapter1}} -y {{args.adapter2}} -Q {{args.minq}} -l {{args.minlen}} {{args.gz | lambda x: '-z' if x else ''}} -o "{{job.outdir}}/tmp" "{{fq1}}" "{{fq2}}"'
		runcmd (cmd)
		outfq1 = "{{job.outdir}}/tmp-trimmed-pair1.fastq"
		outfq2 = "{{job.outdir}}/tmp-trimmed-pair2.fastq"
		if {{args.gz}}:
			outfq1 += ".gz"
			outfq2 += ".gz"
		move (outfq1, "{{outfq1}}")
		move (outfq2, "{{outfq2}}")
			
	else:
		raise Exception ('Tool %s not supported' % tool)
except Exception as ex:
	stderr.write ("Job failed: %s" % str(ex))
	raise
"""

"""
@name:
	pFastqSETrim
@description:
	Trim single-end FASTQ reads
@input:
	`fq:file`:  The input fastq file
@output:
	`outfq:file`: The trimmed fastq file
@args:
	`tool`        : The tools used for trimming. Default: trimmomatic (cutadapt|skewer)
	`cutadapt`    : The path of seqtk. Default: cutadapt
	`skewer`      : The path of fastx toolkit trimmer. Default: skewer
	`trimmomatic` : The path of trimmomatic. Default: trimmomatic
	`params`      : Other params for `tool`. Default: ""
	`nthread`     : Number of threads to be used. Default: 1
	- Not for cutadapt
	`gz`          : Whether gzip output files. Default: True
	`mem`         : The memory to be used. Default: 4G
	- Only for trimmomatic
	`minlen`      : Discard trimmed reads that are shorter than `minlen`. Default: 18
	- For trimmomatic, the number will be `minlen`*2 for MINLEN, as it filters before trimming
	`minq`        : Minimal mean qulity for 4-base window or leading/tailing reads. Default: 3
	`cut5`        : Remove the 5'end reads if they are below qulity. Default: 3
	`cut3`        : Remove the 3'end reads if they are below qulity. Default: 3
	- Not for skewer
	`adapter`     : The adapter for sequence. Default: AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
@requires:
	[`cutadapt`](http://cutadapt.readthedocs.io/en/stable/guide.html)
	[`skewer`](https://github.com/relipmoc/skewer)
	[`trimmomatic`](https://github.com/timflutre/trimmomatic)
"""
pFastqSETrim                  = proc (desc = 'Trim single-end reads in fastq file.')
pFastqSETrim.input            = "fq:file"
pFastqSETrim.output           = "outfq:file:{{fq | fn | lambda x: x.rpartition('.')[0] if x.endswith('fastq') or x.endswith('fq') else x}}.fq{{args.gz | lambda x: '.gz' if x else ''}}"
pFastqSETrim.lang             = "python"
pFastqSETrim.args.tool        = 'skewer'
pFastqSETrim.args.cutadapt    = 'cutadapt'
pFastqSETrim.args.skewer      = 'skewer'
pFastqSETrim.args.trimmomatic = 'trimmomatic'
pFastqSETrim.args.params      = ''
pFastqSETrim.args.nthread     = 1
pFastqSETrim.args.gz          = True
pFastqSETrim.args.mem         = '4G'
pFastqSETrim.args.minlen      = 18
pFastqSETrim.args.minq        = 3
pFastqSETrim.args.cut5        = 3
pFastqSETrim.args.cut3        = 3
pFastqSETrim.args.adapter     = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'
pFastqSETrim.args._memtoJava  = mem.toJava.python
pFastqSETrim.args._runcmd     = runcmd.python
pFastqSETrim.script           = """
from sys import stderr, exit
from shutil import move

{{args._runcmd}}
{{args._memtoJava}}

tool = {{args.tool | quote}}
try:
	if tool == 'trimmomatic':		
		mem    = memtoJava ({{args.mem | quote}})
		minlen = str({{args.minlen}} * 2)
		adfile = "{{job.outdir}}/adapters.fa"
		with open (adfile, "w") as ad:
			ad.write (">TruSeq3_IndexedAdapter\\n")
			ad.write ({{args.adapter | quote}} + "\\n")
		cmd    = '{{args.trimmomatic}} %s SE -threads {{args.nthread}} "{{fq}}" "{{outfq}}" ILLUMINACLIP:%s:2:30:10 LEADING:{{args.cut5}} TRAILING:{{args.cut3}} SLIDINGWINDOW:4:{{args.minq}} MINLEN:%s {{args.params}} ' % (mem, adfile, minlen)
		runcmd (cmd)
	elif tool == 'cutadapt':
		cmd = '{{args.cutadapt}} -a {{args.adapter}} -u {{args.cut5}} -u -{{args.cut3}} -m {{args.minlen}} -q {{args.minq}},{{args.minq}} -o "{{outfq}}" "{{fq}}"'
		runcmd (cmd)
	elif tool == 'skewer':
		cmd = '{{args.skewer}} -m any -t {{args.nthread}} -x {{args.adapter}} -Q {{args.minq}} -l {{args.minlen}} {{args.gz | lambda x: '-z' if x else ''}} -o "{{job.outdir}}/tmp" "{{fq}}"'
		runcmd (cmd)
		outfq = "{{job.outdir}}/tmp-trimmed.fastq"
		if {{args.gz}}:
			outfq += ".gz"
		move (outfq, "{{outfq}}")
			
	else:
		raise Exception ('Tool %s not supported' % tool)
except Exception as ex:
	stderr.write ("Job failed: %s" % str(ex))
	raise
"""


"""
@name:
	pFastqSE2Sam
@description:
	Cleaned paired fastq (.fq, .fq.gz, .fastq, .fastq.gz file to mapped sam/bam file
@args:
	`tool`:   The tool used for alignment. Default: bwa (bowtie2|ngm)
	`bwa`:    Path of bwa, default: bwa
	`ngm`:    Path of ngm, default: ngm
	`bowtie2`:Path of bowtie2, default: bowtie2
	`rg`:     The read group. Default: {'id': '', 'pl': 'Illumina', 'pu': 'unit1', 'lb': 'lib1', 'sm': ''}
	- `id` will be parsed from filename with "_LX_" in it if not given
	- `sm` will be parsed from filename
	`ref`:    Path of reference file
	`params`: Other params for tool, default: ''
"""
pFastqSE2Sam                     = proc (desc = 'Map cleaned single-end fastq file to reference genome.')
pFastqSE2Sam.input               = "fq:file"
pFastqSE2Sam.output              = "outfile:file:{{fq | bn | lambda x: __import__('re').sub(r'(\\.clean)?(\\.fq|\\.fastq)(\\.gz)?$', '', x)}}.{{args.outformat}}"
pFastqSE2Sam.args.outformat      = "sam"
pFastqSE2Sam.args.tool           = 'bwa'
pFastqSE2Sam.args.bwa            = 'bwa'
pFastqSE2Sam.args.ngm            = 'ngm'
pFastqSE2Sam.args.bowtie2        = 'bowtie2'
pFastqSE2Sam.args.bowtie2_build  = 'bowtie2-build'
pFastqSE2Sam.args.rg             = {'id': '', 'pl': 'Illumina', 'pu': 'unit1', 'lb': 'lib1', 'sm': ''}
pFastqSE2Sam.args.ref            = ''
pFastqSE2Sam.args.nthread        = 1
pFastqSE2Sam.args.params         = ''
pFastqSE2Sam.args._runcmd        = runcmd.python
pFastqSE2Sam.args._buildArgIndex = buildArgIndex.python
pFastqSE2Sam.beforeCmd           = checkArgsRef.bash
pFastqSE2Sam.lang                = 'python'
pFastqSE2Sam.script              = """
import re
from os import path, symlink, remove
from sys import stdout, stderr, exit
from time import sleep

{{ args._runcmd }}
{{ args._buildArgIndex }}
		
# determine whether the reference file is specified
ref = "{{args.ref}}"

# detemine default read group
rg = {{ args.rg | json }}
rg = {key.upper():val for key, val in rg.items()}
if not rg['ID']:
	g = re.search (r'[^a-zA-Z0-9]+(L\\d+)[^a-zA-Z0-9]+', "{{outfile | fn}}")
	rg['ID'] = g.group(1) if g else "{{outfile | fn}}.L{{#}}"
if not rg['SM']:
	rg['SM'] = "{{outfile | fn}}"
	
if {{args.tool | quote}} == 'bowtie2':
	# check reference index files
	ref = buildArgIndex (
		{{#}}, 
		ref, 
		["{{ args.ref }}.1.bt2", "{{ args.ref }}.2.bt2", "{{ args.ref }}.3.bt2", "{{ args.ref }}.4.bt2", "{{ args.ref }}.rev.1.bt2", "{{ args.ref }}.rev.2.bt2"], 
		'{{args.bowtie2_build}} --thread {{args.nthread}} "{{args.ref}}" "{{args.ref}}"',
		"{{ proc.workdir }}",
		'{{args.bowtie2_build}} --thread {{args.nthread}} "{{ job.outdir }}/{{ args.ref | bn }}" "{{ job.outdir }}/{{ args.ref | bn }}"')
	
	# do mapping
	cmd = '{{args.bowtie2}} --threads {{args.nthread}} --rg-id %s %s -x "%s" -U "{{fq}}" -S "{{outfile}}"' % (rg['ID'], " ".join(["--rg " + k + ":" + v for k,v in rg.items() if k!='ID']), ref)
	runcmd (cmd)

elif {{args.tool | quote}} == 'bwa':
	# check reference index
	ref = buildArgIndex (
		{{#}}, 
		ref, 
		["{{ args.ref }}.bwt", "{{ args.ref }}.amb", "{{ args.ref }}.ann", "{{ args.ref }}.pac", "{{ args.ref }}.sa"], 
		'{{args.bwa}} index "{{args.ref}}"',
		"{{ proc.workdir }}",
		'{{args.bwa}} index "{{ job.outdir }}/{{ args.ref | bn }}"')
			
	# do mapping
	cmd = '{{args.bwa}} mem -t {{args.nthread}} -R "@RG\\\\tID:%s\\\\t%s" "%s" "{{fq}}" > "{{outfile}}"' % (rg['ID'], "\\\\t".join([k+':'+v for k,v in rg.items() if k!='ID']), ref)
	runcmd (cmd)

elif {{args.tool | quote}} == 'ngm':
	# check reference index
	ref = buildArgIndex (
		{{#}}, 
		ref, 
		["{{ args.ref }}-enc.2.ngm", "{{ args.ref }}-ht-13-2.3.ngm"], 
		'{{args.ngm}} -r "{{args.ref}}" -t {{args.nthread}}',
		"{{ proc.workdir }}",
		'{{args.ngm}} -r "{{ job.outdir }}/{{ args.ref | bn }}" -t {{args.nthread}}')
	
	# do mapping
	b = ""
	if {{args.outformat | quote}} == 'bam':
		b = "-b"
	cmd = '{{args.ngm}} -q "{{fq}}" -r "%s" -o "{{outfile}}" %s --rg-id %s --rg-sm %s --rg-lb %s --rg-pl %s --rg-pu %s -t {{args.nthread}}' % (ref, b, rg['ID'], rg['SM'], rg['LB'], rg['PL'], rg['PU'])
	runcmd (cmd)
		
"""

"""
@name:
	pFastqPE2Sam
@description:
	Cleaned paired fastq (.fq, .fq.gz, .fastq, .fastq.gz file to mapped sam/bam file
@args:
	`tool`   : The tool used for alignment. Default: bwa (bowtie2, ngm, star)
	`bwa`    : Path of bwa, default: bwa
	`ngm`    : Path of ngm, default: ngm
	`star`   : Path of ngm, default: STAR
	`bowtie2`: Path of bowtie2, default: bowtie2
	`rg`:     The read group. Default: {'id': '', 'pl': 'Illumina', 'pu': 'unit1', 'lb': 'lib1', 'sm': ''}
	- `id` will be parsed from filename with "_LX_" in it if not given
	- `sm` will be parsed from filename
	`ref`    : Path of reference file
	`refgene`: The GTF file for STAR to build index. It's not neccessary if index is already been built. Default: ''
	`params` : Other params for tool, default: ''
"""
pFastqPE2Sam                     = proc (desc = 'Map cleaned paired fastq file to reference genome.')
pFastqPE2Sam.input               = "fq1:file, fq2:file"
pFastqPE2Sam.output              = "outfile:file:{{fq1 | bn | lambda x: __import__('re').sub(r'[^a-zA-Z0-9]*1(\\.clean)?(\\.fq|\\.fastq)(\\.gz)?$', '', x)}}.{{args.outformat}}"
pFastqPE2Sam.args.outformat      = "sam"
pFastqPE2Sam.args.tool           = 'bwa'
pFastqPE2Sam.args.bwa            = 'bwa'
pFastqPE2Sam.args.ngm            = 'ngm'
pFastqPE2Sam.args.star           = 'STAR'
pFastqPE2Sam.args.bowtie2        = 'bowtie2'
pFastqPE2Sam.args.bowtie2_build  = 'bowtie2-build'
pFastqPE2Sam.args.rg             = {'id': '', 'pl': 'Illumina', 'pu': 'unit1', 'lb': 'lib1', 'sm': ''}
pFastqPE2Sam.args.ref            = ''
pFastqPE2Sam.args.refgene        = ''
pFastqPE2Sam.args.nthread        = 1
pFastqPE2Sam.args.params         = ''
pFastqPE2Sam.args._runcmd        = runcmd.python
pFastqPE2Sam.args._buildArgIndex = buildArgIndex.python
pFastqPE2Sam.beforeCmd           = checkArgsRef.bash
pFastqPE2Sam.lang                = 'python'
pFastqPE2Sam.script              = """
import re
from os import path, symlink
from sys import stdout, stderr, exit
from time import sleep

{{ args._runcmd }}
{{ args._buildArgIndex }}
		
# determine whether the reference file is specified
ref = "{{args.ref}}"
if not ref or not path.exists (ref):
	stderr.write ('No reference specified or found.')
	exit (1)

# detemine default read group
rg = {{ args.rg | json }}
rg = {key.upper():val for key, val in rg.items()}
if not rg['ID']:
	g = re.search (r'[^a-zA-Z0-9]+(L\\d+)[^a-zA-Z0-9]+', "{{outfile | fn}}")
	rg['ID'] = g.group(1) if g else "{{outfile | fn}}.L{{#}}"
if not rg['SM']:
	rg['SM'] = "{{outfile | fn}}"
	
if {{args.tool | quote}} == 'bowtie2':
	# check reference index files
	ref = buildArgIndex (
		{{#}}, 
		ref, 
		["{{ args.ref }}.1.bt2", "{{ args.ref }}.2.bt2", "{{ args.ref }}.3.bt2", "{{ args.ref }}.4.bt2", "{{ args.ref }}.rev.1.bt2", "{{ args.ref }}.rev.2.bt2"], 
		'{{args.bowtie2_build}} --thread {{args.nthread}} "{{args.ref}}" "{{args.ref}}"',
		"{{ proc.workdir }}",
		'{{args.bowtie2_build}} --thread {{args.nthread}} "{{ job.outdir }}/{{ args.ref | bn }}" "{{ job.outdir }}/{{ args.ref | bn }}"')
	
	# do mapping
	cmd = '{{args.bowtie2}} {{args.params}} --threads {{args.nthread}} --rg-id %s %s -x "%s" -1 "{{fq1}}" -2 "{{fq2}}" -S "{{outfile}}"' % (rg['ID'], " ".join(["--rg " + k + ":" + v for k,v in rg.items() if k!='ID']), ref)
	runcmd (cmd)

elif {{args.tool | quote}} == 'bwa':
	# check reference index
	ref = buildArgIndex (
		{{#}}, 
		ref, 
		["{{ args.ref }}.bwt", "{{ args.ref }}.amb", "{{ args.ref }}.ann", "{{ args.ref }}.pac", "{{ args.ref }}.sa"], 
		'{{args.bwa}} index "{{args.ref}}"',
		"{{ proc.workdir }}",
		'{{args.bwa}} index "{{ job.outdir }}/{{ args.ref | bn }}"')
			
	# do mapping
	cmd = '{{args.bwa}} mem {{args.params}} -t {{args.nthread}} -R "@RG\\\\tID:%s\\\\t%s" "%s" "{{fq1}}" "{{fq2}}" > "{{outfile}}"' % (rg['ID'], "\\\\t".join([k+':'+v for k,v in rg.items() if k!='ID']), ref)
	runcmd (cmd)

elif {{args.tool | quote}} == 'ngm':
	# check reference index
	ref = buildArgIndex (
		{{#}}, 
		ref, 
		["{{ args.ref }}-enc.2.ngm", "{{ args.ref }}-ht-13-2.3.ngm"], 
		'{{args.ngm}} -r "{{args.ref}}" -t {{args.nthread}}',
		"{{ proc.workdir }}",
		'{{args.ngm}} -r "{{ job.outdir }}/{{ args.ref | bn }}" -t {{args.nthread}}')
		
	# do mapping
	b = ""
	if {{args.outformat | quote}} == 'bam':
		b = "-b"
	cmd = '{{args.ngm}} {{args.params}} -1 "{{fq1}}" -2 "{{fq2}}" -r "%s" -o "{{outfile}}" %s --rg-id %s --rg-sm %s --rg-lb %s --rg-pl %s --rg-pu %s -t {{args.nthread}}' % (ref, b, rg['ID'], rg['SM'], rg['LB'], rg['PL'], rg['PU'])
	runcmd (cmd)

elif {{args.tool | quote}} == 'star':
	# check reference index
	ref = buildArgIndex (
		{{#}}, 
		ref, 
		map(lambda x: path.join("{{args.ref | prefix}}.star", x), ["chrLength.txt", "chrNameLength.txt", "chrName.txt", "chrStart.txt", "exonGeTrInfo.tab", "exonInfo.tab", "geneInfo.tab", "Genome", "genomeParameters.txt", "SA", "SAindex", "sjdbInfo.txt", "sjdbList.fromGTF.out.tab", "sjdbList.out.tab", "transcriptInfo.tab"]), 
		'mkdir -p "{{args.ref | prefix}}.star"; {{args.star}} --runMode genomeGenerate --genomeDir "{{args.ref | prefix}}.star" --genomeFastaFiles "{{args.ref}}" --sjdbOverhang 100 --sjdbGTFfile "{{args.refgene}}" --runThreadN {{args.nthread}}',
		"{{ proc.workdir }}",
		'mkdir -p "{{job.outdir}}/{{args.ref | fn}}.star"; {{args.star}} --runMode genomeGenerate --genomeDir "{{job.outdir}}/{{args.ref | fn}}.star" --genomeFastaFiles "{{job.outdir}}/{{args.ref | bn}}" --sjdbOverhang 100 --sjdbGTFfile "{{args.refgene}}" --runThreadN {{args.nthread}}')
	refDir = "{{args.ref | prefix}}.star" if ref == "{{args.ref}}" else "{{job.outdir}}/{{args.ref | fn}}.star"
	rfcmd  = "zcat" if "{{fq1}}".endswith('.gz') else 'bzcat' if "{{fq1}}".endswith('.bz2') else 'cat'
	
	cmd = '{{args.star}} {{args.params}} --genomeDir "%s" --readFilesIn "{{fq1}}" "{{fq2}}" --readFilesCommand %s --readNameSeparator . --outFileNamePrefix "{{job.outdir}}/" --outSAMtype {{args.outformat | .upper()}}' % (refDir, rfcmd)
	runcmd (cmd)
"""