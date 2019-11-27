"""
A set of processes to generate/process fastq/fasta files
"""

from os import path
from pyppl import Proc, Box
from . import params, delefactory, procfactory
from modkit import Modkit
Modkit().delegate(delefactory())

def _getFastqFn(fastq):
	ret = path.basename(fastq)
	if ret.endswith('.gz')   : ret = ret[:-3]
	if ret.endswith('.fastq'): ret = ret[:-6]
	if ret.endswith('.fq')   : ret = ret[:-3]
	if ret.endswith('.clean'): ret = ret[:-6]
	if ret.endswith('_clean'): ret = ret[:-6]
	return ret

def _getCommonName(f1, f2):
	ret = path.basename(path.commonprefix([f1, f2]))
	while not ret[-1].isalnum(): ret = ret[:-1]
	return ret

@procfactory
def _pFastq2Expr():
	"""
	@input:
		fqfile1: The fastq file1.
		fqfile2: The fastq file2.
	@output:
		outfile: The expression file
		outdir : Output direcotry with expression and other output files
	@args:
		params   (Box) : Other parameters for `kallisto quant`.
		idxfile  (path): The kallisto index file.
		kallisto (path): The path to `kallisto`.
		nthread  (int) : # threads to use.
	"""
	return Box(
		desc   = 'Call gene expression from pair-end fastq files using kallisto.',
		lang   = params.python.value,
		input  = 'fqfile1:file, fqfile2:file',
		output = [
			'outfile:file:{{	i.fqfile1, i.fqfile2 | \
							__import__("os").path.commonprefix | \
							.rstrip("._[]") | bn }}.kallisto/{{i.fqfile1, i.fqfile2 | \
							__import__("os").path.commonprefix | \
							.rstrip("._[]") | bn }}.expr.txt',
			'outdir:dir:{{	i.fqfile1, i.fqfile2 | \
							__import__("os").path.commonprefix | \
							.rstrip("._[]") | bn }}.kallisto'
		],
		preCmd = """
			{{"reference.bash" | bashimport}}
			reference kallisto {{args.idxfile | quote}}
		""",
		args = Box(
			kallisto = params.kallisto.value,
			params   = Box(),
			nthread  = 1,
			idxfile  = params.kallistoIdx.value,
		)
	)

@procfactory
def _pFastqSim():
	"""
	@input:
		`seed`: The seed to generate simulation file
			- None: use current timestamp.
	@output:
		fq1: The first pair read file
		fq2: The second pair read file
	@args:
		tool (str) : The tool used for simulation. Could be one of:
			- wgsim/dwgsim/art_illumina
		len1         (int) : The length of first pair read.
		len2         (int) : The length of second pair read.
		num          (int) : The number of read PAIRs.
		gz           (bool): Whether generate gzipped read file.
		wgsim        (path): The path of wgsim.
		dwgsim       (path): The path of wgsim.
		art_illumina (path): The path of art_illumina.
		ref          (path): The reference genome. Required
		params       (Box) : Other params for the tool
	@requires:
		[wgsim](https://github.com/lh3/wgsim)
		[dwgsim](https://github.com/nh13/DWGSIM)
		[art](https://www.niehs.nih.gov/research/resources/software/biostatistics/art/): `conda install -c bioconda art`
	"""
	return Box(
		desc   = 'Simulation of pair-end reads',
		lang   = params.python.value,
		input  = 'seed',
		output = [	'fq1:file:read_seed{{i.seed}}_1.fq{{args.gz | ? | =:".gz" | !:""}}',
					'fq2:file:read_seed{{i.seed}}_2.fq{{args.gz | ? | =:".gz" | !:""}}'],
		args = Box(
			tool         = 'dwgsim',
			wgsim        = params.wgsim.value,
			dwgsim       = params.dwgsim.value,
			art_illumina = params.art_illumina.value,
			len1         = 150,
			len2         = 150,
			num          = 1000000,
			gz           = False,
			params       = Box(),
			ref          = params.ref.value,
		)
	)

@procfactory
def _pFastQC():
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
	pFastQC                 = Proc(desc = 'QC report for fastq file.')
	pFastQC.input           = "fq:file"
	pFastQC.output          = "outdir:dir:{{i.fq | getFastqFn}}.fastqc"
	pFastQC.args.tool       = 'fastqc'
	pFastQC.args.fastqc     = params.fastqc.value
	pFastQC.args.nthread    = 1
	pFastQC.args.params     = Box()
	pFastQC.lang            = params.python.value
	pFastQC.envs.getFastqFn = _getFastqFn
	pFastQC.script          = "file:scripts/fastx/pFastQC.py"
	return pFastQC

@procfactory
def _pFastMC():
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
	pFastMC              = Proc(desc = 'Multi-QC based on pFastQC.')
	pFastMC.input        = "qcdir:file"
	pFastMC.output       = "outdir:dir:{{i.qcdir | fn}}.multiqc"
	pFastMC.args.tool    = 'multiqc'
	pFastMC.args.multiqc = params.multiqc.value
	pFastMC.args.params  = Box()
	pFastMC.lang         = params.python.value
	pFastMC.script       = "file:scripts/fastx/pFastMC.py"
	return pFastMC

@procfactory
def _pFastqTrim():
	"""
	@input:
		fq1: The input fastq file 1
		fq2: The input fastq file 2
	@output:
		outfq1: The trimmed fastq file 1
		outfq2: The trimmed fastq file 2
	@args:
		tool       : The tools used for trimming. Available: trimmomatic, cutadapt or skewer
		cutadapt   : The path of seqtk.
		skewer     : The path of fastx toolkit trimmer.
		trimmomatic: The path of trimmomatic.
		params     : Other params for `tool`.
		nthread    : Number of threads to be used.
			- Not for cutadapt
		gz : Whether gzip output files.
		mem: The memory to be used.
			- Only for trimmomatic
		minlen: Discard trimmed reads that are shorter than `minlen`.
			- For trimmomatic, the number will be `minlen`*2 for MINLEN, as it filters before trimming
		minq: Minimal mean qulity for 4-base window or leading/tailing reads.
		cut5: Remove the 5'end reads if they are below qulity.
		cut3: Remove the 3'end reads if they are below qulity.
			- Not for skewer
		adapter1: The adapter for sequence.
		adapter2: The adapter for pair-end sequence.
	@requires:
		[`cutadapt`](http://cutadapt.readthedocs.io/en/stable/guide.html)
		[`skewer`](https://github.com/relipmoc/skewer)
		[`trimmomatic`](https://github.com/timflutre/trimmomatic)
	"""
	return Box(
		desc   = 'Trim pair-end reads in fastq file.',
		lang   = params.python.value,
		envs   = Box(getFastqFn = _getFastqFn),
		input  = "fq1:file, fq2:file",
		output = [
			"outfq1:file:{{i.fq1 | getFastqFn}}.fastq{% if args.gz %}.gz{% endif %}",
			"outfq2:file:{{i.fq2 | getFastqFn}}.fastq{% if args.gz %}.gz{% endif %}"
		],
		args = Box(
			tool        = 'skewer',
			cutadapt    = params.cutadapt.value,
			skewer      = params.skewer.value,
			trimmomatic = params.trimmomatic.value,
			params      = Box(),
			nthread     = 1,
			gz          = False,
			mem         = params.mem4G.value,
			minlen      = 18,
			minq        = 3,
			cut5        = 3,
			cut3        = 3,
			adapter1    = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC',
			adapter2    = 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA'))

@procfactory
def _pFastqSETrim():
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
	pFastqSETrim                  = Proc(desc = 'Trim single-end reads in fastq file.')
	pFastqSETrim.input            = "fq:file"
	pFastqSETrim.output           = "outfq:file:{{i.fq | getFastqFn }}.fastq{% if args.gz %}.gz{% endif %}"
	pFastqSETrim.lang             = params.python.value
	pFastqSETrim.args.tool        = 'skewer'
	pFastqSETrim.args.cutadapt    = params.cutadapt.value
	pFastqSETrim.args.skewer      = params.skewer.value
	pFastqSETrim.args.trimmomatic = params.trimmomatic.value
	pFastqSETrim.args.params      = Box()
	pFastqSETrim.args.nthread     = 1
	pFastqSETrim.args.gz          = False
	pFastqSETrim.args.mem         = params.mem4G.value
	pFastqSETrim.args.minlen      = 18
	pFastqSETrim.args.minq        = 3
	pFastqSETrim.args.cut5        = 3
	pFastqSETrim.args.cut3        = 3
	pFastqSETrim.args.adapter     = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'
	pFastqSETrim.envs.getFastqFn  = _getFastqFn
	pFastqSETrim.script           = "file:scripts/fastx/pFastqSETrim.py"
	return pFastqSETrim


@procfactory
def _pFastqSE2Sam():
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
	pFastqSE2Sam                    = Proc(desc = 'Map cleaned single-end fastq file to reference genome.')
	pFastqSE2Sam.input              = "fq:file"
	pFastqSE2Sam.output             = "outfile:file:{{i.fq | getFastqFn }}.{{args.outfmt}}"
	pFastqSE2Sam.args.outfmt        = "sam"
	pFastqSE2Sam.args.tool          = 'bwa'
	pFastqSE2Sam.args.bwa           = params.bwa.value
	pFastqSE2Sam.args.ngm           = params.ngm.value
	pFastqSE2Sam.args.star          = params.star.value
	pFastqSE2Sam.args.samtools      = params.samtools.value
	pFastqSE2Sam.args.bowtie2       = params.bowtie2.value
	pFastqSE2Sam.args.bowtie2_build = params.bowtie2.value + '-build'
	pFastqSE2Sam.args.rg            = Box(id = '', pl = 'Illumina', pu = 'unit1', lb = 'lib1', sm = '')
	pFastqSE2Sam.args.ref           = params.ref.value
	pFastqSE2Sam.args.refgene       = params.refgene.value
	pFastqSE2Sam.args.nthread       = 1
	pFastqSE2Sam.args.params        = Box()
	pFastqSE2Sam.envs.getFastqFn    = _getFastqFn
	pFastqSE2Sam.lang               = params.python.value
	pFastqSE2Sam.script             = "file:scripts/fastx/pFastqSE2Sam.py"
	return pFastqSE2Sam

@procfactory
def _pFastq2Sam():
	"""
	@description:
		Cleaned paired fastq (.fq, .fq.gz, .fastq, .fastq.gz file to mapped sam/bam file
	@args:
		`tool`   : The tool used for alignment. Available: bowtie2, ngm or star.
		`bwa`    : Path of bwa
		`ngm`    : Path of ngm
		`star`   : Path of STAR
			- The genomeDir should be at `/path/to/hg19.star` if `args.ref` is `/path/to/hg19.fa`
		`bowtie2`: Path of bowtie2
		`rg`:     The read group.
			- `id` will be parsed from filename with "_LX_" in it if not given
			- `sm` will be parsed from filename
		`ref`    : Path of reference file
		`refexon`: The GTF file for STAR to build index. It's not neccessary if index is already been built.
		`params` : Other params for tool
	"""
	return Box(
		desc   = 'Map cleaned paired fastq file to reference genome.',
		lang   = params.python.value,
		input  = "fq1:file, fq2:file",
		output = "outfile:file:{{i.fq1, i.fq2 | path.commonprefix | bn | .rstrip: '_. ,[]' }}.{{args.outfmt}}",
		envs   = Box(path = path),
		preCmd = """
			{{"reference.bash" | bashimport}}
			export bwa={{args.bwa | squote}}
			export ngm={{args.ngm | squote}}
			export star={{args.star | squote}}
			export samtools={{args.samtools | squote}}
			export bowtie2={{args.bowtie2 | squote}}
			export nthread={{args.nthread}}
			export refexon={{args.refexon | squote}}
			reference {{args.tool | squote}} {{args.ref | squote}}""",
		args = Box(
			tool     = 'bwa',
			outfmt   = 'sam',
			bwa      = params.bwa.value,
			ngm      = params.ngm.value,
			star     = params.star.value,
			samtools = params.samtools.value,
			bowtie2  = params.bowtie2.value,
			rg       = Box(id = '', pl = 'Illumina', pu = 'unit1', lb = 'lib1', sm = ''),
			ref      = params.ref.value,
			refexon  = params.refexon.value,
			nthread  = 1,
			params   = Box()
		)
	)

