"""
A set of processes to generate/process fastq/fasta files
"""

from os import path
from pyppl import Proc, Box
from .utils import mem2, runcmd, buildref, checkref, helpers, txt
from . import params

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

pFastq2Expr        = Proc(desc = 'Use Kallisto to get gene expression from pair-end fastq files.')
pFastq2Expr.input  = "fqfile1:file, fqfile2:file"
pFastq2Expr.output = [
	"outfile:file:{{in.fqfile1, in.fqfile2 | getCommonName}}/{{in.fqfile1, in.fqfile2 | getCommonName}}.expr", 
	"outdir:dir:{{in.fqfile1, in.fqfile2 | getCommonName}}"
]
pFastq2Expr.args.params            = Box()
pFastq2Expr.args.ref               = '' # don't give the whole genome sequence! takes long time!
pFastq2Expr.args.idxfile           = params.kallistoIdx.value
pFastq2Expr.args.kallisto          = params.kallisto.value
pFastq2Expr.args.nthread           = 1
pFastq2Expr.tplenvs.getCommonName  = _getCommonName
pFastq2Expr.tplenvs.runcmd         = runcmd.py
pFastq2Expr.tplenvs.params2CmdArgs = helpers.params2CmdArgs.py
pFastq2Expr.tplenvs.txtFilter = txt.filter.py
pFastq2Expr.beforeCmd              = """
if [[ -z "{{args.idxfile}}" ]]; then
	echo 'Index file (args.idxfile) is required.' 1>&2
	exit 1
fi
if [[ ! -e "{{args.idxfile}}" ]]; then
	if [[ -z "{{args.ref}}" ]]; then
		echo 'Reference file (args.ref) is required to generate index file for kallisto.' 1>&2
		exit 1
	fi
	cmd='{{args.kallisto}} index -i "{{args.idxfile}}" "{{args.ref}}"' 
	echo 'Generating kallisto index ...'
	eval $cmd
fi
"""
pFastq2Expr.lang                   = params.python.value
pFastq2Expr.script                 = "file:scripts/fastx/pFastq2Expr.py"

"""
@name:
	pFastqSim
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
pFastqSim                        = Proc(desc = 'Simulate pair-end reads.')
pFastqSim.input                  = "in"
pFastqSim.output                 = "fq1:file:read{{in.in}}_1.fastq{% if args.gz %}.gz{% endif %}, fq2:file:read{{in.in}}_2.fastq{% if args.gz %}.gz{% endif %}"
pFastqSim.args.tool              = 'wgsim'
pFastqSim.args.wgsim             = params.wgsim.value
pFastqSim.args.dwgsim            = params.dwgsim.value
pFastqSim.args.len1              = 100
pFastqSim.args.len2              = 100
pFastqSim.args.num               = 1000000
pFastqSim.args.gz                = False
pFastqSim.args.seed              = None
pFastqSim.args.params            = Box()
pFastqSim.args.ref               = params.ref.value
pFastqSim.beforeCmd              = checkref.fa.bash
pFastqSim.lang                   = params.python.value
pFastqSim.tplenvs.runcmd         = runcmd.py
pFastqSim.tplenvs.params2CmdArgs = helpers.params2CmdArgs.py
pFastqSim.script                 = "file:scripts/fastx/pFastqSim.py"

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
pFastQC                        = Proc(desc = 'QC report for fastq file.')
pFastQC.input                  = "fq:file"
pFastQC.output                 = "outdir:dir:{{in.fq | getFastqFn}}"
pFastQC.args.tool              = 'fastqc'
pFastQC.args.fastqc            = params.fastqc.value
pFastQC.args.nthread           = 1
pFastQC.args.params            = Box()
pFastQC.lang                   = params.python.value
pFastQC.tplenvs.params2CmdArgs = helpers.params2CmdArgs.py
pFastQC.tplenvs.getFastqFn     = _getFastqFn
pFastQC.tplenvs.runcmd         = runcmd.py
pFastQC.script                 = "file:scripts/fastx/pFastQC.py"

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
pFastMC                        = Proc(desc = 'Multi-QC based on pFastQC.')
pFastMC.input                  = "qcdir:file"
pFastMC.output                 = "outdir:dir:{{in.qcdir | fn}}_multiqc_{{job.index}}"
pFastMC.args.tool              = 'multiqc'
pFastMC.args.multiqc           = params.multiqc.value
pFastMC.args.params            = Box()
pFastMC.lang                   = params.python.value
pFastMC.tplenvs.params2CmdArgs = helpers.params2CmdArgs.py
pFastMC.tplenvs.runcmd         = runcmd.py
pFastMC.script                 = "file:scripts/fastx/pFastMC.py"

"""
@name:
	pFastqTrim
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
pFastqTrim                  = Proc(desc = 'Trim pair-end reads in fastq file.')
pFastqTrim.input            = "fq1:file, fq2:file"
pFastqTrim.output           = [
	"outfq1:file:{{in.fq1 | getFastqFn}}.fastq{% if args.gz %}.gz{% endif %}",
	"outfq2:file:{{in.fq2 | getFastqFn}}.fastq{% if args.gz %}.gz{% endif %}"
]
pFastqTrim.lang                   = params.python.value
pFastqTrim.args.tool              = 'skewer'
pFastqTrim.args.cutadapt          = params.cutadapt.value
pFastqTrim.args.skewer            = params.skewer.value
pFastqTrim.args.trimmomatic       = params.trimmomatic.value
pFastqTrim.args.params            = Box()
pFastqTrim.args.nthread           = 1
pFastqTrim.args.gz                = False
pFastqTrim.args.mem               = params.mem4G.value
pFastqTrim.args.minlen            = 18
pFastqTrim.args.minq              = 3
pFastqTrim.args.cut5              = 3
pFastqTrim.args.cut3              = 3
pFastqTrim.args.adapter1          = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'
pFastqTrim.args.adapter2          = 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA'
pFastqTrim.tplenvs.mem2           = mem2.py
pFastqTrim.tplenvs.runcmd         = runcmd.py
pFastqTrim.tplenvs.getFastqFn     = _getFastqFn
pFastqTrim.tplenvs.params2CmdArgs = helpers.params2CmdArgs.py
pFastqTrim.script                 = "file:scripts/fastx/pFastqTrim.py"

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
pFastqSETrim                        = Proc(desc = 'Trim single-end reads in fastq file.')
pFastqSETrim.input                  = "fq:file"
pFastqSETrim.output                 = "outfq:file:{{in.fq | getFastqFn }}.fastq{% if args.gz %}.gz{% endif %}"
pFastqSETrim.lang                   = params.python.value
pFastqSETrim.args.tool              = 'skewer'
pFastqSETrim.args.cutadapt          = params.cutadapt.value
pFastqSETrim.args.skewer            = params.skewer.value
pFastqSETrim.args.trimmomatic       = params.trimmomatic.value
pFastqSETrim.args.params            = Box()
pFastqSETrim.args.nthread           = 1
pFastqSETrim.args.gz                = False
pFastqSETrim.args.mem               = params.mem4G.value
pFastqSETrim.args.minlen            = 18
pFastqSETrim.args.minq              = 3
pFastqSETrim.args.cut5              = 3
pFastqSETrim.args.cut3              = 3
pFastqSETrim.args.adapter           = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'
pFastqSETrim.tplenvs.mem2           = mem2.py
pFastqSETrim.tplenvs.runcmd         = runcmd.py
pFastqSETrim.tplenvs.getFastqFn     = _getFastqFn
pFastqSETrim.tplenvs.params2CmdArgs = helpers.params2CmdArgs.py
pFastqSETrim.script                 = "file:scripts/fastx/pFastqSETrim.py"


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
pFastqSE2Sam                        = Proc(desc = 'Map cleaned single-end fastq file to reference genome.')
pFastqSE2Sam.input                  = "fq:file"
pFastqSE2Sam.output                 = "outfile:file:{{in.fq | getFastqFn }}.{{args.outformat}}"
pFastqSE2Sam.args.outformat         = "sam"
pFastqSE2Sam.args.tool              = 'bwa'
pFastqSE2Sam.args.bwa               = params.bwa.value
pFastqSE2Sam.args.ngm               = params.ngm.value
pFastqSE2Sam.args.bowtie2           = params.bowtie2.value
pFastqSE2Sam.args.bowtie2_build     = params.bowtie2.value + '-build'
pFastqSE2Sam.args.rg                = Box({'id': '', 'pl': 'Illumina', 'pu': 'unit1', 'lb': 'lib1', 'sm': ''})
pFastqSE2Sam.args.ref               = params.ref.value
pFastqSE2Sam.args.nthread           = 1
pFastqSE2Sam.args.params            = Box()
pFastqSE2Sam.tplenvs.runcmd         = runcmd.py
pFastqSE2Sam.tplenvs.buildrefIndex  = buildref.index.py
pFastqSE2Sam.tplenvs.getFastqFn     = _getFastqFn
pFastqSE2Sam.tplenvs.params2CmdArgs = helpers.params2CmdArgs.py
pFastqSE2Sam.beforeCmd              = checkref.fa.bash
pFastqSE2Sam.lang                   = params.python.value
pFastqSE2Sam.script                 = "file:scripts/fastx/pFastqSE2Sam.py"

"""
@name:
	pFastq2Sam
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
pFastq2Sam                        = Proc(desc = 'Map cleaned paired fastq file to reference genome.')
pFastq2Sam.input                  = "fq1:file, fq2:file"
pFastq2Sam.output                 = "outfile:file:{{in.fq1 | getFastqFn | lambda x: x[:-2] if x.endswith('_1') or x.endswith('_2') else x }}.{{args.outformat}}"
pFastq2Sam.args.outformat         = "sam"
pFastq2Sam.args.tool              = 'bwa'
pFastq2Sam.args.bwa               = params.bwa.value
pFastq2Sam.args.ngm               = params.ngm.value
pFastq2Sam.args.star              = params.star.value
pFastq2Sam.args.bowtie2           = params.bowtie2.value
pFastq2Sam.args.bowtie2_build     = params.bowtie2.value + '-build'
pFastq2Sam.args.rg                = Box({'id': '', 'pl': 'Illumina', 'pu': 'unit1', 'lb': 'lib1', 'sm': ''})
pFastq2Sam.args.ref               = params.ref.value
pFastq2Sam.args.refgene           = params.refgene.value
pFastq2Sam.args.nthread           = 1
pFastq2Sam.args.params            = Box()
pFastq2Sam.tplenvs.runcmd         = runcmd.py
pFastq2Sam.tplenvs.buildrefIndex  = buildref.index.py
pFastq2Sam.tplenvs.getFastqFn     = _getFastqFn
pFastq2Sam.tplenvs.params2CmdArgs = helpers.params2CmdArgs.py
pFastq2Sam.beforeCmd              = checkref.fa.bash
pFastq2Sam.lang                   = 'python'
pFastq2Sam.script                 = "file:scripts/fastx/pFastq2Sam.py"