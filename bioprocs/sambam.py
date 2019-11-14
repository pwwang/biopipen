"""
A set of processes to generate/process sam/bam files
"""
from pyppl import Proc, Box
#from .utils import mem2, runcmd, buildref, checkref, polling, helpers, plot, dirnameFiles
from . import params, bashimport, rimport
from . import delefactory, procfactory
from modkit import Modkit
Modkit().delegate(delefactory())

@procfactory
def _pSam2Bam():
	"""
	@input:
		infile: The input file
	@output:
		outfile: The output bam file
	@args:
		tool             : The tool used to do the sort. Available: sambamba|picard|sambamba|biobambam|samtools
		sambamba         : The path of the sambamba.
		picard           : The path of the picard.
		biobambam_bamsort: The path of the biobambam's bamsort.
		samtools: The path of the samtools.
		sort    : Do sorting?
			- If input is sam, tool is biobambam, this should be True
		index  : Do indexing?
		markdup: Do duplicates marking?
			- `rmdup` for samtools will be called
		rmdup  : Do duplicates removing?
		tmpdir : The tmp dir used to store tmp files.
		sortby : Sort by coordinate or queryname.
		nthread: Number of threads to use.
		infmt  : The format of input file. Available: sam|bam
		params : Other parameters for `tool`.
		mem    : The max memory to use.
			- Unit could be G/g/M/m
			- Will be converted to -Xmx4G, and -Xms will be 1/8 of it
	@requires:
		[sambamba](https://lomereiter.github.io/sambamba/docs/sambamba-view.html) if `args.tool` == samtools or reference used but not indexed.
		[picard](https://broadinstitute.github.io/picard/command-line-overview.html)
		[biobambam](https://github.com/gt1/biobambam2)
		[samtools](https://github.com/samtools/samtools)
	"""
	return Box(
		desc   = 'Deal with mapped sam/bam files, including sort, markdup, rmdup, and/or index.',
		lang   = params.python.value,
		input  = "infile:file",
		output = "outfile:file:{{i.infile | fn}}.bam, outidx:file:{{i.infile | fn}}.bam.bai",
		errhow = 'retry',
		preCmd = """
			{{bashimport}} reference.bash
			export elprep={{args.elprep | quote}}
			if [[ {{args.tool | quote}} == "elprep" ]]; then
				reference elprep {{args.ref | squote}}
			fi""",
		args = Box(
			tool       = "elprep",
			sambamba   = params.sambamba.value,
			picard     = params.picard.value,
			biobambam  = params.biobambam_bamsort.value,
			samtools   = params.samtools.value,
			elprep     = params.elprep.value,
			steps      = Box(sort=True, index=True, markdup=True, rmdup=True, recal=True),
			tmpdir     = params.tmpdir.value,
			sortby     = "coordinate",
			nthread    = 1,
			params     = Box(),
			mem        = params.mem16G.value,
			ref        = params.ref.value,
			knownSites = ''))

@procfactory
def _pBamMarkdup():
	"""
	@name:
		pBamMarkdup
	@description:
		Mark/remove duplicates for bam files
	@input:
		`infile:file`: The input file
	@output:
		`outfile:file`: The output bam file
	@args:
		`tool`             : The tool used to do the sort. Default: sambamba (picard|sambamba|biobambam|samtools|bamutil)
		`sambamba`         : The path of sambamba. Default: sambamba
		`picard`           : The path of picard. Default: picard
		`biobambam_bamsort`: The path of biobambam's bamsort. Default: bamsort
		`samtools`         : The path of samtools. Default: samtools
		`bamutil`          : The path of bamutil. Default: bam
		`rmdup`            : Do duplicates removing? Default: False
		- Samtools will anyway remove the duplicates
		`tmpdir`           : The tmp dir used to store tmp files. Default: <system default tmpdir>
		`nthread`          : Default: 1
		- Not available for samtools and picard
		`params`           : Other parameters for `tool`. Defaut: ""
		`mem`              : The max memory to use. Default: "16G"
		- Unit could be G/g/M/m
		- Will be converted to -Xmx4G, and -Xms will be 1/8 of it
	@requires:
		[sambamba](https://lomereiter.github.io/sambamba/docs/sambamba-view.html)
		[picard](https://broadinstitute.github.io/picard/command-line-overview.html)
		[biobambam](https://github.com/gt1/biobambam2)
		[samtools](https://github.com/samtools/samtools)
		[bamutil](http://genome.sph.umich.edu/wiki/BamUtil#Programs)
	"""
	pBamMarkdup                        = Proc(desc = 'Mark/remove duplicates for bam files.')
	pBamMarkdup.input                  = "infile:file"
	pBamMarkdup.output                 = "outfile:file:{{i.infile | fn}}.bam"
	pBamMarkdup.args.tool              = "biobambam"
	pBamMarkdup.args.sambamba          = params.sambamba.value
	pBamMarkdup.args.picard            = params.picard.value
	pBamMarkdup.args.biobambam_bamsort = params.biobambam_bamsort.value
	pBamMarkdup.args.samtools          = params.samtools.value
	pBamMarkdup.args.bamutil           = params.bamutil.value
	pBamMarkdup.args.rmdup             = False
	pBamMarkdup.args.tmpdir            = params.tmpdir.value
	pBamMarkdup.args.nthread           = 1
	pBamMarkdup.args.params            = Box()
	pBamMarkdup.args.mem               = params.mem16G.value
	#pBamMarkdup.envs.mem2              = mem2.py
	#pBamMarkdup.envs.runcmd            = runcmd.py
	#pBamMarkdup.envs.params2CmdArgs    = helpers.params2CmdArgs.py
	pBamMarkdup.lang                   = params.python.value
	pBamMarkdup.script                 = "file:scripts/sambam/pBamMarkdup.py"
	return pBamMarkdup

@procfactory
def _pBamRecal():
	"""
	@name:
		pBamRecal
	@description:
		Recalibrate a bam file
	@input:
		`infile:file`: The bam file
	@output:
		`outfile:file`: The output bam file
	@args:
		`tool`    : The tool used to recalibrate the bam file. Default: `gatk` (gatk|bamutil)
		`gatk`    : The path of gatk, including java path. Default: `gatk`
		`samtools`: The path of samtools. Default: `samtools`
		`bamutil` : The path of bamutil. Default: `bam`
		`picard`  : The path of picard. Default: `picard`
		`params`  : Other parameters for `bam recab`. Default         : ""
			`RealignerTargetCreator` : Other parameters for `gatk RealignerTargetCreator`. Defaut: ""
			`IndelRealigner`         : Other parameters for `gatk IndelRealigner`. Defaut: ""
			`BaseRecalibrator`       : Other parameters for `gatk BaseRecalibrator`. Defaut: ""
			`PrintReads`             : Other parameters for `gatk PrintReads`. Defaut: ""
		`mem`: The max memory to use. Default: "32G"
		`knownSites`: The known polymorphic sites to mask out. Default: "" (Required for GATK)
		`ref`: The reference file. Required.
			- Will be converted to -Xmx4G, and -Xms will be 1/8 of it
	@requires:
		[gatk](https://software.broadinstitute.org/gatk)
		[samtools](https://github.com/samtools/samtools) if `args.ref` is not indexed, or bamutil is used for bam index file generation.
		[picard](https://broadinstitute.github.io/picard/command-line-overview.html) if `args.ref is not dicted.`
	"""
	pBamRecal               = Proc(desc = 'Recalibrate a bam file.')
	pBamRecal.input         = "infile:file"
	pBamRecal.output        = "outfile:file:{{i.infile | fn}}.bam, outidx:file:{{i.infile | fn}}.bam.bai"
	pBamRecal.args.tool     = "bamutil"
	pBamRecal.args.gatk     = params.gatk.value
	pBamRecal.args.samtools = params.samtools.value
	pBamRecal.args.picard   = params.picard.value
	pBamRecal.args.bamutil  = params.bamutil.value
	pBamRecal.args.params   = Box(
		# For GATK
		# RealignerTargetCreator = Box(),
		# IndelRealigner         = Box(),
		# BaseRecalibrator       = Box(),
		# PrintReads             = Box()
	)
	pBamRecal.args.ref        = params.ref.value
	pBamRecal.args.tmpdir     = params.tmpdir.value
	pBamRecal.args.knownSites = params.dbsnp.value
	pBamRecal.args.nthread    = 1
	pBamRecal.args.mem        = params.mem32G.value
	pBamRecal.envs.bashimport = bashimport
	pBamRecal.preCmd          = """
	{{bashimport}} reference.bash
	export samtools={{args.samtools | squote}}
	reference fasta {{args.ref | squote}}
	"""
	pBamRecal.lang   = params.python.value
	pBamRecal.script = "file:scripts/sambam/pBamRecal.py"
	return pBamRecal

@procfactory
def _pBamReadGroup():
	"""
	@name:
		pBamReadGroup
	@description:
		Add or replace read groups of a bam file
	@input:
		`infile:file`: The bam file
	@output:
		`outfile:file`: The output bam file
	@args:
		`tool`                         : The tool used. Default: `picard` (picard|bamutil)
		`picard`                       : The path of picard. Default: `picard`
		`bamutil`                      : The path of bamutil. Default: `bam`
		`rg`                           : The read group. Default: {'id': '', 'pl': 'Illumina', 'pu': 'unit1', 'lb': 'lib1', 'sm': ''}
		- `id` will be parsed from filename with "_LX_" in it if not given
		- `sm` will be parsed from filename
		`params`                       : Other parameters for `tool`. Defaut: ""
		`mem`                          : The max memory to use. Default: "4G"
		- Will be converted to -Xmx4G, and -Xms will be 1/8 of it
		`tmpdir`                       : The temporary directory. Default: <system tmpdir>
	@requires:
		[gatk](https://lomereiter.github.io/sambamba/docs/sambamba-view.html)
		[samtools](https://github.com/samtools/samtools) if `args.ref` is not indexed.
		[picard](https://broadinstitute.github.io/picard/command-line-overview.html) if `args.ref is not dicted.`
	"""
	pBamReadGroup                     = Proc(desc = 'Add or replace read groups of a bam file.')
	pBamReadGroup.input               = "infile:file"
	pBamReadGroup.output              = "outfile:file:{{i.infile | bn}}"
	pBamReadGroup.args.tool           = "bamutil"
	pBamReadGroup.args.picard         = params.picard.value
	pBamReadGroup.args.bamutil        = params.bamutil.value
	pBamReadGroup.args.rg             = Box({'id': '', 'pl': 'Illumina', 'pu': 'unit1', 'lb': 'lib1', 'sm': ''})
	pBamReadGroup.args.params         = Box()
	pBamReadGroup.args.tmpdir         = params.tmpdir.value
	pBamReadGroup.args.mem            = params.mem4G.value
	#pBamReadGroup.envs.params2CmdArgs = helpers.params2CmdArgs.py
	pBamReadGroup.lang                = params.python.value
	pBamReadGroup.script              = "file:scripts/sambam/pBamReadGroup.py"
	return pBamReadGroup


@procfactory
def _pBamReorder():
	"""
	@name:
		pBamReorder
	@description:
		Reorder a sam/bam file by a given reference file using `picard ReorderSam`
	@input:
		`infile:file`: The sam/bam file
	@output:
		`outfile:file`: The output bam file
	@args:
		`picard`                       : The path of picard. Default: `picard`
		`ref`                          : The reference file. Required
		`params`                       : Other parameters for `picard ReorderSam`. Defaut: ""
		`mem`                          : The max memory to use. Default: "4G"
		- Will be converted to -Xmx4G, and -Xms will be 1/8 of it
		`tmpdir`                       : The temporary directory. Default: <system tmpdir>
	@requires:
		[picard](https://broadinstitute.github.io/picard/command-line-overview.html)
	"""
	pBamReorder                     = Proc(desc = 'Reorder a sam/bam file by a given reference.')
	pBamReorder.input               = "infile:file"
	pBamReorder.output              = "outfile:file:{{i.infile | bn}}"
	pBamReorder.args.picard         = params.picard.value
	pBamReorder.args.params         = Box()
	pBamReorder.args.tmpdir         = params.tmpdir.value
	pBamReorder.args.mem            = params.mem4G.value
	pBamReorder.args.ref            = params.ref.value
	pBamReorder.lang                = params.python.value
	pBamReorder.script              = "file:scripts/sambam/pBamReorder.py"
	return pBamReorder


@procfactory
def _pBamMerge():
	"""
	@name:
		pBamMerge
	@description:
		Merges multiple SAM and/or BAM files (must be sorted by coordinate) into a single file.
	@input:
		`infiles:files`: Input sam/bam files to be merged
	@output:
		`outfile:file`: The merged bam file
	@args:
		`tool`     : The tool used to merge. Default: bamutil (picard|samtools|sambamba)
		`picard`   : The path of picard. Default: `picard`
		`bamutil`  : The path of bamutil. Default: `bam`
		`samtools` : The path of samtools. Default: `samtools`
		`sambamba` : The path of sambamba. Default: `sambamba`
		`params`   : Other parameters for `tool`. Defaut: ""
		`mem`      : The max memory to use. Default: "4G"
		- Will be converted to -Xmx4G, and -Xms will be 1/8 of it, just for picard
		`tmpdir`   : The temporary directory. Default: <system tmpdir>
		`nthread`  : # threads to use. Default: 1
		- For picard, if nthread>1, USE_THREADING=true, otherwise USE_THREADING=false
	@requires:
		[picard](https://broadinstitute.github.io/picard/command-line-overview.html)
	"""
	pBamMerge                     = Proc(desc = 'Merges multiple SAM and/or BAM sorted files into a single file.')
	pBamMerge.input               = "infiles:files"
	pBamMerge.output              = "outfile:file:{{i.infiles[0] | fn2}}.etc.bam"
	pBamMerge.args.tool           = "picard"
	pBamMerge.args.picard         = params.picard.value
	pBamMerge.args.bamutil        = params.bamutil.value
	pBamMerge.args.samtools       = params.samtools.value
	pBamMerge.args.sambamba       = params.sambamba.value
	pBamMerge.args.params         = Box()
	pBamMerge.args.tmpdir         = params.tmpdir.value
	pBamMerge.args.nthread        = 1
	pBamMerge.args.mem            = params.mem4G.value
	#pBamMerge.envs.fsDirname      = dirnameFiles
	pBamMerge.lang                = params.python.value
	pBamMerge.script              = "file:scripts/sambam/pBamMerge.py"
	return pBamMerge


@procfactory
def _pBam2Gmut():
	"""
	@name:
		pBam2Gmut
	@description:
		Call germline (snps and indels) from a call-ready bam file.
	@input:
		`infile:file`: The input bam file
	@output:
		`outfile:file`: The vcf file containing the mutations
	@args:
		`tool`      : The tool used to call mutations. Default: gatk (vardict, snvsniffer, platypus, strelka)
		`gatk`      : The path of gatk. Default: gatk
		`vardict`   : The path of vardict. Default: vardict
		`snvsniffer`: The path of snvsniffer. Default: SNVSniffer
		`samtools`  : The path of samtools. Default: samtools (used to generate reference index)
		`platypus`  : The path of platypus. Default: platypus
		`strelka`   : The path of strelka. Default: configureStrelkaGermlineWorkflow.py
		`cfgParams` : The params for `strelka` configuration. Default: ""
		`picard`    : The path of picard. Default: picard
		`mem`       : The memory to be used. Default: 32G
			- will be converted to -Xms4G -Xmx32G for java programs
		`ref`:          The reference file. Required.
		`gz`:           Gzip output file? Default: False
		`tmpdir`:       The temporary directory. Default: <system tmpdir>
		`params`:       Other params for `tool`. Default: ""
	@requires:
		[gatk](https://lomereiter.github.io/sambamba/docs/sambamba-view.html)
		[samtools](https://github.com/samtools/samtools) if `args.ref` is not indexed.
		[picard](https://broadinstitute.github.io/picard/command-line-overview.html) if `args.ref is not dicted.`
		[vardict](https://github.com/AstraZeneca-NGS/VarDict)
		[snvsniffer](http://snvsniffer.sourceforge.net/homepage.htm#latest)
		[platypus](http://www.well.ox.ac.uk/platypus)
		[strelka@2.7.1+](https://github.com/Illumina/strelka)
	"""
	pBam2Gmut                 = Proc(desc = 'Call germline (snps and indels) from a call-ready bam file.')
	pBam2Gmut.input           = "infile:file"
	pBam2Gmut.output          = "outfile:file:{{i.infile | fn}}.vcf{{args.gz | lambda x: '.gz' if x else ''}}"
	pBam2Gmut.args.tool       = "strelka"
	pBam2Gmut.args.gatk       = params.gatk.value
	pBam2Gmut.args.vardict    = params.vardict.value
	pBam2Gmut.args.snvsniffer = params.snvsniffer.value
	pBam2Gmut.args.samtools   = params.samtools.value # required by SNVSniffer to generate a bam header file
	pBam2Gmut.args.platypus   = params.platypus.value
	pBam2Gmut.args.picard     = params.picard.value
	pBam2Gmut.args.strelka    = params.strelka_germ.value
	pBam2Gmut.args.mem        = params.mem24G.value
	pBam2Gmut.args.ref        = params.ref.value
	pBam2Gmut.args.tmpdir     = params.tmpdir.value
	pBam2Gmut.args.cfgParams  = Box() # only for strelka
	pBam2Gmut.args.params     = Box()
	pBam2Gmut.args.gz         = False
	pBam2Gmut.args.nthread    = 1 # for gatk and platypus
	pBam2Gmut.envs.bashimport = bashimport
	pBam2Gmut.preCmd          = """
	{{bashimport}} reference.bash
	export samtools={{args.samtools | squote}}
	export picard={{args.picard | squote}}
	reference fasta {{args.ref | squote}}
	reference picard {{args.ref | squote}}
	"""
	pBam2Gmut.lang   = params.python.value
	pBam2Gmut.script = "file:scripts/sambam/pBam2Gmut.py"
	return pBam2Gmut

@procfactory
def _pBamPair2Smut():
	"""
	@name:
		pBamPair2Smut
	@description:
		Call somatic mutations from tumor-normal bam pair.
	@input:
		`tumor:file`: The tumor bam file
		`normal:file`: The normal bam file
	@output:
		`outfile:file`: The vcf file
	@args:
		`tool`: The tool used to call mutations. Default: gatk (somaticsniper, strelka, snvsniffer, virmid, varidct)
		`gatk`: The path to gatk. Default: gatk
		`somaticsniper`: The path to gatk. Default: bam-somaticsniper
		`strelka`: The path to gatk. Default: configureStrelkaSomaticWorkflow.py
		`snvsniffer`: The path to gatk. Default: SNVSniffer
		`virmid`: The path to gatk. Default: virmid
		`vardict`: The path to gatk. Default: vardict
		`samtools`: The path to gatk. Default: samtools
		`picard`: The path to gatk. Default: picard
		`configParams`: The configuration parameters for `configureStrelkaSomaticWorkflow.py`. Default: `{}`
		`params`: The parameters for main programs. Default: `{}`
		`meme`: The memory. Default: 24G
		`ref`: The reference genom. Default: `params.ref.value`
		`gz`: Whether gzip the output vcf file. Default: False
		`nthread`: The number of threads to use. Default: 1
		`tmpdir`: The temporary directory. Default: `params.tmpdir.value`
	@requires:
		[gatk](https://lomereiter.github.io/sambamba/docs/sambamba-view.html)
		[samtools](https://github.com/samtools/samtools) if `args.ref` is not indexed.
		[picard](https://broadinstitute.github.io/picard/command-line-overview.html) if `args.ref is not dicted.`
		[vardict](https://github.com/AstraZeneca-NGS/VarDict)
		[snvsniffer](http://snvsniffer.sourceforge.net/homepage.htm#latest)
		[platypus](http://www.well.ox.ac.uk/platypus)
		[strelka@2.7.1+](https://github.com/Illumina/strelka)
	"""
	pBamPair2Smut                    = Proc(desc = 'Call somatic mutations from tumor-normal bam pair.')
	pBamPair2Smut.input              = "tumor:file, normal:file"
	pBamPair2Smut.output             = "outfile:file:{{i.tumor | fn}}-{{i.normal | fn}}.vcf{{args.gz | lambda x: '.gz' if x else ''}}"
	pBamPair2Smut.args.tool          = 'strelka'
	pBamPair2Smut.args.gatk          = params.gatk.value # required for strelka
	pBamPair2Smut.args.somaticsniper = params.somaticsniper.value
	pBamPair2Smut.args.strelka       = params.strelka_soma.value # @2.7.1
	pBamPair2Smut.args.snvsniffer    = params.snvsniffer.value
	pBamPair2Smut.args.virmid        = params.virmid.value
	pBamPair2Smut.args.vardict       = params.vardict.value
	pBamPair2Smut.args.samtools      = params.samtools.value
	pBamPair2Smut.args.picard        = params.picard.value
	pBamPair2Smut.args.configParams  = Box() # only for strelka
	pBamPair2Smut.args.params        = Box()
	pBamPair2Smut.args.mem           = params.mem24G.value
	pBamPair2Smut.args.ref           = params.ref.value
	pBamPair2Smut.args.gz            = False
	pBamPair2Smut.args.nthread       = 1
	pBamPair2Smut.args.tmpdir        = params.tmpdir.value
	pBamPair2Smut.envs.bashimport    = bashimport
	pBamPair2Smut.preCmd             = """
	{{bashimport}} reference.bash
	export samtools={{args.samtools | squote}}
	export picard={{args.picard | squote}}
	reference fasta {{args.ref | squote}}
	reference picard {{args.ref | squote}}
	"""
	pBamPair2Smut.lang   = params.python.value
	pBamPair2Smut.script = "file:scripts/sambam/pBamPair2Smut.py"
	return pBamPair2Smut

@procfactory
def _pBamStats():
	"""
	@name:
		pBamStats
	@description:
		Get read depth from bam files.
	@input:
		`infile:file`: The input bam file
	@output:
		`outfile:file`: The output statistic file
		`outdir:dir`:   The directory containing result files and figures.
	@args:
		`tool`: The tool used to do the job. Default: bamstats
		`bamstats`: The path to bamstats. Default: bamstats
		`params`: Other params to main program. Default: `{}`
		`mem`: The memory to be used. Default: 16G
		`plot`: Whether plot the result. Default: True
	"""
	pBamStats                  = Proc(desc = 'Get read depth from bam files.')
	pBamStats.input            = 'infile:file'
	pBamStats.output           = 'outfile:file:{{i.infile | fn}}/{{i.infile | fn}}.stat.txt, outdir:dir:{{i.infile | fn}}'
	pBamStats.args.tool        = 'bamstats'
	pBamStats.args.bamstats    = params.bamstats.value
	pBamStats.args.params      = Box()
	pBamStats.args.mem         = params.mem16G.value
	pBamStats.args.plot        = True
	pBamStats.args.histplotggs = Box(xlab = {0: 'Read depth'}, ylab = {0: '# samples'})
	pBamStats.args.boxplotggs  = Box(xlab = {0: 'Counts'})
	pBamStats.args.devpars     = Box(res = 300, width = 2000, height = 2000)
	pBamStats.args.cap         = 500
	pBamStats.args.cutoff      = 0
	pBamStats.args.nfeats      = 40
	pBamStats.args.feature     = 'wgs'
	pBamStats.envs.rimport     = rimport
	pBamStats.preCmd        = """
	if [[ "{{args.plot | R}}" == "TRUE" && {{proc.forks}} -lt 2 ]]; then
		echo "Plots can only be done with proc.forks >= 2." 1>&2
		exit 1
	fi
	"""
	pBamStats.lang   = params.Rscript.value
	pBamStats.script = "file:scripts/sambam/pBamStats.r"
	return pBamStats

@procfactory
def _pBam2Fastq():
	"""
	@input:
		infile: The sam/bam file.
			- Sam files only available for biobambam, picard
	@output:
		fqfile1: The 1st match of paired reads
		fqfile2: The 2nd match of paired reads
	@args:
		`tool`     : The tool to use. Default: biobambam (bedtools, samtools, picard)
		`biobambam`: The path of bamtofastq of biobambam. Default: bamtofastq
		`bedtools` : The path of bedtools. Default: bedtools
		`samtools` : The path of samtools. Default: samtools
		`picard`   : The path of picard. Default: picard
		`mem`      : The memory to be used by picard. Default: 8G
		`gz`       : Whether gzip the output files. Default: True
		`params`:  : Other params for `tool`. Default: ''
		`tmpdir`   : The tmpdir. Default: `__import__('tempfile').gettempdir()`
	@requires:
		[picard](https://broadinstitute.github.io/picard/command-line-overview.html)
		[biobambam](https://github.com/gt1/biobambam2)
		[samtools](https://github.com/samtools/samtools)
		[bedtools](http://bedtools.readthedocs.io/en/latest/content/bedtools-suite.html)
	"""
	return Box(
		desc   = 'Convert sam/bam files to pair-end fastq files.',
		input  = "infile:file",
		output = [
			"fqfile1:file:{{ i.infile | fn }}_1.fastq{% if args.gz %}.gz{% endif %}",
			"fqfile2:file:{{ i.infile | fn }}_2.fastq{% if args.gz %}.gz{% endif %}"
		],
		lang = params.python.value,
		args = Box(
			tool      = 'biobambam',
			biobambam = params.biobambam_bamtofastq.value,
			bedtools  = params.bedtools.value,
			samtools  = params.samtools.value,
			picard    = params.picard.value,
			mem       = params.mem8G.value, # only for picard
			gz        = False,
			params    = Box(),
			tmpdir    = params.tmpdir.value))

@procfactory
def _pBam2FastqSE():
	"""
	@name:
		pBam2FastqSE
	@description:
		Convert sam/bam files to single-end fastq files.
	@input:
		`infile:file`: The sam/bam file.
			- Sam files only available for biobambam, picard
	@output:
		`fqfile:file`: The fastq file
	@args:
		`tool`     : The tool to use. Default: biobambam (bedtools, samtools, picard)
		`biobambam`: The path of bamtofastq of biobambam. Default: bamtofastq
		`bedtools` : The path of bedtools. Default: bedtools
		`samtools` : The path of samtools. Default: samtools
		`picard`   : The path of picard. Default: picard
		`mem`      : The memory to be used by picard. Default: 8G
		`gz`       : Whether gzip the output files. Default: True
		`params`:  : Other params for `tool`. Default: ''
		`tmpdir`   : The tmpdir. Default: `__import__('tempfile').gettempdir()`
	@requires:
		[picard](https://broadinstitute.github.io/picard/command-line-overview.html)
		[biobambam](https://github.com/gt1/biobambam2)
		[samtools](https://github.com/samtools/samtools)
		[bedtools](http://bedtools.readthedocs.io/en/latest/content/bedtools-suite.html)
	"""
	pBam2FastqSE                     = Proc(desc = 'Convert bam files to single-end fastq files.')
	pBam2FastqSE.input               = "infile:file"
	pBam2FastqSE.output              = "fqfile:file:{{i.infile | fn }}.fq{% if args.gz %}.gz{% endif %}"
	pBam2FastqSE.args.tool           = 'biobambam'
	pBam2FastqSE.args.biobambam      = params.biobambam_bamtofastq.value
	pBam2FastqSE.args.bedtools       = params.bedtools.value
	pBam2FastqSE.args.samtools       = params.samtools.value
	pBam2FastqSE.args.picard         = params.picard.value
	pBam2FastqSE.args.mem            = params.mem8G.value # only for picard
	pBam2FastqSE.args.gz             = False
	pBam2FastqSE.args.params         = Box()
	pBam2FastqSE.args.tmpdir         = params.tmpdir.value
	#pBam2FastqSE.envs.runcmd         = runcmd.py
	#pBam2FastqSE.envs.mem2           = mem2.py
	#pBam2FastqSE.envs.params2CmdArgs = helpers.params2CmdArgs.py
	pBam2FastqSE.lang                = params.python.value
	pBam2FastqSE.script              = "file:scripts/sambam/pBam2FastqSE.py"
	return pBam2FastqSE

@procfactory
def _pBam2Counts():
	"""
	@name:
		pBam2Counts
	@description:
		Extract read counts from RNA-seq bam files.
	@input:
		`infile:file`: The input bam files
	@outfile:
		`outfile:file`: The count file
	@args:
		`tool`: The tool used to extract counts. Default: ht-seq
		`htseq`: The path of htseq-count.
		`params`: Other params for main program.
		`refgene`: The reference gene in GTF format.
	@requires:
		[`htseq`](https://htseq.readthedocs.io/)
	"""
	pBam2Counts                     = Proc(desc = 'Extract read counts from RNA-seq bam files.')
	pBam2Counts.input               = 'infile:file'
	pBam2Counts.output              = 'outfile:file:{{i.infile | fn}}.counts'
	pBam2Counts.args.tool           = 'htseq'
	pBam2Counts.args.htseq          = params.htseq_count.value
	pBam2Counts.args.params         = Box()
	pBam2Counts.args.refgene        = params.refexon.value
	pBam2Counts.lang                = params.python.value
	pBam2Counts.script              = "file:scripts/sambam/pBam2Counts.py"
	return pBam2Counts

@procfactory
def _pBamIndex():
	"""
	@name:
		pBamIndex
	@description:
		Index bam files.
	@input:
		`infile:file`: The input bam file
	@output:
		`outfile:file`: The symbolic link to the input file
		`outidx:file` : The index file
	@args:
		`samtools`: Path to samtools. Default: `params.samtools`
		`params`  : Other parameters for samtools. Default: `Box(b = True)`
		`nthread` : # threads to use. Default: `1`
	"""
	pBamIndex               = Proc(desc = 'Index bam files')
	pBamIndex.input         = 'infile:file'
	pBamIndex.output        = 'outfile:file:{{i.infile | bn}}, outidx:file:{{i.infile | bn}}.bai'
	pBamIndex.args.tool     = 'samtools'
	pBamIndex.args.samtools = params.samtools.value
	pBamIndex.args.params   = Box(b = True)
	pBamIndex.args.nthread  = 1
	pBamIndex.lang          = params.python.value
	pBamIndex.script        = "file:scripts/sambam/pBamIndex.py"
	return pBamIndex

