"""
A set of processes to generate/process sam/bam files
"""
from pyppl import Proc
from diot import Diot
from . import params, proc_factory

pSam2Bam = proc_factory(
	desc   = 'Deal with mapped sam/bam files, including sort, markdup, rmdup, and/or index.',
	config = Diot(annotate = """
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
	""", runcmd_pre = """
		{{"reference.bash" | bashimport}}
		export elprep={{args.elprep | quote}}
		if [[ {{args.tool | quote}} == "elprep" ]]; then
			reference elprep {{args.ref | squote}}
		fi
	"""),
	lang   = params.python.value,
	input  = "infile:file",
	output = "outfile:file:{{i.infile | fn}}.bam, outidx:file:{{i.infile | fn}}.bam.bai",
	errhow = 'retry',
	args   = Diot(
		tool       = "elprep",
		sambamba   = params.sambamba.value,
		picard     = params.picard.value,
		biobambam  = params.biobambam_bamsort.value,
		samtools   = params.samtools.value,
		elprep     = params.elprep.value,
		steps      = Diot(sort=True, index=True, markdup=True, rmdup=True, recal=True),
		tmpdir     = params.tmpdir.value,
		sortby     = "coordinate",
		nthread    = 1,
		params     = Diot(),
		mem        = params.mem16G.value,
		ref        = params.ref.value,
		knownSites = ''))

pBamMarkdup = proc_factory(
	desc = 'Mark/remove duplicates for bam files.',
	config = Diot(annotate = """
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
	"""))
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
pBamMarkdup.args.params            = Diot()
pBamMarkdup.args.mem               = params.mem16G.value
#pBamMarkdup.envs.mem2              = mem2.py
#pBamMarkdup.envs.runcmd            = runcmd.py
#pBamMarkdup.envs.params2CmdArgs    = helpers.params2CmdArgs.py
pBamMarkdup.lang                   = params.python.value

pBamRecal = proc_factory(
	desc = 'Recalibrate a bam file.',
	config = Diot(annotate = """
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
	""", runcmd_pre = """
		{{"reference.bash" | bashimport}}
		export samtools={{args.samtools | squote}}
		reference fasta {{args.ref | squote}}
	"""))
pBamRecal.input         = "infile:file"
pBamRecal.output        = "outfile:file:{{i.infile | fn}}.bam, outidx:file:{{i.infile | fn}}.bam.bai"
pBamRecal.args.tool     = "bamutil"
pBamRecal.args.gatk     = params.gatk.value
pBamRecal.args.samtools = params.samtools.value
pBamRecal.args.picard   = params.picard.value
pBamRecal.args.bamutil  = params.bamutil.value
pBamRecal.args.params   = Diot(
	# For GATK
	# RealignerTargetCreator = Diot(),
	# IndelRealigner         = Diot(),
	# BaseRecalibrator       = Diot(),
	# PrintReads             = Diot()
)
pBamRecal.args.ref        = params.ref.value
pBamRecal.args.tmpdir     = params.tmpdir.value
pBamRecal.args.knownSites = params.dbsnp.value
pBamRecal.args.nthread    = 1
pBamRecal.args.mem        = params.mem32G.value
pBamRecal.lang            = params.python.value

pBamReadGroup = proc_factory(
	desc = 'Add or replace read groups of a bam file.',
	config = Diot(annotate = """
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
	"""))
pBamReadGroup.input               = "infile:file"
pBamReadGroup.output              = "outfile:file:{{i.infile | bn}}"
pBamReadGroup.args.tool           = "bamutil"
pBamReadGroup.args.picard         = params.picard.value
pBamReadGroup.args.bamutil        = params.bamutil.value
pBamReadGroup.args.rg             = Diot({'id': '', 'pl': 'Illumina', 'pu': 'unit1', 'lb': 'lib1', 'sm': ''})
pBamReadGroup.args.params         = Diot()
pBamReadGroup.args.tmpdir         = params.tmpdir.value
pBamReadGroup.args.mem            = params.mem4G.value
#pBamReadGroup.envs.params2CmdArgs = helpers.params2CmdArgs.py
pBamReadGroup.lang                = params.python.value

pBamReorder = proc_factory(
	desc = 'Reorder a sam/bam file by a given reference.',
	config = Diot(annotate = """
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
	"""))
pBamReorder.input               = "infile:file"
pBamReorder.output              = "outfile:file:{{i.infile | bn}}"
pBamReorder.args.picard         = params.picard.value
pBamReorder.args.params         = Diot()
pBamReorder.args.tmpdir         = params.tmpdir.value
pBamReorder.args.mem            = params.mem4G.value
pBamReorder.args.ref            = params.ref.value

pBamMerge = proc_factory(
	desc = 'Merges multiple SAM and/or BAM sorted files into a single file.',
	config = Diot(annotate = """
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
	"""))
pBamMerge.input               = "infiles:files"
pBamMerge.output              = "outfile:file:{{i.infiles[0] | fn2}}.etc.bam"
pBamMerge.args.tool           = "picard"
pBamMerge.args.picard         = params.picard.value
pBamMerge.args.bamutil        = params.bamutil.value
pBamMerge.args.samtools       = params.samtools.value
pBamMerge.args.sambamba       = params.sambamba.value
pBamMerge.args.params         = Diot()
pBamMerge.args.tmpdir         = params.tmpdir.value
pBamMerge.args.nthread        = 1
pBamMerge.args.mem            = params.mem4G.value
#pBamMerge.envs.fsDirname      = dirnameFiles
pBamMerge.lang                = params.python.value

pBam2Gmut = proc_factory(
	desc   = 'Call germline (snps and indels) from a call-ready bam file.',
	lang   = params.python.value,
	config = Diot(annotate = """
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
	"""))
pBam2Gmut.input             = "infile:file"
pBam2Gmut.output            = "outfile:file:{{i.infile | fn}}.vcf{{args.gz | lambda x: '.gz' if x else ''}}"
pBam2Gmut.args.tool         = "strelka"
pBam2Gmut.args.gatk         = params.gatk.value
pBam2Gmut.args.vardict      = params.vardict.value
pBam2Gmut.args.snvsniffer   = params.snvsniffer.value
pBam2Gmut.args.samtools     = params.samtools.value # required by SNVSniffer to generate a bam header file
pBam2Gmut.args.platypus     = params.platypus.value
pBam2Gmut.args.picard       = params.picard.value
pBam2Gmut.args.strelka      = params.strelka_germ.value
pBam2Gmut.args.mem          = params.mem24G.value
pBam2Gmut.args.ref          = params.ref.value
pBam2Gmut.args.tmpdir       = params.tmpdir.value
pBam2Gmut.args.cfgParams    = Diot() # only for strelka
pBam2Gmut.args.params       = Diot()
pBam2Gmut.args.gz           = False
pBam2Gmut.args.nthread      = 1 # for gatk and platypus
pBam2Gmut.config.runcmd_pre = """
{{"reference.bash" | bashimport}}
export samtools={{args.samtools | squote}}
export picard={{args.picard | squote}}
reference fasta {{args.ref | squote}}
reference picard {{args.ref | squote}}
"""

pBamPair2Smut = proc_factory(
	desc   = 'Call somatic mutations from tumor-normal bam pair.',
	config = Diot(annotate = """
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
	""", runcmd_pre = """
		{{"reference.bash" | bashimport}}
		export samtools={{args.samtools | squote}}
		export picard={{args.picard | squote}}
		reference fasta {{args.ref | squote}}
		reference picard {{args.ref | squote}}
	"""),
	lang   = params.python.value,
	input  = "tumor:file, normal:file",
	output = "outfile:file:{{i.tumor | fn}}-{{i.normal | fn}}.vcf{{args.gz | lambda x: '.gz' if x else ''}}",
	args   = Diot(
		tool          = 'strelka',
		gatk          = params.gatk.value, # required for strelka,
		somaticsniper = params.somaticsniper.value,
		strelka       = params.strelka_soma.value, # @2.7.1,
		snvsniffer    = params.snvsniffer.value,
		virmid        = params.virmid.value,
		vardict       = params.vardict.value,
		samtools      = params.samtools.value,
		picard        = params.picard.value,
		configParams  = Diot(), # only for strelka,
		params        = Diot(),
		mem           = params.mem24G.value,
		ref           = params.ref.value,
		gz            = False,
		nthread       = 1,
		tmpdir        = params.tmpdir.value,
	)
)

pBamStats = proc_factory(
	desc = 'Get read depth from bam files.',
	config = Diot(annotate = """
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
	"""))
pBamStats.input             = 'infile:file'
pBamStats.output            = 'outfile:file:{{i.infile | fn}}/{{i.infile | fn}}.stat.txt, outdir:dir:{{i.infile | fn}}'
pBamStats.args.tool         = 'bamstats'
pBamStats.args.bamstats     = params.bamstats.value
pBamStats.args.params       = Diot()
pBamStats.args.mem          = params.mem16G.value
pBamStats.args.plot         = True
pBamStats.args.histplotggs  = Diot(xlab = {0: 'Read depth'}, ylab = {0: '# samples'})
pBamStats.args.boxplotggs   = Diot(xlab = {0: 'Counts'})
pBamStats.args.devpars      = Diot(res = 300, width = 2000, height = 2000)
pBamStats.args.cap          = 500
pBamStats.args.cutoff       = 0
pBamStats.args.nfeats       = 40
pBamStats.args.feature      = 'wgs'
pBamStats.config.runcmd_pre = """
if [[ "{{args.plot | R}}" == "TRUE" && {{proc.forks}} -lt 2 ]]; then
	echo "Plots can only be done with proc.forks >= 2." 1>&2
	exit 1
fi
"""
pBamStats.lang   = params.Rscript.value

pBam2Fastq = proc_factory(
	desc   = 'Convert sam/bam files to pair-end fastq files.',
	config = Diot(annotate = """
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
	"""),
	input  = "infile:file",
	output = [
		"fqfile1:file:{{ i.infile | fn }}_1.fastq{% if args.gz %}.gz{% endif %}",
		"fqfile2:file:{{ i.infile | fn }}_2.fastq{% if args.gz %}.gz{% endif %}"
	],
	lang = params.python.value,
	args = Diot(
		tool      = 'biobambam',
		biobambam = params.biobambam_bamtofastq.value,
		bedtools  = params.bedtools.value,
		samtools  = params.samtools.value,
		picard    = params.picard.value,
		mem       = params.mem8G.value, # only for picard
		gz        = False,
		params    = Diot(),
		tmpdir    = params.tmpdir.value))

pBam2FastqSE = proc_factory(
	desc = 'Convert bam files to single-end fastq files.',
	config = Diot(annotate = """
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
	"""))
pBam2FastqSE.input               = "infile:file"
pBam2FastqSE.output              = "fqfile:file:{{i.infile | fn }}.fq{% if args.gz %}.gz{% endif %}"
pBam2FastqSE.args.tool           = 'biobambam'
pBam2FastqSE.args.biobambam      = params.biobambam_bamtofastq.value
pBam2FastqSE.args.bedtools       = params.bedtools.value
pBam2FastqSE.args.samtools       = params.samtools.value
pBam2FastqSE.args.picard         = params.picard.value
pBam2FastqSE.args.mem            = params.mem8G.value # only for picard
pBam2FastqSE.args.gz             = False
pBam2FastqSE.args.params         = Diot()
pBam2FastqSE.args.tmpdir         = params.tmpdir.value
#pBam2FastqSE.envs.runcmd         = runcmd.py
#pBam2FastqSE.envs.mem2           = mem2.py
#pBam2FastqSE.envs.params2CmdArgs = helpers.params2CmdArgs.py
pBam2FastqSE.lang                = params.python.value

pBam2Counts = proc_factory(
	desc = 'Extract read counts from RNA-seq bam files.',
	config = Diot(annotate = """
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
	"""))
pBam2Counts.input        = 'infile:file'
pBam2Counts.output       = 'outfile:file:{{i.infile | fn}}.counts'
pBam2Counts.args.tool    = 'htseq'
pBam2Counts.args.htseq   = params.htseq_count.value
pBam2Counts.args.params  = Diot()
pBam2Counts.args.refgene = params.refexon.value
pBam2Counts.lang         = params.python.value

pBamIndex = proc_factory(
	desc = 'Index bam files',
	config = Diot(annotate = """
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
		`params`  : Other parameters for samtools. Default: `Diot(b = True)`
		`nthread` : # threads to use. Default: `1`
	"""))
pBamIndex.input         = 'infile:file'
pBamIndex.output        = 'outfile:file:{{i.infile | bn}}, outidx:file:{{i.infile | bn}}.bai'
pBamIndex.args.tool     = 'samtools'
pBamIndex.args.samtools = params.samtools.value
pBamIndex.args.params   = Diot(b = True)
pBamIndex.args.nthread  = 1
pBamIndex.lang          = params.python.value
