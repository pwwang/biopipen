"""
A set of processes to generate/process sam/bam files
"""
from pyppl import Proc, Box
from .utils import mem2, runcmd, buildref, checkref, polling, helpers
from . import params

"""
@name:
	pSam2Bam
@description:
	Deal with mapped sam/bam files, including sort, markdup, and/or index
@input:
	`infile:file`: The input file
@output:
	`outfile:file`: The output bam file
	`idxfile:file`: The index of the output bam file
	- If args.index == False, it'll a link to outfile and should be never used
@args:
	`tool`             : The tool used to do the sort. Default: sambamba (picard|sambamba|biobambam|samtools)
	`sambamba`         : The path of the sambamba. Default: sambamba 
	`picard`           : The path of the picard. Default: picard 
	`biobambam_bamsort`: The path of the biobambam's bamsort. Default: bamsort 
	`samtools`         : The path of the samtools. Default: samtools 
	`sort`             : Do sorting? Default: True 
	- If input is sam, tool is biobambam, this should be True
	`index`            : Do indexing? Default: True
	`markdup`          : Do duplicates marking? Default: False
	- `rmdup` for samtools will be called
	`rmdup`            : Do duplicates removing? Default: False
	`tmpdir`           : The tmp dir used to store tmp files. Default: <system default tmpdir>
	`sortby`           : Sort by coordinate or queryname. Default: coordinate
	`nthread`          : Default: 1
	`informat`         : The format of input file. Default: <detect from extension> (sam|bam)
	`params`           : Other parameters for `tool`. Defaut: ""
	`mem`              : The max memory to use. Default: "16G"
	- Unit could be G/g/M/m
	- Will be converted to -Xmx4G, and -Xms will be 1/8 of it
@requires:
	[sambamba](https://lomereiter.github.io/sambamba/docs/sambamba-view.html) if `args.tool` == samtools or reference used but not indexed.
	[picard](https://broadinstitute.github.io/picard/command-line-overview.html)
	[biobambam](https://github.com/gt1/biobambam2)
	[samtools](https://github.com/samtools/samtools)
"""
pSam2Bam                        = Proc(desc = 'Deal with mapped sam/bam files, including sort, markdup, rmdup, and/or index.')
pSam2Bam.input                  = "infile:file"
pSam2Bam.output                 = "outfile:file:{{in.infile | fn}}.bam, idxfile:file:{{in.infile | fn}}.bam.bai"
pSam2Bam.args.tool              = "biobambam"
pSam2Bam.args.sambamba          = params.sambamba.value
pSam2Bam.args.picard            = params.picard.value
pSam2Bam.args.biobambam_bamsort = params.biobambam_bamsort.value
pSam2Bam.args.samtools          = params.samtools.value
pSam2Bam.args.sort              = True
pSam2Bam.args.index             = True
pSam2Bam.args.markdup           = False
pSam2Bam.args.rmdup             = False
pSam2Bam.args.tmpdir            = params.tmpdir.value
pSam2Bam.args.sortby            = "coordinate"
pSam2Bam.args.nthread           = 1
pSam2Bam.args.informat          = ""
pSam2Bam.args.params            = Box()
pSam2Bam.args.mem               = params.mem16G.value
pSam2Bam.tplenvs.mem2           = mem2.py
pSam2Bam.tplenvs.runcmd         = runcmd.py
pSam2Bam.tplenvs.params2CmdArgs = helpers.params2CmdArgs.py
pSam2Bam.lang                   = params.python.value
pSam2Bam.script                 = "file:scripts/sambam/pSam2Bam.py"

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
pBamMarkdup.output                 = "outfile:file:{{in.infile | fn}}.bam"
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
pBamMarkdup.tplenvs.mem2           = mem2.py
pBamMarkdup.tplenvs.runcmd         = runcmd.py
pBamMarkdup.tplenvs.params2CmdArgs = helpers.params2CmdArgs.py
pBamMarkdup.lang                   = params.python.value
pBamMarkdup.script                 = "file:scripts/sambam/pBamMarkdup.py"

"""
@name:
	pBamRecal
@description:
	Recalibrate a bam file
@input:
	`infile:file`: The bam file
@brings:
	`infile`: {{in.infile | bn}}.bai, the index file of bam
@output:
	`outfile:file`: The output bam file
@args:
	`tool`                         : The tool used to recalibrate the bam file. Default: `gatk` (gatk|bamutil)
	`gatk`                         : The path of gatk, including java path. Default: `gatk`
	`samtools`                     : The path of samtools. Default: `samtools`
	`bamutil`                      : The path of bamutil. Default: `bam`
	`picard`                       : The path of picard. Default: `picard`
	`paramsRealignerTargetCreator` : Other parameters for `gatk RealignerTargetCreator`. Defaut: ""
	`paramsIndelRealigner`         : Other parameters for `gatk IndelRealigner`. Defaut: ""
	`paramsBaseRecalibrator`       : Other parameters for `gatk BaseRecalibrator`. Defaut: ""
	`paramsPrintReads`             : Other parameters for `gatk PrintReads`. Defaut: ""
	`params`                       : Other parameters for `bam recab`. Default: ""
	`mem`                          : The max memory to use. Default: "32G"
	`knownSites`                   : The known polymorphic sites to mask out. Default: "" (Required for GATK)
	`ref`                          : The reference file. Required.
	- Will be converted to -Xmx4G, and -Xms will be 1/8 of it
@requires:
	[gatk](https://software.broadinstitute.org/gatk)
	[samtools](https://github.com/samtools/samtools) if `args.ref` is not indexed, or bamutil is used for bam index file generation.
	[picard](https://broadinstitute.github.io/picard/command-line-overview.html) if `args.ref is not dicted.`
"""
pBamRecal                                   = Proc(desc = 'Recalibrate a bam file.')
pBamRecal.input                             = "infile:file"
pBamRecal.brings                            = {"infile": "{{in.infile | bn}}.bai"}
pBamRecal.output                            = "outfile:file:{{in.infile | fn}}.bam, idxfile:file:{{in.infile | fn}}.bam.bai"
pBamRecal.args.tool                         = "bamutil"
pBamRecal.args.gatk                         = params.gatk.value
pBamRecal.args.samtools                     = params.samtools.value
pBamRecal.args.picard                       = params.picard.value
pBamRecal.args.bamutil                      = params.bamutil.value
pBamRecal.args.paramsRealignerTargetCreator = Box()
pBamRecal.args.paramsIndelRealigner         = Box()
pBamRecal.args.paramsBaseRecalibrator       = Box()
pBamRecal.args.paramsPrintReads             = Box()
pBamRecal.args.params                       = Box()
pBamRecal.args.ref                          = params.ref.value
pBamRecal.args.tmpdir                       = params.tmpdir.value
pBamRecal.args.knownSites                   = ""
pBamRecal.args.mem                          = params.mem32G.value
pBamRecal.tplenvs.mem2                      = mem2.py
pBamRecal.tplenvs.runcmd                    = runcmd.py
pBamRecal.tplenvs.buildrefIndex             = buildref.index.py
pBamRecal.tplenvs.params2CmdArgs            = helpers.params2CmdArgs.py
pBamRecal.beforeCmd                         = checkref.fa.bash + buildref.fai.bash + buildref.dict.bash
pBamRecal.lang                              = params.python.value
pBamRecal.script                            = "file:scripts/sambam/pBamRecal.py"

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
pBamReadGroup                        = Proc(desc = 'Add or replace read groups of a bam file.')
pBamReadGroup.input                  = "infile:file"
pBamReadGroup.output                 = "outfile:file:{{in.infile | bn}}"
pBamReadGroup.args.tool              = "bamutil"
pBamReadGroup.args.picard            = params.picard.value
pBamReadGroup.args.bamutil           = params.bamutil.value
pBamReadGroup.args.rg                = Box({'id': '', 'pl': 'Illumina', 'pu': 'unit1', 'lb': 'lib1', 'sm': ''})
pBamReadGroup.args.params            = Box()
pBamReadGroup.args.tmpdir            = params.tmpdir.value
pBamReadGroup.args.mem               = params.mem4G.value
pBamReadGroup.tplenvs.mem2           = mem2.py
pBamReadGroup.tplenvs.runcmd         = runcmd.py
pBamReadGroup.tplenvs.params2CmdArgs = helpers.params2CmdArgs.py
pBamReadGroup.lang                   = params.python.value
pBamReadGroup.script                 = "file:scripts/sambam/pBamReadGroup.py"


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
pBamReorder                        = Proc(desc = 'Reorder a sam/bam file by a given reference.')
pBamReorder.input                  = "infile:file"
pBamReorder.output                 = "outfile:file:{{in.infile | bn}}"
pBamReorder.args.picard            = params.picard.value
pBamReorder.args.params            = Box()
pBamReorder.args.tmpdir            = params.tmpdir.value
pBamReorder.args.mem               = params.mem4G.value
pBamReorder.args.ref               = params.ref.value
pBamReorder.tplenvs.mem2           = mem2.py
pBamReorder.tplenvs.runcmd         = runcmd.py
pBamReorder.tplenvs.params2CmdArgs = helpers.params2CmdArgs.py
pBamReorder.beforeCmd              = checkref.fa.bash + buildref.dict.bash
pBamReorder.lang                   = params.python.value
pBamReorder.script                 = "file:scripts/sambam/pBamReorder.py"


"""
@name:
	pBamMerge
@description:
	Merges multiple SAM and/or BAM files (must be sorted by coordinate) into a single file.
@input:
	`inlist:file`: The directory containing sam/bam files to be merged
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
pBamMerge                        = Proc(desc = 'Merges multiple SAM and/or BAM sorted files into a single file.')
pBamMerge.input                  = "inlist:file"
pBamMerge.output                 = "outfile:file:{{in.inlist | fn | lambda x: x + '_merged'}}.bam"
pBamMerge.args.tool              = "picard"
pBamMerge.args.picard            = params.picard.value
pBamMerge.args.bamutil           = params.bamutil.value
pBamMerge.args.samtools          = params.samtools.value
pBamMerge.args.sambamba          = params.sambamba.value
pBamMerge.args.params            = Box()
pBamMerge.args.tmpdir            = params.tmpdir.value
pBamMerge.args.nthread           = 1
pBamMerge.args.mem               = params.mem4G.value
pBamMerge.tplenvs.mem2           = mem2.py
pBamMerge.tplenvs.runcmd         = runcmd.py
pBamMerge.tplenvs.params2CmdArgs = helpers.params2CmdArgs.py
pBamMerge.lang                   = params.python.value
pBamMerge.script                 = "file:scripts/sambam/pBamMerge.py"


"""
@name:
	pBam2Gmut
@description:
	Call germline (snps and indels) from a call-ready bam file.
@input:
	`infile:file`: The input bam file
@brings:
	`infile`: `{{in.infile | bn}}.bai`, the bam index file
@output:
	`outfile:file`: The vcf file containing the mutations
@args:
	`tool`:         The tool used to call mutations. Default: gatk (vardict, snvsniffer, platypus, strelka)
	`gatk`:         The path of gatk. Default: gatk
	`vardict`:      The path of vardict. Default: vardict
	`snvsniffer`:   The path of snvsniffer. Default: SNVSniffer
	`samtools`:     The path of samtools. Default: samtools (used to generate reference index)
	`platypus`:     The path of platypus. Default: platypus
	`strelka`:      The path of strelka. Default: configureStrelkaGermlineWorkflow.py
	`configParams`: The params for `strelka` configuration. Default: ""
	`picard`:       The path of picard. Default: picard
	`mem`:          The memory to be used. Default: 32G
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
pBam2Gmut                        = Proc(desc = 'Call germline (snps and indels) from a call-ready bam file.')
pBam2Gmut.input                  = "infile:file"
pBam2Gmut.brings                 = {"infile": "{{in.infile | bn}}.bai"}
pBam2Gmut.output                 = "outfile:file:{{in.infile | fn}}.vcf{{args.gz | lambda x: '.gz' if x else ''}}"
pBam2Gmut.lang                   = params.python.value
pBam2Gmut.args.tool              = "gatk"
pBam2Gmut.args.gatk              = params.gatk.value
pBam2Gmut.args.vardict           = params.vardict.value
pBam2Gmut.args.snvsniffer        = params.snvsniffer.value
pBam2Gmut.args.samtools          = params.samtools.value # required by SNVSniffer to generate a bam header file
pBam2Gmut.args.platypus          = params.platypus.value
pBam2Gmut.args.picard            = params.picard.value
pBam2Gmut.args.strelka           = params.strelka_germ.value
pBam2Gmut.args.mem               = params.mem24G.value
pBam2Gmut.args.ref               = params.ref.value
pBam2Gmut.args.tmpdir            = params.tmpdir.value
pBam2Gmut.args.configParams      = Box() # only for strelka
pBam2Gmut.args.params            = Box()
pBam2Gmut.args.gz                = False
pBam2Gmut.args.nthread           = 1 # for gatk and platypus
pBam2Gmut.tplenvs.mem2           = mem2.py
pBam2Gmut.tplenvs.runcmd         = runcmd.py
pBam2Gmut.tplenvs.params2CmdArgs = helpers.params2CmdArgs.py
pBam2Gmut.beforeCmd              = checkref.fa.bash + buildref.fai.bash + buildref.dict.bash
pBam2Gmut.script                 = "file:scripts/sambam/pBam2Gmut.py"

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
pBamPair2Smut                        = Proc(desc = 'Call somatic mutations from tumor-normal bam pair.')
pBamPair2Smut.input                  = "tumor:file, normal:file"
pBamPair2Smut.brings                 = {"tumor": "{{in.tumor | bn}}.bai", "normal": "{{in.normal | bn}}.bai"}
pBamPair2Smut.output                 = "outfile:file:{{in.tumor | fn | fn}}-{{in.normal | fn | fn}}.vcf{{args.gz | lambda x: '.gz' if x else ''}}"
pBamPair2Smut.args.tool              = 'gatk'
pBamPair2Smut.args.gatk              = params.gatk.value # required for strelka
pBamPair2Smut.args.somaticsniper     = params.somaticsniper.value
pBamPair2Smut.args.strelka           = params.strelka_soma.value # @2.7.1
pBamPair2Smut.args.snvsniffer        = params.snvsniffer.value
pBamPair2Smut.args.virmid            = params.virmid.value
pBamPair2Smut.args.vardict           = params.vardict.value
pBamPair2Smut.args.samtools          = params.samtools.value
pBamPair2Smut.args.picard            = params.picard.value
pBamPair2Smut.args.configParams      = Box() # only for strelka
pBamPair2Smut.args.params            = Box()
pBamPair2Smut.args.mem               = params.mem24G.value
pBamPair2Smut.args.ref               = params.ref.value
pBamPair2Smut.args.gz                = False
pBamPair2Smut.args.nthread           = 1
pBamPair2Smut.args.tmpdir            = params.tmpdir.value
pBamPair2Smut.tplenvs.mem2           = mem2.py
pBamPair2Smut.tplenvs.runcmd         = runcmd.py
pBamPair2Smut.tplenvs.params2CmdArgs = helpers.params2CmdArgs.py
pBamPair2Smut.beforeCmd              = checkref.fa.bash + buildref.fai.bash + buildref.dict.bash
pBamPair2Smut.lang                   = params.python.value
pBamPair2Smut.script                 = "file:scripts/sambam/pBamPair2Smut.py"

"""
@name:
	pBam2Cnv
@description:
	Detect copy number variation from bam files.
@input:
	`input:file`: The bam file
@brings:
	`infile`: "{{in.infile | bn}}.bai" The bam index file
@output:
	`outfile:file`: The output vcf file
	`outdir`: The output directory containing other result files
@args:
	`gz`                    : Whether to gzip the output vcf file. Default: False
	`tool`                  : The tool used to call cnv. Default: 'cnvkit'
	`cnvnator`              : The path of cnvnator. Default: 'cnvnator'
	`cnvnator2vcf`          : The path of cnvnator2VCF. Default: 'cnvnator2VCF.pl'
	`cnvkit`                : The path of cnvkit. Default: 'cnvkit.py'
	`wandy`                 : Tha path of Wandy. Default: 'Wandy'. A `tool.info` file should be with the executable file.
	`ref`                   : The reference file. Required by cnvkit to generate access file. Default: ''
	`cnvkitAccessParams`    : The params for cnvkit access command. Default: '-s 5000'
	`cnvkitTargetParams`    : The params for cnvkit target command. Default: '--split --short-names'
	`cnvkitCoverageParams`  : The params for cnvkit coverage command. Default: ''
	`cnvkitReferenceParams` : The params for cnvkit reference command. Default: '--no-edge'
	`cnvkitFixParams`       : The params for cnvkit fix command. Default: '--no-edge'
	`cnvkitSegmentParams`   : The params for cnvkit segment command. Default: ''
	`cnvkitCallParams`      : The params for cnvkit call command. Default: ''
	`cnvkitPlotParams`      : The params for cnvkit plot command. Default: ''
	`cnvkitBreaksParams`    : The params for cnvkit breaks command. Default: ''
	`cnvkitGainlossParams`  : The params for cnvkit gainloss command. Default: ''
	`cnvkitMetricsParams`   : The params for cnvkit metrics command. Default: ''
	`cnvkitSegmetricsParams`: The params for cnvkit segmetrics command. Default: '--iqr'
	`cnvkitExportParams`    : The params for cnvkit export command. Default: ''
	`cnvkitScatterParams`   : The params for cnvkit scatter command. Default: [''] # multiple scatter plots
	`cnvkitHeatmapParams`   : The params for cnvkit heatmap command. Default: [''] # multiple heatmap plots
	`cnvkitDiagramParams`   : The params for cnvkit diagram command. Default: ''
	`cnvkitReport`          : Generate cnvkit reports? Default: True
	`cnvkitPlot`            : Generate cnvkit plots? Default: True
	`cnvnatorBinsize`       : Bin size for cnvnator. Default: 100
	`cnvnatorGenome`        : Genome for cnvnator. Default: 'hg19'. (NCBI36, hg18, GRCh37, hg19)
	`params`                : The params for `tool`. Default: '-t 1' # wandy 1:hg19 solid cell/blood, 2:hg19 cell free/plamsa, 3:hg38 solid cell/blood, 4:hg38 cell free/plamsa
	`mem`                   : The memory used. Default: '20G' # only for wandy
	`nthread`               : The # threads to use. Default: 1	 # only for cnvkit
@requires:
	[`cnvkit`](http://cnvkit.readthedocs.io/en/stable/index.html)
	[`cnvnator`](https://github.com/abyzovlab/CNVnator)
	`wandy`: Inside cnv caller
"""
pBam2Cnv                             = Proc(desc = 'Detect copy number variation from bam files.')
pBam2Cnv.input                       = 'infile:file'
pBam2Cnv.brings                      = {"infile": "{{in.infile | bn}}.bai"}
pBam2Cnv.output                      = [
	"outfile:file:{{in.infile | fn}}.{{args.tool}}/{{in.infile | fn}}.{{args.tool}}.vcf{% if args.gz %}.gz{% endif %}", 
	"outdir:dir:{{in.infile | fn}}.{{args.tool}}"
]
pBam2Cnv.args.gz                     = False
pBam2Cnv.args.tool                   = 'cnvkit'
pBam2Cnv.args.cnvnator               = params.cnvnator.value
pBam2Cnv.args.cnvnator2vcf           = params.cnvnator2vcf.value
pBam2Cnv.args.cnvkit                 = params.cnvkit.value
pBam2Cnv.args.wandy                  = params.wandy.value
pBam2Cnv.args.ref                    = params.ref.value
pBam2Cnv.args.cnvkitAccessParams     = Box({'s': 5000})
pBam2Cnv.args.cnvkitTargetParams     = Box({'split': True, 'short-names': True})
pBam2Cnv.args.cnvkitCoverageParams   = Box()
pBam2Cnv.args.cnvkitReferenceParams  = Box({'no-edge': True})
pBam2Cnv.args.cnvkitFixParams        = Box({'no-edge': True})
pBam2Cnv.args.cnvkitSegmentParams    = Box()
pBam2Cnv.args.cnvkitCallParams       = Box()
pBam2Cnv.args.cnvkitPlotParams       = Box()
pBam2Cnv.args.cnvkitBreaksParams     = Box()
pBam2Cnv.args.cnvkitGainlossParams   = Box()
pBam2Cnv.args.cnvkitMetricsParams    = Box()
pBam2Cnv.args.cnvkitSegmetricsParams = Box({'iqr': True})
pBam2Cnv.args.cnvkitExportParams     = Box()
pBam2Cnv.args.cnvkitScatterParams    = [Box()] # multiple scatter plots
pBam2Cnv.args.cnvkitHeatmapParams    = [Box()] # multiple heatmap plots
pBam2Cnv.args.cnvkitDiagramParams    = Box()
pBam2Cnv.args.cnvkitReport           = True
pBam2Cnv.args.cnvkitPlot             = True
pBam2Cnv.args.cnvnatorBinsize        = 100
pBam2Cnv.args.cnvnatorGenome         = 'hg19'
pBam2Cnv.args.wandyType              = 1  # wandy 1:hg19 solid cell/blood, 2:hg19 cell free/plamsa, 3:hg38 solid cell/blood, 4:hg38 cell free/plamsa
pBam2Cnv.args.params                 = Box()
pBam2Cnv.args.mem                    = params.mem24G.value # for wandy
pBam2Cnv.args.nthread                = 1 # for cnvkit
pBam2Cnv.tplenvs.pollingFirst        = polling.first.py
pBam2Cnv.tplenvs.pollingAll          = polling.all.py
pBam2Cnv.tplenvs.runcmd              = runcmd.py
pBam2Cnv.tplenvs.params2CmdArgs      = helpers.params2CmdArgs.py
pBam2Cnv.beforeCmd                   = """
if [[ ! -e "{{args.ref}}" && "{{args.tool}}" == "cnvkit" ]]; then
	echo "No reference file specified." 1>&2
	exit 1
fi
if [[ "{{args.tool}}" == "cnvkit" && {{proc.forks}} -lt {{proc.size}} ]]; then
	echo "Cnvkit requires all jobs run simultaneously (proc.forks >= # jobs)." 1>&2
	echo "Because it needs all coverage files to generate reference coverage file." 1>&2
	exit 1
fi
"""
pBam2Cnv.lang   = params.python.value
pBam2Cnv.script = "file:scripts/sambam/pBam2Cnv.py"

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
pBamStats                    = Proc(desc = 'Get read depth from bam files.')
pBamStats.input              = 'infile:file'
pBamStats.output             = 'outfile:file:{{in.infile | fn}}/{{in.infile | fn}}.stat.txt, outdir:dir:{{in.infile | fn}}'
pBamStats.args.tool          = 'bamstats'
pBamStats.args.bamstats      = params.bamstats.value
pBamStats.args.params        = Box()
pBamStats.args.mem           = params.mem16G.value
pBamStats.args.plot          = True
pBamStats.tplenvs.runcmd     = runcmd.r
pBamStats.tplenvs.mem2       = mem2.r
pBamStats.tplenvs.pollingAll = polling.all.r
pBamStats.beforeCmd          = """
if [[ "{{args.plot | R}}" == "TRUE" && {{proc.forks}} -lt {{proc.size}} ]]; then
	echo "Plots can only be done with all jobs run simultaneously (proc.forks >= # jobs)." 1>&2
	exit 1
fi
"""
pBamStats.lang   = params.Rscript.value
pBamStats.script = "file:scripts/sambam/pBamStats.r"

"""
@name:
	pBam2Fastq
@description:
	Convert sam/bam files to pair-end fastq files.
@input:
	`infile:file`: The sam/bam file. 
		- Sam files only available for biobambam, picard
@output:
	`fqfile1:file`: The 1st match of paired reads
	`fqfile2:file`: The 2nd match of paired reads
@args:
	`tool`                : The tool to use. Default: biobambam (bedtools, samtools, picard)
	`biobambam_bamtofastq`: The path of bamtofastq of biobambam. Default: bamtofastq
	`bedtools`            : The path of bedtools. Default: bedtools
	`samtools`            : The path of samtools. Default: samtools
	`picard`              : The path of picard. Default: picard
	`mem`                 : The memory to be used by picard. Default: 8G
	`gz`                  : Whether gzip the output files. Default: True
	`params`:             : Other params for `tool`. Default: ''
	`tmpdir`              : The tmpdir. Default: `__import__('tempfile').gettempdir()`
@requires:
	[picard](https://broadinstitute.github.io/picard/command-line-overview.html)
	[biobambam](https://github.com/gt1/biobambam2)
	[samtools](https://github.com/samtools/samtools)
	[bedtools](http://bedtools.readthedocs.io/en/latest/content/bedtools-suite.html)
"""
pBam2Fastq                           = Proc(desc = 'Convert bam files to pair-end fastq files.')
pBam2Fastq.input                     = "infile:file"
pBam2Fastq.output                    = [
	"fqfile1:file:{{ in.infile | fn | fn | fn }}_1.fastq{% if args.gz %}.gz{% endif %}", 
	"fqfile2:file:{{ in.infile | fn | fn | fn }}_2.fastq{% if args.gz %}.gz{% endif %}"
]
pBam2Fastq.args.tool                 = 'biobambam'
pBam2Fastq.args.biobambam_bamtofastq = params.biobambam_bamtofastq.value
pBam2Fastq.args.bedtools             = params.bedtools.value
pBam2Fastq.args.samtools             = params.samtools.value
pBam2Fastq.args.picard               = params.picard.value
pBam2Fastq.args.mem                  = params.mem8G.value # only for picard
pBam2Fastq.args.gz                   = False
pBam2Fastq.args.params               = Box()
pBam2Fastq.args.tmpdir               = params.tmpdir.value
pBam2Fastq.tplenvs.runcmd            = runcmd.py
pBam2Fastq.tplenvs.params2CmdArgs    = helpers.params2CmdArgs.py
pBam2Fastq.tplenvs.mem2              = mem2.py
pBam2Fastq.lang                      = params.python.value
pBam2Fastq.script                    = "file:./scripts/sambam/pBam2Fastq.py"

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
	`tool`                : The tool to use. Default: biobambam (bedtools, samtools, picard)
	`biobambam_bamtofastq`: The path of bamtofastq of biobambam. Default: bamtofastq
	`bedtools`            : The path of bedtools. Default: bedtools
	`samtools`            : The path of samtools. Default: samtools
	`picard`              : The path of picard. Default: picard
	`mem`                 : The memory to be used by picard. Default: 8G
	`gz`                  : Whether gzip the output files. Default: True
	`params`:             : Other params for `tool`. Default: ''
	`tmpdir`              : The tmpdir. Default: `__import__('tempfile').gettempdir()`
@requires:
	[picard](https://broadinstitute.github.io/picard/command-line-overview.html)
	[biobambam](https://github.com/gt1/biobambam2)
	[samtools](https://github.com/samtools/samtools)
	[bedtools](http://bedtools.readthedocs.io/en/latest/content/bedtools-suite.html)
"""
pBam2FastqSE                           = Proc(desc = 'Convert bam files to single-end fastq files.')
pBam2FastqSE.input                     = "infile:file"
pBam2FastqSE.output                    = "fqfile:file:{{in.infile | fn | fn | fn }}.fq{% if args.gz %}.gz{% endif %}"
pBam2FastqSE.args.tool                 = 'biobambam'
pBam2FastqSE.args.biobambam_bamtofastq = params.biobambam_bamtofastq.value
pBam2FastqSE.args.bedtools             = params.bedtools.value
pBam2FastqSE.args.samtools             = params.samtools.value
pBam2FastqSE.args.picard               = params.picard.value
pBam2FastqSE.args.mem                  = params.mem8G.value # only for picard
pBam2FastqSE.args.gz                   = False
pBam2FastqSE.args.params               = Box()
pBam2FastqSE.args.tmpdir               = params.tmpdir.value
pBam2FastqSE.tplenvs.runcmd            = runcmd.py
pBam2FastqSE.tplenvs.mem2              = mem2.py
pBam2FastqSE.tplenvs.params2CmdArgs    = helpers.params2CmdArgs.py
pBam2FastqSE.lang                      = params.python.value
pBam2FastqSE.script                    = "file:scripts/sambam/pBam2FastqSE.py"

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
pBam2Counts                        = Proc(desc = 'Extract read counts from RNA-seq bam files.')
pBam2Counts.input                  = 'infile:file'
pBam2Counts.output                 = 'outfile:file:{{in.infile | fn}}.counts'
pBam2Counts.args.tool              = 'htseq'
pBam2Counts.args.htseq             = params.htseq_count.value
pBam2Counts.args.params            = Box()
pBam2Counts.args.refgene           = params.refgene.value
pBam2Counts.tplenvs.runcmd         = runcmd.py
pBam2Counts.tplenvs.params2CmdArgs = helpers.params2CmdArgs.py
pBam2Counts.beforeCmd              = checkref.gene.bash
pBam2Counts.lang                   = params.python.value
pBam2Counts.script                 = "file:scripts/sambam/pBam2Counts.py"

