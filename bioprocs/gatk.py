"""GATK utilities"""
from pyppl import Proc
from . import delefactory, procfactory
from modkit import Modkit
Modkit().delegate(delefactory())

#############################
# GATK utilities            #
#############################

_sReferenceFai = """
# make sure reference index is generated
if [[ ! -e "{{brings.reffile#fai}}" ]]; then
	echo "[BioProcs] Cannot find index file for reference: {{reffile.orig}}, generating ..." 1>&2
	({{args.samtools}} faidx "{{reffile | realpath}}" ]] && \\
		ln -s "{{reffile | realpath}}.fai" "{{reffile}}.fai" && \\
		echo "[BioProcs] Generated: {{reffile | realpath}}.fai" 1>&2) || \\
	({{args.samtools}} faidx "{{reffile.orig}}" ]] && \\
		ln -s "{{reffile.orig}}.fai" "{{reffile}}.fai" && \\
		echo "[BioProcs] Generated: {{reffile.orig}}.fai" 1>&2) || \\
	({{args.samtools}} faidx "{{reffile}}" && \\
		echo "[BioProcs] Generated: {{reffile}}.fai" 1>&2)
fi
"""

_sReferenceDict = """
# make sure reference dict is generated
if [[ ! -e "{{brings.reffile#dict}}" ]]; then
	echo "[BioProcs] Cannot find dict file for reference: {{reffile.orig}}, generating ..." 1>&2
	({{args.picard}} CreateSequenceDictionary R="{{reffile}}" O="{{reffile | realpath | prefix}}.dict" && \\
		ln -s "{{reffile | realpath | prefix}}.dict" "{{reffile | prefix}}.dict" && \\
		echo "[BioProcs] Generated: {{reffile | realpath | prefix}}.dict" 1>&2) || \\
	({{args.picard}} CreateSequenceDictionary R="{{reffile}}" O="{{reffile.orig | prefix}}.dict" && \\
		ln -s "{{reffile.orig | prefix}}.dict" "{{reffile | prefix}}.dict" && \\
		echo "[BioProcs] Generated: {{reffile.orig | prefix}}.dict" 1>&2) || \\
	({{args.picard}} CreateSequenceDictionary R="{{reffile}}" O="{{reffile | prefix}}.dict" && \\
		echo "[BioProcs] Generated: {{reffile | prefix}}.dict" 1>&2)
fi
"""

_sBamBai = """
# get index file of bam file
if [[ ! -e "{{brings.bamfile}}" ]]; then
	echo "[BioProcs] Cannot find index file for input file: {{bamfile.orig}}, generating ..." 1>&2
	({{args.samtools}} index "{{bamfile | realpath}}" && \\
		ln -s "{{bamfile | realpath}}.bai" "{{bamfile}}.bai" && \\
		echo "[BioProcs] Generated: {{bamfile | realpath}}.bai" 1>&2) || \\
	({{args.samtools}} index "{{bamfile.orig}}" && \\
		ln -s "{{bamfile.orig}}.bai" "{{bamfile}}.bai" && \\
		echo "[BioProcs] Generated: {{bamfile.orig}}.bai" 1>&2) || \\
	({{args.samtools}} index "{{bamfile}}" && \\
		echo "[BioProcs] Generated: {{bamfile}}.bai" 1>&2)
fi
"""

_sTumorBai = """
if [[ ! -e "{{brings.tumor}}" ]]; then
	echo "[BioProcs] Cannot find index file for tumor: {{tumor.orig}}, generating ..." 1>&2
	({{args.samtools}} index "{{tumor | realpath}}" && \\
		ln -s "{{tumor | realpath}}.bai" "{{tumor}}.bai" && \\
		echo "[BioProcs] Generated: {{tumor | realpath}}.bai" 1>&2) || \\
	({{args.samtools}} index "{{tumor.orig}}" && \\
		ln -s "{{tumor.orig}}.bai" "{{tumor}}.bai" && \\
		echo "[BioProcs] Generated: {{tumor.orig}}.bai" 1>&2) || \\
	({{args.samtools}} index "{{tumor}}" && \\
		echo "[BioProcs] Generated: {{tumor}}.bai" 1>&2)
fi
"""

_sNormalBai = """
if [[ ! -e "{{brings.normal}}" ]]; then
	echo "[BioProcs] Cannot find index file for tumor: {{tumor.orig}}, generating ..." 1>&2
	({{args.samtools}} index "{{normal | realpath}}" && \\
		ln -s "{{normal | realpath}}.bai" "{{normal}}.bai" && \\
		echo "[BioProcs] Generated: {{normal | realpath}}.bai" 1>&2) || \\
	({{args.samtools}} index "{{normal.orig}}" && \\
		ln -s "{{normal.orig}}.bai" "{{normal}}.bai" &&
		echo "[BioProcs] Generated: {{normal.orig}}.bai" 1>&2) || \\
	({{args.samtools}} index "{{normal}}" && \\
		echo "[BioProcs] Generated: {{normal}}.bai" 1>&2)
fi
"""

@procfactory
def _pRealignerTargetCreator():
	"""
	@name:
		pRealignerTargetCreator
	@description:
		The local realignment process is designed to consume one or more BAM files and to locally realign reads such that the number of mismatching bases is minimized across all the reads. In general, a large percent of regions requiring local realignment are due to the presence of an insertion or deletion (indels) in the individual's genome with respect to the reference genome. Such alignment artifacts result in many bases mismatching the reference near the misalignment, which are easily mistaken as SNPs. Moreover, since read mapping algorithms operate on each read independently, it is impossible to place reads on the reference genome such that mismatches are minimized across all reads. Consequently, even when some reads are correctly mapped with indels, reads covering the indel near just the start or end of the read are often incorrectly mapped with respect the true indel, also requiring realignment. Local realignment serves to transform regions with misalignments due to indels into clean reads containing a consensus indel suitable for standard variant discovery approaches.
		Note that indel realignment is no longer necessary for variant discovery if you plan to use a variant caller that performs a haplotype assembly step, such as HaplotypeCaller or MuTect2. However it is still required when using legacy callers such as UnifiedGenotyper or the original MuTect. There are 2 steps to the realignment process:
		- Determining (small) suspicious intervals which are likely in need of realignment (RealignerTargetCreator)
		- Running the realigner over those intervals (see the IndelRealigner tool)
		For more details, see [the indel realignment method documentation](http://www.broadinstitute.org/gatk/guide/article?id=38).
	@input:
		`bamfile:file`:  The aligned bam file
		`reffile`: The reference file
	@brings:
		`bamfile`: `{{bamfile | bn}}.bai` The index file of input bam file
		`reffile#fai`: `{{reffile | bn}}.fai`
		`reffile#dict`: `{{reffile | bn}}.dict`
	@output:
		`outfile:file`: A list of target intervals to pass to the IndelRealigner.
	@args:
		`gatk`:     The gatk executable, default: "gatk"
		`picard`:   The picard executable, default: "picard"
		`params`:   Other parameters for RealignerTargetCreator, default: ""
		`samtools`: The samtools executable, default: "samtools"
		`tmpdir`:  The tmpdir to use. Default: /tmp
		`javamem`: The memory for java vm. Default: "-Xms1g -Xmx8g"
	@requires:
		[GATK](https://software.broadinstitute.org/gatk)
		[samtools](http://www.htslib.org/) if `reffile` is not indexed or `bamfile` is not indexed.
		[picard](https://broadinstitute.github.io/picard/) if `reffile` is not dicted.
	"""
	pRealignerTargetCreator = Proc()
	pRealignerTargetCreator.input  = "bamfile:file, reffile:file"
	# pRealignerTargetCreator.brings = {"bamfile": "{{bamfile | bn}}.bai", "reffile#dict": "{{reffile | fn}}.dict", "reffile#fai": "{{reffile | bn}}.fai"}
	pRealignerTargetCreator.output = "outfile:file:{{bamfile | fn}}.indelRealigner.intervals"
	pRealignerTargetCreator.args   = { "gatk": "gatk", "picard": "picard", "samtools": "samtools", "javamem": "-Xms1g -Xmx8g", "tmpdir": "/tmp", "params": "" }
	pRealignerTargetCreator.script = """
	%s%s%s
	# run cmd
	{{args.gatk}} -T RealignerTargetCreator {{args.javamem}} -Djava.io.tmpdir="{{args.tmpdir}}" -R "{{reffile}}" -I "{{bamfile}}" -o "{{outfile}}" {{args.params}}
	""" % (_sReferenceFai, _sReferenceDict, _sBamBai)
	return pRealignerTargetCreator

@procfactory
def _pIndelRealigner ():
	"""
	@name:
		pIndelRealigner
	@description:
		The local realignment process is designed to consume one or more BAM files and to locally realign reads such that the number of mismatching bases is minimized across all the reads. In general, a large percent of regions requiring local realignment are due to the presence of an insertion or deletion (indels) in the individual's genome with respect to the reference genome. Such alignment artifacts result in many bases mismatching the reference near the misalignment, which are easily mistaken as SNPs. Moreover, since read mapping algorithms operate on each read independently, it is impossible to place reads on the reference genome such at mismatches are minimized across all reads. Consequently, even when some reads are correctly mapped with indels, reads covering the indel near just the start or end of the read are often incorrectly mapped with respect the true indel, also requiring realignment. Local realignment serves to transform regions with misalignments due to indels into clean reads containing a consensus indel suitable for standard variant discovery approaches.
		Note that indel realignment is no longer necessary for variant discovery if you plan to use a variant caller that performs a haplotype assembly step, such as HaplotypeCaller or MuTect2. However it is still required when using legacy callers such as UnifiedGenotyper or the original MuTect.
		There are 2 steps to the realignment process:
		- Determining (small) suspicious intervals which are likely in need of realignment (see the RealignerTargetCreator tool)
		- Running the realigner over those intervals (IndelRealigner)
		For more details, see [the indel realignment method documentation](http://www.broadinstitute.org/gatk/guide/article?id=38).
	@input:
		`bamfile:file`: The aligned bam file
		`intfile:file`: Intervals file output from RealignerTargetCreator
		`reffile:file`: The reference file
	@brings:
		`bamfile`: `{{bamfile | bn}}.bai` The index file of input bam file
		`reffile#fai`: `{{reffile | bn}}.fai`
		`reffile#dict`: `{{reffile | bn}}.dict`
	@output:
		`outfile:file`: A realigned version of input BAM file.
	@args:
		`gatk`:     The gatk executable, default: "gatk"
		`picard`:   The picard executable, default: "picard"
		`params`:  Other parameters for IndelRealigner, default: ""
		`samtools`: The samtools executable, default: samtools
		`tmpdir`:  The tmpdir to use. Default: /tmp
		`javamem`: The memory for java vm. Default: "-Xms1g -Xmx8g"
	@requires:
		[GATK](https://software.broadinstitute.org/gatk)
		[samtools](http://www.htslib.org/) if `reffile` is not indexed or `bamfile` is not indexed.
		[picard](https://broadinstitute.github.io/picard/) if `reffile` is not dicted.
	"""
	pIndelRealigner = Proc()
	pIndelRealigner.input  = "bamfile:file, intfile:file, reffile:file"
	# pIndelRealigner.brings = {"bamfile": "{{bamfile | bn}}.bai", "reffile#dict": "{{reffile | fn}}.dict", "reffile#fai": "{{reffile | bn}}.fai"}
	pIndelRealigner.output = "outfile:file:{{bamfile | fn | lambda x: __import__('re').sub(r'(\\.dedup)?$', '', x)}}.realigned.bam"
	pIndelRealigner.args   = { "gatk": "gatk", "samtools": "samtools", "picard": "picard", "javamem": "-Xms1g -Xmx8g", "tmpdir": "/tmp", "params": "" }
	pIndelRealigner.script = """
	%s%s%s
	return pIndelRealigner
	{{args.gatk}} -T IndelRealigner {{args.javamem}} -Djava.io.tmpdir="{{args.tmpdir}}" -R "{{reffile}}" -I "{{bamfile}}" -o "{{outfile}}" -targetIntervals "{{intfile}}" {{args.params}}
	""" % (_sReferenceFai, _sReferenceDict, _sBamBai)
	return pIndelRealigner

@procfactory
def _pBaseRecalibrator():
	"""
	@name:
		pBaseRecalibrator
	@description:
		Variant calling algorithms rely heavily on the quality scores assigned to the individual base calls in each sequence read. These scores are per-base estimates of error emitted by the sequencing machines. Unfortunately the scores produced by the machines are subject to various sources of systematic technical error, leading to over- or under-estimated base quality scores in the data. Base quality score recalibration (BQSR) is a process in which we apply machine learning to model these errors empirically and adjust the quality scores accordingly. This allows us to get more accurate base qualities, which in turn improves the accuracy of our variant calls. The base recalibration process involves two key steps: first the program builds a model of covariation based on the data and a set of known variants (which you can bootstrap if there is none available for your organism), then it adjusts the base quality scores in the data based on the model. There is an optional but highly recommended step that involves building a second model and generating before/after plots to visualize the effects of the recalibration process. This is useful for quality control purposes. This tool performs the first step described above: it builds the model of covariation and produces the recalibration table. It operates only at sites that are not in dbSNP; we assume that all reference mismatches we see are therefore errors and indicative of poor base quality. This tool generates tables based on various user-specified covariates (such as read group, reported quality score, cycle, and context). Assuming we are working with a large amount of data, we can then calculate an empirical probability of error given the particular covariates seen at this site, where p(error) = num mismatches / num observations. The output file is a table (of the several covariate values, number of observations, number of mismatches, empirical quality score).
	@input:
		`bamfile:file`: A BAM file containing data that needs to be recalibrated.
		`reffile:file`: The reference file
	@brings:
		`bamfile`: `{{bamfile | bn}}.bai` The index file of input bam file
		`reffile#fai`: `{{reffile | bn}}.fai`
		`reffile#dict`: `{{reffile | bn}}.dict`
	@output:
		`outfile:file`: A GATKReport file with many tables:
			- The list of arguments
			- The quantized qualities table
			- The recalibration table by read group
			- The recalibration table by quality score
			- The recalibration table for all the optional covariates
	@args:
		`gatk`:     The gatk executable, default: "gatk"
		`params`:  Other parameters for BaseRecalibrator, default: ""
		`knownSites`: The known polymorphic sites to mask out, required
		`samtools`: The samtools executable, default: samtools
		`picard`:   The picard executable, default: "picard"
		`tmpdir`:  The tmpdir to use. Default: /tmp
		`javamem`: The memory for java vm. Default: "-Xms1g -Xmx8g"
	@requires:
		[GATK](https://software.broadinstitute.org/gatk)
		[samtools](http://www.htslib.org/) if `reffile` is not indexed or `bamfile` is not indexed.
		[picard](https://broadinstitute.github.io/picard/) if `reffile` is not dicted.
	"""
	pBaseRecalibrator = Proc()
	pBaseRecalibrator.input  = "bamfile:file, reffile:file"
	# pBaseRecalibrator.brings = {"bamfile": "{{bamfile | bn}}.bai", "reffile#dict": "{{reffile | fn}}.dict", "reffile#fai": "{{reffile | bn}}.fai"}
	pBaseRecalibrator.output = "outfile:file:{{bamfile | fn}}.recal.table"
	pBaseRecalibrator.args   = { "gatk": "gatk", "samtools": "samtools", "picard": "picard", "javamem": "-Xms1g -Xmx8g", "tmpdir": "/tmp", "knownSites": "", "params": "" }
	pBaseRecalibrator.script = """
	if [[ -z "{{args.knownSites}}" ]]; then
		echo "knownSites file not specified!" 1>&2
		exit 1
	fi
	%s%s%s
	{{args.gatk}} -T BaseRecalibrator {{args.javamem}} -Djava.io.tmpdir="{{args.tmpdir}}" -R "{{reffile}}" -I "{{bamfile}}" -o "{{outfile}}" -knownSites "{{args.knownSites}}" {{args.params}}
	""" % (_sReferenceFai, _sReferenceDict, _sBamBai)
	return pBaseRecalibrator

@procfactory
def _pPrintReads():
	"""
	@name:
		pPrintReads
	@description:
		PrintReads is a generic utility tool for manipulating sequencing data in SAM/BAM format. It can dynamically merge the contents of multiple input BAM files, resulting in merged output sorted in coordinate order. It can also optionally filter reads based on various read properties such as read group tags using the `--read_filter/-rf` command line argument (see documentation on read filters for more information).
		Note that when PrintReads is used as part of the Base Quality Score Recalibration workflow, it takes the `--BQSR` engine argument, which is listed under Inherited Arguments > CommandLineGATK below.
	@input:
		`bamfile:file`:    A BAM file.
		`recaltable:file`: The GATKReport file
		`reffile:file`:    The reference file
	@brings:
		`bamfile`: `{{bamfile | bn}}.bai` The index file of input bam file
		`reffile#fai`: `{{reffile | bn}}.fai`
		`reffile#dict`: `{{reffile | bn}}.dict`
	@output:
		`outfile:file`: A single processed bam file.
	@args:
		`gatk`:     The gatk executable, default: "gatk"
		`params`:  Other parameters for PrintReads, default: ""
		`samtools`: The samtools executable, default: samtools
		`picard`:   The picard executable, default: "picard"
		`tmpdir`:  The tmpdir to use. Default: /tmp
		`javamem`: The memory for java vm. Default: "-Xms1g -Xmx8g"
	@requires:
		[GATK](https://software.broadinstitute.org/gatk)
		[samtools](http://www.htslib.org/) if `reffile` is not indexed or `infile` is not indexed.
		[picard](https://broadinstitute.github.io/picard/) if `reffile` is not dicted.
	"""
	pPrintReads = Proc()
	pPrintReads.input  = "bamfile:file, recaltable:file, reffile:file"
	# pPrintReads.brings = {"bamfile": "{{bamfile | bn}}.bai", "reffile#dict": "{{reffile | fn}}.dict", "reffile#fai": "{{reffile | bn}}.fai"}
	pPrintReads.output = "outfile:file:{{bamfile | fn | lambda x: __import__('re').sub(r'(\\.realigned)?$', '', x)}}.bam"
	pPrintReads.args   = { "gatk": "gatk", "samtools": "samtools", "picard": "picard", "javamem": "-Xms1g -Xmx8g", "tmpdir": "/tmp", "params": "" }
	pPrintReads.script = """
	%s%s%s
	{{args.gatk}} -T PrintReads {{args.javamem}} -Djava.io.tmpdir="{{args.tmpdir}}" -R "{{reffile}}" -I "{{bamfile}}" -o "{{outfile}}" -BQSR "{{recaltable}}" {{args.params}}
	""" % (_sReferenceFai, _sReferenceDict, _sBamBai)
	return pPrintReads

@procfactory
def _pHaplotypeCaller():
	"""
	@name:
		pHaplotypeCaller
	@description:
		PrintReads is a generic utility tool for manipulating sequencing data in SAM/BAM format. It can dynamically merge the contents of multiple input BAM files, resulting in merged output sorted in coordinate order. It can also optionally filter reads based on various read properties such as read group tags using the `--read_filter/-rf` command line argument (see documentation on read filters for more information).
		Note that when PrintReads is used as part of the Base Quality Score Recalibration workflow, it takes the `--BQSR` engine argument, which is listed under Inherited Arguments > CommandLineGATK below.
	@input:
		`bamfile:file`: A BAM file.
		`reffile:file`: The reference file
	@brings:
		`bamfile`: `{{bamfile | bn}}.bai` The index file of input bam file
		`reffile#fai`: `{{reffile | bn}}.fai`
		`reffile#dict`: `{{reffile | bn}}.dict`
	@output:
		`outfile:file`: Either a VCF or gVCF file with raw, unfiltered SNP and indel calls.
	@args:
		`gatk`    : The gatk executable, default: "gatk"
		`params`  : Other parameters for HaplotypeCaller, default: ""
		`samtools`: The samtools executable, default: samtools
		`picard`:   The picard executable, default: "picard"
		`tmpdir`:  The tmpdir to use. Default: /tmp
		`javamem`: The memory for java vm. Default: "-Xms1g -Xmx8g"
		`nthread`: Corresponding to -nct option
	@requires:
		[GATK](https://software.broadinstitute.org/gatk)
		[samtools](http://www.htslib.org/) if `reffile` is not indexed or `infile` is not indexed.
		[picard](https://broadinstitute.github.io/picard/) if `reffile` is not dicted.
	"""
	pHaplotypeCaller = Proc()
	pHaplotypeCaller.input  = "bamfile:file, reffile:file"
	# pHaplotypeCaller.brings = {"bamfile": "{{bamfile | bn}}.bai", "reffile#dict": "{{reffile | fn}}.dict", "reffile#fai": "{{reffile | bn}}.fai"}
	pHaplotypeCaller.output = "outfile:file:{{bamfile | fn}}.vcf"
	pHaplotypeCaller.args   = { "gatk": "gatk", "samtools": "samtools", "picard": "picard", "javamem": "-Xms1g -Xmx8g", "tmpdir": "/tmp", "nthread": 1, "params": "" }
	pHaplotypeCaller.script = """
	%s%s%s
	{{args.gatk}} -T HaplotypeCaller {{args.javamem}} -Djava.io.tmpdir="{{args.tmpdir}}" -R "{{reffile}}" -I "{{bamfile}}" -o "{{outfile}}" -nct {{args.nthread}} {{args.params}}
	""" % (_sReferenceFai, _sReferenceDict, _sBamBai)
	return pHaplotypeCaller

@procfactory
def _pSelectVariants():
	"""
	@name:
		pSelectVariants
	@description:
		Often, a VCF containing many samples and/or variants will need to be subset in order to facilitate certain analyses (e.g. comparing and contrasting cases vs. controls; extracting variant or non-variant loci that meet certain requirements, displaying just a few samples in a browser like IGV, etc.). SelectVariants can be used for this purpose.
		There are many different options for selecting subsets of variants from a larger callset:
		- Extract one or more samples from a callset based on either a complete sample name or a pattern match.
		- Specify criteria for inclusion that place thresholds on annotation values, e.g. "DP > 1000" (depth of coverage greater than 1000x), "AF < 0.25" (sites with allele frequency less than 0.25). These - criteria are written as "JEXL expressions", which are documented in the article about using JEXL expressions.
		- Provide concordance or discordance tracks in order to include or exclude variants that are also present in other given callsets.
		- Select variants based on criteria like their type (e.g. INDELs only), evidence of mendelian violation, filtering status, allelicity, and so on.
		There are also several options for recording the original values of certain annotations that are recalculated when a subsetting the new callset, trimming alleles, and so on.
	@input:
		`vcffile:file`: A variant call set from which to select a subset.
		`reffile:file`: The reference file
	@brings:
		`reffile#fai`: `{{reffile | bn}}.fai`
		`reffile#dict`: `{{reffile | bn}}.dict`
	@output:
		`outfile:file`: A new VCF file containing the selected subset of variants.
	@args:
		`gatk`:     The gatk executable, default: "gatk"
		`params`:  Other parameters for SelectVariants, default: ""
		`samtools`: The samtools executable, default: samtools
		`picard`:   The picard executable, default: "picard"
		`tmpdir`:  The tmpdir to use. Default: /tmp
		`javamem`: The memory for java vm. Default: "-Xms1g -Xmx8g"
	@requires:
		[GATK](https://software.broadinstitute.org/gatk)
		[samtools](http://www.htslib.org/) if `reffile` is not indexed or `infile` is not indexed.
		[picard](https://broadinstitute.github.io/picard/) if `reffile` is not dicted.
	"""
	pSelectVariants = Proc()
	pSelectVariants.input  = "vcffile:file, reffile:file"
	# pSelectVariants.brings = {"reffile#fai": "{{reffile | bn}}.fai", "reffile#dict": "{{reffile | fn}}.dict"}
	pSelectVariants.output = "outfile:file:{{vcffile | fn}}.selected.vcf"
	pSelectVariants.args   = {"gatk": "gatk", "samtools": "samtools", "picard": "picard", "javamem": "-Xms1g -Xmx8g", "tmpdir": "/tmp", "params": ""}
	pSelectVariants.script = """
	%s%s
	{{args.gatk}} -T SelectVariants {{args.javamem}} -Djava.io.tmpdir="{{args.tmpdir}}" -R "{{reffile}}" -V "{{vcffile}}" -o "{{outfile}}" {{args.params}}
	""" % (_sReferenceFai, _sReferenceDict)
	return pSelectVariants

@procfactory
def _pVariantFiltration():
	"""
	@name:
		pVariantFiltration
	@description:
		This tool is designed for hard-filtering variant calls based on certain criteria. Records are hard-filtered by changing the value in the FILTER field to something other than PASS. Filtered records will be preserved in the output unless their removal is requested in the command line.
		The most common way of specifying filtering criteria is by using JEXL queries. See the article on JEXL expressions in the documentation Guide for detailed information and examples.
	@input:
		`vcffile:file`: A variant call set from which to select a subset.
		`reffile:file`: The reference file
	@brings:
		`reffile#fai`: `{{reffile | bn}}.fai`
		`reffile#dict`: `{{reffile | bn}}.dict`
	@output:
		`outfile:file`: A filtered VCF.
	@args:
		`gatk`:     The gatk executable, default: "gatk -T VariantFiltration"
		`params`:  Other parameters for VariantFiltration, default: ""
		`samtools`: The samtools executable, default: samtools
		`picard`:   The picard executable, default: "picard"
		`tmpdir`:  The tmpdir to use. Default: /tmp
		`javamem`: The memory for java vm. Default: "-Xms1g -Xmx8g"
	@requires:
		[GATK](https://software.broadinstitute.org/gatk)
		[samtools](http://www.htslib.org/) if `reffile` is not indexed or `infile` is not indexed.
		[picard](https://broadinstitute.github.io/picard/) if `reffile` is not dicted.
	"""
	pVariantFiltration = Proc()
	pVariantFiltration.input  = "vcffile:file, reffile:file"
	# pVariantFiltration.brings = {"reffile#fai": "{{reffile | bn}}.fai", "reffile#dict": "{{reffile | fn}}.dict"}
	pVariantFiltration.output = "outfile:file:{{vcffile | fn | lambda x: __import__('re').sub(r'(\\.selected)?$', '', x)}}.filtered.vcf"
	pVariantFiltration.args   = { "gatk": "gatk", "samtools": "samtools", "picard": "picard", "javamem": "-Xms1g -Xmx8g", "tmpdir": "/tmp", "params": "" }
	pVariantFiltration.script = """
	%s%s
	{{args.gatk}} -T VariantFiltration {{args.javamem}} -Djava.io.tmpdir="{{args.tmpdir}}" -R "{{reffile}}" -V "{{vcffile}}" -o "{{outfile}}" {{args.params}}
	""" % (_sReferenceFai, _sReferenceDict)
	return pVariantFiltration

@procfactory
def _pMuTect2():
	"""
	@name:
		pMuTect2
	@description:
		MuTect2 is a somatic SNP and indel caller that combines the DREAM challenge-winning somatic genotyping engine of the original MuTect ([Cibulskis et al., 2013](http://www.nature.com/nbt/journal/v31/n3/full/nbt.2514.html)) with the assembly-based machinery of HaplotypeCaller. The basic operation of MuTect2 proceeds similarly to that of the HaplotypeCaller.
		NOTE: only Tumor/Normal variant calling implemented in bioprocs
	@input:
		`tumor:file`:   the tumor bam file
		`normal:file`:  the normal bam file
		`reffile:file`: the reference file
	@brings:
		`tumor`:  `{{tumor | bn}}.bai` the index file of tumor
		`normal`: `{{normal | bn}}.bai` the index file of normal
		`reffile#fai`: `{{reffile | bn}}.fai`
		`reffile#dict`: `{{reffile | bn}}.dict`
	@output:
		`outfile:file`: The vcf file containing somatic mutations
	@args:
		`gatk`:     The gatk executable, default: "gatk"
		`samtools`: The samtools executable, default: samtools
		`params`:   Other parameters for MuTect2, default: ""
		`picard`:   The picard executable, default: "picard"
		`tmpdir`:  The tmpdir to use. Default: /tmp
		`javamem`: The memory for java vm. Default: "-Xms1g -Xmx8g"
	@requires:
		[GATK](https://software.broadinstitute.org/gatk)
		[samtools](http://www.htslib.org/) if index files of input files are not found
		[picard](https://broadinstitute.github.io/picard/) if `reffile` is not dicted.
	"""
	pMuTect2 = Proc()
	pMuTect2.input  = "tumor:file, normal:file, reffile:file"
	#pMuTect2.brings = {"tumor": "{{tumor | bn}}.bai", "normal": "{{normal | bn}}.bai", "reffile#fai": "{{reffile | bn}}.fai", "reffile#dict": "{{reffile | fn}}.dict"}
	pMuTect2.output = "outfile:file:{{tumor | fn}}-{{normal | fn}}.vcf"
	pMuTect2.args   = { "gatk": "gatk", "samtools": "samtools", "picard": "picard", "javamem": "-Xms1g -Xmx8g", "tmpdir": "/tmp", "params": "" }
	pMuTect2.script = """
	%s%s%s%s
	{{args.gatk}} -T MuTect2 {{args.javamem}} -Djava.io.tmpdir="{{args.tmpdir}}" -R "{{reffile}}" -I:tumor "{{tumor}}" -I:normal "{{normal}}" -o "{{outfile}}" {{args.params}}
	""" % (_sNormalBai, _sTumorBai, _sReferenceFai, _sReferenceDict)
	return pMuTect2

@procfactory
def _pMuTect2Interval():
	"""
	@name:
		pMuTect2Interval
	@description:
		Use interval file model of MuTect2
	@input:
		`tumor:file`:   the tumor bam file
		`normal:file`:  the normal bam file
		`reffile:file`: the reference file
	@brings:
		`tumor`:  `{{tumor | bn}}.bai` the index file of tumor
		`normal`: `{{normal | bn}}.bai` the index file of normal
		`reffile#fai`: `{{reffile | bn}}.fai`
		`reffile#dict`: `{{reffile | bn}}.dict`
	@output:
		`outfile:file`: The vcf file containing somatic mutations
	@args:
		`gatk`:     The gatk executable, default: "gatk"
		`samtools`: The samtools executable, default: samtools
		`params`:   Other parameters for MuTect2, default: ""
		`picard`:   The picard executable, default: "picard"
		`tmpdir`:  The tmpdir to use. Default: /tmp
		`javamem`: The memory for java vm. Default: "-Xms1g -Xmx8g"
	@requires:
		[GATK](https://software.broadinstitute.org/gatk)
		[samtools](http://www.htslib.org/) if index files of input files are not found
		[picard](https://broadinstitute.github.io/picard/) if `reffile` is not dicted.
	"""
	pMuTect2Interval = Proc()
	pMuTect2Interval.input  = "tumor:file, normal:file, reffile:file"
	#pMuTect2Interval.brings = {"tumor": "{{tumor | bn}}.bai", "normal": "{{normal | bn}}.bai", "reffile#fai": "{{reffile | bn}}.fai", "reffile#dict": "{{reffile | fn}}.dict"}
	pMuTect2Interval.output = "outfile:file:{{tumor | fn}}-{{normal | fn}}.vcf"
	pMuTect2Interval.args   = { "gatk": "gatk", "picard": "picard", "samtools": "samtools", "javamem": "-Xms1g -Xmx8g", "tmpdir": "/tmp", "params": "" }
	pMuTect2Interval.script = """
	%s%s%s%s
	# generate interval file
	intvfile="{{job.outdir}}/interval.list"
	{{args.samtools}} idxstats "{{tumor}}" | head -n -1 | cut -f1 > $intvfile
	{{args.gatk}} -T MuTect2 {{args.javamem}} -Djava.io.tmpdir="{{args.tmpdir}}" -R "{{reffile}}" -I:tumor "{{tumor}}" -I:normal "{{normal}}" -o "{{outfile}}" -L $intvfile {{args.params}}
	""" % (_sNormalBai, _sTumorBai, _sReferenceFai, _sReferenceDict)
	return pMuTect2Interval

