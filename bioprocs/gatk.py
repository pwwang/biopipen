from pyppl import proc

#############################
# GATK utilities            #
#############################

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
	`infile:file`:  The aligned bam file 
@output:
	`outfile:file`: A list of target intervals to pass to the IndelRealigner.
@args:
	`bin`:     The gatk executable, default: "gatk -T RealignerTargetCreator"
	`params`:  Other parameters for RealignerTargetCreator, default: ""
	`reffile`: The reference file
@requires:
	[GATK](https://software.broadinstitute.org/gatk)
"""
pRealignerTargetCreator = proc()
pRealignerTargetCreator.input  = "infile:file"
pRealignerTargetCreator.output = "outfile:file:{{infile.fn}}.indelRealigner.intervals"
pRealignerTargetCreator.args   = { "bin": "gatk -T RealignerTargetCreator", "params": "", "reffile": "" }
pRealignerTargetCreator.script = """
if [[ -z "{{proc.args.reffile}}" ]]; then
	echo "Reference file not specified!" 1>&2
	exit 1
fi
{{proc.args.bin}} -R "{{proc.args.reffile}}" -I "{{infile}}" -o "{{outfile}}" {{params}}
"""

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
@output:
	`outfile:file`: A realigned version of input BAM file.
@args:
	`bin`:     The gatk executable, default: "gatk -T IndelRealigner"
	`params`:  Other parameters for IndelRealigner, default: ""
	`reffile`: The reference file
@requires:
	[GATK](https://software.broadinstitute.org/gatk)
"""
pIndelRealigner = proc()
pIndelRealigner.input  = "bamfile:file, intfile:file"
pIndelRealigner.output = "outfile:file:{{bamfile.fn}}.realigned.bam"
pIndelRealigner.args   = { "bin": "gatk -T IndelRealigner", "params": "", "reffile": "" }
pIndelRealigner.script = """
if [[ -z "{{proc.args.reffile}}" ]]; then
	echo "Reference file not specified!" 1>&2
	exit 1
fi
{{proc.args.bin}} -R "{{proc.args.reffile}}" -I "{{bamfile}}" -o "{{outfile}}" -targetIntervals "{{intfile}}" {{params}}
"""

"""
@name:
	pBaseRecalibrator  
@description:
	Variant calling algorithms rely heavily on the quality scores assigned to the individual base calls in each sequence read. These scores are per-base estimates of error emitted by the sequencing machines. Unfortunately the scores produced by the machines are subject to various sources of systematic technical error, leading to over- or under-estimated base quality scores in the data. Base quality score recalibration (BQSR) is a process in which we apply machine learning to model these errors empirically and adjust the quality scores accordingly. This allows us to get more accurate base qualities, which in turn improves the accuracy of our variant calls. The base recalibration process involves two key steps: first the program builds a model of covariation based on the data and a set of known variants (which you can bootstrap if there is none available for your organism), then it adjusts the base quality scores in the data based on the model. There is an optional but highly recommended step that involves building a second model and generating before/after plots to visualize the effects of the recalibration process. This is useful for quality control purposes. This tool performs the first step described above: it builds the model of covariation and produces the recalibration table. It operates only at sites that are not in dbSNP; we assume that all reference mismatches we see are therefore errors and indicative of poor base quality. This tool generates tables based on various user-specified covariates (such as read group, reported quality score, cycle, and context). Assuming we are working with a large amount of data, we can then calculate an empirical probability of error given the particular covariates seen at this site, where p(error) = num mismatches / num observations. The output file is a table (of the several covariate values, number of observations, number of mismatches, empirical quality score).
@input:
	`bamfile:file`: A BAM file containing data that needs to be recalibrated.
@output:
	`outfile:file`: A GATKReport file with many tables:
		- The list of arguments
		- The quantized qualities table
		- The recalibration table by read group
		- The recalibration table by quality score
		- The recalibration table for all the optional covariates
@args:
	`bin`:     The gatk executable, default: "gatk -T BaseRecalibrator"
	`params`:  Other parameters for BaseRecalibrator, default: ""
	`reffile`: The reference file
@requires:
	[GATK](https://software.broadinstitute.org/gatk)
"""
pBaseRecalibrator = proc()
pBaseRecalibrator.input  = "bamfile:file"
pBaseRecalibrator.output = "outfile:file:{{bamfile.fn}}.recal.table"
pBaseRecalibrator.args   = { "bin": "gatk -T BaseRecalibrator", "params": "", "reffile": "" }
pBaseRecalibrator.script = """
if [[ -z "{{proc.args.reffile}}" ]]; then
	echo "Reference file not specified!" 1>&2
	exit 1
fi
{{proc.args.bin}} -R "{{proc.args.reffile}}" -I "{{bamfile}}" -o "{{outfile}}" {{params}}
"""

"""
@name:
	pPrintReads   
@description:
	PrintReads is a generic utility tool for manipulating sequencing data in SAM/BAM format. It can dynamically merge the contents of multiple input BAM files, resulting in merged output sorted in coordinate order. It can also optionally filter reads based on various read properties such as read group tags using the `--read_filter/-rf` command line argument (see documentation on read filters for more information).
	Note that when PrintReads is used as part of the Base Quality Score Recalibration workflow, it takes the `--BQSR` engine argument, which is listed under Inherited Arguments > CommandLineGATK below.
@input:
	`bamfile:file`: A BAM file.
	`recaltable:file`: The GATKReport file
@output:
	`outfile:file`: A single processed bam file.
@args:
	`bin`:     The gatk executable, default: "gatk -T PrintReads"
	`params`:  Other parameters for PrintReads, default: ""
	`reffile`: The reference file
@requires:
	[GATK](https://software.broadinstitute.org/gatk)
"""
pPrintReads = proc()
pPrintReads.input  = "bamfile:file, recaltable:file"
pPrintReads.output = "outfile:file:{{bamfile.fn}}.recal.bam"
pPrintReads.args   = { "bin": "gatk -T PrintReads", "params": "", "reffile": "" }
pPrintReads.script = """
if [[ -z "{{proc.args.reffile}}" ]]; then
	echo "Reference file not specified!" 1>&2
	exit 1
fi
{{proc.args.bin}} -R "{{proc.args.reffile}}" -I "{{bamfile}}" -o "{{outfile}}" -BQSR "{{recaltable}}" {{params}}
"""

"""
@name:
	pHaplotypeCaller 
@description:
	PrintReads is a generic utility tool for manipulating sequencing data in SAM/BAM format. It can dynamically merge the contents of multiple input BAM files, resulting in merged output sorted in coordinate order. It can also optionally filter reads based on various read properties such as read group tags using the `--read_filter/-rf` command line argument (see documentation on read filters for more information).
	Note that when PrintReads is used as part of the Base Quality Score Recalibration workflow, it takes the `--BQSR` engine argument, which is listed under Inherited Arguments > CommandLineGATK below.
@input:
	`bamfile:file`: A BAM file.
@output:
	`outfile:file`: Either a VCF or gVCF file with raw, unfiltered SNP and indel calls.
@args:
	`bin`:     The gatk executable, default: "gatk -T HaplotypeCaller"
	`params`:  Other parameters for HaplotypeCaller, default: ""
	`reffile`: The reference file
@requires:
	[GATK](https://software.broadinstitute.org/gatk)
"""
pHaplotypeCaller = proc()
pHaplotypeCaller.input  = "bamfile:file"
pHaplotypeCaller.output = "outfile:file:{{bamfile.fn}}.vcf"
pHaplotypeCaller.args   = { "bin": "gatk -T HaplotypeCaller", "params": "", "reffile": "" }
pHaplotypeCaller.script = """
if [[ -z "{{proc.args.reffile}}" ]]; then
	echo "Reference file not specified!" 1>&2
	exit 1
fi
{{proc.args.bin}} -R "{{proc.args.reffile}}" -I "{{bamfile}}" -o "{{outfile}}" {{params}}
"""

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
@output:
	`outfile:file`: A new VCF file containing the selected subset of variants.
@args:
	`bin`:     The gatk executable, default: "gatk -T SelectVariants"
	`params`:  Other parameters for SelectVariants, default: ""
	`reffile`: The reference file
@requires:
	[GATK](https://software.broadinstitute.org/gatk)
"""
pSelectVariants = proc()
pSelectVariants.input  = "vcffile:file"
pSelectVariants.output = "outfile:file:{{vcffile.fn}}.selected.vcf"
pSelectVariants.args   = { "bin": "gatk -T SelectVariants", "params": "", "reffile": "" }
pSelectVariants.script = """
if [[ -z "{{proc.args.reffile}}" ]]; then
	echo "Reference file not specified!" 1>&2
	exit 1
fi
{{proc.args.bin}} -R "{{proc.args.reffile}}" -V "{{vcffile}}" -o "{{outfile}}" {{params}}
"""

"""
@name:
	pVariantFiltration
@description:
	This tool is designed for hard-filtering variant calls based on certain criteria. Records are hard-filtered by changing the value in the FILTER field to something other than PASS. Filtered records will be preserved in the output unless their removal is requested in the command line.
	The most common way of specifying filtering criteria is by using JEXL queries. See the article on JEXL expressions in the documentation Guide for detailed information and examples.
@input:
	`vcffile:file`: A variant set to filter.
@output:
	`outfile:file`: A filtered VCF.
@args:
	`bin`:     The gatk executable, default: "gatk -T VariantFiltration"
	`params`:  Other parameters for VariantFiltration, default: ""
	`reffile`: The reference file
@requires:
	[GATK](https://software.broadinstitute.org/gatk)
"""
pVariantFiltration = proc()
pVariantFiltration.input  = "vcffile:file"
pVariantFiltration.output = "outfile:file:{{vcffile.fn}}.filtered.vcf"
pVariantFiltration.args   = { "bin": "gatk -T VariantFiltration", "params": "", "reffile": "" }
pVariantFiltration.script = """
if [[ -z "{{proc.args.reffile}}" ]]; then
	echo "Reference file not specified!" 1>&2
	exit 1
fi
{{proc.args.bin}} -R "{{proc.args.reffile}}" -V "{{vcffile}}" -o "{{outfile}}" {{params}}
"""
