from pyppl import proc

#############################
# picard utilities          #
#############################

"""
@name:
	pMarkDuplicates
@description:
	Identifies duplicate reads.

	This tool locates and tags duplicate reads in a BAM or SAM file, where duplicate reads are defined as originating from a single fragment of DNA. Duplicates can arise during sample preparation e.g. library construction using PCR. See also EstimateLibraryComplexity for additional notes on PCR duplication artifacts. Duplicate reads can also result from a single amplification cluster, incorrectly detected as multiple clusters by the optical sensor of the sequencing instrument. These duplication artifacts are referred to as optical duplicates.

	The MarkDuplicates tool works by comparing sequences in the 5 prime positions of both reads and read-pairs in a SAM/BAM file. An BARCODE_TAG option is available to facilitate duplicate marking using molecular barcodes. After duplicate reads are collected, the tool differentiates the primary and duplicate reads using an algorithm that ranks reads by the sums of their base-quality scores (default method).

	The tool's main output is a new SAM or BAM file, in which duplicates have been identified in the SAM flags field for each read. Duplicates are marked with the hexadecimal value of 0x0400, which corresponds to a decimal value of 1024. If you are not familiar with this type of annotation, please see the following [blog post](https://www.broadinstitute.org/gatk/blog?id=7019) for additional information.

	Although the bitwise flag annotation indicates whether a read was marked as a duplicate, it does not identify the type of duplicate. To do this, a new tag called the duplicate type (DT) tag was recently added as an optional output in the 'optional field' section of a SAM/BAM file. Invoking the TAGGING_POLICY option, you can instruct the program to mark all the duplicates (All), only the optical duplicates (OpticalOnly), or no duplicates (DontTag). The records within the output of a SAM/BAM file will have values for the 'DT' tag (depending on the invoked TAGGING_POLICY), as either library/PCR-generated duplicates (LB), or sequencing-platform artifact duplicates (SQ). This tool uses the READ_NAME_REGEX and the OPTICAL_DUPLICATE_PIXEL_DISTANCE options as the primary methods to identify and differentiate duplicate types. Set READ_NAME_REGEX to null to skip optical duplicate detection, e.g. for RNA-seq or other data where duplicate sets are extremely large and estimating library complexity is not an aim. Note that without optical duplicate counts, library size estimation will be inaccurate.

	MarkDuplicates also produces a metrics file indicating the numbers of duplicates for both single- and paired-end reads.

	The program can take either coordinate-sorted or query-sorted inputs, however the behavior is slightly different. When the input is coordinate-sorted, unmapped mates of mapped records and supplementary/secondary alignments are not marked as duplicates. However, when the input is query-sorted (actually query-grouped), then unmapped mates and secondary/supplementary reads are not excluded from the duplication test and can be marked as duplicate reads.

	If desired, duplicates can be removed using the REMOVE_DUPLICATE and REMOVE_SEQUENCING_DUPLICATES options.
@input:
	`infile:file`:  The bam file 
@output:
	`outfile:file`: The marked bam file
@args:
	`bin`:     The picard executable, default: "picard MarkDuplicates"
	`params`:  Other parameters for picard MarkDuplicates, default: ""
	`tmpdir`:  The tmpdir to use. Default: /tmp
@requires:
	[picard](https://broadinstitute.github.io/picard/)
"""
pMarkDuplicates = proc()
pMarkDuplicates.input  = "infile:file"
pMarkDuplicates.output = "outfile:file:{{infile | fn | lambda x: __import__('re').sub(r'(\\.sort|\\.sorted)?$', '', x)}}.dedup.bam"
pMarkDuplicates.args   = { "bin": "picard MarkDuplicates", "params": "" }
pMarkDuplicates.script = """
tmpdir="{{proc.args.tmpdir}}/{{proc.id}}_{{#}}_{{infile | fn}}"
mkdir -p "$tmpdir"
mfile="{{job.outdir}}/{{infile | fn}}.metrics.txt"
{{proc.args.bin}} -Djava.io.tmpdir="$tmpdir" TMP_DIR="$tmpdir" I="{{infile}}" O="{{outfile}}" M="$mfile" {{proc.args.params}}
rm -rf "$tmpdir"
"""

"""
@name:
	pAddOrReplaceReadGroups
@description:
	Replace read groups in a BAM file.This tool enables the user to replace all read groups in the INPUT file with a single new read group and assign all reads to this read group in the OUTPUT BAM file.

	For more information about read groups, see the [GATK Dictionary entry](https://www.broadinstitute.org/gatk/guide/article?id=6472). 
	
	This tool accepts INPUT BAM and SAM files or URLs from the Global Alliance for Genomics and Health (GA4GH) (see http://ga4gh.org/#/documentation).
@input:
	`infile:file`:  The bam file
	`rg`:           The read group information. For example:
		- "RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=20"
@output:
	`outfile:file`: The bam file with read group added
@args:
	`bin`:     The picard executable, default: "picard AddOrReplaceReadGroups"
	`params`:  Other parameters for picard AddOrReplaceReadGroups, default: ""
@requires:
	[picard](https://broadinstitute.github.io/picard/)
"""
pAddOrReplaceReadGroups = proc()
pAddOrReplaceReadGroups.input  = "infile:file, rg"
pAddOrReplaceReadGroups.output = "outfile:file:{{infile | fn}}.rg.bam"
pAddOrReplaceReadGroups.args   = { "bin": "picard AddOrReplaceReadGroups", "params": "" }
pAddOrReplaceReadGroups.script = """
rg="{{rg}}"
if [[ "$rg" != *"RGPL="* ]]; then rg="$rg RGPL=illumina"; fi
if [[ "$rg" != *"RGPU="* ]]; then rg="$rg RGPU=unit1"; fi
if [[ "$rg" != *"RGLB="* ]]; then rg="$rg RGLB=lib1"; fi
{{proc.args.bin}} I="{{infile}}" O="{{outfile}}" $rg {{proc.args.params}}
"""

"""
@name:
	pCreateSequenceDictionary
@description:
	Creates a sequence dictionary for a reference sequence. This tool creates a sequence dictionary file (with ".dict" extension) from a reference sequence provided in FASTA format, which is required by many processing and analysis tools. The output file contains a header but no SAMRecords, and the header contains only sequence records.

	The reference sequence can be gzipped (both .fasta and .fasta.gz are supported).
@input:
	`infile:file`:  The fasta file 
@output:
	`outfile:file`: The same fasta file, but with dict file created
@args:
	`bin`:     The picard executable, default: "picard CreateSequenceDictionary"
	`params`:  Other parameters for picard CreateSequenceDictionary, default: ""
@requires:
	[picard](https://broadinstitute.github.io/picard/)
"""
pCreateSequenceDictionary = proc()
pCreateSequenceDictionary.input  = "infile:file"
pCreateSequenceDictionary.output = "outfile:file:{{infile | bn}}"
pCreateSequenceDictionary.args   = { "bin": "picard CreateSequenceDictionary", "params": "" }
pCreateSequenceDictionary.script = """
link="{{job.outdir}}/{{infile | bn}}"
ln -s "{{infile}}" "$link"
{{proc.args.bin}} R="$link" {{proc.args.params}}
"""

"""
@name:
	pCollectWgsMetrics
@description:
	Collect metrics about coverage and performance of whole genome sequencing (WGS) experiments.

	This tool collects metrics about the fractions of reads that pass base- and mapping-quality filters as well as coverage (read-depth) levels for WGS analyses. Both minimum base- and mapping-quality values as well as the maximum read depths (coverage cap) are user defined.
	
	Note: Metrics labeled as percentages are actually expressed as fractions!
@input:
	`infile:file`:  The bam file 
@output:
	`outfile:file`: The metrics file
@args:
	`bin`:     The picard executable, default: "picard CollectWgsMetrics"
	`params`:  Other parameters for `picard CollectWgsMetrics`, default: ""
	`reffile`: The reference file, default: ""
@requires:
	[picard](https://broadinstitute.github.io/picard/)
"""
pCollectWgsMetrics = proc()
pCollectWgsMetrics.input  = "infile:file"
pCollectWgsMetrics.output = "outfile:file:{{infile | bn}}.metrics.txt"
pCollectWgsMetrics.args   = { "bin": "picard CollectWgsMetrics", "params": "", "reffile": "" }
pCollectWgsMetrics.script = """
if [[ -z "{{proc.args.reffile}}" ]]; then
	echo "Reference file not specified!" 1>&2
	exit 1
fi
{{proc.args.bin}} I="{{infile}}" O="{{outfile}}" R="{{proc.args.reffile}}" {{proc.args.params}}
"""

"""
@name:
	pSortSam
@description:
	Use `picard SortSam` to sort sam or bam file
@input:
	`infile:file`:  The sam or bam file to be sorted
@output:
	`outfile:file`: The sorted sam or bam file
@args:
	`bin`:     The picard executable, default: "picard SortSam"
	`order`:   The sort order, default: coordinate. Possible: unsorted, queryname, coordinate, duplicate
	`outtype`: The type of output file, sam or bam. Default: bam
	`params`:  Other parameters for `picard SortSame`, default: "-Xms1g -Xmx8g"
	`tmpdir`:  The tmpdir to use. Default: /tmp
@requires:
	[picard](http://broadinstitute.github.io/picard/command-line-overview.html)
"""
pSortSam = proc()
pSortSam.input  = "infile:file"
pSortSam.output = "outfile:file:{{infile | fn}}.sorted.{{proc.args.outtype}}"
pSortSam.args   = { "bin": "picard SortSam", "order": "coordinate", "outtype": "bam", "params": "", "tmpdir": "/tmp" }
pSortSam.script = """
tmpdir="{{proc.args.tmpdir}}/{{proc.id}}_{{#}}_{{infile | fn}}"
mkdir -p "$tmpdir"
{{proc.args.bin}} -Djava.io.tmpdir="$tmpdir" TMP_DIR="$tmpdir" {{proc.args.params}} I="{{infile}}" O="{{outfile}}" SORT_ORDER={{proc.args.order}}
rm -rf "$tmpdir"
"""

"""
@name:
	pIndexBam
@description:
	Use `picard BuildBamIndex` to index bam file
@input:
	`infile:file`:  The bam file 
@output:
	`outfile:file`: The same bam file (link) but with .bai file in `proc.outdir`
@args:
	`bin`:    The picard executable, default: "picard BuildBamIndex"
	`params`:  Other parameters for picard , default: "-Xms1g -Xmx8g"
@requires:
	[picard](http://broadinstitute.github.io/picard/command-line-overview.html)
"""
pIndexBam = proc()
pIndexBam.input  = "infile:file"
pIndexBam.output = "outfile:file:{{infile | bn}}"
pIndexBam.args   = { "bin": "picard BuildBamIndex", "params": "-Xms1g -Xmx8g" }
pIndexBam.script = """
ln -s "{{infile}}" "{{outfile}}"
baifile="{{job.outdir}}/{{infile | fn}}.bai"
{{proc.args.bin}} I="{{outfile}}" O="$baifile" {{proc.args.params}}
# make sure .bai file is generated:
if [[ ! -e $baifile ]]; then
	echo "Index file $baifile not generated!" 1>&2
	exit 1
fi
"""
