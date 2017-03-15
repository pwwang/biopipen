from pyppl import proc

"""
@name:
	pTrimmomaticPE
@description:
	Trimming Illumina NGS paired-end data
@input:
	`fqfile1:file`: The 1st fastq file (could be in .gz format)
	`fqfile2:file`: The 2nd fastq file
@output:
	`outfile1:file`: The 1st output file
	`outfile2:file`: The 2nd output file
@args:
	`bin`:    The trimmomatic executable, default: "trimmomatic"
	`phred`:  "phred33" (default) or "phred64"
	`params`: Other params for trimmomatric, default: "ILLUMINACLIP:{adapter}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"
		- have to replace `{adapter}` with the path of the adapter file
	`nthread`: 1
@requires:
	[trimmomatic](http://www.usadellab.org/cms/index.php?page=trimmomatic)
"""
pTrimmomaticPE = proc ()
pTrimmomaticPE.input   = "fqfile1:file, fqfile2:file"
pTrimmomaticPE.output  = "outfile1:file:{{fqfile1.fn}}.clean{{fqfile1.ext}}, outfile2:file:{{fqfile2.fn}}.clean{{fqfile2.ext}}"
pTrimmomaticPE.args    = {
	#"bin": "java -jar trimmomatic.jar"
	"bin":    "trimmomatic",
	"phred":  "phred33",
	"params": "ILLUMINACLIP:{adapter}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36",
	"nthread": 1
}
pTrimmomaticPE.script  = """
unpaired1="{{proc.outdir}}/{{fqfile1}}.unpaired{{fqfile1.ext}}"
unpaired2="{{proc.outdir}}/{{fqfile2}}.unpaired{{fqfile2.ext}}"
logfile="{{proc.outdir}}/{{fqfile1.fn}}-{{fqfile2.fn}}.log"
{{proc.args.bin}} PE -{{phred}} -threads {{proc.args.nthread}} -trimlog $logfile "{{fqfile1}}" "{{fqfile2}}" "{{outfile1}}" $unpaired1 "{{outfile2}}" $unpaired2 {{proc.args.params}}
"""

"""
@name:
	pTrimmomaticSE
@description:
	Trimming Illumina NGS single-end data
@input:
	`fqfile:file`: The fastq file (could be in .gz format)
@output:
	`outfile:file`: The output file
@args:
	`bin`:    The trimmomatic executable, default: "trimmomatic"
	`phred`:  "phred33" (default) or "phred64"
	`params`: Other params for trimmomatric, default: "ILLUMINACLIP:{adapter}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"
		- have to replace `{adapter}` with the path of the adapter file
	`nthread`: 1
@requires:
	[trimmomatic](http://www.usadellab.org/cms/index.php?page=trimmomatic)
"""
pTrimmomaticSE = proc ()
pTrimmomaticSE.input   = "fqfile:file"
pTrimmomaticSE.output  = "outfile:file:{{fqfile1.fn}}.clean{{fqfile1.ext}}"
pTrimmomaticSE.args    = {
	#"bin": "java -jar trimmomatic.jar"
	"bin":    "trimmomatic",
	"phred":  "phred33",
	"params": "ILLUMINACLIP:{adapter}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36",
	"nthread": 1
}
pTrimmomaticSE.script  = """
logfile="{{proc.outdir}}/{{fqfile1.fn}}-{{fqfile2.fn}}.log"
{{proc.args.bin}} SE -{{phred}} -threads {{proc.args.nthread}} -trimlog $logfile "{{fqfile}}" "{{outfile}}" {{proc.args.params}}
"""

"""
@name:
	pAlignPEByBWA
@description:
	Align paired-end reads to reference genome using bwa mem
@input:
	`infile1:file`: read file 1 (fastq, or fastq gzipped)
	`infile2:file`: read file 2 (fastq, or fastq gzipped)
@output:
	`outfile:file`: The output sam file
@args:
	`bin`:    The bwa executable, default: bwa
	`params`: Other params for bwa mem, default: "-M"
	`nthread`: 1
	`reffile`: The reference file
@requires:
	[bwa](https://github.com/lh3/bwa)
"""
pAlignPEByBWA = proc ()
pAlignPEByBWA.input   = "infile1:file, infile2:file"
pAlignPEByBWA.output  = "outfile:file:{{infile1.fn}}.sam"
pAlignPEByBWA.args    = {
	#"bin": "java -jar trimmomatic.jar"
	"bin":     "bwa",
	"params":  "-M",
	"nthread": 1,
	"reffile": ""
}
pAlignPEByBWA.script  = """
if [[ -z "{{proc.args.reffile}}" ]]; then
	echo "Reference file not specified!" 1>&2
	exit 1
fi
{{proc.args.bin}} mem {{proc.args.params}} -t {{proc.args.nthread}} "{{proc.args.reffile}}" "{{infile1}}" "{{infile2}}" > "{{outfile}}"
"""

"""
@name:
	pAlignSEByBWA
@description:
	Align paired-end reads to reference genome using bwa mem
@input:
	`infile:file`:  read file (fastq, or fastq gzipped)
@output:
	`outfile:file`: The output sam file
@args:
	`bin`:    The bwa executable, default: bwa
	`params`: Other params for bwa mem, default: "-M"
	`nthread`: 1
	`reffile`: The reference file, required
@requires:
	[bwa](https://github.com/lh3/bwa)
"""
pAlignSEByBWA = proc ()
pAlignSEByBWA.input   = "infile:file"
pAlignSEByBWA.output  = "outfile:file:{{infile.fn}}.sam"
pAlignSEByBWA.args    = {
	#"bin": "java -jar trimmomatic.jar"
	"bin":     "bwa",
	"params":  "-M",
	"nthread": 1,
	"reffile": ""
}
pAlignSEByBWA.script  = """
if [[ -z "{{proc.args.reffile}}" ]]; then
	echo "Reference file not specified!" 1>&2
	exit 1
fi
{{proc.args.bin}} mem {{proc.args.params}} -t {{proc.args.nthread}} "{{reffile}}" "{{infile}}" > "{{outfile}}"
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
	`bin`:    The picard executable, default: "picard SortSam"
	`order`:  The sort order, default: coordinate. Possible: unsorted, queryname, coordinate, duplicate
	`outtype`:The type of output file, sam or bam. Default: bam
@requires:
	[picard](http://broadinstitute.github.io/picard/command-line-overview.html)
"""
pSortSam = proc()
pSortSam.input  = "infile:file"
pSortSam.output = "outfile:file:{{infile.fn}}.sorted.{{proc.args.outtype}}"
pSortSam.args   = { "bin": "picard SortSam", "order": "coordinate", "outtype": "bam" }
pSortSam.script = """
{{proc.args.bin}} I="{{infile}}" O="{{outfile}}" SORT_ORDER={{proc.args.order}}
"""

"""
@name:
	pMarkDup
@description:
	Use `picard MarkDuplicates` to  mark duplicates for bam file
@input:
	`infile:file`:  The bam file 
@output:
	`outfile:file`: The marked bam file
@args:
	`bin`:    The picard executable, default: "picard MarkDuplicates"
	`params`:  Other parameters for picard MarkDuplicates, default: ""
@requires:
	[picard](http://broadinstitute.github.io/picard/command-line-overview.html)
"""
pMarkDup = proc()
pMarkDup.input  = "infile:file"
pMarkDup.output = "outfile:file:{{infile.fn}}.dedup.{{proc.args.outtype}}"
pMarkDup.args   = { "bin": "picard MarkDuplicates", "params": "" }
pMarkDup.script = """
metrics="{{proc.outdir}}/{{infile.fn}}.metrics.txt"
{{proc.args.bin}} I="{{infile}}" O="{{outfile}}" M="$metrics" {{proc.args.params}}
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
	`params`:  Other parameters for picard , default: ""
@requires:
	[picard](http://broadinstitute.github.io/picard/command-line-overview.html)
"""
pIndexBam = proc()
pIndexBam.input  = "infile:file"
pIndexBam.output = "outfile:file:{{infile.bn}}"
pIndexBam.args   = { "bin": "picard BuildBamIndex", "params": "" }
pIndexBam.script = """
ln -s "{{infile}}" "{{outfile}}"
{{proc.args.bin}} I="{{outfile}}" {{proc.args.params}}
# make sure .bai file is generated:
baifile="{{proc.outdir}}/{{infile.fn}}.bai"
if [[ ! -e $baifile ]]; then
	echo "Index file $baifile not generated!" 1>&2
	exit 1
fi
"""

"""
@name:
	pCNVnator
@description:
	Use `CNVnator` to call CNVs from bam file
@input:
	`infile:file`:  The bam file 
@output:
	`outfile:file`: The vcf file
@args:
	`bin`:      The CNVnator executable, default: "cnvnator"
	`bin-vcf`:  The converter executable to convert CNVnator results to vcf, default: "cnvnator2VCF.pl"
	`binsize`:  The bin_size, default: 100
	`genome`:   The genome: default: hg19
	`chrom`:    Chromosome names, default: "" (all chromosomes)
	`chrdir`:   The dir contains reference sequence of chromosomes, default: "" (don't specify)
	
@requires:
	[CNVnator](https://github.com/abyzovlab/CNVnator)
"""
pCNVnator = proc()
pCNVnator.input  = "infile:file"
pCNVnator.output = "outfile:file:{{infile.bn}}.vcf"
pCNVnator.args   = { "bin": "cnvnator", "bin-vcf": "cnvnator2VCF.pl", "binsize": 100, "genome": "hg19", "chrom": "", "chrdir": "" }
pCNVnator.script = """
rootfile="{{outfile}}.root"
[[ -z "{{proc.args.genome}}" ]] && genomeparam="" || genomeparam="-genome {{proc.args.genome}}"
[[ -z "{{proc.args.chrom}}" ]]  && chromparam=""  || chromparam="-chrom {{proc.args.chrom}}"
[[ -z "{{proc.args.chrdir}}" ]] && chrdirparam="" || chrdirparam="-d '{{proc.args.chrdir}}'"
# EXTRACTING READ MAPPING FROM BAM/SAM FILES
{{proc.args.bin}} $genomeparam -root "$rootfile" $chromparam -tree "{{infile}}"
# GENERATING HISTOGRAM
{{proc.args.bin}} $genomeparam -root "$rootfile" $chromparam -his {{proc.args.binsize}} $chrdirparam
# CALCULATING STATISTICS
{{proc.args.bin}} -root "$rootfile" $chromparam -stat {{proc.args.binsize}}
# RD SIGNAL PARTITIONING
{{proc.args.bin}} -root "$rootfile" $chromparam  -partition {{proc.args.binsize}}
# CNV CALLING
{{proc.args.bin}} -root "$rootfile" $chromparam -call {{proc.args.binsize}} > "{{outfile}}.cnvnator"
# Convert results to VCF:
{{proc.args.bin-vcf}} "{{outfile}}.cnvnator" > "{{outfile}}"
"""





