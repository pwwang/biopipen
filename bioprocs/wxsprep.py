from pyppl import proc

"""
Prepare WXS data, including alignment, QC,...
"""

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
pTrimmomaticPE.output  = [
	"outfile1:file:{{fqfile1.fn | (lambda x: __import__('re').sub(r'(\\.fq|\\.fastq)$', '', x))(_)}}.clean.fq.gz",
	"outfile2:file:{{fqfile2.fn | (lambda x: __import__('re').sub(r'(\\.fq|\\.fastq)$', '', x))(_)}}.clean.fq.gz"
]
pTrimmomaticPE.args    = {
	#"bin": "java -jar trimmomatic.jar"
	"bin":    "trimmomatic",
	"phred":  "phred33",
	"params": "ILLUMINACLIP:{adapter}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36",
	"nthread": 1
}
pTrimmomaticPE.script  = """
#determine default adapter
params="{{proc.args.params}}"
bin="{{proc.args.bin}}"
if [[ "$params" == *"{adapter}"* ]]; then
	exe=($bin)
	exe="${exe[${#exe[@]}-1]}"
	trimdir=$(readlink -f $(which $exe))
	trimdir=$(dirname $trimdir)
	r="{adapter}"
	t="$trimdir/adapters/TruSeq3-PE.fa"
	params=${params//$r/$t}
fi
unpaired1="{{proc.outdir}}/{{fqfile1.fn | (lambda x: __import__('re').sub(r'(\\.fq|\\.fastq)$', '', x))(_)}}.unpaired.fq.gz"
unpaired2="{{proc.outdir}}/{{fqfile2.fn | (lambda x: __import__('re').sub(r'(\\.fq|\\.fastq)$', '', x))(_)}}.unpaired.fq.gz"
logfile="{{proc.outdir}}/{{fqfile1.fn | (lambda x: __import__('re').sub(r'(\\.fq|\\.fastq)$', '', x))(_)}}-{{fqfile2.fn | (lambda x: __import__('re').sub(r'(\\.fq|\\.fastq)$', '', x))(_)}}.log"
$bin PE -{{proc.args.phred}} -threads {{proc.args.nthread}} -trimlog $logfile "{{fqfile1}}" "{{fqfile2}}" "{{outfile1}}" $unpaired1 "{{outfile2}}" $unpaired2 $params
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
pTrimmomaticSE.output  = "outfile:file:{{fqfile.fn | (lambda x: __import__('re').sub(r'(\\.fq|\\.fastq)$', '', x))(_)}}.clean.fq.gz"
pTrimmomaticSE.args    = {
	#"bin": "java -jar trimmomatic.jar"
	"bin":    "trimmomatic",
	"phred":  "phred33",
	"params": "ILLUMINACLIP:{adapter}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36",
	"nthread": 1
}
pTrimmomaticSE.script  = """
#determine default adapter
params="{{proc.args.params}}"
bin="{{proc.args.bin}}"
if [[ "$params" == *"{adapter}"* ]]; then
	exe=($bin)
	exe="${exe[${#exe[@]}-1]}"
	trimdir=$(readlink -f $(which $exe))
	trimdir=$(dirname $trimdir)
	r="{adapter}"
	t="$trimdir/adapters/TruSeq3-SE.fa"
	params=${params//$r/$t}
fi

logfile="{{proc.outdir}}/{{fqfile.fn | (lambda x: __import__('re').sub(r'(\\.fq|\\.fastq)$', '', x))(_)}}.log"
$bin SE -{{proc.args.phred}} -threads {{proc.args.nthread}} -trimlog $logfile "{{fqfile}}" "{{outfile}}"  $params
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
pAlignPEByBWA.output  = "outfile:file:{{infile1.fn | (lambda x: __import__('re').sub(r'[^\\w]?1(\\.clean)?(\\.fq|\\.fastq)?$', '', x))(_)}}.sam"
pAlignPEByBWA.args    = {
	#"bin": "java -jar trimmomatic.jar"
	"bin":     "bwa",
	"params":  "-M -R '@RG\\tID:bwa\\tSM:sample\\tPL:ILLUMINA'", # -R makes mapping very slow
	#"params":  "-M", # -R makes mapping very slow, use picard AddOrReplaceReadGroups
	"nthread": 1,
	"reffile": ""
}
pAlignPEByBWA.script  = """
if [[ -z "{{proc.args.reffile}}" ]]; then
	echo "Reference file not specified!" 1>&2
	exit 1
fi
if [[ ! -e "{{proc.args.reffile}}.bwt" ]]; then
	{{proc.args.bin}} index "{{proc.args.reffile}}"
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
pAlignSEByBWA.output  = "outfile:file:{{infile.fn | (lambda x: __import__('re').sub(r'(\\.clean)?(\\.fq|\\.fastq)?$', '', x))(_)}}.sam"
pAlignSEByBWA.args    = {
	#"bin": "java -jar trimmomatic.jar"
	"bin":     "bwa",
	"params":  "-M -R '@RG\\tID:bwa\\tSM:sample\\tPL:ILLUMINA'", # -R makes mapping very slow
	#"params":  "-M", # -R makes mapping very slow, use picard AddOrReplaceReadGroups
	"nthread": 1,
	"reffile": ""
}
pAlignSEByBWA.script  = """
if [[ -z "{{proc.args.reffile}}" ]]; then
	echo "Reference file not specified!" 1>&2
	exit 1
fi
if [[ ! -e "{{proc.args.reffile}}.bwt" ]]; then
	{{proc.args.bin}} index "{{proc.args.reffile}}"
fi
{{proc.args.bin}} mem {{proc.args.params}} -t {{proc.args.nthread}} "{{proc.args.reffile}}" "{{infile}}" > "{{outfile}}"
"""

"""
@name:
	pAlignPEByNGM
@description:
	Align paired-end reads to reference genome using NextGenMap
@input:
	`infile1:file`: read file 1 (fastq, or fastq gzipped)
	`infile2:file`: read file 2 (fastq, or fastq gzipped)
@output:
	`outfile:file`: The output sam/bam file
@args:
	`bin`:    The NextGenMap executable, default: ngm
	`params`: Other params for ngm, default: "--rg-id ngm --rg-sm sample"
	`nthread`: 1
	`reffile`: The reference file
	`outtype`: sam or bam, default: bam
@requires:
	[NextGenMap](https://github.com/Cibiv/NextGenMap/wiki)
"""
pAlignPEByNGM = proc ()
pAlignPEByNGM.input   = "infile1:file, infile2:file"
pAlignPEByNGM.output  = "outfile:file:{{infile1.fn | (lambda x: __import__('re').sub(r'[^\\w]?1(\\.clean)?(\\.fq|\\.fastq)?$', '', x))(_)}}.{{proc.args.outtype}}"
pAlignPEByNGM.args    = {
	"bin":     "ngm",
	"params":  "--rg-id ngm --rg-sm sample",
	"nthread": 1,
	"reffile": "",
	"outtype": "bam"
}
pAlignPEByNGM.script  = """
if [[ -z "{{proc.args.reffile}}" ]]; then
	echo "Reference file not specified!" 1>&2
	exit 1
fi
bam=""
if [[ "{{proc.args.outtype}}" == "bam" ]]; then
	bam="--bam"
fi
{{proc.args.bin}} {{proc.args.params}} -r "{{proc.args.reffile}}" -t {{proc.args.nthread}} -1 "{{infile1}}" -2 "{{infile2}}" -o "{{outfile}}" $bam
"""

"""
@name:
	pAlignSEByNGM
@description:
	Align single-end reads to reference genome using NextGenMap
@input:
	`infile:file`: read file (fastq, or fastq gzipped)
@output:
	`outfile:file`: The output sam/bam file
@args:
	`bin`:    The NextGenMap executable, default: ngm
	`params`: Other params for ngm, default: "--rg-id ngm --rg-sm sample"
	`nthread`: 1
	`reffile`: The reference file
	`outtype`: sam or bam, default: bam
@requires:
	[NextGenMap](https://github.com/Cibiv/NextGenMap/wiki)
"""
pAlignSEByNGM = proc ()
pAlignSEByNGM.input   = "infile:file"
pAlignSEByNGM.output  = "outfile:file:{{infile.fn}}.{{proc.args.outtype}}"
pAlignSEByNGM.args    = {
	"bin":     "ngm",
	"params":  "--rg-id ngm --rg-sm sample",
	"nthread": 1,
	"reffile": "",
	"outtype": "bam"
}
pAlignSEByNGM.script  = """
if [[ -z "{{proc.args.reffile}}" ]]; then
	echo "Reference file not specified!" 1>&2
	exit 1
fi
bam=""
if [[ "{{proc.args.outtype}}" == "bam" ]]; then
	bam="--bam"
fi
{{proc.args.bin}} {{proc.args.params}} -r "{{proc.args.reffile}}" -t {{proc.args.nthread}} -q "{{infile}}" -o "{{outfile}}" $bam
"""

"""
@name:
	pMergeBams
@description:
	Merge bam files
@input:
	`sname`:        the sample name
	`bams:files`:   the bam files to be merged
@output:
	`outfile:file`: the merged bam file
@args:
	`bin-samtools`: the executable path of samtools, default: "samtools"
	`nthread`:      Number of BAM/CRAM compression threads
	`params`:       Other parameters for `samtools merge`, default: ""
@requires:
	[samtools](http://www.htslib.org/)
"""
pMergeBams = proc ()
pMergeBams.input   = "sname, bams:files"
pMergeBams.output  = "outfile:file:{{sname}}.bam"
pMergeBams.args    = {"bind-samtools": "samtools", "nthread": 1, "params": ""}
pMergeBams.script  = """
{{proc.args.bin-samtools}} merge -@ {{proc.args.nthread}} {{proc.args.params}} {{outfile}} "{{bams | '" "'.join(_)}}"
"""

