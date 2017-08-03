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
	`trimmomatic`:    The trimmomatic executable, default: "trimmomatic"
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
	"outfile1:file:{{fqfile1 | fn | lambda x: __import__('re').sub(r'(\\.fq|\\.fastq)$', '', x)}}.clean.fq.gz",
	"outfile2:file:{{fqfile2 | fn | lambda x: __import__('re').sub(r'(\\.fq|\\.fastq)$', '', x)}}.clean.fq.gz"
]
pTrimmomaticPE.args    = {
	#"trimmomatic": "java -jar trimmomatic.jar"
	"trimmomatic":    "trimmomatic",
	"phred":  "phred33",
	"params": "ILLUMINACLIP:{adapter}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36",
	"nthread": 1
}
pTrimmomaticPE.script  = """
#determine default adapter
params="{{args.params}}"
trimmomatic="{{args.trimmomatic}}"
if [[ "$params" == *"{adapter}"* ]]; then
	exe=($trimmomatic)
	exe="${exe[${#exe[@]}-1]}"
	trimdir=$(readlink -f $(which $exe))
	trimdir=$(dirname $trimdir)
	r="{adapter}"
	t="$trimdir/adapters/TruSeq3-PE.fa"
	params=${params//$r/$t}
fi
unpaired1="{{job.outdir}}/{{fqfile1 | fn | lambda x: __import__('re').sub(r'(\\.fq|\\.fastq)$', '', x)}}.unpaired.fq.gz"
unpaired2="{{job.outdir}}/{{fqfile2 | fn | lambda x: __import__('re').sub(r'(\\.fq|\\.fastq)$', '', x)}}.unpaired.fq.gz"
logfile="{{job.outdir}}/{{fqfile1 | fn | lambda x: __import__('re').sub(r'(\\.fq|\\.fastq)$', '', x)}}-{{fqfile2 | fn | lambda x: __import__('re').sub(r'(\\.fq|\\.fastq)$', '', x)}}.log"
$trimmomatic PE -{{args.phred}} -threads {{args.nthread}} -trimlog $logfile "{{fqfile1}}" "{{fqfile2}}" "{{outfile1}}" $unpaired1 "{{outfile2}}" $unpaired2 $params
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
	`trimmomatic`:    The trimmomatic executable, default: "trimmomatic"
	`phred`:  "phred33" (default) or "phred64"
	`params`: Other params for trimmomatric, default: "ILLUMINACLIP:{adapter}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"
		- have to replace `{adapter}` with the path of the adapter file
	`nthread`: 1
@requires:
	[trimmomatic](http://www.usadellab.org/cms/index.php?page=trimmomatic)
"""
pTrimmomaticSE = proc ()
pTrimmomaticSE.input   = "fqfile:file"
pTrimmomaticSE.output  = "outfile:file:{{fqfile | fn | lambda x: __import__('re').sub(r'(\\.fq|\\.fastq)$', '', x)}}.clean.fq.gz"
pTrimmomaticSE.args    = {
	#"trimmomatic": "java -jar trimmomatic.jar"
	"trimmomatic":    "trimmomatic",
	"phred":  "phred33",
	"params": "ILLUMINACLIP:{adapter}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36",
	"nthread": 1
}
pTrimmomaticSE.script  = """
#determine default adapter
params="{{args.params}}"
trimmomatic="{{args.trimmomatic}}"
if [[ "$params" == *"{adapter}"* ]]; then
	exe=($trimmomatic)
	exe="${exe[${#exe[@]}-1]}"
	trimdir=$(readlink -f $(which $exe))
	trimdir=$(dirname $trimdir)
	r="{adapter}"
	t="$trimdir/adapters/TruSeq3-SE.fa"
	params=${params//$r/$t}
fi

logfile="{{job.outdir}}/{{fqfile | fn | lambda x: __import__('re').sub(r'(\\.fq|\\.fastq)$', '', x)}}.log"
$trimmomatic SE -{{args.phred}} -threads {{args.nthread}} -trimlog $logfile "{{fqfile}}" "{{outfile}}"  $params
"""

"""
@name:
	pAlignPEByBWA
@description:
	Align paired-end reads to reference genome using bwa mem
@input:
	`infile1:file`: read file 1 (fastq, or fastq gzipped)
	`infile2:file`: read file 2 (fastq, or fastq gzipped)
	`reffile:file`: The reference file
@output:
	`outfile:file`: The output sam file
@args:
	`bwa`:    The bwa executable, default: bwa
	`params`: Other params for bwa mem, default: "-M"
	`nthread`: 1
@requires:
	[bwa](https://github.com/lh3/bwa)
"""
pAlignPEByBWA = proc ()
pAlignPEByBWA.input   = "infile1:file, infile2:file, reffile:file"
pAlignPEByBWA.output  = "outfile:file:{{infile1 | bn | lambda x: __import__('re').sub(r'[^a-zA-Z0-9]*1(\\.clean)?(\\.fq|\\.fastq)(\\.gz)?$', '', x)}}.sam"
pAlignPEByBWA.args    = {
	"bwa":     "bwa",
	"params":  "", 
	"nthread": 1
}
pAlignPEByBWA.script  = """
lane="L{{#}}"
sample="{{infile1 | fn | .split('_')[0]}}"
rg="-R @RG\\tID:$lane\\tSM:$sample\\tPL:ILLUMINA"
params="{{args.params}}"
if [[ "$params" == *"@RG"* ]]; then
	rg=""
	params=${params/\{lane\}/$lane}
	params=${params/\{sample\}/$sample}
fi
{{args.bwa}} mem $rg {{args.params}} -t {{args.nthread}} "{{reffile}}" "{{infile1}}" "{{infile2}}" > "{{outfile}}"
"""

"""
@name:
	pAlignSEByBWA
@description:
	Align paired-end reads to reference genome using bwa mem
@input:
	`infile:file`:  read file (fastq, or fastq gzipped)
	`reffile:file`: The reference file
@brings:
	`reffile#bwt`: "{{reffile | bn}}.bwt", 
	`reffile#sa`:  "{{reffile | bn}}.sa",
	`reffile#ann`: "{{reffile | bn}}.ann",
	`reffile#amb`: "{{reffile | bn}}.amb",
	`reffile#pac`: "{{reffile | bn}}.pac"
@output:
	`outfile:file`: The output sam file
@args:
	`bwa`:    The bwa executable, default: bwa
	`params`: Other params for bwa mem, default: "-M"
	`nthread`: 1
	`reffile`: The reference file, required
@requires:
	[bwa](https://github.com/lh3/bwa)
"""
pAlignSEByBWA = proc ()
pAlignSEByBWA.input   = "infile:file, reffile:file"
pAlignSEByBWA.output  = "outfile:file:{{infile | fn | lambda x: __import__('re').sub(r'(\\.clean)?(\\.fq|\\.fastq)(\\.gz)?$', '', x)}}.sam"
pAlignSEByBWA.args    = {
	"bwa":     "bwa",
	"params":  "-M -R '@RG\\tID:bwa\\tSM:sample\\tPL:ILLUMINA'", # -R makes mapping very slow
	"nthread": 1
}
pAlignSEByBWA.script  = """
lane="L{{#}}"
sample="{{infile1 | fnnodot | .split('_')[0]}}"
rg="-R @RG\\tID:$lane\\tSM:$sample\\tPL:ILLUMINA"
params="{{args.params}}"
if [[ "$params" == *"@RG"* ]]; then
	rg=""
	params=${params/\{lane\}/$lane}
	params=${params/\{sample\}/$sample}
fi
{{args.bwa}} mem $rg {{args.params}} -t {{args.nthread}} "{{reffile}}" "{{infile}}" > "{{outfile}}"
"""

"""
@name:
	pAlignPEByNGM
@description:
	Align paired-end reads to reference genome using NextGenMap
@input:
	`infile1:file`: read file 1 (fastq, or fastq gzipped)
	`infile2:file`: read file 2 (fastq, or fastq gzipped)
	`reffile:file`: The reference file
@output:
	`outfile:file`: The output sam/bam file
@args:
	`ngm`:    The NextGenMap executable, default: ngm
	`nthread`: 1
	`outtype`: sam or bam, default: sam (only sam for now, due to bug of ngm 0.5.3 (fixed in 0.5.4))
	`params`: Other params for ngm, default: "--rg-id ngm --rg-sm sample"
@requires:
	[NextGenMap](https://github.com/Cibiv/NextGenMap/wiki)
"""
pAlignPEByNGM = proc ()
pAlignPEByNGM.input   = "infile1:file, infile2:file, reffile:file"
pAlignPEByNGM.output  = "outfile:file:{{infile1 | fn | lambda x: __import__('re').sub(r'[^a-zA-Z0-9]?1(\\.clean)?(\\.fq|\\.fastq)?$', '', x)}}.{{args.outtype}}"
pAlignPEByNGM.args    = {
	"ngm":     "ngm",
	"params":  "--rg-id ngm --rg-sm sample",
	"nthread": 1,
	"outtype": "sam"
}
pAlignPEByNGM.script  = """
bam=""
if [[ "{{args.outtype}}" == "bam" ]]; then
	bam="--bam"
fi
{{args.ngm}} {{args.params}} -r "{{reffile}}" -t {{args.nthread}} -1 "{{infile1}}" -2 "{{infile2}}" -o "{{outfile}}" $bam
"""

"""
@name:
	pAlignSEByNGM
@description:
	Align single-end reads to reference genome using NextGenMap
@input:
	`infile1:file`: read file 1 (fastq, or fastq gzipped)
	`infile2:file`: read file 2 (fastq, or fastq gzipped)
	`reffile:file`: The reference file
@output:
	`outfile:file`: The output sam/bam file
@args:
	`ngm`:    The NextGenMap executable, default: ngm
	`nthread`: 1
	`outtype`: sam or bam, default: sam (only sam for now, due to bug of ngm 0.5.3 (fixed in 0.5.4))
	`params`: Other params for ngm, default: "--rg-id ngm --rg-sm sample"
@requires:
	[NextGenMap](https://github.com/Cibiv/NextGenMap/wiki)
"""
pAlignSEByNGM = proc ()
pAlignSEByNGM.input   = "infile:file, reffile:file"
pAlignSEByNGM.output  = "outfile:file:{{infile | fn}}.{{args.outtype}}"
pAlignSEByNGM.args    = {
	"ngm":     "ngm",
	"params":  "--rg-id ngm --rg-sm sample",
	"nthread": 1,
	"outtype": "sam"
}
pAlignSEByNGM.script  = """
bam=""
if [[ "{{args.outtype}}" == "bam" ]]; then
	bam="--bam"
fi
{{args.ngm}} {{args.params}} -r "{{reffile}}" -t {{args.nthread}} -q "{{infile}}" -o "{{outfile}}" $bam
"""

"""
@name:
	pMergeBams
@description:
	Merge bam files
@input:
	`bamdir:dir`:   the dir containing bam files 
@output:
	`outfile:file`: the merged bam file
@args:
	`samtools`: the executable path of samtools, default: "samtools"
	`nthread`:      Number of BAM/CRAM compression threads
	`params`:       Other parameters for `samtools merge`, default: ""
@requires:
	[samtools](http://www.htslib.org/)
"""
pMergeBams = proc ()
pMergeBams.input   = "bamdir:dir"
pMergeBams.output  = "outfile:file:{{bamdir | fn}}.merged.bam"
pMergeBams.args    = {"samtools": "samtools", "nthread": 1, "params": ""}
pMergeBams.script  = """
{{args.samtools}} merge -@ {{args.nthread}} {{args.params}} {{outfile}} {{bamdir}}/*.bam
"""

