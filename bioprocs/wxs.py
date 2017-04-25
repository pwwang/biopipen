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
pCNVnator.output = "outfile:file:{{infile.fn}}.cnv.vcf"
pCNVnator.args   = { "bin": "cnvnator", "bin-vcf": "cnvnator2VCF.pl", "binsize": 100, "genome": "hg19", "chrom": "", "chrdir": "" }
pCNVnator.script = """
rootfile="{{outfile}}.root"
[[ -z "{{proc.args.genome}}" ]] && genomeparam="" || genomeparam="-genome {{proc.args.genome}}"
[[ -z "{{proc.args.chrom}}" ]]  && chromparam=""  || chromparam="-chrom {{proc.args.chrom}}"
[[ -z "{{proc.args.chrdir}}" ]] && chrdirparam="" || chrdirparam="-d {{proc.args.chrdir}}"
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

"""
@name:
	pVcf2Stat
@description:
	Convert vcf to stat files for pVcfStats
@input:
	`vcffile:file`: The vcf file
@output:
	`outfile:file`: The stat file
@args:
	`chroms`: SNPs on chromosomes to consider, default: "*" (all chroms)
	- use "chr1-22, chrX, chrY" for chr1 to chr22, chrX and chrY
@requires:
	[`pyvcf`](https://github.com/jamescasbon/PyVCF)
"""
pVcf2Stat = proc()
pVcf2Stat.input   = "vcffile:file"
pVcf2Stat.output  = "outfile:file:{{vcffile.fn}}.stat"
pVcf2Stat.lang    = "python"
pVcf2Stat.args    = {"chroms": "*"}
pVcf2Stat.script  = """
import vcf
# parse chroms:
chroms = "{{proc.args.chroms}}"
mychrs = []
if chroms != "*":
	chroms = chroms.replace(" ", "").split(",")
	for chrom in chroms:
		if "-" not in chrom:
			mychrs.append (chrom)
			continue
		parts = chrom.split("-")
		i = len (parts[0]) - 1
		while parts[0][i].isdigit() and i >= 0: i -= 1
		i += 1
		prefix = parts[0][:i]
		num1   = int (parts[0][i:])
		num2   = int (parts[1][i:]) if parts[1].startswith (prefix) else int(parts[1])
		for i in range(num1, num2+1):
			mychrs.append (prefix + str(i))

def isMySnp (chr):
	if not mychrs: return True
	return chr in mychrs
	
vcf_reader = vcf.Reader (open("{{vcffile}}"))
ti    = tv   = 0
hetro = homo = 0
snps  = []
for record in vcf_reader:
	chr = record.CHROM
	if not isMySnp: continue
	
	gt  = record.samples[0]['GT']
	if gt == "0/1" or gt == "0|1": hetro += 1
	else: homo += 1
	
	alt = ''.join (map(str, record.ALT))
	r   = record.REF + alt
	if (r=="AG" or r=="GA" or r == "CT" or r == "TC") :
		ti += 1
	else: tv += 1
	
	snps.append (chr + ':' + str(record.POS) + '_' + record.REF + '/' + alt)

with open ("{{outfile}}", "w") as fout:
	fout.write ("Heterozygosity\\t%.3f\\n" % (float(hetro) / (hetro + homo)))
	fout.write ("TiTv\\t%.3f\\n" % (float(ti)/tv))
	for snp in snps:
		fout.write ("%s\\t1\\n" % snp)

"""

"""
@name:
	pVcfStat2Mat
@description:
	Convert vcf stat files to matrix for clustering, rownames are sample names
	First 2 lines TiTv and Heterozygosity will be ignored.
@input:
	`statdir:file`: The input directory containing vcf stat files
@output:
	`outfile:file`: The matrix file
"""
pVcfStat2Mat = proc()
pVcfStat2Mat.input   = "statdir:file"
pVcfStat2Mat.output  = "outfile:file:{{statdir.fn}}.mat"
pVcfStat2Mat.lang    = "Rscript"
pVcfStat2Mat.script  = """
setwd ("{{statdir}}")
statfiles = list.files()
samples   = noquote (unlist(lapply(statfiles, tools::file_path_sans_ext)))
nsample   = length (samples)

mat       = matrix (nrow = nsample)
rownames(mat) = samples
for (statfile in statfiles) {
	print (paste("Reading", statfile, "...", sep=" "))
	sample   = tools::file_path_sans_ext (statfile)
	data     = read.table (statfile, sep = "\\t", row.names=1, check.names=F, strip.white=T, stringsAsFactors=F, skip=2)
	
	snps     = rownames(data)
	matsnps  = colnames(mat)
	ovsnps   = intersect (snps, matsnps)
	diffsnps = setdiff (snps, matsnps)
	rmat     = matrix (0, nrow=1, ncol=ncol(mat))
	rownames (rmat) = c(sample)
	colnames (rmat) = colnames(mat)
	rmat [sample, ovsnps] = 1
	mat      = rbind (mat, rmat)
	cmat     = matrix (0, nrow=nrow(mat), ncol=length(diffsnps))
	colnames(cmat) = diffsnps
	rownames(cmat) = rownames(mat)
	cmat[sample, ] = 1
	mat      = cbind (mat, cmat)
}

write.table (mat, "{{outfile}}", quote=F, sep="\\t")
"""

"""
@name:
	pVcfStats
@description:
	Calculate sample/snp call rate and heterozygosity, TiTv ratio from single-sample vcfs using result from pVcf2Stat
@input:
	`statdir:file`:   The directory containing vcf stat files
@output:
	`outsample:file`: The report of call rate for each sample
	`figsample:file`: The bar chat of sample call rates
	`outsnp:file`:    The report of call rate for each snp
	`figsnp:file`:    The bar chat of snp call rates
	`outhetero:file`: The report of heterozygosity of each sample
	`fighetero:file`: The bar chat of snp call rates
	`outtitv:file`:   The report of transition/transversion ratio of each sample
	`figtitv:file`:   The bar chat of transition/transversion ratio
"""
pVcfStats = proc()
pVcfStats.input     = "statdir:file"
pVcfStats.output    = [
	"outsample:file:{{statdir.fn}}.sampleCallRate.txt",
	"figsample:file:{{statdir.fn}}.sampleCallRate.png",
	"outsnp:file:{{statdir.fn}}.snpCallRate.txt",
	"figsnp:file:{{statdir.fn}}.snpCallRate.png",
	"outhetero:file:{{statdir.fn}}.heterozygosity.txt",
	"fighetero:file:{{statdir.fn}}.heterozygosity.png",
	"outtitv:file:{{statdir.fn}}.titv.txt",
	"figtitv:file:{{statdir.fn}}.titv.png"
]
pVcfStats.lang      = "Rscript"
pVcfStats.script    = """
library ('plyr')
setwd ("{{statdir}}")
statfiles = list.files()
samples   = noquote (unlist(lapply(statfiles, tools::file_path_sans_ext)))
nsample   = length (samples)

hetreo    = matrix (nrow = nsample, ncol = 1)
titv      = matrix (nrow = nsample, ncol = 1)
samplecr  = matrix (nrow = nsample, ncol = 1) # snpcalled
snpcr     = {} # samplescalled
rownames(hetreo)   = samples
rownames(titv)     = samples
rownames(samplecr) = samples
for (statfile in statfiles) {
	print (paste("Reading", statfile, "...", sep=" "))
	sample   = tools::file_path_sans_ext (statfile)
	data     = read.table (statfile, sep = "\\t", row.names=1, check.names=F, strip.white=T, stringsAsFactors=F)
	hetreo [sample, 1]  = data ["Heterozygosity", 1]
	titv   [sample, 1]  = data ["TiTv", 1]
	data     = data[!rownames(data) %in% c("Heterozygosity", "TiTv"),, drop=F]
	samplecr[sample, 1] = nrow(data)
	snps     = rownames(data)
	crsnps   = rownames(snpcr)
	ovsnps   = intersect (snps, crsnps)
	diffsnps = setdiff (snps, crsnps)
	snpcr[ovsnps, ]  = snpcr[ovsnps,] + 1
	snpmat   = matrix (1, ncol=1, nrow=length(diffsnps))
	rownames(snpmat) = diffsnps
	snpcr    = rbind (snpcr, snpmat)
}
totalsnp = nrow(snpcr)
samplecr = samplecr / totalsnp
snpcr    = snpcr / nsample

# plot frequency
plotFreq = function (obj, figure, xlab, ylab="Frequency") {
	png (file=figure)
	h = hist (obj, freq=T, xlab=xlab, ylab=ylab, col="gray", main=paste(xlab, "distribution", sep=" "), axes=F)
	minb = min(h$breaks)
	maxb = max(h$breaks)
	maxc = max(h$counts)
	lenb = length(h$breaks)
	stpb = (maxb-minb)/(lenb-1)
	axis(1, pos=0, labels=T, at=seq(minb,maxb,stpb))
	lab0 = floor(log10(maxc))
	stpc = ceiling(maxc/(10**lab0)) * (10 ** (lab0-1))
	axis(2, pos=minb, labels=T, at=seq(0, maxc, stpc))
	dev.off()
}

# heterozygosity
write.table (hetreo, "{{outhetero}}", quote=F, sep="\\t", col.names=F)
plotFreq (hetreo, "{{fighetero}}", "Heterozygosity", "# samples")

# titv
write.table (titv, "{{outtitv}}", quote=F, sep="\\t", col.names=F)
plotFreq (titv, "{{figtitv}}", "TiTv ratio", "# samples")

# sample/snp call rate
write.table (samplecr, "{{outsample}}", quote=F, sep="\\t", col.names=F)
plotFreq (samplecr, "{{figsample}}", xlab="Sample call rate", "# samples")
write.table (snpcr, "{{outsnp}}", quote=F, sep="\\t", col.names=F)
plotFreq (snpcr, "{{figsnp}}", xlab="SNP call rate", "# snps")
"""

"""
@name:
	pCoverageByBamstats
@description:
	Use `bamstats` to calculate coverage for bam file
@input:
	`infile:file`:  The bam file
@output:
	`outfile:file`:    The report of coverage for the bam file
@args:
	`bin`: The `bamstats` executable, default: "bamstats"
	`params`: Other parameters for `bamstats`, default: ""
@requires:
	[bamstats](http://bamstats.sourceforge.net/)
"""
pCoverageByBamstats = proc()
pCoverageByBamstats.input     = "infile:file"
pCoverageByBamstats.output    = "outfile:file:{{infile.fn}}.bamstats.txt"
pCoverageByBamstats.args      = { "bin": "bamstats", "params": "" }
pCoverageByBamstats.script    = """
{{proc.args.bin}} -i "{{infile}}" -o "{{outfile}}" {{proc.args.params}}
"""

"""
@name:
	pPlotBamstats
@description:
	Plot coverage use output files generated by `bamstats` or `wxs.pCoverageByBamstats`
@input:
	`indir:file`: The directory containing bamstats output files
@args:
	`chroms`: Chromosomes to plot. Default: "" (all chroms)
	- Note: Whether to have "chr" prefix or not depends on your reference when mapping.
	- You can do a scope assignment: "chr1-chr22, chrX, chrY"
@output:
	`outdir:file`: The directory containing output figures
"""
pPlotBamstats = proc()
pPlotBamstats.input     = "indir:file"
pPlotBamstats.output    = "outdir:file:{{indir.bn}}-covplots"
pPlotBamstats.args      = {"chroms": ""}
pPlotBamstats.script    = """
#!/usr/bin/env Rscript
dir.create("{{outdir}}", showWarnings = F, recursive = T)
setwd ("{{indir}}")
chroms  = "{{proc.args.chroms}}"
# remove all spaces
chroms  = gsub("[[:blank:]]", "", chroms)
chroms  = noquote(unlist(strsplit(chroms, ",")))
chrs2   = vector(mode="character")
for (chrom in chroms) {
	chrparts = noquote(unlist(strsplit(chrom, "-")))
	if (length(chrparts) < 2) {
		chrs2 = c (chrs2, chrom)
		next
	}
	chr1   = gsub("[^0-9]", "", chrparts[1])
	chr2   = gsub("[^0-9]", "", chrparts[2])
	pos    = unlist(gregexpr (chr1, chrparts[1]))
	pos    = pos[length(pos)]
	prefix = substr(chrparts[1], 1, pos-1)
	for (chr in chr1:chr2) {
		chrs2 = c (chrs2, paste(prefix, chr, sep=""))
	}
}

bsfiles = list.files()
means   = matrix(ncol=1, nrow=length(bsfiles))
chrs    = vector()
rownames(means) = bsfiles
colnames(means) = c("Average coverage")
for (bsfile in bsfiles) {
	print (paste("Reading", bsfile, "...", sep=" "))
	sample = basename (bsfile)
	sample = tools::file_path_sans_ext(sample)
	stat   = read.table (bsfile, sep="", header=T, check.names=F, row.names=1)
	stat   = stat[chrs2, ]
	stat[, "N"] = as.numeric(gsub(",", "", stat[, "N"]))
	means[bsfile, 1] = sum(stat[, "N"] * stat[, "mean"])/sum(stat[, "N"])
	col2in = stat[, "mean", drop=F]
	colnames(col2in) = sample
	if (is.vector(chrs)) {
		chrs = col2in
	} else {
		chrs = cbind(chrs, col2in)
	}
}
# plot average coverage
plotFreq = function (obj, figure, xlab, ylab="Frequency") {
	png (file=figure)
	h = hist (obj, freq=T, xlab=xlab, ylab=ylab, col="gray", main=paste(xlab, "distribution", sep=" "), axes=F)
	minb = min(h$breaks)
	maxb = max(h$breaks)
	maxc = max(h$counts)
	lenb = length(h$breaks)
	stpb = (maxb-minb)/(lenb-1)
	axis(1, pos=0, labels=T, at=seq(minb,maxb,stpb))
	lab0 = floor(log10(maxc))
	stpc = ceiling(maxc/(10**lab0)) * (10 ** (lab0-1))
	axis(2, pos=minb, labels=T, at=seq(0, maxc, stpc))
	dev.off()
}
print ("Plotting average coverages ...")
write.table (means, "{{outdir}}/avgCoverage.txt", quote=F, sep="\\t")
plotFreq (means, "{{outdir}}/avgCoverage.png", xlab="Average coverage")

# plot chromosomes
print ("Plotting chromosome coverages ...")
png ("{{outdir}}/chrCoverage.png")
#colnames(chrs) = NULL # just show index
boxplot(t(chrs), ylab="Coverage")
dev.off()
"""


"""
@name:
	pMutSig
@description:
	MutSig stands for "Mutation Significance".  MutSig analyzes lists of mutations discovered in DNA sequencing, to identify genes that were mutated more often than expected by chance given background mutation processes.

	For more information, see Lawrence, M. et al. Mutational heterogeneity in cancer and the search for new cancer-associated genes. Nature 499, 214-218 (2013).
	
	See [dcumentation](http://archive.broadinstitute.org/cancer/cga/mutsig_run)
@input:
	`maffile:file`: mutation table
	`cvgfile:file`: coverage table
	`cvrfile:file`: covariates table
	`mutdict:file`: mutation_type_dictionary_file
	`chrdir:file`:  chr_files_hg18 or chr_files_hg19 
@output:
	`outfile:file`: The output file
@args:
	`bin`: The path to `run_MutSigCV.sh`, default: 'mutsig'
	`mcr`: The Matlab MCR path
@requires:
	[MutSing](http://archive.broadinstitute.org/cancer/cga/mutsig_download)
"""
pMutSig = proc()
pMutSig.input     = "maffile:file, cvgfile:file, cvrfile:file, mutdict:file, chrdir:file"
pMutSig.output    = "outfile:file:{{maffile.fn}}.mutsig.txt"
pMutSig.args      = {"bin": "mutsig", "mcr": ""}
pMutSig.script    = """
{{proc.args.bin}} "{{proc.args.mcr}}" "{{maffile}}" "{{cvgfile}}" "{{cvrfile}}" "{{outfile}}" "{{mutdict}}" "{{chrdir}}"
"""

"""
@name:
	pSnpEff2Maf
@description:
	Convert a snpEff-annotated somatic mutation vcf file (with normal and tumor samples) to [maf](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+(MAF)+Specification) file
@input:
	`vcffile:file`: snpEff-annotated vcf file
@output:
	`outfile:file`: The maf file
@args:
	`params`: A dict to specify a constant for columns. For example: `proc.args.params['NCBI_Build'] = 'hg19'` will set column 'NCBI_Build' as 'hg19' for all records. Default: {}
	- Keys could be one of these: ['Hugo_Symbol', 'Entrez_Gene_Id', 'Center', 'NCBI_Build', 'Chromosome', 'Start_Position', 'End_Position', 'Strand', 'Variant_Classification', 'Variant_Type', 'Reference_Allele', 'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2', 'dbSNP_RS', 'dbSNP_Val_Status', 'Tumor_Sample_Barcode', 'Matched_Norm_Sample_Barcode', 'Match_Norm_Seq_Allele1', 'Match_Norm_Seq_Allele2', 'Tumor_Validation_Allele1', 'Tumor_Validation_Allele2', 'Match_Norm_Validation_Allele1', 'Match_Norm_Validation_Allele2', 'Verification_Status', 'Validation_Status', 'Mutation_Status', 'Sequencing_Phase', 'Sequence_Source', 'Validation_Method', 'Score', 'BAM_File', 'Tumor_Sample_UUID', 'Matched_Norm_Sample_UUID']
	`normal`: The normal sample position in `record.samples` from `pyvcf`, 0-based. Default: 0
	`tumor`:  The tumor sample position in `record.samples` from `pyvcf`, 0-based. Default: 1
@requires:
	[pyvcf](https://github.com/jamescasbon/PyVCF)
"""
pSnpEff2Maf = proc()
pSnpEff2Maf.input     = "vcffile:file"
pSnpEff2Maf.output    = "outfile:file:{{vcffile.fn}}.maf"
pSnpEff2Maf.args      = {"params": {}, "normal": 0, "tumor": 1}
pSnpEff2Maf.script    = """
#!/usr/bin/env python
import vcf, shlex, mygene, json, re, os
mg = mygene.MyGeneInfo()
#Variant_Classification: Frame_Shift_Del, Frame_Shift_Ins, In_Frame_Del, In_Frame_Ins, Missense_Mutation, Nonsense_Mutation, Silent, Splice_Site, Translation_Start_Site, Nonstop_Mutation, 3'UTR, 3'Flank, 5'UTR, 5'Flank, IGR1 , Intron, RNA, Targeted_Region
vclass = {
	"frameshift_variant": ["Frame_Shift_Del", "Frame_Shift_Ins"],
	"inframe_deletion": "In_Frame_Del",
	"inframe_insertion": "In_Frame_Ins",
	"missense_variant": "Missense_Mutation",
	"initiator_codon_variant": "Missense_Mutation",
	"stop_retained_variant": "Missense_Mutation",
	"rare_amino_acid_variant": "Missense_Mutation",
	"stop_gained": "Nonsense_Mutation",
	"synonymous_variant": "Silent",
	"splice_acceptor_variant": "Splice_Site",
	"splice_donor_variant": "Splice_Site",
	"start_lost": "Translation_Start_Site",
	"start_retained": "Translation_Start_Site",
	"stop_lost": "Nonstop_Mutation",
	"3_prime_UTR_variant": "3'UTR",
	"3_prime_UTR_truncation": "3'UTR",
	"5_prime_UTR_variant": "5'UTR",
	"5_prime_UTR_truncation": "5'UTR",
	"intergenic_region": "IGR",
	"conserved_intergenic_variant": "IGR",
	"intron_variant": "Intron",
	"conserved_intron_variant": "Intron",
	"transcript_variant": "RNA",
	"miRNA": "RNA",
}

header = ['Hugo_Symbol', 'Entrez_Gene_Id', 'Center', 'NCBI_Build', 'Chromosome', 'Start_Position', 'End_Position', 'Strand', 'Variant_Classification', 'Variant_Type', 'Reference_Allele', 'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2', 'dbSNP_RS', 'dbSNP_Val_Status', 'Tumor_Sample_Barcode', 'Matched_Norm_Sample_Barcode', 'Match_Norm_Seq_Allele1', 'Match_Norm_Seq_Allele2', 'Tumor_Validation_Allele1', 'Tumor_Validation_Allele2', 'Match_Norm_Validation_Allele1', 'Match_Norm_Validation_Allele2', 'Verification_Status', 'Validation_Status', 'Mutation_Status', 'Sequencing_Phase', 'Sequence_Source', 'Validation_Method', 'Score', 'BAM_File', 'Tumor_Sample_UUID', 'Matched_Norm_Sample_UUID']

params = json.loads("{{proc.args.params | __import__('json').dumps(_)}}")
def getValue (key, value):
	return params[key] if params.has_key(key) else value

fout = open("{{outfile}}", "w")
fout.write ("\\t".join(header) + "\\n")
vr = vcf.Reader(open('{{vcffile}}'))
# Get NCBI_Build
NCBI_Build = shlex.split(vr.metadata['SnpEffCmd'][0][1:-1])[-2]
NCBI_Build = getValue ('NCBI_Build', NCBI_Build)
Center     = getValue ('Center', '.') # fake ones

for r in vr:
	# filter non coding region variants
	if not r.INFO.has_key('ANN') or len(r.INFO['ANN']) < 1:
		continue
	for ann in r.INFO['ANN']:
		tmp            = ann.split('|')
		annotations    = tmp[1].split('&')
		found          = False
		for annotation in annotations:
			if vclass.has_key(annotation):
				found                      = True
				Hugo_Symbol                = tmp[3]
				Hugo_Symbol                = getValue('Hugo_Symbol', Hugo_Symbol)
				gene                       = tmp[3]
				rgene                      = mg.query(gene, scope="symbol, alias", fields='entrezgene', size=1, species=NCBI_Build)
				if not rgene.has_key('hits') or len(rgene['hits']) < 1 or not rgene['hits'][0].has_key('entrezgene'): continue
				Entrez_Gene_Id             = rgene['hits'][0]['entrezgene']
				Entrez_Gene_Id             = getValue('Entrez_Gene_Id', Entrez_Gene_Id)
				Chromosome                 = getValue('Chromosome', r.CHROM)
				Start_Position             = getValue('Start_Position', r.POS)
				End_Position               = getValue('End_Position', r.POS)
				Strand                     = getValue('Strand', '+')
				Variant_Classification     = vclass[annotation]
				ref                        = r.REF
				alt                        = ','.join(map(str, r.ALT))
				reflen                     = len(ref)
				altlen                     = max(map(len, r.ALT))
				if annotation == 'frameshift_variant':
					if reflen >= altlen:
						Variant_Classification = getValue ('Variant_Classification', Variant_Classification[0])
					else:
						Variant_Classification = getValue ('Variant_Classification', Variant_Classification[1])
				vartype                    = 'SNP'
				if reflen == altlen:
					if ref == 2: vartype   = 'DNP'
					elif ref == 3: vartype = 'TNP'
					elif ref > 3: vartype  = 'ONP'
				elif reflen > altlen:
					vartype                = 'DEL'
				else:
					vartype                = 'INS'
				Variant_Type               = getValue('Variant_Type', vartype)
				Reference_Allele           = getValue('Reference_Allele', ref)
				
				normal                     = r.samples[{{proc.args.normal}}]
				tumor                      = r.samples[{{proc.args.tumor}}]
				tumor_genotype             = tumor['GT']
				tumor_genotypes            = re.split(r'[|\/]', tumor_genotype)
				if tumor_genotypes[0] == '0':
					Tumor_Seq_Allele1      = getValue ('Tumor_Seq_Allele1', Reference_Allele)
				else:
					Tumor_Seq_Allele1      = getValue ('Tumor_Seq_Allele1', alt)
				if tumor_genotypes[1] == '0':
					Tumor_Seq_Allele2      = getValue ('Tumor_Seq_Allele2', Reference_Allele)
				else:
					Tumor_Seq_Allele2      = getValue ('Tumor_Seq_Allele2', alt)
				
				if r.ID and r.ID.startswith('rs'):
					dbSNP_RS               = r.ID
				else:
					dbSNP_RS               = 'novel'
				dbSNP_RS                   = getValue ('dbSNP_RS', dbSNP_RS)
				
				dbSNP_Val_Status           = getValue ('dbSNP_Val_Status', '')
				Tumor_Sample_Barcode       = getValue ('Tumor_Sample_Barcode', tumor.sample)
				Matched_Norm_Sample_Barcode= getValue ('Matched_Norm_Sample_Barcode' ,normal.sample)
				
				normal_genotype            = normal['GT']
				normal_genotypes           = re.split(r'[|\/]', normal_genotype)
				if normal_genotypes[0] == '0':
					Match_Norm_Seq_Allele1 = getValue ('Match_Norm_Seq_Allele1', Reference_Allele)
				else:
					Match_Norm_Seq_Allele1 = getValue ('Match_Norm_Seq_Allele1', alt)
				if normal_genotypes[1] == '0':
					Match_Norm_Seq_Allele2 = getValue ('Match_Norm_Seq_Allele2', Reference_Allele)
				else:
					Match_Norm_Seq_Allele2 = getValue ('Match_Norm_Seq_Allele2', alt)
					
				Tumor_Validation_Allele1   = getValue ('Tumor_Validation_Allele1', '')
				Tumor_Validation_Allele2   = getValue ('Tumor_Validation_Allele2', '')
				Match_Norm_Validation_Allele1   = getValue ('Match_Norm_Validation_Allele1', '')
				Match_Norm_Validation_Allele2   = getValue ('Match_Norm_Validation_Allele2', '')
				Verification_Status        = getValue ('Verification_Status', '')
				Validation_Status          = getValue ('Validation_Status', 'Untested')				
				Mutation_Status            = getValue ('Mutation_Status', 'Unknown')
				Sequencing_Phase           = getValue ('Sequencing_Phase', '')
				Sequence_Source            = getValue ('Sequence_Source', 'WGS')
				Validation_Method          = getValue ('Validation_Method', 'none')
				Score                      = getValue ('Score', '')
				BAM_File                   = getValue ('BAM_File', '')
				Sequencer                  = getValue ('Sequencer', 'Illumina HiSeq')
				Tumor_Sample_UUID          = getValue ('Tumor_Sample_UUID', tumor.sample)
				Matched_Norm_Sample_UUID   = getValue ('Matched_Norm_Sample_UUID', normal.sample)
				
				outs = map (lambda x: str(eval(x)), header)
				fout.write ("\\t".join(outs) + "\\n")
				break
		if found: break
		
fout.close()	
"""
