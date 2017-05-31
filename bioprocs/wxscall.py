from pyppl import proc

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


