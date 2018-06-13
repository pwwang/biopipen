from pyppl import Proc

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
	`cnvnator`:      The CNVnator executable, default: "cnvnator"
	`cnv2vcf`:  The converter executable to convert CNVnator results to vcf, default: "cnvnator2VCF.pl"
	`binsize`:  The bin_size, default: 100
	`genome`:   The genome: default: hg19
	`chrom`:    Chromosome names, default: "" (all chromosomes)
	`chrdir`:   The dir contains reference sequence of chromosomes, default: "" (don't specify)
	
@requires:
	[CNVnator](https://github.com/abyzovlab/CNVnator)
"""
pCNVnator = Proc()
pCNVnator.input  = "infile:file"
pCNVnator.output = "outfile:file:{{infile | fn}}.cnv.vcf"
pCNVnator.args   = { "cnvnator": "cnvnator", "cnv2vcf": "cnvnator2VCF.pl", "binsize": 100, "genome": "hg19", "chrom": "", "chrdir": "" }
pCNVnator.script = """
rootfile="{{outfile}}.root"
[[ -z "{{args.genome}}" ]] && genomeparam="" || genomeparam="-genome {{args.genome}}"
[[ -z "{{args.chrom}}" ]]  && chromparam=""  || chromparam="-chrom {{args.chrom}}"
[[ -z "{{args.chrdir}}" ]] && chrdirparam="" || chrdirparam="-d {{args.chrdir}}"
# EXTRACTING READ MAPPING FROM BAM/SAM FILES
{{args.cnvnator}} $genomeparam -root "$rootfile" $chromparam -tree "{{infile}}"
# GENERATING HISTOGRAM
{{args.cnvnator}} $genomeparam -root "$rootfile" $chromparam -his {{args.binsize}} $chrdirparam
# CALCULATING STATISTICS
{{args.cnvnator}} -root "$rootfile" $chromparam -stat {{args.binsize}}
# RD SIGNAL PARTITIONING
{{args.cnvnator}} -root "$rootfile" $chromparam  -partition {{args.binsize}}
# CNV CALLING
{{args.cnvnator}} -root "$rootfile" $chromparam -call {{args.binsize}} > "{{outfile}}.cnvnator"
# Convert results to VCF:
{{args.cnv2vcf}} "{{outfile}}.cnvnator" > "{{outfile}}"
"""

