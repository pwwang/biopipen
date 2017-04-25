from pyppl import proc

#############################
# snpEff utilities          #
#############################

"""
@name:
	pAnn
@description:
	This is the default command. It is used for annotating variant filed (e.g. VCF files).
@input:
	`infile:file`:  The input file 
@output:
	`outdir:file`: The directory containing output anntated file, snpEff_genes.txt and snpEff_summary.html
@args:
	`bin`:       The snpEff executable, default: "snpEff"
	`params`:    Other parameters for `snpEff`, default: "-Xms1g -Xmx4g -v"
	`genome`:    The genome used for annotation, default: "hg19"
	`informat`:  The format of input file [vcf or bed], default: "vcf"
	`outformat`: The format of output file [vcf, gatk, bed, bedAnn], default: "vcf"
	`csvStats`:  Whether to generate csv stats file, default: "{{infile.fn}}.stats.csv", set False to disable.
	`stats`:     The name of the summary file, default: "{{infile.fn}}.summary.html", set False to disable.
@requires:
	[snpEff](http://snpeff.sourceforge.net/SnpEff_manual.html)
"""
pAnn = proc()
pAnn.input  = "infile:file"
pAnn.output = "outdir:file:{{infile.fn}}_snpeff"
pAnn.args   = { "bin": "snpEff", "params": "-Xms1g -Xmx4g -v", "genome": "hg19", "informat": "vcf", "outformat": "vcf", "csvStats": True, "stats": True }
pAnn.script = """
mkdir -p "{{outdir}}"
cd "{{outdir}}"
csvStats=""
if [[ "{{proc.args.csvStats}}" == "True" ]]; then
	csvStats='-csvStats "{{infile.fn}}.stats.csv"'
elif [[ "{{proc.args.csvStats}}" != "False" ]]; then
	csvStats='-csvStats "{{proc.args.csvStats}}"'
fi
stats="-noStats"
if [[ "{{proc.args.stats}}" == "True" ]]; then
	stats='-stats "{{infile.fn}}.summary.html"'
elif [[ "{{proc.args.stats}}" != "False" ]]; then
	stats='-stats "{{proc.args.stats}}'
fi
{{proc.args.bin}} -i {{proc.args.informat}} -o {{proc.args.outformat}} $csvStats $stats {{proc.args.params}} {{proc.args.genome}} "{{infile}}" > "{{infile.fn}}.{{proc.args.outformat}}"
"""

"""
@name:
	pCount
@description:
	Count how many intervals (from a BAM, BED or VCF file) overlap with each genomic interval.
@input:
	`infiles:files`:  The input files
@output:
	`outdir:file`: The directory containg summary (html) file and output file
@args:
	`bin`:     The snpEff executable, default: "snpEff"
	`params`:  Other parameters for `snpEff`, default: "-Xms1g -Xmx4g -v"
	`genome`:  The genome used for counting, default: "hg19"
@requires:
	[snpEff](http://snpeff.sourceforge.net/SnpEff_manual.html)
"""
pCount = proc()
pCount.input  = "infiles:files"
pCount.output = "outdir:file:{{infiles.fn | [0]}}_etc_snpEff_count"
pCount.args   = { "bin": "snpEff", "params": "-v", "genome": "hg19" }
pCount.script = """
mkdir -p "{{outdir}}"
cd "{{outdir}}"
{{proc.args.bin}} count {{proc.args.params}} {{proc.args.genome}} "{{infiles | '" "'.join(_)}}" > "{{outdir}}/{{infiles.fn | [0]}}.etc.snpEff.count.txt"
"""