from pyppl import proc

"""
@name:
	pSnpEff
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
	`csvStats`:  Whether to generate csv stats file, default: True.
	`htmlStats`: Whether to generate the html summary file, default: False.
@requires:
	[snpEff](http://snpeff.sourceforge.net/SnpEff_manual.html)
"""
pSnpEff = proc()
pSnpEff.input  = "infile:file"
pSnpEff.output = "outdir:dir:{{infile | fn}}.snpeff"
pSnpEff.args   = { "bin": "snpEff", "params": "-Xms1g -Xmx4g -v", "genome": "hg19", "informat": "vcf", "outformat": "vcf", "csvStats": True, "htmlStats": False }
pSnpEff.script = """
csvfile="{{outdir}}/{{infile | fn}}.csvstat"
sumfile="{{outdir}}/{{infile | fn}}.html"
outfile="{{outdir}}/{{infile | fn}}.snpEff.vcf"
csvStats=""
if [[ "{{proc.args.csvStats}}" == "True" ]]; then
	csvStats="-csvStats \\"$csvfile\\""
fi
stats=""
if [[ "{{proc.args.htmlStats}}" == "True" ]]; then
	stats="-stats \\"$sumfile\\""
fi
{{proc.args.bin}} -i {{proc.args.informat}} -o {{proc.args.outformat}} $csvStats $stats {{proc.args.params}} {{proc.args.genome}} "{{infile}}" > "$outfile"
"""
