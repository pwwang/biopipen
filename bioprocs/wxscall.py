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
pCNVnator.output = "outfile:file:{{infile | fn}}.cnv.vcf"
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
	pCNVkitTarget
@description:
	Generate targets file for CNVkit using access file and annotate file (`cnvkit.py target`)
@input:
	`acfile:file`: The access file
	`anfile:file`: The annotate file
@output:
	`outfile:file`: The targets file
@args:
	`bin-cnvkit`:  The executable of cnvkit. Default: 'cnvkit.py'
@requires:
	[CNVkit](http://cnvkit.readthedocs.io/)
"""
pCNVkitTarget = proc ()
pCNVkitTarget.input  = "acfile:file, anfile:file"
pCNVkitTarget.output = "outfile:file:{{acfile | fn}}-{{anfile | fn}}.targets"
pCNVkitTarget.args   = {"bin-cnvkit": "cnvkit.py"}
pCNVkitTarget.script = """
{{proc.args.bin-cnvkit}} target "{{acfile}}" --split --short-names --annotate "{{anfile}}" -o "{{outfile}}"
"""

"""
@name:
	pCNVkitCov
@description:
	Calculate coverage in the given regions from BAM read depths.
@input:
	`infile:file`: The bam file
@output:
	`outfile:file`: The output cnn file
@args:
	`tgfile:file`: The target file
	`bin-cnvkit`:  The executable of cnvkit. Default: 'cnvkit.py'
	`nthread`:     The number of threads to use. Default: 1
@requires:
	[CNVkit](http://cnvkit.readthedocs.io/)
"""
pCNVkitCov = proc ()
pCNVkitCov.input  = "infile:file"
pCNVkitCov.output = "outfile:file:{{infile | fn}}.cov.cnn"
pCNVkitCov.args   = {"bin-cnvkit": "cnvkit.py", "nthread": 1, "tgfile": ""}
pCNVkitCov.script = """
if [[ -z "{{proc.args.tgfile}}" ]]; then
	echo "Missing target file" 1>&2
	exit 1
fi
{{proc.args.bin-cnvkit}} coverage "{{infile}}" "{{proc.args.tgfile}}" -p {{proc.args.nthread}} -o "{{outfile}}"
"""

"""
@name:
	pCNVkitRef
@description:
	Compile a copy-number reference from the given files or directory (containing normal samples). If given a reference genome (-f option), also calculate the GC content and repeat-masked proportion of each region.
@input:
	`indir:file`:  The input directory containing the cnn files
@output:
	`outfile:file`: The output reference cnn file
@args:
	`bin-cnvkit`:  The executable of cnvkit. Default: 'cnvkit.py'
@requires:
	[CNVkit](http://cnvkit.readthedocs.io/)
"""
pCNVkitRef = proc ()
pCNVkitRef.input  = "indir:file"
pCNVkitRef.output = "outfile:file:{{indir | fn}}.ref.cnn"
pCNVkitRef.args   = {"bin-cnvkit": "cnvkit.py"}
pCNVkitRef.script = """
{{proc.args.bin-cnvkit}} reference "{{indir}}"/*.cov.cnn --no-edge -o "{{outfile}}"
"""

"""
@name:
	pCNVkitFix
@description:
	Combine the uncorrected target and antitarget coverage tables (.cnn) and correct for biases in regional coverage and GC content, according to the given reference. Output a table of copy number ratios (.cnr)
@input:
	`infile:file`:  The cnn file to be fixed
	`rcfile:file`:  The reference cnn file
@output:
	`outfile:file`: The fixed cnn file
@args:
	`bin-cnvkit`:  The executable of cnvkit. Default: 'cnvkit.py'
@requires:
	[CNVkit](http://cnvkit.readthedocs.io/)
"""
pCNVkitFix = proc ()
pCNVkitFix.input  = "infile:file, rcfile:file"
pCNVkitFix.output = "outfile:file:{{infile | fn}}.cnr"
pCNVkitFix.args   = {"bin-cnvkit": "cnvkit.py"}
pCNVkitFix.script = """
{{proc.args.bin-cnvkit}} fix "{{infile}}" <(echo "") "{{rcfile}}" --no-edge -o "{{outfile}}"
"""

"""
@name:
	pCNVkitSeg
@description:
	Infer discrete copy number segments from the given coverage table
@input:
	`infile:file`:  The cnn file to be fixed
	`rcfile:file`:  The reference cnn file
@output:
	`outfile:file`: The fixed cnn file
@args:
	`bin-cnvkit`:  The executable of cnvkit. Default: 'cnvkit.py'
	`nthread`:     The number of threads to use. Default: 1
@requires:
	[CNVkit](http://cnvkit.readthedocs.io/)
"""
pCNVkitSeg = proc ()
pCNVkitSeg.input  = "infile:file"
pCNVkitSeg.output = "outfile:file:{{infile | fn}}.cns"
pCNVkitSeg.args   = {"bin-cnvkit": "cnvkit.py", "nthread": 1}
pCNVkitSeg.script = """
{{proc.args.bin-cnvkit}} fix "{{infile}}" <(echo "") "{{rcfile}}" -p {{proc.args.nthread}} --no-edge -o "{{outfile}}"
"""

"""
@name:
	pCNVkitPlot
@description:
	Plot CNVkit results
@input:
	`cnrdir:file`:  The directory containing copy number ratio files
	`cnsdir:file`:  The directory containing copy number segment files
@output:
	`outdir:dir`:   The output directory
@args:
	`bin-cnvkit`:   The executable of cnvkit. Default: 'cnvkit.py'
	`region`:       The region for zoom-in plots. Default: '' (don't plot zoom-in view)
	`gene`:         The genes to be highlighted. Default: ''
	`scatter`:      Whether to generate the scatter plot. Default: True
	`diagram`:      Whether to generate the diagram plot. Default: True
	`heatmap`:      Whether to generate the heatmap plot. Default: True
@requires:
	[CNVkit](http://cnvkit.readthedocs.io/)
"""
pCNVkitPlot = proc ()
pCNVkitPlot.input  = "cnrdir:file, cnrdir:file"
pCNVkitPlot.output = "outdir:dir:{{cnrdir | fn}}.cnvkit.plots"
pCNVkitPlot.args   = {"bin-cnvkit": "cnvkit.py", "region": "", "gene": "", "scatter": True, "diagram": True, "heatmap": True}
pCNVkitPlot.script = """
region="{{proc.args.region}}"
gene="{{proc.args.gene}}"
fn="$region$gene"
fn=${fn//:/-}
fn=${fn//,/-}
regnargs=""
geneargs=""
if [[ ! -z "$region" ]]; then regnargs="-c $region"; fi
if [[ ! -z "$gene"   ]]; then geneargs="-g $gene"; fi
for cnrfile in "{{cnrdir}}"/*.cnr; do
	cnrname=$(basename $cnrfile .cnr)
	cnsfile="{{cnsdir}}/$cnrname.cns"
	if [[ "{{proc.args.scatter}}" == "True" ]]; then
		scatterout1="{{outdir}}/$cnrname.scatter.pdf"
		{{proc.args.bin-cnvkit}} scatter "$cnrfile" -s "$cnsfile" -o "$scatterout1"
		if [[ ! -z "$region" ]] || [[ ! -z "$gene" ]]; then
			scatterout2="{{outdir}}/${cnrname}_$fn.scatter.pdf"
			{{proc.args.bin-cnvkit}} scatter "$cnrfile" -s "$cnsfile" $regnargs $geneargs -o "$scatterout2"
		fi
	fi
done

if [[ "{{proc.args.diagram}}" == "True" ]]; then
	diagramout="{{outdir}}/{{cnrdir | fn}}.diagram.pdf"
	{{proc.args.bin-cnvkit}} diagram -s "$cnsfile" "$cnrfile" -o "$diagramout"
fi

if [[ "{{proc.args.heatmap}}" == "True" ]]; then
		heatmapout1="{{outdir}}/{{cnrdir | fn}}.heatmap.pdf"
		{{proc.args.bin-cnvkit}} heatmap "{{cnsdir}}"/*.cns -d -o "$heatmapout1"
		if [[ ! -z "$region" ]]; then
			heatmapout2="{{outdir}}/{{cnrdir | fn}}_$fn.heatmap.pdf"
			{{proc.args.bin-cnvkit}} heatmap "{{cnsdir}}"/*.cns -d $regnargs -o "$heatmapout2"
		fi
	fi
"""


"""
@name:
	pCNVkitRpt
@description:
	Report CNVkit results
@input:
	`cnrdir:file`:  The directory containing copy number ratio files
	`cnsdir:file`:  The directory containing copy number segment files
@output:
	`outdir:dir`:   The output directory
@args:
	`bin-cnvkit`:   The executable of cnvkit. Default: 'cnvkit.py'
	`breaks`:       Whether to report breakpoints. Default: True
	`gainloss`:     Whether to report gainloss. Default: True
	`metrics`:      Whether to report metrics. Default: True
	`segmetrics`:   Whether to report segmetrics. Default: True
@requires:
	[CNVkit](http://cnvkit.readthedocs.io/)
"""
pCNVkitRpt = proc ()
pCNVkitRpt.input  = "cnrdir:file, cnrdir:file"
pCNVkitRpt.output = "outdir:dir:{{cnrdir | fn}}.cnvkit.reports"
pCNVkitRpt.args   = {"bin-cnvkit": "cnvkit.py", "breaks": True, "gainloss": True, "metrics": True, "segmetrics": True}
pCNVkitRpt.script = """
for cnrfile in "{{cnrdir}}"/*.cnr; do
	cnrname=$(basename $cnrfile .cnr)
	cnsfile="{{cnsdir}}/$cnrname.cns"
	if [[ "{{proc.args.breaks}}" == "True" ]]; then
		breaks="{{outdir}}/$cnrname.breaks.txt"
		{{proc.args.bin-cnvkit}} breaks "$cnrfile" "$cnsfile" -o "$breaks"
	fi
	
	if [[ "{{proc.args.gainloss}}" == "True" ]]; then
		gainloss="{{outdir}}/$cnrname.gainloss.txt"
		{{proc.args.bin-cnvkit}} gainloss "$cnrfile" -s "$cnsfile" -o "$gainloss"
	fi
	
	if [[ "{{proc.args.metrics}}" == "True" ]]; then
		metrics="{{outdir}}/$cnrname.metrics.txt"
		{{proc.args.bin-cnvkit}} metrics "$cnrfile" -s "$cnsfile" -o "$metrics"
	fi
	
	if [[ "{{proc.args.segmetrics}}" == "True" ]]; then
		segmetrics="{{outdir}}/$cnrname.segmetrics.txt"
		{{proc.args.bin-cnvkit}} segmetrics "$cnrfile" -s "$cnsfile" --iqr -o "$segmetrics"
	fi
done

"""
