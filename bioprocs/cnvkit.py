from pyppl import proc

"""
@name:
	pCNVkitAccess
@description:
	Calculate the sequence-accessible coordinates in chromosomes from the given reference genome, output as a BED file.
@input:
	`fafile:file`: The fasta file
@output:
	`outfile:file`: The output file
@args:
	`params`: Other parameters for `cnvkit.py access`
	`cnvkit`: The executable of cnvkit. Default: 'cnvkit.py'
@requires:
	[CNVkit](http://cnvkit.readthedocs.io/)
"""
pCNVkitAccess = proc ()
pCNVkitAccess.input  = "fafile:file"
pCNVkitAccess.args   = {"params": "", "cnvkit": "cnvkit.py"}
pCNVkitAccess.output = "outfile:file:{{fafile | fn}}.access.bed"
pCNVkitAccess.script = """
{{args.cnvkit}} access -o "{{outfile}}" {{args.params}} "{{fafile}}"
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
	`cnvkit`:  The executable of cnvkit. Default: 'cnvkit.py'
	`params`: Other parameters for `cnvkit.py target`
@requires:
	[CNVkit](http://cnvkit.readthedocs.io/)
"""
pCNVkitTarget = proc ()
pCNVkitTarget.input  = "acfile:file, anfile:file"
pCNVkitTarget.output = "outfile:file:{{acfile | fn}}-{{anfile | fn}}.targets"
pCNVkitTarget.args   = {"cnvkit": "cnvkit.py", "params": ""}
pCNVkitTarget.script = """
{{args.cnvkit}} target "{{acfile}}" --split --short-names --annotate "{{anfile}}" -o "{{outfile}}" {{args.params}}
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
	`tgfile`:  The target file
	`cnvkit`:  The executable of cnvkit. Default: 'cnvkit.py'
	`nthread`: The number of threads to use. Default: 1
	`params`:  Other parameters for `cnvkit.py coverage`
@requires:
	[CNVkit](http://cnvkit.readthedocs.io/)
"""
pCNVkitCov = proc ()
pCNVkitCov.input  = "infile:file"
pCNVkitCov.output = "outfile:file:{{infile | fn}}.covcnn"
pCNVkitCov.args   = {"cnvkit": "cnvkit.py", "nthread": 1, "tgfile": "", "params": ""}
pCNVkitCov.script = """
if [[ -z "{{args.tgfile}}" ]]; then
	echo "Missing target file" 1>&2
	exit 1
fi
{{args.cnvkit}} coverage "{{infile}}" "{{args.tgfile}}" -p {{args.nthread}} -o "{{outfile}}" {{args.params}}
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
	`cnvkit`:  The executable of cnvkit. Default: 'cnvkit.py'
	`params`:  Other parameters for `cnvkit.py reference`, default: " --no-edge "
@requires:
	[CNVkit](http://cnvkit.readthedocs.io/)
"""
pCNVkitRef = proc ()
pCNVkitRef.input  = "indir:file"
pCNVkitRef.output = "outfile:file:{{indir | fn}}.refcnn"
pCNVkitRef.args   = {"cnvkit": "cnvkit.py", "params": " --no-edge "}
pCNVkitRef.script = """
{{args.cnvkit}} reference "{{indir}}"/*.covcnn -o "{{outfile}}" {{args.params}}
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
	`outfile:file`: The cnr file
@args:
	`cnvkit`:  The executable of cnvkit. Default: 'cnvkit.py'
	`params`:  Other parameters for `cnvkit.py fix`, default: " --no-edge "
@requires:
	[CNVkit](http://cnvkit.readthedocs.io/)
"""
pCNVkitFix = proc ()
pCNVkitFix.input  = "infile:file, rcfile:file"
pCNVkitFix.output = "outfile:file:{{infile | fn}}.cnr"
pCNVkitFix.args   = {"cnvkit": "cnvkit.py", "params": " --no-edge "}
pCNVkitFix.script = """
{{args.cnvkit}} fix "{{infile}}" <(echo "") "{{rcfile}}" -o "{{outfile}}" {{args.params}}
"""

"""
@name:
	pCNVkitSeg
@description:
	Infer discrete copy number segments from the given coverage table
@input:
	`infile:file`:  The cnr file 
@output:
	`outfile:file`: The cns file
@args:
	`cnvkit`:  The executable of cnvkit. Default: 'cnvkit.py'
	`nthread`: The number of threads to use. Default: 1
	`params`:  Other parameters for `cnvkit.py segment`, default: ""
@requires:
	[CNVkit](http://cnvkit.readthedocs.io/)
"""
pCNVkitSeg = proc ()
pCNVkitSeg.input  = "infile:file"
pCNVkitSeg.output = "outfile:file:{{infile | fn}}.cns"
pCNVkitSeg.args   = {"cnvkit": "cnvkit.py", "nthread": 1, "params": ""}
pCNVkitSeg.script = """
{{args.cnvkit}} segment -o "{{outfile}}" -p {{args.nthread}} {{args.params}} "{{infile}}"
"""

"""
@name:
	pCNVkitCall
@description:
	Given segmented log2 ratio estimates (.cns), derive each segment's absolute integer copy number 
@input:
	`infile:file`:  The cns file 
@output:
	`outfile:file`: The callcns file
@args:
	`cnvkit`:  The executable of cnvkit. Default: 'cnvkit.py'
	`params`:  Other parameters for `cnvkit.py segment`, default: ""
@requires:
	[CNVkit](http://cnvkit.readthedocs.io/)
"""
pCNVkitCall = proc (desc="Given segmented log2 ratio estimates (.cns), derive each segment's absolute integer copy number")
pCNVkitCall.input  = "infile:file"
pCNVkitCall.output = "outfile:file:{{infile | fn}}.callcns"
pCNVkitCall.args   = {"cnvkit": "cnvkit.py", "params": ""}
pCNVkitCall.script = """
{{args.cnvkit}} call -o "{{outfile}}" {{args.params}} "{{infile}}"
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
	`cnvkit`:   The executable of cnvkit. Default: 'cnvkit.py'
	`region`:       The region for zoom-in plots. Default: '' (don't plot zoom-in view)
	`gene`:         The genes to be highlighted. Default: ''
	`scatter`:      Whether to generate the scatter plot. Default: True
	`diagram`:      Whether to generate the diagram plot. Default: True
	`heatmap`:      Whether to generate the heatmap plot. Default: True
@requires:
	[CNVkit](http://cnvkit.readthedocs.io/)
"""
pCNVkitPlot = proc ()
pCNVkitPlot.input  = "cnrdir:file, cnsdir:file"
pCNVkitPlot.output = "outdir:dir:{{cnrdir | fn}}.cnvkit.plots"
pCNVkitPlot.args   = {"cnvkit": "cnvkit.py", "region": "", "gene": "", "scatter": True, "diagram": True, "heatmap": True}
pCNVkitPlot.script = """
region="{{args.region}}"
gene="{{args.gene}}"
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
	if [[ "{{args.scatter}}" == "True" ]]; then
		scatterout1="{{outdir}}/$cnrname.scatter.pdf"
		{{args.cnvkit}} scatter "$cnrfile" -s "$cnsfile" -o "$scatterout1"
		if [[ ! -z "$region" ]] || [[ ! -z "$gene" ]]; then
			scatterout2="{{outdir}}/${cnrname}_$fn.scatter.pdf"
			{{args.cnvkit}} scatter "$cnrfile" -s "$cnsfile" $regnargs $geneargs -o "$scatterout2"
		fi
	fi
done

if [[ "{{args.diagram}}" == "True" ]]; then
	diagramout="{{outdir}}/{{cnrdir | fn}}.diagram.pdf"
	{{args.cnvkit}} diagram -s "$cnsfile" "$cnrfile" -o "$diagramout"
fi

if [[ "{{args.heatmap}}" == "True" ]]; then
		heatmapout1="{{outdir}}/{{cnrdir | fn}}.heatmap.pdf"
		{{args.cnvkit}} heatmap "{{cnsdir}}"/*.cns -d -o "$heatmapout1"
		if [[ ! -z "$region" ]]; then
			heatmapout2="{{outdir}}/{{cnrdir | fn}}_$fn.heatmap.pdf"
			{{args.cnvkit}} heatmap "{{cnsdir}}"/*.cns -d $regnargs -o "$heatmapout2"
		fi
	fi
"""


"""
@name:
	pCNVkitRpt
@description:
	Report CNVkit results
@input:
	`cnrfile:file`:  The file containing copy number ratio
	`cnsfile:file`:  The file containing copy number segment
@output:
	`outdir:dir`:   The output directory
@args:
	`cnvkit`:   The executable of cnvkit. Default: 'cnvkit.py'
	`breaks`:       Whether to report breakpoints. Default: True
	`gainloss`:     Whether to report gainloss. Default: True
	`metrics`:      Whether to report metrics. Default: True
	`segmetrics`:   Whether to report segmetrics. Default: True
@requires:
	[CNVkit](http://cnvkit.readthedocs.io/)
"""
pCNVkitRpt = proc ()
pCNVkitRpt.input  = "cnrfile:file, cnsfile:file"
pCNVkitRpt.output = "outdir:dir:{{cnrfile | fn}}.cnvkit.reports"
pCNVkitRpt.args   = {"cnvkit": "cnvkit.py", "breaks": True, "gainloss": True, "metrics": True, "segmetrics": True}
pCNVkitRpt.script = """
bn="{{cnrfile | fn}}"
if [[ "{{args.breaks}}" == "True" ]]; then
	{{args.cnvkit}} breaks "{{cnrfile}}" "{{cnsfile}}" -o "{{outdir}}/$bn.breaks.txt"
fi

if [[ "{{args.gainloss}}" == "True" ]]; then
	{{args.cnvkit}} gainloss "{{cnrfile}}" -s "{{cnsfile}}" -o "{{outdir}}/$bn.gainloss.txt"
fi

if [[ "{{args.metrics}}" == "True" ]]; then
	{{args.cnvkit}} metrics "{{cnrfile}}" -s "{{cnsfile}}" -o "{{outdir}}/$bn.metrics.txt"
fi

if [[ "{{args.segmetrics}}" == "True" ]]; then
	{{args.cnvkit}} segmetrics "{{cnrfile}}" -s "{{cnsfile}}" --iqr -o "{{outdir}}/$bn.segmetrics.txt"
fi
"""

"""
@name:
	pCNVkit2Vcf
@description:
	Output vcf file for cnvkit results
@input:
	`cnsfile:file`: The cns file
@output:
	`outfile:file`: The vcf file
@args:
	`cnvkit`:   The executable of cnvkit. Default: 'cnvkit.py'
	`params`:   Other params for `cnvkit.py export`
@requires:
	[CNVkit](http://cnvkit.readthedocs.io/)
"""
pCNVkit2Vcf = proc ()
pCNVkit2Vcf.input  = "cnsfile:file"
pCNVkit2Vcf.output = "outfile:file:{{cnsfile | fn}}.cnvkit.vcf"
pCNVkit2Vcf.args   = {"cnvkit": "cnvkit.py", "params": ""}
pCNVkit2Vcf.script = """
cnvkit.py export vcf "{{cnsfile}}" -o "{{outfile}}" {{args.params}}
"""