"""Utilities of cnvkit"""
from pyppl import Proc, Diot
from . import params
from .utils import fs2name
from . import delefactory, procfactory
from modkit import Modkit
Modkit().delegate(delefactory())

pCNVkitPrepare = proc_factory(
	desc = 'Generate target files for cnvkit.',
	config = Diot(annotate = """
	@name:
		pCNVkitPrepare
	@description:
		Generate target files for cnvkit, using probably cnvkit's access, target, autobin commands.
	@input:
		`infiles:files`: The bam files. Indexes are necessary.
			- Hint: if indexes are not with the input files, you probably need `pCNVkitPrepare.infile = 'origin'`
	@output:
		`target:file`:     The (autobinned) target file
		`antitarget:file`: The (autobinned) target file
	@args:
		`cnvkit`:  The executable of cnvkit. Default: 'cnvkit.py'
		`baits` : The bait file for the regions you captured in the experiment.
			- See https://github.com/AstraZeneca-NGS/reference_data/tree/master/hg19/bed
		`accfile`: Directly use the access file. Default: generating from the reference file.
			- See https://github.com/etal/cnvkit/tree/master/data
		`nthread`: The number of threads to use. Default: 1
		`ref`    : The reference genome.
		`params` : The extra parameters for cnvkit's `access`, `target` and `autobin` command. Default:
			```python
			Diot(
				target  = Diot({'short-name': True, 'split': True}),
				access  = Diot(s = '5000'),
				autobin = Diot()
			)
			```
	@requires:
		[CNVkit](http://cnvkit.readthedocs.io/)
	"""))
pCNVkitPrepare.input        = 'infiles:files'
pCNVkitPrepare.output       = 'target:file:{{i.infiles | fs2name}}.target.bed, antitarget:file:{{i.infiles | fs2name}}.antitarget.bed'
pCNVkitPrepare.args.cnvkit  = params.cnvkit.value
pCNVkitPrepare.args.baits   = params.refexon.value
pCNVkitPrepare.args.accfile = ''
pCNVkitPrepare.args.nthread = 1
pCNVkitPrepare.args.ref     = params.ref.value
pCNVkitPrepare.args.params  = Diot(
	target  = Diot({'short-name': True, 'split': True}),
	access  = Diot(s = '5000'),
	autobin = Diot()
)
pCNVkitPrepare.envs.fs2name = fs2name
pCNVkitPrepare.lang         = params.python.value

pCNVkitCov = proc_factory(
	desc = 'Calculate coverage in the given regions from BAM read depths.',
	config = Diot(annotate = """
	@name:
		pCNVkitCov
	@description:
		Calculate coverage in the given regions from BAM read depths.
	@input:
		`infile:file`:  The bam file
		`tgfile:file`:  The target file
		`atgfile:file`: The antitarget file
	@output:
		`outfile:file`: The output cnn file
	@args:
		`cnvkit`:  The executable of cnvkit. Default: 'cnvkit.py'
		`nthread`: The number of threads to use. Default: 1
		`params`:  Other parameters for `cnvkit.py coverage`
	@requires:
		[CNVkit](http://cnvkit.readthedocs.io/)
	"""))
pCNVkitCov.input        = "infile:file, tgfile:file, atgfile:file"
pCNVkitCov.output       = "outfile:file:{{i.infile | fn}}.target.cnn, antifile:file:{{i.infile | fn}}.antitarget.cnn"
pCNVkitCov.args.cnvkit  = params.cnvkit.value
pCNVkitCov.args.nthread = 1
pCNVkitCov.args.params  = Diot()
pCNVkitCov.lang         = params.python.value

pCNVkitRef = proc_factory(
	desc = 'Compile a copy-number reference from the given files or directory',
	config = Diot(annotate = """
	@name:
		pCNVkitRef
	@description:
		Compile a copy-number reference from the given files or directory (containing normal samples). If given a reference genome (-f option), also calculate the GC content and repeat-masked proportion of each region.
	@input:
		`infiles:files`:  The input reference coverage files
	@output:
		`outfile:file`: The output reference cnn file
	@args:
		`ref`   :  The reference file.
		`cnvkit`:  The executable of cnvkit. Default: 'cnvkit.py'
		`nthread`: The number of threads to use. Default: 1
		`params`:  Other parameters for `cnvkit.py reference`, default: " --no-edge "
	@requires:
		[CNVkit](http://cnvkit.readthedocs.io/)
	"""))
pCNVkitRef.input        = "infiles:files"
pCNVkitRef.output       = "outfile:file:{{i.infiles | fs2name}}.reference.cnn"
pCNVkitRef.args.cnvkit  = params.cnvkit.value
pCNVkitRef.args.ref     = params.ref.value
pCNVkitRef.args.nthread = 1
pCNVkitRef.args.params  = Diot({'no-edge': True})
pCNVkitRef.envs.fs2name = fs2name
pCNVkitRef.lang         = params.python.value

pCNVkitFlatRef = proc_factory(
	desc = 'Compile a copy-number flat reference without normal samples',
	config = Diot(annotate = """
	@name:
		pCNVkitFlatRef
	@description:
		Generate reference coverage if there are no normal samples.
	@input:
		`tgfile:file`:  The target file
		`atgfile:file`: The antitarget file
	@output:
		`outfile:file`: The output reference cnn file
	@args:
		`ref`   :  The reference file.
		`cnvkit`:  The executable of cnvkit. Default: 'cnvkit.py'
		`params`:  Other parameters for `cnvkit.py reference`, default: `{}`
	@requires:
		[CNVkit](http://cnvkit.readthedocs.io/)
	"""))
pCNVkitFlatRef.input        = "tgfile:file, atgfile:file"
pCNVkitFlatRef.output       = "outfile:file:{{i.tgfile | fn}}.reference.cnn"
pCNVkitFlatRef.args.cnvkit  = params.cnvkit.value
pCNVkitFlatRef.args.ref     = params.ref.value
pCNVkitFlatRef.args.params  = Diot()
pCNVkitFlatRef.envs.fs2name = fs2name
pCNVkitFlatRef.lang         = params.python.value

pCNVkitFix = proc_factory(
	desc = 'Combine the uncorrected target and antitarget coverage tables and correct them.',
	config = Diot(annotate = """
	@name:
		pCNVkitFix
	@description:
		Combine the uncorrected target and antitarget coverage tables (.cnn) and correct for biases in regional coverage and GC content, according to the given reference. Output a table of copy number ratios (.cnr)
	@input:
		`tgfile:file`:  The target coverage file
		`atgfile:file`: The antitarget coverage file
		`rcfile:file`:  The reference cnn file
	@output:
		`outfile:file`: The cnr file
	@args:
		`nthread`: The number of threads to use. Default: 1
		`cnvkit`:  The executable of cnvkit. Default: 'cnvkit.py'
		`params`:  Other parameters for `cnvkit.py fix`, default: " --no-edge "
	@requires:
		[CNVkit](http://cnvkit.readthedocs.io/)
	"""))
pCNVkitFix.input        = "tgfile:file, atgfile:file, rcfile:file"
pCNVkitFix.output       = "outfile:file:{{i.tgfile | fn}}.cnr"
pCNVkitFix.args.cnvkit  = params.cnvkit.value
pCNVkitFix.args.params  = Diot({'no-edge': True})
pCNVkitFix.args.nthread = 1
pCNVkitFix.lang         = params.python.value

pCNVkitSeg = proc_factory(
	desc = 'Infer discrete copy number segments from the given coverage table.',
	config = Diot(annotate = """
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
	"""))
pCNVkitSeg.input        = "infile:file"
pCNVkitSeg.output       = "outfile:file:{{i.infile | fn}}.cns"
pCNVkitSeg.args.cnvkit  = params.cnvkit.value
pCNVkitSeg.args.nthread = 1
pCNVkitSeg.args.params  = Diot()
pCNVkitSeg.lang         = params.python.value

pCNVkitCall = proc_factory(
	desc="Given segmented log2 ratio estimates (.cns), derive each segment's absolute integer copy number",
	config = Diot(annotate = """
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
	"""))
pCNVkitCall.input       = "infile:file"
pCNVkitCall.output      = "outfile:file:{{i.infile | fn}}.callcns"
pCNVkitCall.args.cnvkit = params.cnvkit.value
pCNVkitCall.args.params = Diot()
pCNVkitCall.lang        = params.python.value

pCNVkitScatter = proc_factory(
	desc = 'Generate scatter plot for CNVkit results.',
	config = Diot(annotate = """
	@name:
		pCNVkitScatter
	@description:
		Generate scatter plot for CNVkit results.
	@input:
		`cnrfile:file`: The cnr file
		`cnsfile:file`: The cns file from call
	@output:
		`outdir:dir`: The output directory
	@args:
		`cnvkit`:  The executable of cnvkit. Default: 'cnvkit.py'
		`nthread`: The number of threads to use. Default: 1
		`params`:  Other parameters for `cnvkit.py scatter`
		`regions`: The regoins to plot. Default: `['']`
			- You can have extra specific regions, format:
			- `chr5:100-50000000:TERT` or `chr7:BRAF,MET` (genes are used to highlight)
	"""))
pCNVkitScatter.input        = 'cnrfile:file, cnsfile:file'
pCNVkitScatter.output       = 'outdir:dir:{{i.cnrfile | fn}}.scatters'
pCNVkitScatter.args.cnvkit  = params.cnvkit.value
pCNVkitScatter.args.nthread = 1
pCNVkitScatter.args.params  = Diot()
pCNVkitScatter.args.regions = [
	'', # plot whole genome
	# extra regions, format: chr5:100-50000000:TERT
	# or chr7:BRAF,MET
]
pCNVkitScatter.lang   = params.python.value

pCNVkitDiagram = proc_factory(
	desc = 'Generate diagram plot for CNVkit results.',
	config = Diot(annotate = """
	@name:
		pCNVkitDiagram
	@description:
		Generate diagram plot for CNVkit results.
	@input:
		`cnrfile:file`: The cnr file
		`cnsfile:file`: The cns file from call
	@output:
		`outfile:file`: The output file
	@args:
		`cnvkit`:  The executable of cnvkit. Default: 'cnvkit.py'
		`nthread`: The number of threads to use. Default: 1
		`params`:  Other parameters for `cnvkit.py scatter`
	"""))
pCNVkitDiagram.input       = 'cnrfile:file, cnsfile:file'
pCNVkitDiagram.output      = 'outfile:file:{{i.cnrfile | fn}}.diagram.pdf'
pCNVkitDiagram.args.cnvkit = params.cnvkit.value
pCNVkitDiagram.args.nthread = 1
pCNVkitDiagram.args.params = Diot()
pCNVkitDiagram.lang        = params.python.value

pCNVkitHeatmap = proc_factory(
	desc = 'Generate heatmap plot for CNVkit results.',
	config = Diot(annotate = """
	@name:
		pCNVkitHeatmap
	@description:
		Generate heatmap plot for CNVkit results.
	@input:
		`cnfiles:files`: The cnr or cns files.
	@output:
		`outdir:dir`: The output directory
	@args:
		`cnvkit`:  The executable of cnvkit. Default: 'cnvkit.py'
		`params`:  Other parameters for `cnvkit.py scatter`
		`regions`: The regoins to plot. Default: `['']`
			- You can have extra specific regions, format:
			- `chr5:100-50000000` or `chr7` (genes are used to highlight)
	"""))
pCNVkitHeatmap.input        = 'cnfiles:files'
pCNVkitHeatmap.output       = 'outdir:dir:{{i.cnfiles | fs2name}}.heatmaps'
pCNVkitHeatmap.args.cnvkit  = params.cnvkit.value
pCNVkitHeatmap.args.params  = Diot()
pCNVkitHeatmap.args.nthread = 1
pCNVkitHeatmap.args.regions = [
	'', # plot whole genome
	# extra regions, format: chr5:100-50000000
	# or chr7
]
pCNVkitHeatmap.envs.fs2name = fs2name
pCNVkitHeatmap.lang         = params.python.value

pCNVkitReport = proc_factory(
	desc = 'Report CNVkit results',
	config = Diot(annotate = """
	@name:
		pCNVkitReport
	@description:
		Report CNVkit results
	@input:
		`cnrfile:file`:  The file containing copy number ratio
		`cnsfile:file`:  The file containing copy number segment
	@output:
		`outdir:dir`:   The output directory
	@args:
		`cnvkit` : The executable of cnvkit. Default: 'cnvkit.py'
		`nthread`: The number of threads to use. Default: 1
		`params` : Extra parameters to the commands.
			- `breaks`:       Whether to report breakpoints. Default: True
			- `gainloss`:     Whether to report gainloss. Default: True
			- `metrics`:      Whether to report metrics. Default: True
			- `segmetrics`:   Whether to report segmetrics. Default: True
	@requires:
		[CNVkit](http://cnvkit.readthedocs.io/)
	"""))
pCNVkitReport.input        = "cnrfile:file, cnsfile:file"
pCNVkitReport.output       = "outdir:dir:{{i.cnrfile | fn}}.cnvkit.reports"
pCNVkitReport.args.cnvkit  = params.cnvkit.value
pCNVkitReport.args.nthread = 1
pCNVkitReport.args.params  = Diot(
	breaks     = True, # None to disable
	gainloss   = True,
	metrics    = True,
	segmetrics = Diot(iqr = True)
)
pCNVkitReport.lang   = params.python.value

pCNVkit2Vcf = proc_factory(
	desc = 'Output vcf file for cnvkit results',
	config = Diot(annotate = """
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
		`nthread`: The number of threads to use. Default: 1
		`params`:   Other params for `cnvkit.py export`
	@requires:
		[CNVkit](http://cnvkit.readthedocs.io/)
	"""))
pCNVkit2Vcf.input        = "cnsfile:file"
pCNVkit2Vcf.output       = "outfile:file:{{i.cnsfile | fn}}.cnvkit.vcf"
pCNVkit2Vcf.args.cnvkit  = params.cnvkit.value
pCNVkit2Vcf.args.nthread = 1
pCNVkit2Vcf.args.params  = Diot()
pCNVkit2Vcf.lang         = params.python.value

pCNVkit2Theta = proc_factory(
	desc = 'Convert the results to THetA2 interval input.',
	config = Diot(annotate = """
	@name:
		pCNVkit2Theta
	@description:
		Convert the results to THetA2 interval input.
	@input:
		`cnsfile:file`: The cns file
		`cnnfile:file`: The reference cnn file or the cnr file for paired Normal sample. Could be empty.
	@output:
		`outfile:file`: The interval file for THetA2
	@args:
		`nthread` : Number threads to use. Default: `1`
		`cnvkit`  : The executable of cnvkit. Default: `cnvkit.py`
		`params`  : Other params for `cnvkit.py export theta`
	"""))
pCNVkit2Theta.input        = 'cnsfile:file, cnnfile:file'
pCNVkit2Theta.output       = 'outfile:file:{{i.cnsfile | fn2}}.interval.txt'
pCNVkit2Theta.args.nthread = 1
pCNVkit2Theta.args.params  = Diot()
pCNVkit2Theta.args.cnvkit  = params.cnvkit.value
pCNVkit2Theta.lang         = params.python.value
