from pyppl import Proc, Box
#from .utils import runcmd, plot, helpers
from . import params, rimport
from .utils import fs2name


"""
@name:
	pVcfStatsPlot
@description:
	Convert csvstat file from snpEff to R-readable matrix and plot them.
@input:
	`indir:file`: The directory containing the csv stat files from `snpEff ann`
@output:
	`outdir:dir`: The output directory
@args:
	`chroms`:     The chromsome filter. Default: "" (all chroms)
	- Note: snpEff csvstat file has no "chr" prefix
"""
pVcfStatsPlot                  = Proc(desc = 'Convert csvstat file from snpEff to R-readable matrix and plot them.')
pVcfStatsPlot.input            = "indir:file"
pVcfStatsPlot.output           = "outdir:dir:{{i.indir | fn}}-{{job.index}}.statplots"
pVcfStatsPlot.args.chroms      = ""
pVcfStatsPlot.args.Rscript     = params.Rscript.value
pVcfStatsPlot.args.devpars     = Box({'res':300, 'width':2000, 'height':2000})
pVcfStatsPlot.args.histplotggs = []
pVcfStatsPlot.args.boxplotggs  = []
pVcfStatsPlot.lang             = params.python.value
pVcfStatsPlot.script           = "file:scripts/vcfnext/pVcfStatsPlot.py"

"""
@name:
	pCallRate
@description:
	Calculate sample/snp call rate from single sample vcfs
@input:
	`indir:file`:     The dir containing the vcfs
@output:
	`outsample:file`: The report of call rate for each sample
	`figsample:file`: The bar chat of sample call rates
	`outsnp:file`:    The report of call rate for each snp
	`figsnp:file`:    The bar chat of snp call rates
"""
pCallRate                  = Proc()
pCallRate.input            = "indir:file"
pCallRate.output           = "outdir:dir:{{ i.indir | fn }}.callrate"
pCallRate.args.histplotggs = []
pCallRate.args.devpars     = Box({'res':300, 'width':2000, 'height':2000})
pCallRate.lang             = params.Rscript.value
pCallRate.script           = "file:scripts/vcfnext/pCallRate.r"

"""
@name:
	pCepip
@description:
	Run CEPIP.
@input:
	`infile:file`: The input file (vcf or avinput)
@output:
	`outfile:file`: The cepip result file
@args:
	`cepip`:    The path of cepip
	`cell` :    The related cell line
	`params`:   Other params for cepip
@requires:
	[`cepip`](http://jjwanglab.org/cepip/)
"""
pCepip                     = Proc(desc = 'Run cepip for input mutations.')
pCepip.input               = "infile:file"
pCepip.output              = "outfile:file:{{i.infile | fn}}.cepip.txt"
pCepip.args.cepip          = params.cepip.value
pCepip.args.cell           = ""
pCepip.args.params         = Box()
pCepip.lang                = params.python.value
pCepip.script              = "file:scripts/vcfnext/pCepip.py"

"""
@name:
	pMutSig
@description:
	MutSig stands for "Mutation Significance".  MutSig analyzes lists of mutations discovered in DNA sequencing, to identify genes that were mutated more often than expected by chance given background mutation processes.
	For more information, see Lawrence, M. et al. Mutational heterogeneity in cancer and the search for new cancer-associated genes. Nature 499, 214-218 (2013).
	
	See [dcumentation](http://archive.broadinstitute.org/cancer/cga/mutsig_run)
@input:
	`infile:file`: mutation table
@output:
	`outdir:dir`: The output directory
@args:
	`mutsig` : The path to `run_MutSigCV.sh`, default: 'mutsig'
	`mcr`    : The Matlab MCR path
	`cvrg`   : coverage table
	`cvrt`   : covariates table
	`mutdict`: mutation_type_dictionary_file
	`chrdir` : chr_files_hg18 or chr_files_hg19
@requires:
	[MutSig](http://archive.broadinstitute.org/cancer/cga/mutsig_download)
"""
pMutSig              = Proc(desc = 'Run MutSig.')
pMutSig.input        = 'infile:file'
pMutSig.output       = "outdir:dir:{{i.infile | fn}}.mutsig"
pMutSig.args.cvrg    = params.mutsig_cvrg.value
pMutSig.args.cvrt    = params.mutsig_cvrt.value
pMutSig.args.mutdict = params.mutsig_mutdict.value
pMutSig.args.chrdir  = params.mutsig_chrdir.value
pMutSig.args.mutsig  = params.mutsig.value
pMutSig.args.mcr     = params.mcr.value
pMutSig.script       = "file:scripts/vcfnext/pMutSig.bash"

"""
@name:
	pMafLiftover
@description:
	Liftover maf file from one assembly to another
@input:
	`infile:file`: The input maf file
@output:
	`outfile:file`: The output maf file
@args:
	`liftover`: The liftOver program.
	`lochain`:  The liftOver chain file.
	`genome`:   The target genome.
@requires:
	liftOver from UCSC
"""
pMafLiftover               = Proc(desc = 'Liftover a maf file from one assembly to another')
pMafLiftover.input         = 'infile:file'
pMafLiftover.output        = 'outfile:file:{{i.infile | fn | lambda x: x if x.endswith(".maf") else x + ".maf"}}'
pMafLiftover.args.liftover = params.liftover.value
pMafLiftover.args.lochain  = params.lochai.value
pMafLiftover.args.genome   = params.genome.value
pMafLiftover.lang          = params.python.value
pMafLiftover.script        = "file:scripts/vcfnext/pMafLiftOver.py"


"""
@name:
	pMafMerge
@description:
	Merge maf files.
@input:
	`infiles:files`: The maf files
@output:
	`outfile:file`: The merged maf file
@args:
	`excols`: How to deal with extra columns other than 34 standard columns from TCGA.
		- merge(default): Merge the columns, if one not exists, fill with an empty string.
		- discard: Just discard the extra columns, with only 34 columns left. So you can also put just one maf file in the indir with some columns missed to fill it with standard columns.
"""
pMafMerge              = Proc(desc = 'Merge maf files.')
pMafMerge.input        = 'infiles:files'
pMafMerge.output       = 'outfile:file:{{i.infiles | fs2name}}.maf'
pMafMerge.args.excols  = 'merge' # discard
pMafMerge.envs.fs2name = fs2name
pMafMerge.lang         = params.python.value
pMafMerge.script       = "file:scripts/vcfnext/pMafMerge.py"

"""
@name:
	pMaf2Mat
@description:
	Convert maf file to a gene(row)-sample(column) matrix
@input:
	`infile:file`: The input file
@output:
	`outfile:file`: The output matrix
@args:
	`mutypes`: Provide manual list of variant classifications to be counted, only effective when `args.binary = False`. Default: `None` (all counted)
	`binary` : Just generate a binary matrix instead of a count matrix. Default: `False`
	`na`: What value to use for no mutations reported on a gene. Default: `0`
	`samfn`  : A function (in r) to transform the sample names. Default: `function(sample) sample`
"""
pMaf2Mat              = Proc(desc = 'Convert maf file to a gene-based mutation matrix')
pMaf2Mat.input        = 'infile:file'
pMaf2Mat.output       = 'outfile:file:{{i.infile | fn}}.mat.txt'
pMaf2Mat.args.binary  = False
pMaf2Mat.args.mutypes = None
pMaf2Mat.args.na      = 0
pMaf2Mat.args.samfn   = 'function(sample) sample'
pMaf2Mat.lang         = params.Rscript.value
pMaf2Mat.script       = "file:scripts/vcfnext/pMaf2Mat.r"

"""
@name:
	pMaftools
@description:
	Use maftools to draw plots.
@input:
	`indir:dir`: The input directory. Could contain:
		- `*.maf` or `*.maf.gz` file (required)
		- `*.annot.tsv` or `*.annot.txt` file (see: https://github.com/PoisonAlien/maftools/blob/master/inst/extdata/tcga_laml_annot.tsv)
		- `all_lesions.conf_*.txt`: Gistic cnv data
		- `amp_genes.conf_*.txt`: Gistic cnv data
		- `del_genes.conf_*.txt`: Gistic cnv data
		- `scores.gistic`: Gistic cnv data
		- `*.seg.txt`: CBS segments data
		- `*sig_genes.txt` or `*sig_genes.txt.gz`: Mutsig results, to do pancancer somparison.
@output:
	`outdir:dir`: The output directory
@args:
	`ngenes` : Top number of genes to plot for some plots. Default: `10`
	`mutypes`: Provide manual list of variant classifications to be considered as non-synonymous. Rest will be considered as silent variants. Default: `["Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site","Nonsense_Mutation", "Nonstop_Mutation", "In_Frame_Del","In_Frame_Ins", "Missense_Mutation"]`
	`isTCGA`:  If the maf file is from TCGA? Default: `False`
	`ref`   :  The reference file for signature plot.
	`plot`  :  Which plots to plot. 
		- Default:
		  ```python
		  Box(
		  	summary        = True,
		  	oncoplot       = True,
		  	oncostrip      = True,
		  	titv           = True,
		  	lollipop       = True,
		  	cbsseg         = True,
		  	rainfall       = True,
		  	tcgacomp       = True,
		  	vaf            = True,
		  	genecloud      = True,
		  	gisticGenome   = True,
		  	gisticBubble   = True,
		  	gisticOncoplot = True,
		  	somInteraction = True,
		  	oncodrive      = True,
		  	pfam           = True,
		  	pancan         = True,
		  	survival       = True,
		  	heterogeneity  = True,
		  	signature      = True,
		  )
		  ```
	`params`:  The extra parameters for each plot function.
		- Default:
		  ```python
		  Box(
		  	summary        = Box(rmOutlier = True, addStat = 'median', dashboard = True),
		  	oncoplot       = Box(),
		  	oncostrip      = Box(),
		  	titv           = Box(),
		  	lollipop       = Box(AACol = 'Protein_Change'),
		  	cbsseg         = Box(labelAll = True),
		  	rainfall       = Box(detectChangePoints = True),
		  	tcgacomp       = Box(),
		  	vaf            = Box(flip = True),
		  	genecloud      = Box(minMut = 3),
		  	gisticGenome   = Box(markBands = 'all'),
		  	gisticBubble   = Box(),
		  	gisticOncoplot = Box(),
		  	somInteraction = Box(),
		  	oncodrive      = Box(AACol = 'Protein_Change', minMut = 5, pvalMethod = 'zscore', fdrCutOff = 0.1, useFraction = True),
		  	pfam           = Box(AACol = 'Protein_Change'),
		  	pancan         = Box(qval = 0.1, label = 1, normSampleSize = True),
		  	survival       = Box(),
		  	heterogeneity  = Box(),
		  	signature      = Box(nTry = 6, plotBestFitRes = False),
		  )
		  ```
	`devpars`:  The parameters for plot device. Default: `Box(res = 300, height = 2000, width = 2000)`
	`nthread`:  Number of threads used for multiple plot of one type. Default: `1`
@requires:
	[Maftools](https://bioconductor.org/packages/devel/bioc/vignettes/maftools/inst/doc/maftools.html)
"""
pMaftools              = Proc(desc = 'Use maftools to draw plots.')
pMaftools.input        = 'indir:dir'
pMaftools.output       = 'outdir:dir:{{i.indir | fn}}.maftools'
pMaftools.args.ngenes  = 10
pMaftools.args.isTCGA  = False
pMaftools.args.ref     = params.ref.value # for signature
pMaftools.args.mutypes = ["Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site","Nonsense_Mutation", "Nonstop_Mutation", "In_Frame_Del","In_Frame_Ins", "Missense_Mutation"]
pMaftools.args.plot    = Box(
	summary        = True,
	oncoplot       = True,
	oncostrip      = True,
	titv           = True,
	lollipop       = True,
	cbsseg         = True,
	rainfall       = True,
	tcgacomp       = True,
	vaf            = True,
	genecloud      = True,
	gisticGenome   = True,
	gisticBubble   = True,
	gisticOncoplot = True,
	somInteraction = True,
	oncodrive      = True,
	pfam           = True,
	pancan         = True,
	survival       = True,
	heterogeneity  = True,
	signature      = True,
)
pMaftools.args.params  = Box(
	summary        = Box(rmOutlier = True, addStat = 'median', dashboard = True),
	oncoplot       = Box(),
	oncostrip      = Box(),
	titv           = Box(),
	lollipop       = Box(),
	cbsseg         = Box(labelAll = True),
	rainfall       = Box(detectChangePoints = True),
	tcgacomp       = Box(),
	vaf            = Box(flip = True),
	genecloud      = Box(minMut = 3),
	gisticGenome   = Box(markBands = 'all'),
	gisticBubble   = Box(),
	gisticOncoplot = Box(),
	somInteraction = Box(),
	oncodrive      = Box(minMut = 5, pvalMethod = 'zscore', fdrCutOff = 0.1, useFraction = True),
	pfam           = Box(),
	pancan         = Box(qval = 0.1, label = 1, normSampleSize = True),
	survival       = Box(),
	heterogeneity  = Box(),
	signature      = Box(nTry = 6, plotBestFitRes = False),
)
pMaftools.args.devpars = Box(res = 300, height = 2000, width = 2000)
pMaftools.args.nthread = 1
pMaftools.envs.rimport = rimport
pMaftools.lang         = params.Rscript.value
pMaftools.script       = "file:scripts/vcfnext/pMaftools.r"

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
	`snpEff`:       The snpEff executable, default: "snpEff"
	`params`:    Other parameters for `snpEff`, default: "-Xms1g -Xmx4g -v"
	`genome`:    The genome used for annotation, default: "hg19"
	`informat`:  The format of input file [vcf or bed], default: "vcf"
	`outformat`: The format of output file [vcf, gatk, bed, bedAnn], default: "vcf"
	`csvStats`:  Whether to generate csv stats file, default: True.
	`htmlStats`: Whether to generate the html summary file, default: False.
	`javamem`:   The memory to use. Default: '-Xms1g -Xmx8g'
@requires:
	[snpEff](http://snpeff.sourceforge.net/SnpEff_manual.html)
"""
pSnpEff = Proc()
pSnpEff.input  = "infile:file"
pSnpEff.output = "outdir:dir:{{infile | fn}}.snpeff"
pSnpEff.args   = { "snpEff": "snpEff", "javamem": "-Xms1g -Xmx8g", "genome": "hg19", "informat": "vcf", "outformat": "vcf", "csvStats": True, "htmlStats": False, "params": "" }
pSnpEff.script = """
csvfile="{{outdir}}/{{infile | fn}}.csvstat"
sumfile="{{outdir}}/{{infile | fn}}.html"
outfile="{{outdir}}/{{infile | fn}}.snpEff.vcf"
csvStats=""
if [[ "{{args.csvStats}}" == "True" ]]; then
	csvStats="-csvStats \\"$csvfile\\""
fi
stats=""
if [[ "{{args.htmlStats}}" == "True" ]]; then
	stats="-stats \\"$sumfile\\""
fi
echo {{args.snpEff}} {{args.javamem}} -i {{args.informat}} -o {{args.outformat}} $csvStats $stats {{args.params}} {{args.genome}} "{{infile}}"
{{args.snpEff}} {{args.javamem}} -i {{args.informat}} -o {{args.outformat}} $csvStats $stats {{args.params}} {{args.genome}} "{{infile}}" > "$outfile"
"""
