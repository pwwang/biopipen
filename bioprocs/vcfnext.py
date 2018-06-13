from pyppl import Proc, Box
#from .utils import runcmd, plot, helpers
from . import params


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
pVcfStatsPlot.output           = "outdir:dir:{{in.indir | fn}}-{{job.index}}.statplots"
pVcfStatsPlot.args.chroms      = ""
pVcfStatsPlot.args.Rscript     = params.Rscript.value
pVcfStatsPlot.args.devpars     = Box({'res':300, 'width':2000, 'height':2000})
pVcfStatsPlot.args.histplotggs = []
pVcfStatsPlot.args.boxplotggs  = []
#pVcfStatsPlot.envs.runcmd      = runcmd.py
#pVcfStatsPlot.envs.plotHist    = plot.hist.r
#pVcfStatsPlot.envs.plotBoxplot = plot.boxplot.r
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
pCallRate.output           = "outdir:dir:{{ in.indir | fn }}.callrate"
pCallRate.args.histplotggs = []
pCallRate.args.devpars     = Box({'res':300, 'width':2000, 'height':2000})
#pCallRate.envs.runcmd      = runcmd.r
#pCallRate.envs.cbindfill   = helpers.cbindfill.r
#pCallRate.envs.plotHist    = plot.hist.r
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
pCepip.output              = "outfile:file:{{in.infile | fn}}.cepip.txt"
pCepip.args.cepip          = params.cepip.value
pCepip.args.cell           = ""
pCepip.args.params         = Box()
#pCepip.envs.runcmd         = runcmd.py
#pCepip.envs.params2CmdArgs = helpers.params2CmdArgs.py
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
pMutSig.output       = "outdir:dir:{{in.infile | fn}}.mutsig"
pMutSig.args.cvrg    = params.mutsig_cvrg.value
pMutSig.args.cvrt    = params.mutsig_cvrt.value
pMutSig.args.mutdict = params.mutsig_mutdict.value
pMutSig.args.chrdir  = params.mutsig_chrdir.value
pMutSig.args.mutsig  = params.mutsig.value
pMutSig.args.mcr     = params.mcr.value
pMutSig.script       = "file:scripts/vcfnext/pMutSig.bash"

"""
@name:
	pMafMerge
@description:
	Merge maf files.
@input:
	`indir:dir`: The directory containing the maf files
@output:
	`outfile:file`: The merged maf file
@args:
	`excols`: How to deal with extra columns other than 34 standard columns from TCGA.
		- merge(default): Merge the columns, if one not exists, fill with an empty string.
		- discard: Just discard the extra columns, with only 34 columns left. So you can also put just one maf file in the indir with some columns missed to fill it with standard columns.
"""
pMafMerge             = Proc(desc = 'Merge maf files.')
pMafMerge.input       = 'indir:dir'
pMafMerge.output      = 'outfile:file:{{in.indir | fn}}.maf'
pMafMerge.args.excols = 'merge' # discard
pMafMerge.lang        = params.python.value
pMafMerge.script      = "file:scripts/vcfnext/pMafMerge.py"

"""
@name:
	pMaftools
@description:
	Use maftools to draw plots.
@args:
	`ngenes`: 
@requires:
	[``]
"""
pMaftools              = Proc(desc = 'Use maftools to draw plots.')
pMaftools.input        = 'indir:dir'
pMaftools.output       = 'outdir:dir:{{in.indir | fn}}.maftools'
pMaftools.args.ngenes  = 10
pMaftools.args.isTCGA  = False
pMaftools.args.ref     = params.ref.value # for signature
pMaftools.args.mutypes = ["Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site","Nonsense_Mutation", "Nonstop_Mutation", "In_Frame_Del","In_Frame_Ins", "Missense_Mutation"]
pMaftools.args.plots   = Box(
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
pMaftools.args.devpars = Box(res = 300, height = 2000, width = 2000)
pMaftools.args.nthread = 1
#pMaftools.envs.heatmap = plot.heatmap.r
pMaftools.lang         = params.Rscript.value
pMaftools.script       = "file:scripts/vcfnext/pMaftools.r"
