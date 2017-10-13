from pyppl import Proc, Box
from .utils import runcmd, plot, helpers
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
pVcfStatsPlot                     = Proc(desc = 'Convert csvstat file from snpEff to R-readable matrix and plot them.')
pVcfStatsPlot.input               = "indir:file"
pVcfStatsPlot.output              = "outdir:dir:{{in.indir | fn}}-{{job.index}}.statplots"
pVcfStatsPlot.args.chroms         = ""
pVcfStatsPlot.args.Rscript        = params.Rscript.value
pVcfStatsPlot.args.devpars        = Box({'res':300, 'width':2000, 'height':2000})
pVcfStatsPlot.args.histplotggs    = []
pVcfStatsPlot.args.boxplotggs     = []
pVcfStatsPlot.tplenvs.runcmd      = runcmd.py
pVcfStatsPlot.tplenvs.plotHist    = plot.hist.r
pVcfStatsPlot.tplenvs.plotBoxplot = plot.boxplot.r
pVcfStatsPlot.lang                = params.python.value
pVcfStatsPlot.script              = "file:scripts/vcfnext/pVcfStatsPlot.py"

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
pCallRate.tplenvs.runcmd    = runcmd.r
pCallRate.tplenvs.cbindfill = helpers.cbindfill.r
pCallRate.tplenvs.plotHist  = plot.hist.r
pCallRate.lang              = params.Rscript.value
pCallRate.script            = "file:scripts/vcfnext/pCallRate.r"

