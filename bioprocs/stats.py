from pyppl import Proc, Box
from . import params
from .utils import helpers

"""
@name:
	pMetaPval
@description:
	Combine p-values in the files from input directory
@input:
	`indir:dir`: The directory containing the input files
@output:
	`outfile:file`: The output file containing the meta-pvalues
@args:
	`args.pattern`: The pattern used to filter the input files. Default: '*'
	`args.header`: Whether the input files contains a header. Default: True
		- Could be a list to specify it for each file.
		- The order should be concordant with the file names
	`args.pcol`: Which column is the p-value. Default: -1 (last column)
	`args.poutonly`: Only output pvalues. Default: False (output all possible information)
	`args.outheader`: Whether output the header. Default: True
	`args.method`: The method used to calculate the meta-pvalue. Default: sumlog (Fisher's method)
		- Other available methods: logitp, sumz, votep, sump, meanp and wilkinsonp
		- See: https://www.rdocumentation.org/packages/metap/versions/0.8
@requires:
	[`r-matep`](https://www.rdocumentation.org/packages/metap/)
"""
pMetaPval                   = Proc(desc = "Combine p-values.")
pMetaPval.input             = "indir:dir"
pMetaPval.output            = "outfile:file:{{in.indir | fn}}.metapval.txt"
pMetaPval.args.pattern      = "*"
pMetaPval.args.header       = True
pMetaPval.args.pcol         = -1
pMetaPval.args.poutonly     = False
pMetaPval.args.outheader    = True
pMetaPval.args.method       = 'sumlog' # fisher's method
pMetaPval.tplenvs.cbindfill = helpers.cbindfill.r
pMetaPval.lang              = params.Rscript.value
pMetaPval.script            = "file:scripts/stats/pMetaPval.r"

"""
@name:
	pMetaPval1
@description:
	Combine p-values in a single file by rows.
@input:
	`infile:file`: The input file
@output:
	`outfile:file`: The output file containing the meta-pvalues
@args:
	`args.header`: Whether the input files contains a header. Default: True
	`args.pcol`: Which column is the p-value. Default: -1 (last column)
	`args.poutonly`: Only output pvalues. Default: False (output all possible information)
	`args.outheader`: Whether output the header. Default: True
	`args.method`: The method used to calculate the meta-pvalue. Default: sumlog (Fisher's method)
		- Other available methods: logitp, sumz, votep, sump, meanp and wilkinsonp
		- See: https://www.rdocumentation.org/packages/metap/versions/0.8
@requires:
	[`r-matep`](https://www.rdocumentation.org/packages/metap/)
"""
pMetaPval1 = Proc(desc = "Combine p-values in a single file by rows.")
pMetaPval1.input             = "infile:file"
pMetaPval1.output            = "outfile:file:{{in.infile | fn}}.metapval.txt"
pMetaPval1.args.header       = True
pMetaPval1.args.pcol         = -1
pMetaPval1.args.poutonly     = False
pMetaPval1.args.outheader    = True
pMetaPval1.args.method       = 'sumlog' # fisher's method
pMetaPval1.lang              = params.Rscript.value
pMetaPval1.script            = "file:scripts/stats/pMetaPval1.r"

"""
@name:
	pSurvival
@description:
	Survival analysis
@input:
	`infile:file`: The input file (header is required).
		- col1: rownames if args.rnames = True
		- col2: the survival time
		- col3: the status. 0/1 for alive/dead or 1/2 for alive dead
		- col4: group1.
		- ... other groups
@output:
	`outdir:dir`: The output directory containing the pval files and plots
@args:
	`inunit`    : The time unit in input file. Default: days
	`outunit`   : The output unit for plots. Default: days
	`nthread`   : Number of threads used to perform analysis for groups. Default: 1
	`rnames`    : Whether input file has row names. Default: True
	`combine`   : Whether combine groups in the same plot. Default: True
	`devpars`   : The device parameters for png. Default: `{res:300, height:2000, width:2000}`
		- The height and width are for each survival plot. If args.combine is True, the width and height will be multiplied by `max(gridParams.ncol, gridParams.nrow)`
	`plotParams`: The parameters for `ggsurvplot`. Default: `{risk.table: True, conf.int = True}`
	`gridParams`: The parameters for `arrange_ggsurvplots`.
	`pval`      : Whether print pvalue on the plot. Default: True
@requires:
	[`r-survival`](https://rdrr.io/cran/survival/)
	[`r-survminer`](https://rdrr.io/cran/survminer/)
"""
pSurvival                 = Proc(desc = "Survival analysis.")
pSurvival.input           = 'infile:file'
pSurvival.output          = 'outdir:dir:{{in.infile | fn}}.survival'
pSurvival.args.inunit     = 'days' # months, weeks, years
pSurvival.args.outunit    = 'days'
pSurvival.args.nthread    = 1
pSurvival.args.rnames     = True
pSurvival.args.combine    = True
pSurvival.args.devpars    = Box(res = 300, height = 2000, width = 2000)
pSurvival.args.plotParams = Box({'risk.table': True, 'conf.int': True})
pSurvival.args.gridParams = Box() # ncol, nrow
pSurvival.args.pval       = True
pSurvival.lang            = params.Rscript.value
pSurvival.script          = "file:scripts/stats/pSurvival.r"
