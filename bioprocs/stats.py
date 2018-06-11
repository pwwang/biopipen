from pyppl import Proc, Box
from . import params, rimport

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
pMetaPval.args.inopts       = Box(
	cnames = True,
	pcol   = -1
)
pMetaPval.args.outopts      = Box(
	ponly  = False,
	head   = True
)
pMetaPval.args.method       = 'sumlog' # fisher's method
pMetaPval.envs.rimport      = rimport
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
pMetaPval1.args.inopts       = Box(
	cnames = True,
	pcol   = -1
)
pMetaPval1.args.outopts      = Box(
	ponly  = False,
	head   = True
)
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
		- col1: rownames if args.inopts.rnames = True
		- col2: the survival time
		- col3: the status. 0/1 for alive/dead or 1/2 for alive dead
		- col4: var1.
		- ... other variables
@output:
	`outfile:file`: The outfile containing the pvalues
	`outdir:dir`  : The output directory containing the pval files and plots
@args:
	`inunit`    : The time unit in input file. Default: days
	`outunit`   : The output unit for plots. Default: days
	`nthread`   : Number of threads used to perform analysis for groups. Default: 1
	`inopts`    : The options for input file
		- `rnames`: Whether input file has row names. Default: True
	`combine`   : Whether combine groups in the same plot. Default: True
	`devpars`   : The device parameters for png. Default: `{res:300, height:2000, width:2000}`
		- The height and width are for each survival plot. If args.combine is True, the width and height will be multiplied by `max(arrange.ncol, arrange.nrow)`
	`covfile`   : The covariant file. Require rownames in both this file and input file.
	`ngroups`   : Number of curves to plot (the continuous number will divided into `ngroups` groups.
	`plot`      : The params for plot.
		- `params` : The params for `ggsurvplot`. Default: `Box({'risk.table': True, 'conf.int': True, 'font.legend': 13, 'pval': '{method}\np = {pval}'})`
			- You may do `ylim.min` to set the min ylim. Or you can set it as 'auto'. Default: 0. 
		- `arrange`: How to arrange multiple survival plots in one if `args.combine = True`.
			- `nrow`: The number of rows. Default: 1
			- `ncol`: The number of cols. Default: 1
	`ggs`       : Extra ggplot2 elements for main plot. `ggs.table` is for the risk table.
	`pval`      : The method to calculate the pvalue shown on the plot. Default: True (logrank)
		- Could also be `waldtest`, `likeratio` (Likelihoold ratio test)
	`method`    : The method to do survival analysis. 
@requires:
	[`r-survival`](https://rdrr.io/cran/survival/)
	[`r-survminer`](https://rdrr.io/cran/survminer/)
"""
pSurvival                 = Proc(desc = "Survival analysis.")
pSurvival.input           = 'infile:file'
pSurvival.output          = [
	'outfile:file:{{in.infile | fn2}}.dir/{{in.infile | fn2}}.survival.txt', 
	'outdir:dir:{{in.infile | fn2}}.dir'
]
pSurvival.args.inunit     = 'days' # months, weeks, years
pSurvival.args.outunit    = 'days'
pSurvival.args.method     = 'cox' # tm or auto 
pSurvival.args.covfile    = None
pSurvival.args.nthread    = 1
pSurvival.args.inopts     = Box(rnames = True)
pSurvival.args.combine    = False
pSurvival.args.devpars    = Box(res = 300, height = 2000, width = 2000)
pSurvival.args.ngroups    = 2 # how many curves to plot, typically 2. The values will divided into <ngroups> groups for the var
pSurvival.args.autogroup  = True # False to use median, else find the best binary split spot, only applicable when args.ngroup = 2
pSurvival.args.plot = Box(
	params  = Box({'font.legend': 13, 'pval': '{method}\np = {pval}', 'risk.table': True}), # params for ggsurvplot
	arrange = Box() # params for arrange_ggsurvplots if args.combine = T. Typically nrow or ncol is set. If args.plot.arrange.ncol = 3, that means {ncol: 3, nrow: 1}. If ncol is not set, then it defaults to 1.
)
pSurvival.args.ggs        = Box(table = Box())
pSurvival.args.pval       = True # 'logrank', 'waldtest', 'likeratio'
pSurvival.envs.rimport    = rimport
pSurvival.lang            = params.Rscript.value
pSurvival.script          = "file:scripts/stats/pSurvival.r"

"""
@name:
	pChiSquare
@description:
	Do chi-square test.
@input:
	`infile:file`: The input file.
@output:
	`outfile:file` : The output file containing Xsquare, df, pval and method
	`obsvfile:file`: The observation matrix
	`exptfile:file`: The expectation matrix
@args:
	`intype`: The type of the input file:
		- `count` (default): The contingency table
		```
		#         | Disease | Healthy |
		# --------+---------+---------+
		#   mut   |   40    |   12    |
		# non-mut |   23    |   98    |
		# --------+---------+---------+
		```
		- `raw`: The raw values:
		```
		# Contingency table rows: Mut, Non
		# Contingency table cols: Disease, Healthy
		#
		#         | S1 | S2 | ... | Sn |
		# --------+----+----+-----+----+
		# Disease | 1  | 0  | ... | 1  |
		# Healthy | 0  | 1  | ... | 0  |
		# --------+----+----+-----+----+
		# Mut     | 1  | 0  | ... | 1  |
		# Non     | 0  | 1  | ... | 0  |
		```
	`ctcols`: The colnames of contingency table if input file is raw values
		- You may also specify them in the head of the input file
"""
pChiSquare = Proc(desc = "Do chi-square test.")
pChiSquare.input = "infile:file"
pChiSquare.output = "outfile:file:{{in.infile | fn2}}.chi2.txt, obsvfile:file:{{in.infile | fn2}}.obsv.txt, exptfile:file:{{in.infile | fn2}}.expt.txt"
pChiSquare.args.intype = 'cont' # raw
pChiSquare.args.ctcols = ''
pChiSquare.lang = params.Rscript.value
pChiSquare.script = "file:scripts/stats/pChiSquare.r"

"""
@name:
	pFisherExact
@description:
	Do fisher exact test.
@input:
	`infile:file`: The input file.
@output:
	`outfile:file` : The output file containing confInt1, confInt2, oddsRatio, pval, alternative and method.
@args:
	`intype`: The type of the input file:
		- `count` (default): The contingency table
		```
		#         | Disease | Healthy |
		# --------+---------+---------+
		#   mut   |   40    |   12    |
		# non-mut |   23    |   98    |
		# --------+---------+---------+
		```
		- `raw`: The raw values:
		```
		# Contingency table rows: Mut, Non
		# Contingency table cols: Disease, Healthy
		#
		#         | S1 | S2 | ... | Sn |
		# --------+----+----+-----+----+
		# Disease | 1  | 0  | ... | 1  |
		# Healthy | 0  | 1  | ... | 0  |
		# --------+----+----+-----+----+
		# Mut     | 1  | 0  | ... | 1  |
		# Non     | 0  | 1  | ... | 0  |
		```
	`ctcols`: The colnames of contingency table if input file is raw values
		- You may also specify them in the head of the input file
"""
pFisherExact = Proc(desc = "Do fisher exact test.")
pFisherExact.input = "infile:file"
pFisherExact.output = "outfile:file:{{in.infile | fn2}}.fexact.txt"
pFisherExact.args.intype = 'cont' # raw
pFisherExact.args.ctcols = ''
pFisherExact.lang = params.Rscript.value
pFisherExact.script = "file:scripts/stats/pFisherExact.r"

"""
@name:
	pPWFisherExact
@description:
	Do pair-wise fisher exact test.
	Commonly used for co-occurrence/mutual-exclusivity analysis.
	P-value indicates if the pairs are significantly co-occurred or mutually exclusive.
	Co-occurrence: Odds ratio > 1
	Mutual-exclusivity: Odds ratio < 1
@input:
	`infile:file`: The input file.
@output:
	`outfile:file` : The output file containing confInt1, confInt2, oddsRatio, pval, qval, alternative and method.
@args:
	`intype`: The type of the input file:
		- `pairs`: The contingency table
		```
		#
		# A+	B+	4
		# A-	B-	175
		# A+	B-	12
		# A-	B+	1
		#
		```
		- `raw` (default): The raw values:
		```
		#
		#         | S1 | S2 | ... | Sn |
		# --------+----+----+-----+----+
		# A       | 1  | 0  | ... | 1  |
		# B       | 0  | 1  | ... | 0  |
		# ...     |           ...      |
		# X       | 0  | 1  | ... | 0  |
		# --------+----+----+-----+----+
		#
		```
	`padj`: The p-value adjustment method, see `p.adjust.methods` in R. Default: `BH`
"""
pPWFisherExact              = Proc(desc = "Do pair-wise fisher exact test.")
pPWFisherExact.input        = "infile:file"
pPWFisherExact.output       = "outfile:file:{{in.infile | fn2}}.pwfexact.txt"
pPWFisherExact.args.intype  = 'raw' # pairs
pPWFisherExact.args.padj    = 'BH'
pPWFisherExact.envs.rimport = rimport
pPWFisherExact.lang         = params.Rscript.value
pPWFisherExact.script       = "file:scripts/stats/pPWFisherExact.r"

"""
@name:
	pMediation
@description:
	Do mediation analysis
@input:
	`infile:file`: The input file (a matrix or data.frame).
@output:
	`outfile:file`: The result file.
@args:
	`inopts`: The options for input file.
		- `cnames`: Whether the input file has column names
		- `rnames`: Whether the input file has row names
	`medopts`: The options for mediation analysis.
		- `modelm`: The model for M ~ X. Default: `lm(M ~ X)`
		- `modely`: The model for Y ~ X + M. Default: `lm(Y ~ X + M)`
		- `mediator`: Tell the model which column is the mediator
		- `treat`: Tell the model which column is the variable
		- `boot`: Use bootstrap?
		- `sims`: How many time simulations?
"""
pMediation = Proc(desc = "Do mediation analysis.")
pMediation.input  = 'infile:file'
pMediation.output = 'outfile:file:{{in.infile | fn2}}.mediation.txt'
pMediation.args.inopts = Box(
	cnames   = True,
	rnames   = True
)
pMediation.args.medopts = Box(
	modelm   = 'lm(M ~ X)',
	modely   = 'lm(Y ~ X + M)',
	mediator = 'M',
	treat    = 'X',
	boot     = True,
	sims     = 500,
)
pMediation.lang = params.Rscript.value
pMediation.script = "file:scripts/stats/pMediation.r"

"""
@name:
	pHypergeom
@description:
	Do hypergeometric test.
@input:
	`infile:file`: The input file, could be raw data (presence (1) and absence (0) of elements) or number of overlapped elements and elements in each category.
		- Set `args.intype` as `raw` if it is raw data. The population size `args.N` is required
		- Set `args.intype` as `numbers` (or any string except `raw`) if it is numbers. You can specified explicit header: `k` = overlapped elements, `m` = size of set 1, `n` = size of set 2 and `N` = the population size. If `N` not included, then `args.N` is required
@output:
	`outfile:file`: The output file
@args:
	`intype`: the type of input file. Default: `raw`. See `infile:file`
	`inopts`: The options for input file.
		- `cnames`: Whether the input file has column names
		- `rnames`: Whether the input file has row names
	`N`: The population size. Default: `None`
"""
pHypergeom             = Proc(desc = "Do hypergeometric test.")
pHypergeom.input       = 'infile:file'
pHypergeom.output      = 'outfile:file:{{in.infile | fn2}}.hypergeom.txt'
pHypergeom.args.intype = 'raw' # numbers
pHypergeom.args.inopts = Box(
	cnames = True,
	rnames = True
)
pHypergeom.args.N       = None
pHypergeom.envs.rimport = rimport
pHypergeom.lang         = params.Rscript.value
pHypergeom.script       = "file:scripts/stats/pHypergeom.r"
