"""Some statistic processes"""
from pyppl import Proc, Box
from . import params, rimport
from . import delefactory, procfactory
from modkit import Modkit
Modkit().delegate(delefactory())

@procfactory
def _pStats():
	"""
	@name:
		pStats
	"""
	return Box(
		desc   = 'Data statistics',
		lang   = params.Rscript.value,
		input  = 'infile:file',
		output = 'outdir:dir:{{i.infile | stem}}.stats',
		args   = Box(inopts  = Box(rnames = True, cnames = True, delimit = "\t"),
			types   = Box(),
			groups  = [],
			ignore  = [],
			devpars = Box(res = 300, height = 2000, width = 2000))
	)

@procfactory
def _pMetaPval():
	"""
	@name:
		pMetaPval
	@description:
		Calculate a meta-pvalue using different methods
	@input:
		`infile:file`: The infile containing multiple pvalues for each entries.
			- Could be two types (see `args.intype`)
			- `matrix`: A matrix with rows as entries (rownames are optional), columns as cases
				```
						Case1   Case2   ...  CaseN
				Entry1  1.33e-2 NA      ...  1.77e-10
				Entry2  2.66e-2 4.22e-5 ...  1.71e-3
				... ...
				EntryM  NA      0.00013 ...  4.11e-3
				```
			- `melt`: Rows are entries from each case
			`args.rnames` should be `False`, but entry names are required in the first column,
			column names are optional (have to be 3 columns)
				```
				Entry   Pvalue   Case
				Entry1  1.33e-2  Case1
				Entry1  1.77e-10 CaseN
				Entry2  2.66e-2  Case1
				Entry2  4.22e-5  Case2
				Entry2  1.71e-3  CaseN
				... ...
				EntryM  0.00013  Case2
				EntryM  4.11e-3  CaseN
				```
	@output:
		`outfile:file`: The output file containing the meta-pvalues. Default: `{{i.infile | fn}}.meta{{i.infile | ext}}`
	@args:
		`intype`: The type of the input file. Default: `matrix` (see `i.infile`)
		`inopts`: The input options to read the input file. Default: `Box(rnames = True, cnames = True)`
		`method`: The method used to calculate the meta-pvalue. Default: sumlog (Fisher's method, aka `fisher`)
			- Other available methods: logitp, sumz, votep, sump, meanp and wilkinsonp
			- See: https://www.rdocumentation.org/packages/metap/versions/0.8
		`na`    : How to deal with `NA` p-values. Default: `skip` (just don't count it)
			- Or a numeric value to replace it with (e.g.: `1`).
	@requires:
		[`r-matep`](https://www.rdocumentation.org/packages/metap/)
	"""
	pMetaPval              = Proc(desc = 'Calculate meta p-values.')
	pMetaPval.input        = 'infile:file'
	pMetaPval.output       = 'outfile:file:{{i.infile | fn}}.meta{{i.infile | ext}}'
	pMetaPval.args.intype  = 'matrix'
	pMetaPval.args.inopts  = Box(rnames = True, cnames = True)
	pMetaPval.args.method  = 'sumlog'
	pMetaPval.args.na      = 'skip' #0, 1
	pMetaPval.envs.rimport = rimport
	pMetaPval.lang         = params.Rscript.value
	pMetaPval.script       = "file:scripts/stats/pMetaPval.r"
	return pMetaPval

@procfactory
def _pSurvival():
	"""
	@input:
		infile: The input file (header is required).
			- col1: rownames if args.inopts.rnames = True
			- col2: the survival time
			- col3: the status. 0/1 for alive/dead or 1/2 for alive dead
			- col4: var1.
			- ... other variables
			- Note that covariales should be put in covfile, each variable in this file will be analyzed separately.
		covfile: The covariate file. If not provided, use `args.covfile`
	@output:
		outfile: The outfile containing the pvalues and other details.
		outdir : The output directory containing the pval files and plots, including:
			- `outfile`
			- cutplot: `<feature>.cut.png`, only available when `args.cutting = 'maxstat'`
			- forest plot: `<feature>.forest.png`, only available using cox
			- data with subpopulations: `<feature>.data.txt`
	@args:
		`inunit`    : The time unit in input file. Default: days
		`outunit`   : The output unit for plots. Default: days
		`inopts`    : The options for input file
			- `rnames`: Whether input file has row names.
		`combine`   : Whether combine variables in the same plot.
			- `nrow`: The number of rows. Default: 1
			- `ncol`: The number of cols. Default: 1
		`devpars`   : The device parameters for png.
			- The height and width are for each survival plot.
			- If args.combine is True, the width and height will be multiplied by `max(combine.ncol, combine.nrow)`
		`covfile`   : The covariant file. Require rownames in both this file and input file.
		`ggs`       : Extra ggplot2 elements for main plot. `ggs$table` is for the risk table.
		`params`    : The params for `ggsurvplot`. Extra arguments:
			- cut: Method to cut the continous value into subpopulations. One of mean, median, q25, q75, asis and maxstat (default).
			- cutlabels: Labels for the subpopulations (from small to large).
		`method`    : The method to do survival analysis. km (Kaplan Merier) or cox.
	@requires:
		[`r-survival`](https://rdrr.io/cran/survival/)
		[`r-survminer`](https://rdrr.io/cran/survminer/)
	"""
	return Box(desc = "Survival analysis",
		lang   = params.Rscript.value,
		input  = 'infile:file, covfile:file',
		output = [
			'outfile:file:{{i.infile | fn2}}.dir/{{i.infile | fn2}}.survival.txt',
			'outdir:dir:{{i.infile | fn2}}.dir'
		],
		args = Box(inunit = 'days', # months, weeks, years
			outunit   = 'days',
			method    = 'cox', # or km
			covfile   = None,
			inopts    = Box(rnames = True),
			# params for arrange_ggsurvplots.
			# Typically nrow or ncol is set.
			# If combine.ncol = 3, that means {ncol: 3, nrow: 1}.
			# If ncol is not set, then it defaults to 1.
			# If empty, the figures will not be combined
			combine   = Box(ncol = 2),
			devpars   = Box(res = 300, height = 2000, width = 2000),
			params    = Box(),
			ggs       = Box(table = Box())))

@procfactory
def _pPostSurvival():
	"""
	@name:
		pPostSurvival
	@description:
		Statistic comparison between groups after survival analysis.
	@input:
		`infile:file`: The result file from `pSurvival`
		`survfile:file`: The survival data. See format of infile of `pSurvival`
	@output:
		`outfile:file`: The output excel file.
	@args:
		`covfile`: The covariant file. Require rownames in both this file and input file.
		`methods`: A list of testing methods
			- `wilcox`: Wilcox rank sum test
			- `t`: t-test
			- `chisq`: chisquare-test
		`inopts`: The input options for `i.survfile`.
			- `rnames`: whether the file has row names. This has to be True if `args.covfile` provided.
	"""
	pPostSurvival              = Proc(desc = "Post survival analysis: statistics on variables in different groups")
	pPostSurvival.input        = 'infile:file, survfile:file'
	pPostSurvival.output       = 'outfile:file:{{i.infile | fn2}}.stats.xlsx'
	pPostSurvival.args.chi2n   = 10
	pPostSurvival.args.inopts  = Box(rnames = True)
	pPostSurvival.args.covfile = None
	pPostSurvival.envs.rimport = rimport
	pPostSurvival.lang         = params.Rscript.value
	pPostSurvival.script       = "file:scripts/stats/pPostSurvival.r"
	return pPostSurvival

@procfactory
def _pBin():
	"""
	@name:
		pBin
	@description:
		Bin the data in columns.
	@input:
		`infile:file`: The input file
	@output:
		`outfile:file`: The output file. Default: `{{i.infile | stem}}.binned{{i.infile | ext}}`
	@args:
		`inopts`: The input options.
			- `delimit`: The delimiter. Default: `\t`
			- `rnames`: Whether input file has row names. Default: `False`
			- `cnames`: Whether input file has column names. Default: `True`
			- Other arguments available for `read.table`
		`binopts`: The default bin options:
			- `nbin`: Number of bins.
			- `step`: The step of binning.
			- `nan`:  What to do if the value is not a number. Default: `skip`
				- `skip/keep`: Keep it
				- `as0`: Treat it as 0
			- `out`: The out value. Default: `step`
				- `step`: Use the step breaks
				- `lower/min`: Use the min value of the records in the bin
				- `upper/max`: Use the max value of the records in the bin
				- `mean`: Use the mean value of the records in the bin
				- `median`: Use the median value of the records in the bin
				- `binno`: Use the bin number (empty bins will be skipped).
		`cols`: The detailed bin options for each column.
			- If not provided (`None`), all columns will use `binopts`.
			- If column specified, only the specified column will be binned.
			- Column indices can be used. It's 1-based.
	"""
	pBin              = Proc(desc = "Bin the data")
	pBin.input        = 'infile:file'
	pBin.output       = 'outfile:file:{{i.infile | stem}}.binned{{i.infile | ext}}'
	pBin.args.inopts  = Box(delimit = "\t", rnames = False, cnames = True)
	pBin.args.binopts = Box(nbin = None, step = None, nan = 'skip', out = 'step') # lower/min, upper/max, mean, median, binno
	pBin.args.cols    = None
	pBin.envs.rimport = rimport
	pBin.lang         = params.Rscript.value
	pBin.script       = "file:scripts/stats/pBin.r"
	return pBin

@procfactory
def _pQuantileNorm():
	"""
	@name:
		pQuantileNorm
	@description:
		Do quantile normalization
	@input:
		`infile:file`: The input matrix
	@output:
		`outfile:file`: The output matrix. Default: `{{i.infile | bn}}`
	"""
	pQuantileNorm              = Proc(desc = 'Do quantile normalization')
	pQuantileNorm.input        = 'infile:file'
	pQuantileNorm.output       = 'outfile:file:{{i.infile | bn}}'
	pQuantileNorm.args.inopts  = Box(rnames = True, cnames = True, delimit = "\t", skip = 0)
	pQuantileNorm.envs.rimport = rimport
	pQuantileNorm.lang         = params.Rscript.value
	pQuantileNorm.script       = "file:scripts/stats/pQuantileNorm.r"
	return pQuantileNorm

@procfactory
def _pChiSquare():
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
	pChiSquare.output = "outfile:file:{{i.infile | fn2}}.chi2.txt, obsvfile:file:{{i.infile | fn2}}.obsv.txt, exptfile:file:{{i.infile | fn2}}.expt.txt"
	pChiSquare.args.intype = 'cont' # raw
	pChiSquare.args.ctcols = ''
	pChiSquare.lang = params.Rscript.value
	pChiSquare.script = "file:scripts/stats/pChiSquare.r"
	return pChiSquare

@procfactory
def _pFisherExact():
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
			# Contingency table rows: Disease, Healthy
			# Contingency table cols: Mut, Non
			#
			#    | Disease Healthy | Mut  Non  |
			# ---+--------+--------+-----+-----+
			# S1 |    1   |    0   |  0  |  1  |
			# S2 |    0   |    1   |  1  |  0  |
			# .. |   ...  |   ...  | ... | ... |
			# Sn |    0   |    1   |  0  |  1  |
			#
			```
		`ctcols`: The colnames of contingency table if input file is raw values
			- You may also specify them in the head of the input file
	"""
	pFisherExact             = Proc(desc = "Do fisher exact test.")
	pFisherExact.input       = "infile:file"
	pFisherExact.output      = "outfile:file:{{i.infile | fn2}}.fexact.txt"
	pFisherExact.args.intype = 'cont' # raw
	pFisherExact.args.ctcols = []
	pFisherExact.lang        = params.Rscript.value
	pFisherExact.script      = "file:scripts/stats/pFisherExact.r"
	return pFisherExact

@procfactory
def _pPWFisherExact():
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
			#    | A | B | ... | X |
			# ---+---+---+-----+---+
			# S1 | 1 | 0 | ... | 1 |
			# S2 | 0 | 1 | ... | 0 |
			# .. | 0 | 0 | ... | 1 |
			# Sn | 0 | 1 | ... | 1 |
			#
			```
		`padj`: The p-value adjustment method, see `p.adjust.methods` in R. Default: `BH`
		`rnames`: If the input file has rownames for `raw` input type.
	"""
	pPWFisherExact              = Proc(desc = "Do pair-wise fisher exact test.")
	pPWFisherExact.input        = "infile:file"
	pPWFisherExact.output       = "outfile:file:{{i.infile | fn2}}.pwfexact.txt"
	pPWFisherExact.args.rnames  = True
	pPWFisherExact.args.intype  = 'raw' # pairs
	pPWFisherExact.args.padj    = 'BH'
	pPWFisherExact.envs.rimport = rimport
	pPWFisherExact.lang         = params.Rscript.value
	pPWFisherExact.script       = "file:scripts/stats/pPWFisherExact.r"
	return pPWFisherExact

@procfactory
def _pMediation():
	"""
	@name:
		pMediation
	@description:
		Do mediation analysis
	@input:
		`infile:file`: The input file (a matrix or data.frame). Example:
			```
				V1   V2   V3
			S1   1    2    3
			S2   4    1    8
			... ...
			Sn   3    3    1
			```
		`casefile:file`: The mediation options. Example:
			```
			Case1   lm      V3~V2|V1
			Case2   lm,glm  V3~V1|V2
			```
			- No column names, but implies `Case`, 'Model' and `Formua`.
			- `\t` as delimiter.
			- This file is optional. If it is not provided, `args.case` is required.
			- If this file is provided, `args.case` is ignored
			- For different models, model for Mediator comes last. For Case2, the models will be:
				- lm(V3 ~ V2 + V1) and glm(V2 ~ V1)
	@output:
		`outfile:file`: The result file.
		`outdir:dir`  : The output directory containing output file and plots.
	@args:
		`inopts`: The options for input file. Default: `Box(cnames = True, rnames = True)`
			- `cnames`: Whether the input file has column names
			- `rnames`: Whether the input file has row names
		`medopts`: The options for mediation analysis.
			- `boot`: Use bootstrap?
			- `sims`: How many time simulations?
		`cov`: The covariate file. Default: ``
		`pval`: The pvalue cutoff. Default: `0.05`
		`fdr` : Method to calculate fdr. Use `False` to disable. Default: `True` (`BH`)
		`plot`: Parameters for `plot.mediate`? Use `False` to disable plotting. Default: `Box()`
			- Only case with pvalue < `args.pval` will be plotted.
			- To plot all cases, use `args.pval = 1`
		`nthread`: Number of threads to use for different cases. Default: `1`
		`devpars`: device parameters for the plot. Default: `Box(res = 300, width = 2000, height = 2000)`
		`case`   : Define cases, each case should have `model` and `fmula`.
			- If you only have one case, then it could be: `Box(model = 'lm', fmula = 'Y~X|M')`
			In this case, `{{i.infile | fn2}}` will be used as case name
			- For multiple cases, this should be a dict of cases:
			`Box(Case1 = Box(model='lm', fmula='Y~X|M'), Case2 = ...)`
	"""
	pMediation = Proc(desc = "Do mediation analysis.")
	pMediation.input  = 'infile:file, casefile:file'
	pMediation.output = [
		'outfile:file:{{i.infile | fn2}}.mediation/{{i.infile | fn2}}.mediation.txt',
		'outdir:dir:{{i.infile | fn2}}.mediation'
	]
	pMediation.args.inopts  = Box(cnames = True, rnames = True)
	pMediation.args.medopts = Box(boot = True, sims = 500)
	pMediation.args.cov     = ''
	pMediation.args.pval    = 0.05
	pMediation.args.fdr     = True # BH, or other methods for p.adjust
	pMediation.args.plot    = Box()
	pMediation.args.case    = Box(model = 'lm', fmula = 'Y~X|M')
	pMediation.args.nthread = 1
	pMediation.args.devpars = Box(res = 300, width = 2000, height = 2000)
	pMediation.envs.rimport = rimport
	pMediation.lang         = params.Rscript.value
	pMediation.script       = "file:scripts/stats/pMediation.r"
	return pMediation

@procfactory
def _pLiquidAssoc():
	"""
	@name:
		pLiquidAssoc
	@description:
		Do liquid association analysis
	@input:
		`infile:file`: The input file with input data, where LA will be done on rows.
			```
				S1   S2 ... ... Sn
			G1   1    2  ... ... 9
			G2   3    1  ... ... 1
			... ...
			Gm   9    2  ... ... 3
			```
		`casefile:file`: Defining the groups (X, Y, Z) and the cases. If case (3rd col) is not provided, all will be treated as one case.
			- Group "Z" is required. You can also specify group "X", then the rest will be group "Y"
			```
			G1   X   Case1
			G2   X   Case1
			Gx   Z   Case1
			Gy   Z   Case1
			```
	@output:
		`outfile:file`: The results of the analysis
		`outdir:dir`  : The output directory containing the result file and plots.
	@args:
		`inopts` : The options for reading input file
		`zcat`   : Whether the group "Z" is categorical. Default: `False`
			- If it is, then `stein lemma` is not suitable, we will calculate LA manually (`E(g'(z)`)
		`pval`   : The pval cutoff. Default: `0.05`
		`fdr`    : The method to calculate FDR. Use `False` to disable. Default: `True` (BH)
		`nthread`: The number of threads to use. Default: `1` (WCGNA requires)
		`plot`   : Whether do plotting or not. Default: `False`
		`devpars`: device parameters for the plot. Default: `Box(res = 300, width = 2000, height = 2000)`
	@requires:
		[r-fastLiquidAssociation](https://github.com/pwwang/fastLiquidAssociation)
	"""
	pLiquidAssoc        = Proc()
	pLiquidAssoc.input  = 'infile:file, casefile:file'
	pLiquidAssoc.output = [
		'outfile:file:{{i.infile | fn2}}.la/{{i.infile | fn2}}.la.txt',
		'outdir:dir:{{i.infile | fn2}}.la'
	]
	pLiquidAssoc.args.inopts = Box(cnames = True, rnames = True)
	pLiquidAssoc.args.zcat    = False
	pLiquidAssoc.args.pval    = 0.05
	pLiquidAssoc.args.fdr     = True # BH, or other methods for p.adjust
	pLiquidAssoc.args.fdrfor  = 'case' # all
	pLiquidAssoc.args.nthread = 1
	pLiquidAssoc.args.plot    = False
	pLiquidAssoc.args.ggs     = Box()
	pLiquidAssoc.args.devpars = Box(res = 300, width = 2000, height = 2000)
	pLiquidAssoc.envs.rimport = rimport
	pLiquidAssoc.lang         = params.Rscript.value
	pLiquidAssoc.script       = "file:scripts/stats/pLiquidAssoc.r"
	return pLiquidAssoc

@procfactory
def _pHypergeom():
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
	pHypergeom              = Proc(desc = "Do hypergeometric test.")
	pHypergeom.input        = 'infile:file'
	pHypergeom.output       = 'outfile:file:{{i.infile | fn2}}.hypergeom.txt'
	pHypergeom.args.intype  = 'raw' # numbers
	pHypergeom.args.inopts  = Box(cnames = True, rnames = True)
	pHypergeom.args.N       = None
	pHypergeom.envs.rimport = rimport
	pHypergeom.lang         = params.Rscript.value
	pHypergeom.script       = "file:scripts/stats/pHypergeom.r"
	return pHypergeom

@procfactory
def _pChow():
	"""
	@name:
		pChow
	@description:
		Do Chow-Test
	@input:
		`infile:file`: The input file for data to do the regressions. Example:
			```
				X1  X2  X3  X4 ... Y
			S1  1   2   1   4  ... 9
			S2  2   3   1   1  ... 3
			... ...
			Sm  3   9   1   7  ... 8
			```
		`groupfile:file`: Specify the groups to compare. You may also specify the cases. The Chow-Test will be done between the group for each case. Example:
			```
					Case1	Case2
			S1	Group1	Group1
			S2	Group1	NA          # Exclude S2 in Case2
			S3	Group2	Group1
			... ...
			Sm	Group2	Group2
			```
			- In such case, the test will be done between Group1 and Group2 for Case1 and Case2, respectively.
			- Instances can be resued (Gt in the example)
			- If cases not provided, all will be treated as one case.
		`casefile:file`: Define the formula (which columns to use for each case). Example:
			```
			Case1	Y ~ X1
			Case2	Y ~ X2
			```
			- This file is optional, then formula `Y ~ .` will be used for all cases.
	@output:
		`outfile:file`: The result file of chow test. Default: `{{i.infile | fn}}.chow/{{i.infile | fn}}.chow.txt`
		`outdir:dir`: The output directory, containing the output file, results of regressions and plots.
	@args:
		`inopts`: The options for input file.
			- `cnames`: Whether the input file has column names. Default: `True`
			- `rnames`: Whether the input file has row names. Default: `True`
		`cov`: The covariate file. `inopts.rnames` required and this file should have row names too. Default: `''`
		`fdr`   : Calculate FDR or not. Use `False` to disable. If `True` will use `BH` method, otherwise, specify the method (see `R`'s `p.adjust`).
		`pval`: The pvalue cutoff. Default: `0.05`
		`plot`: Whether plot the regressions. Default: `False`
		`ggs` : The extra ggs for the plot.
		`devpars`: device parameters for the plot. Default: `Box(res = 300, width = 2000, height = 2000)`
	"""
	pChow              = Proc(desc = "Do Chow-Test")
	pChow.input        = 'infile:file, groupfile:file, casefile:file'
	pChow.output       = 'outfile:file:{{i.infile | fn}}.chow/{{i.infile | fn}}.chow.txt, outdir:dir:{{i.infile | fn}}.chow'
	pChow.args.inopts  = Box(cnames = True, rnames = True)
	pChow.args.cov     = '' # co-variates, inopts.rnames required, and must in same order
	pChow.args.pval    = 0.05
	pChow.args.fdr     = True
	pChow.args.plot    = True
	pChow.args.devpars = Box(res = 300, width = 2000, height = 2000)
	pChow.args.ggs     = Box()
	pChow.envs.rimport = rimport
	pChow.lang         = params.Rscript.value
	pChow.script       = "file:scripts/stats/pChow.r"
	return pChow

@procfactory
def _pAnovaModel():
	"""
	@name:
		pAnovaModel
	"""
	pAnovaModel              = Proc(desc = "Run anova test on model to test the significance with/without terms")
	pAnovaModel.input        = 'infile:file, casefile:file'
	pAnovaModel.output       = 'outfile:file:{{i.infile | fn}}.anova/{{i.infile | fn}}.anova.txt, outdir:dir:{{i.infile | fn}}.anova'
	pAnovaModel.args.model   = 'lm'
	pAnovaModel.args.fmula   = None
	pAnovaModel.args.inopts  = Box(cnames = True, rnames = True)
	pAnovaModel.args.cov     = '' # requires inopts.rnames
	pAnovaModel.envs.rimport = rimport
	pAnovaModel.lang         = params.Rscript.value
	pAnovaModel.script       = "file:scripts/stats/pAnovaModel.r"
	return pAnovaModel

@procfactory
def _pCorr():
	"""
	@input:
		infile: The input file of data to calculate correlations.
			- The columns are instances and rows variables when `args.byrow = True`.
			- Transposed otherwise.
		groupfile: The group file defining the groups of variables, then the correlation will be only calculated between variables in different groups.
			- First column is the variable names, and second is the groups. Only two groups allowed.
			- Overlaps allowed in different groups.
			- If this is not provided, `args.groupfile` will be checked.
			- If none of them provided, correlation will be calculated for each pair of variables.
			- We can also specify variables directly here by 'var1,var2;var3,var4'
	@output:
		outfile: The output file containing the correlation coefficients
		outdir : The output directory containing the outfile, pvalue file and the plot
	@args:
		outfmt: The output format. Could be `matrix` or `pairs`.
		method: The method used to calculate the correlation coefficient.
			- Could be `pearson`, `spearman` or `kendall`
		byrow: Whether the rows are variables or not.
		pval:  Whether output the pvalues as well.
			- Will be generated at `<outdir>/<i.infile | stem>.<args.method>.pval.txt`
		groupfile: The group file. See `in.groupfile`
		inopts: The input options:
			- `cnames`: Required if `args.byrow = False`. Whether the input file has header.
			- `rnames`: Required if `args.byrow = True`. Whether the input file has row names.
			- `delimit`: The separator of columns.
		plot   : Whether output a correlation plot.
		params : The params for `plot.heatmap` in `utils/plot.r`
		ggs    : The extra ggplot2 statements.
		devpars: The parameters for the plot device.
	@requires:
		R packages: `ggplot2` and `reshape`
	"""
	return Box(
		desc   = 'Correlation Coefficients between variables',
		lang   = params.Rscript.value,
		input  = 'infile:file, groupfile:var',
		output = [
			'outfile:file:{{i.infile | stem}}.{{args.method}}/{{i.infile | stem}}.{{args.method}}.txt',
			'outdir:dir:{{i.infile | stem}}.{{args.method}}'
		],
		args = Box(
			outfmt    = 'pairs',
			method    = 'pearson',
			byrow     = True,
			pval      = False,
			groupfile = None,
			inopts    = Box(cnames = True,	rnames = True,	delimit = "\t"),
			plot      = False,
			params    = Box(), # the parameters for plot.heatmap2
			devpars   = Box(height = 2000, width = 2000, res = 300)
		)
	)

@procfactory
def _pCorr2():
	"""
	@name:
		pCorr2
	@description:
		Calculate correlation coefficient between instances of two files
		Don't do it between instances within the same file.
	@input:
		`infile1:file`: The first file. See input of `pCorr`
		`infile2:file`: The second file.
			- must have same number of columns with `infile1`
	@output:
		`outfile:file`: The output file.
		`outdir:dir`  : The output directory containing output file and other files:
			- pvalues/fdr file and plots
	@args:
		`pval`  : Whether output pvalue. Default: `False`
		`fdr`   : Whether output qvalue. Default: `False`
		`outfmt`: The output format. `pairs` (default) or `matrix`
		`plot`  : Whether plot a heatmap or not. Default: `False`
		`params`: The params for `plot.heatmap` in `utils/plot.r`
		`ggs`:    The extra ggplot2 statements.
		`devpars`:The parameters for the plot device. Default: `Box(height = 2000, width = 2000, res = 300)`
	"""
	pCorr2        = Proc(desc = 'Calculate correlation coefficient between instances of two files')
	pCorr2.input  = 'infile1:file, infile2:file'
	pCorr2.output = [
		'outfile:file:{{i.infile1 | fn2}}-{{i.infile2 | fn2}}.corr/{{i.infile1 | fn2}}-{{i.infile2 | fn2}}.corr.txt',
		'outdir:dir:{{i.infile1 | fn2}}-{{i.infile2 | fn2}}.corr'
	]
	pCorr2.args.inopts1 = Box()
	pCorr2.args.inopts2 = Box()
	pCorr2.args.method  = 'pearson' # spearman, kendall
	pCorr2.args.pval    = False
	pCorr2.args.fdr     = False # method: 'BH'
	pCorr2.args.outfmt  = 'pairs' # matrix
	pCorr2.args.plot    = False
	pCorr2.args.params  = Box() # the parameters for plot.heatmap
	pCorr2.args.ggs     = Box() # extra ggplot statements
	pCorr2.args.devpars = Box(height = 2000, width = 2000, res = 300)
	pCorr2.envs.rimport = rimport
	pCorr2.lang         = params.Rscript.value
	pCorr2.script       = "file:scripts/stats/pCorr2.r"
	return pCorr2

@procfactory
def _pDiffCorr():
	"""
	@name:
		pDiffCorr
	@description:
		Test correlation differences using Fisher Z method.
	@input:
		`infile:file`: The entire dataset used to calculate correlations. Rownames and colnames are required. Example:
			```
					S1  S2  S3  S4 ... Sn
				G1  1   2   1   4  ... 9
				G2  2   3   1   1  ... 3
				... ...
				Gm  3   9   1   7  ... 8
			```
		`samfile:file`: The sample groups, between which you want to compare the correlations. You can also specify one sample to multiple groups, and assign with different cases. Example:
			```
				S1	Healthy
				S2	Healthy
				S3	Disease
				... ...
				Sn	Disease
			```
		`casefile:file`: Assign the cases to compare. If not provided, it will do for every possible combination. Example:
			```
				Healthy	Disease
			```
		`groupfile:file`: Specify groups for rows, then the correlation will be only done within the pairs, each of which is from different groups (only 2 allowed). If not provided, it will investigate for every possible row pairs. Example:
			```
				G1	Kinase
				G2	Kinase
				... ...
				Gm	TF
			```
	@output:
		`outfile:file`: The pairs under different cases that their correlations have been changed significantly
		`outdir:dir`: The output directory containing plots and other output files.
	@args:
		`inopts`: The input options for `infile`. `cnames` and `rnames` have to be `True`
		`method`: The method to calculate the correlation. Default: `pearson`
		`pval`  : The pvalue cutoff to define correlation change significance. Default: `0.05`
		`fdr`   : Calculate FDR or not. Use `False` to disable. If `True` will use `BH` method, otherwise, specify the method (see `R`'s `p.adjust`).
		`fdrfor`: Do FDR calculation for each case (`case`) or all instances (`all`). Default: `case`
		`plot`  : Plot the correlation for cases? Default: `False`
		`ggs`   : `ggs` items for the plot.
		`devpars`: The device parameters for the plot.
	"""
	pDiffCorr = Proc(desc = 'Test correlation differences.')
	pDiffCorr.input  = 'infile:file, samfile:file, casefile:file, groupfile:file'
	pDiffCorr.output = [
		'outfile:file:{{i.infile | fn2}}.diffcorr/{{i.infile | fn2}}.diffcorr.txt',
		'outdir:dir:{{i.infile | fn2}}.diffcorr'
	]
	pDiffCorr.args.inopts  = Box(cnames = True, rnames = True)
	pDiffCorr.args.method  = 'pearson' # spearman
	pDiffCorr.args.pval    = 0.05
	pDiffCorr.args.fdr     = True # BH
	pDiffCorr.args.fdrfor  = 'case'
	pDiffCorr.args.plot    = False
	pDiffCorr.args.ggs     = Box() # extra ggplot statements
	pDiffCorr.args.devpars = Box(height = 2000, width = 2000, res = 300)
	pDiffCorr.envs.rimport = rimport
	pDiffCorr.lang         = params.Rscript.value
	pDiffCorr.script       = "file:scripts/stats/pDiffCorr.r"
	return pDiffCorr

@procfactory
def _pBootstrap():
	"""
	@name:
		pBootstrap
	@description:
		Do bootstrapping resampling
	@input:
		`infile:file`: The input data file
	@output:
		`outfile:file`: The output file with the bootstrapped statistics values
			- depends on the `args.stats` function
		`outdir:dir`: The directory to save the outfile and figures.
	@args:
		`inopts`:  The options to read the input file. Default: `Box(cnames = True, rnames = True)`
		`params`:  Other parameters for `boot` function from R `boot` package
		`nthread`: # of threads(cores) to use. Default: `1`
		`n`: Sampling how many times? Default: `1000`
		`stats`: The function to generate statistics for output. Default: `function(x) x`
			- Default to use all data
			- This function can return a multiple statistics in a vector
			- The argument `x` is the data generate for each sampling.
			- Unlink the `statistic` argument from `boot`, to make it convenient, we don't put the `index` here.
		`plot`: Plot the statistics? Default: `all` (plot all statistics)
			- You may also specify indices. For example: `[1, 2]` to plot the 1st and 2nd statistics
			- Use `False` or `None` to disable plotting
		`devpars`: The device parameters for the plot.
	"""
	pBootstrap        = Proc(desc = 'Do bootstrapping')
	pBootstrap.input  = 'infile:file'
	pBootstrap.output = [
		'outfile:file:{{i.infile | fn2}}.boot/{{i.infile | fn2}}.boot.txt',
		'outdir:dir:{{i.infile | fn2}}.boot'
	]
	pBootstrap.args.inopts  = Box(cnames = True, rnames = True)
	pBootstrap.args.params  = Box()
	pBootstrap.args.nthread = 1
	pBootstrap.args.n       = 1000
	pBootstrap.args.stats   = 'function(x) x'
	pBootstrap.args.plot    = 'all'
	pBootstrap.args.devpars = Box(height = 2000, width = 2000, res = 300)
	pBootstrap.envs.rimport = rimport
	pBootstrap.lang         = params.Rscript.value
	pBootstrap.script       = "file:scripts/stats/pBootstrap.r"
	return pBootstrap


@procfactory
def _pPCA():
	"""
	@description:
		Perform PCA analysis using PCAtools.
		See: https://bioconductor.org/packages/release/bioc/vignettes/PCAtools/inst/doc/PCAtools.html
	@input:
		infile: The matrix to do the analysis
			- Columns are features, rows are samples
			- If not, you may use `args.params.transposed = True`
		metafile: The metadata file
			- Columns are features, rows are samples
			- rnames and cnames are required
	@output:
		outfile: The file with the components
		oudir  : The directory containing the output file and plots.
	@args:
		devpars (Box)        : The parameters for plotting device.
		inopts  (Box)        : The options to read the input files.
		na      (bool|number): How to deal with `NA` values.
			- A logistic/boolean value will remove them (use `complete.cases`)
			- Otherwise, it will be replaced by the given value.
		plots (Box): Plotting parameters supported by `PCAtools`. Set `False` to disable a plot.
			- scree: A scree plot. See `?screeplot`
				- Since `PCAtools::screeplot` returns a ggplot object, we have `args.plots.scree.ggs` to extend the plot.
			- bi: A bi-plot. See `?biplot`
			- pairs: A pairs plot. See `?pairsplot`
			- loadings: A loadings plot. See `?plotloadings`
			- eigencor: An eigencor plot. See `?eigencorplot`
		params (Box): Other parameters for `PCAtools::pca`
		npc (int|str): The number of PCs to write to the output file.
			- A fixed number of PCs; or
			- one of `horn` or `elbow` to determine the optimal number of PCs.
	"""
	return Box(
		desc   = 'Perform PCA analysis using PCAtools',
		input  = 'infile:file, metafile:file',
		output = [	'outfile:file:{{i.infile | stem | stem}}.pca/{{i.infile | stem | stem}}.pcs.txt',
					'outdir:dir:{{i.infile | stem | stem}}.pca'],
		lang = params.Rscript.value,
		args = Box(
			devpars = Box(height = 2000, width = 2000, res = 300),
			inopts  = Box(cnames = True, rnames = True),
			params  = Box(),
			na      = 0,
			plots   = Box(scree = True, bi = True, pairs = True, loadings = True, eigencor = True),
			npc     = 'elbow'
		)
	)

