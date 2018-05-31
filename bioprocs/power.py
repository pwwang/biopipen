from pyppl import Proc, Box
from . import params, rimport

"""
@name:
	pSurvivalPower
@description:
	Do power analysis for survival analysis.
	See http://www.sample-size.net/sample-size-survival-analysis/
@input:
	`infile:file`: The input file, could be either:
	 	- detailed suvival data with [`patient`, ]`time`, `status`, `variable1`, `variable2`, ...; or
		- ratios with `variable`, `survrate1`, `survrate2`, `ssratio`, where `survrate1` and
			`survrate2` are survival rates in group1 and group2, respectively,
			and `ssratio` is sample size ratio in group1/group2
@output:
	`outfile:file`: The output file with columns:
		- Variable: the variable (including paired groups)
		- Alpha: the alpha value
		- Beta: the beta value (1-power)
		- SSize1: the sample size for group1
		- SSize2: the sample size for group2
		- Total: the total sample size
"""
pSurvivalPower              = Proc(desc = "Survival analysis.")
pSurvivalPower.input        = 'infile:file'
pSurvivalPower.output       = 'outdir:dir:{{in.infile | fn}}.power'
pSurvivalPower.args.rnames  = True
pSurvivalPower.args.plot    = False
pSurvivalPower.args.ngroup  = 2
pSurvivalPower.args.intype  = 'detailed' # ratios or ratio
pSurvivalPower.args.alphas  = [.005, .01, .05, .1] # two-sided, *2 if you want one-sided
pSurvivalPower.args.betas   = [.05, .1, .2]
pSurvivalPower.envs.rimport = rimport
pSurvivalPower.lang         = params.Rscript.value
pSurvivalPower.script       = "file:scripts/power/pSurvivalPower.r"
