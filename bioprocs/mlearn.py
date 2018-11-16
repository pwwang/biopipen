from pyppl import Proc, Box
from . import params, rimport

"""
@name:
	pLinearRegTrain
@description:
	Train a linear regression model
@input:
	`infile:file`: The input file (Last column as Y)
@output:
	`outmodel:file`: The output model (RData file)
	`outdir:dir`   : The output directory containing model, plots and other files
@args:
	`plot`   : Whether plot the lm probability. Default: `True`
	`formula`: The formula to perform the regression. Default: `None`.
		- If `None`, will use all first N-1 columns as features.
	`inopts` : The input options.
	`yval`   : The type of y values. Default: `categ`
		- `categ`  : categorical values
		- `numeric`: numeric values
"""
pLinearRegTrain = Proc(desc = 'Train a linear regression model')
pLinearRegTrain.input  = 'infile:file'
pLinearRegTrain.output = [
	'outmodel:file:{{i.infile | stem}}.lm/{{i.infile | stem}}.lm.rds', 
	'outdir:dir:{{i.infile | stem}}.lm'
]
pLinearRegTrain.args.plot    = True
pLinearRegTrain.args.formula = None
pLinearRegTrain.args.devpars = Box(res = 300, height = 2000, width = 2000)
pLinearRegTrain.args.ggs     = Box(
	geom_smooth = Box({
		"method"     : "lm",
		"method.args": Box(family = "binomial"),
		"se"         : True
	})
)
pLinearRegTrain.args.inopts  = Box(
	cnames  = True,
	rnames  = True,
	delimit = "\t"
)
pLinearRegTrain.args.yval    = 'categ'
pLinearRegTrain.envs.rimport = rimport
pLinearRegTrain.lang         = params.Rscript.value
pLinearRegTrain.script       = "file:scripts/mlearn/pLinearRegTrain.r"

"""
@name:
	pLinearRegPredict
@description:
	Use a trained linear regression model to predict
@input:
	`infile:file`: The input file 
	`model:file` : The trained model by `pLinearRegTrain`
@output:
	`outdir:dir`: The output directory
@args:
	`inopts` : The input options.
	`outprob`: Also output probabilities? Default: True
"""
pLinearRegPredict = Proc(desc = 'Use a trained linear regression model to predict')
pLinearRegPredict.input        = 'infile:file, model:file'
pLinearRegPredict.output       = 'outdir:dir:{{i.infile | stem}}.pred'
pLinearRegPredict.args.outprob = True
pLinearRegPredict.args.outauc  = True
pLinearRegPredict.args.params  = Box(labels = False, showAUC = True, combine = True)
pLinearRegPredict.args.ggs     = Box({
	'style_roc': {},
	# show legend at bottom right corner
	'theme#auc': {'legend.position': [1, 0], 'legend.justification': [1, 0]} 
})
pLinearRegPredict.args.devpars = Box(res = 300, height = 2000, width = 2000)
pLinearRegPredict.args.inopts  = Box(
	cnames  = True,
	rnames  = True,
	delimit = "\t"
)
pLinearRegPredict.envs.rimport = rimport
pLinearRegPredict.lang         = params.Rscript.value
pLinearRegPredict.script       = "file:scripts/mlearn/pLinearRegPredict.r"

"""
@name:
	pLogitRegTrain
@description:
	Train a linear regression model
@input:
	`infile:file`: The input file (Last column as Y)
@output:
	`outmodel:file`: The output model (RData file)
	`outdir:dir`   : The output directory containing model, plots and other files
@args:
	`plot`   : Whether plot the glm probability. Default: `True`
	`formula`: The formula to perform the regression. Default: `None`.
		- If `None`, will use all first N-1 columns as features.
	`inopts` : The input options.
	`yval`   : The type of y values. Default: `categ`
		- `categ`  : categorical values
		- `prob`   : probabilities
		- `numeric`: numeric values
"""
pLogitRegTrain = Proc(desc = 'Train a linear regression model')
pLogitRegTrain.input  = 'infile:file'
pLogitRegTrain.output = [
	'outmodel:file:{{i.infile | stem}}.glm/{{i.infile | stem}}.glm.rds', 
	'outdir:dir:{{i.infile | stem}}.glm'
]
pLogitRegTrain.args.plot    = True
pLogitRegTrain.args.formula = None
pLogitRegTrain.args.devpars = Box(res = 300, height = 2000, width = 2000)
pLogitRegTrain.args.ggs     = Box(
	geom_smooth = Box({
		"method"     : "glm",
		"method.args": Box(family = "binomial"),
		"se"         : True
	})
)
pLogitRegTrain.args.inopts  = Box(
	cnames  = True,
	rnames  = True,
	delimit = "\t"
)
pLogitRegTrain.args.yval    = 'categ'
pLogitRegTrain.envs.rimport = rimport
pLogitRegTrain.lang         = params.Rscript.value
pLogitRegTrain.script       = "file:scripts/mlearn/pLogitRegTrain.r"

"""
@name:
	pLogitRegPredict
@description:
	Use a trained linear regression model to predict
@input:
	`infile:file`: The input file 
	`model:file` : The trained model by `pLogitRegTrain`
@output:
	`outdir:dir`: The output directory
@args:
	`inopts` : The input options.
	`outprob`: Also output probabilities? Default: True
"""
pLogitRegPredict = Proc(desc = 'Use a trained linear regression model to predict')
pLogitRegPredict.input        = 'infile:file, model:file'
pLogitRegPredict.output       = 'outdir:dir:{{i.infile | stem}}.pred'
pLogitRegPredict.args.outprob = True
pLogitRegPredict.args.outauc  = True
pLogitRegPredict.args.params  = Box(labels = False, showAUC = True, combine = True)
pLogitRegPredict.args.ggs     = Box({
	'style_roc': {},
	# show legend at bottom right corner
	'theme#auc': {'legend.position': [1, 0], 'legend.justification': [1, 0]} 
})
pLogitRegPredict.args.devpars = Box(res = 300, height = 2000, width = 2000)
pLogitRegPredict.args.inopts  = Box(
	cnames  = True,
	rnames  = True,
	delimit = "\t"
)
pLogitRegPredict.envs.rimport = rimport
pLogitRegPredict.lang         = params.Rscript.value
pLogitRegPredict.script       = "file:scripts/mlearn/pLogitRegPredict.r"

"""
@name:
	pRandomForestTrain
@description:
	Train a random forest model
@input:
	`infile:file`: The input file (Last column as Y)
@output:
	`outmodel:file`: The output model (RData file)
	`outdir:dir`   : The output directory containing model, plots and other files
@args:
	`plot`   : Whether plot the feature importance. Default: `True`
	`formula`: The formula to perform the regression. Default: `None`.
		- If `None`, will use all first N-1 columns as features.
	`inopts` : The input options.
	`na`     : Replace NAs with? Default: `0`
	`devpars`: The device parameters for the plot. Default: `Box(res = 300, height = 2000, width = 2000)`
@requires:
	`r-randomForst`
"""
pRandomForestTrain = Proc(desc = "Train a random forest model.")
pRandomForestTrain.input = 'infile:file'
pRandomForestTrain.output = [
	'outmodel:file:{{i.infile | stem}}.rforest/{{i.infile | stem}}.rforest.rds', 
	'outdir:dir:{{i.infile | stem}}.rforest'
]
pRandomForestTrain.args.plot    = True
pRandomForestTrain.args.formula = None
pRandomForestTrain.args.na      = 0
pRandomForestTrain.args.devpars = Box(res = 300, height = 2000, width = 2000)
pRandomForestTrain.args.inopts  = Box(
	cnames  = True,
	rnames  = True,
	delimit = "\t"
)
pRandomForestTrain.envs.rimport = rimport
pRandomForestTrain.lang         = params.Rscript.value
pRandomForestTrain.script       = "file:scripts/mlearn/pRandomForestTrain.r"

"""
@name:
	pDecisionTreeTrain
@description:
	Train a decision tree model
@input:
	`infile:file`: The input file (Last column as Y)
@output:
	`outmodel:file`: The output model (RData file)
	`outdir:dir`   : The output directory containing model, plots and other files
@args:
	`plot`   : Whether plot the feature importance. Default: `True`
	`formula`: The formula to perform the regression. Default: `None`.
		- If `None`, will use all first N-1 columns as features.
	`inopts` : The input options.
	`na`     : Replace NAs with? Default: `0`
	`devpars`: The device parameters for the plot. Default: `Box(res = 300, height = 2000, width = 2000)`
@requires:
	`r-rpart`
"""
pDecisionTreeTrain = Proc(desc = "Train a decision tree model")
pDecisionTreeTrain.input = 'infile:file'
pDecisionTreeTrain.output = [
	'outmodel:file:{{i.infile | stem}}.dtree/{{i.infile | stem}}.dtree.rds', 
	'outdir:dir:{{i.infile | stem}}.dtree'
]
pDecisionTreeTrain.args.plot    = True
pDecisionTreeTrain.args.formula = None
pDecisionTreeTrain.args.na      = 0
pDecisionTreeTrain.args.devpars = Box(res = 300, height = 2000, width = 2000)
pDecisionTreeTrain.args.inopts  = Box(
	cnames  = True,
	rnames  = True,
	delimit = "\t"
)
pDecisionTreeTrain.envs.rimport = rimport
pDecisionTreeTrain.lang         = params.Rscript.value
pDecisionTreeTrain.script       = "file:scripts/mlearn/pDecisionTreeTrain.r"
