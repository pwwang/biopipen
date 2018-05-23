# stats
<!-- toc -->
{% raw %}

## pMetaPval

### description
Combine p-values in the files from input directory

### input
#### `indir:dir`:: The directory containing the input files  

### output
#### `outfile:file`:: The output file containing the meta-pvalues  

### args
#### `args.pattern`:: The pattern used to filter the input files. Default: '*'  
#### `args.header`:: Whether the input files contains a header. Default: True  
	- Could be a list to specify it for each file.
	- The order should be concordant with the file names
#### `args.pcol`:: Which column is the p-value. Default: -1 (last column)  
#### `args.poutonly`:: Only output pvalues. Default: False (output all possible information)  
#### `args.outheader`:: Whether output the header. Default: True  
#### `args.method`:: The method used to calculate the meta-pvalue. Default: sumlog (Fisher's method)  
	- Other available methods: logitp, sumz, votep, sump, meanp and wilkinsonp
	- See: https://www.rdocumentation.org/packages/metap/versions/0.8

### requires
[`r-matep`](https://www.rdocumentation.org/packages/metap/)

## pMetaPval1

### description
Combine p-values in a single file by rows.

### input
#### `infile:file`:: The input file  

### output
#### `outfile:file`:: The output file containing the meta-pvalues  

### args
#### `args.header`:: Whether the input files contains a header. Default: True  
#### `args.pcol`:: Which column is the p-value. Default: -1 (last column)  
#### `args.poutonly`:: Only output pvalues. Default: False (output all possible information)  
#### `args.outheader`:: Whether output the header. Default: True  
#### `args.method`:: The method used to calculate the meta-pvalue. Default: sumlog (Fisher's method)  
	- Other available methods: logitp, sumz, votep, sump, meanp and wilkinsonp
	- See: https://www.rdocumentation.org/packages/metap/versions/0.8

### requires
[`r-matep`](https://www.rdocumentation.org/packages/metap/)

## pSurvival

### description
Survival analysis

### input
#### `infile:file`:: The input file (header is required).  
	- col1: rownames if args.inopts.rnames = True
	- col2: the survival time
	- col3: the status. 0/1 for alive/dead or 1/2 for alive dead
	- col4: var1.
	- ... other variables

### output
#### `outfile:file`:: The outfile containing the pvalues  
#### `outdir:dir`  :: The output directory containing the pval files and plots  

### args
#### `inunit`    :: The time unit in input file. Default: days  
#### `outunit`   :: The output unit for plots. Default: days  
#### `nthread`   :: Number of threads used to perform analysis for groups. Default: 1  
#### `inopts`    :: The options for input file  
	- `rnames`: Whether input file has row names. Default: True
#### `combine`   :: Whether combine groups in the same plot. Default: True  
#### `devpars`   :: The device parameters for png. Default: `{res:300, height:2000, width:2000}`  
	- The height and width are for each survival plot. If args.combine is True, the width and height will be multiplied by `max(arrange.ncol, arrange.nrow)`
#### `covfile`   :: The covariant file. Require rownames in both this file and input file.  
#### `plot`      :: The params for plot.  
	- `ncurves`: Number of curves to plot (the continuous number will divided into `ncurves` groups.
	- `params` : The params for `ggsurvplot`. Default: `Box({'risk.table': True, 'conf.int': True, 'font.legend': 13, 'pval': '{method}\np = {pval}'})`
	- `arrange`: How to arrange multiple survival plots in one if `args.combine = True`.
		- `nrow`: The number of rows. Default: 1
		- `ncol`: The number of cols. Default: 1
#### `ggs`       :: Extra ggplot2 elements for main plot. `ggs.table` is for the risk table.  
#### `pval`      :: The method to calculate the pvalue shown on the plot. Default: True (logrank)  
	- Could also be `waldtest`, `likeratio` (Likelihoold ratio test)

### requires
[`r-survival`](https://rdrr.io/cran/survival/)
[`r-survminer`](https://rdrr.io/cran/survminer/)

## pChiSquare

### description
Do chi-square test.

### input
#### `infile:file`:: The input file.  

### output
#### `outfile:file` :: The output file containing Xsquare, df, pval and method  
#### `obsvfile:file`:: The observation matrix  
#### `exptfile:file`:: The expectation matrix  

### args
#### `intype`:: The type of the input file:  
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
#### `ctcols`:: The colnames of contingency table if input file is raw values  
	- You may also specify them in the head of the input file

## pFisherExact

### description
Do fisher exact test.

### input
#### `infile:file`:: The input file.  

### output
#### `outfile:file` :: The output file containing confInt1, confInt2, oddsRatio, pval, alternative and method.  

### args
#### `intype`:: The type of the input file:  
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
#### `ctcols`:: The colnames of contingency table if input file is raw values  
	- You may also specify them in the head of the input file

## pPWFisherExact

### description
Do pair-wise fisher exact test.
Commonly used for co-occurrence/mutual-exclusivity analysis.
P-value indicates if the pairs are significantly co-occurred or mutually exclusive.
Co-occurrence: Odds ratio > 1
Mutual-exclusivity: Odds ratio < 1

### input
#### `infile:file`:: The input file.  

### output
#### `outfile:file` :: The output file containing confInt1, confInt2, oddsRatio, pval, qval, alternative and method.  

### args
#### `intype`:: The type of the input file:  
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
#### `padj`:: The p-value adjustment method, see `p.adjust.methods` in R. Default: `BH`  

## pMediation

### description
Do mediation analysis

### input
#### `infile:file`:: The input file (a matrix or data.frame).  

### output
#### `outfile:file`:: The result file.  

### args
#### `inopts`:: The options for input file.  
	- `cnames`: Whether the input file has column names
	- `rnames`: Whether the input file has row names
#### `medopts`:: The options for mediation analysis.  
	- `modelm`: The model for M ~ X. Default: `lm(M ~ X)`
	- `modely`: The model for Y ~ X + M. Default: `lm(Y ~ X + M)`
	- `mediator`: Tell the model which column is the mediator
	- `treat`: Tell the model which column is the variable
	- `boot`: Use bootstrap?
	- `sims`: How many time simulations?

## pHypergeom

### description
Do hypergeometric test.

### input
#### `infile:file`:: The input file, could be raw data (presence (1) and absence (0) of elements) or number of overlapped elements and elements in each category.  
	- Set `args.intype` as `raw` if it is raw data. The population size `args.N` is required
	- Set `args.intype` as `numbers` (or any string except `raw`) if it is numbers. You can specified explicit header: `k` = overlapped elements, `m` = size of set 1, `n` = size of set 2 and `N` = the population size. If `N` not included, then `args.N` is required

### output
#### `outfile:file`:: The output file  

### args
#### `intype`:: the type of input file. Default: `raw`. See `infile:file`  
#### `inopts`:: The options for input file.  
	- `cnames`: Whether the input file has column names
	- `rnames`: Whether the input file has row names
#### `N`:: The population size. Default: `None`  
{% endraw %}
