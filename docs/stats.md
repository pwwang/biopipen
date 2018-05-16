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

## pSurvival

### description
Survival analysis

### input
#### `infile:file`:: The input file (header is required).  
	- col1: rownames if args.rnames = True
	- col2: the survival time
	- col3: the status. 0/1 for alive/dead or 1/2 for alive dead
	- col4: group1.
	- ... other groups

### output
#### `outdir:dir`:: The output directory containing the pval files and plots  

### args
#### `inunit`    :: The time unit in input file. Default: days  
#### `outunit`   :: The output unit for plots. Default: days  
#### `nthread`   :: Number of threads used to perform analysis for groups. Default: 1  
#### `rnames`    :: Whether input file has row names. Default: True  
#### `combine`   :: Whether combine groups in the same plot. Default: True  
#### `devpars`   :: The device parameters for png. Default: `{res:300, height:2000, width:2000}`  
	- The height and width are for each survival plot. If args.combine is True, the width and height will be multiplied by `max(gridParams.ncol, gridParams.nrow)`
#### `plotParams`:: The parameters for `ggsurvplot`. Default: `{risk.table: True, conf.int = True}`  
#### `gridParams`:: The parameters for `arrange_ggsurvplots`.  
#### `pval`      :: Whether print pvalue on the plot. Default: True  
#### `noerror`   :: Do not report error if error happens. Generate ouput file anyway.  

## pChiSquare

### description
Do chi-square test.

### input
#### `infile:file`:: The input file.  

### output
#### `outfile:file` :: The output file containing Xsquare, df, pval and method  
#### `obsvfile:file`:: The observation matrix  
#### `exptfile:file`:: The expectation matrix  

## pFisherExact

### description
Do fisher exact test.

### input
#### `infile:file`:: The input file.  

### output
#### `outfile:file` :: The output file containing confInt1, confInt2, oddsRatio, pval, alternative and method.  

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
{% endraw %}
