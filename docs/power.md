# power
<!-- toc -->
{% raw %}

## pSurvivalPower

### description
Do power analysis for survival analysis.
See http://www.sample-size.net/sample-size-survival-analysis/

### input
#### `infile:file`:: The input file, could be either:  
 	- detailed suvival data with [`patient`, ]`time`, `status`, `variable1`, `variable2`, ...; or
	- ratios with `variable`, `survrate1`, `survrate2`, `ssratio`, where `survrate1` and
		`survrate2` are survival rates in group1 and group2, respectively,
		and `ssratio` is sample size ratio in group1/group2
{% endraw %}
