# {{title}}

Principal Component Analysis (PCA) extracts the fundamental structure of the data without the need to build any model to represent it. This 'summary' of the data is arrived at through a process of reduction that can transform the large number of variables into a lesser number that are uncorrelated (i.e. the 'principal components'), whilst at the same time being capable of easy interpretation on the original data.

PCAtools[1] provides functions for data exploration via PCA, and allows the user to generate publication-ready figures.

{% for job in jobs %}
:::: {.tab}

	{%- if forloop.length > 1 %}
## {{job.i.infile | stem}}
	{%- endif %}
	{%- assign hash = forloop.length | ?:_>1 | :'#' | :'' %}

	{%- assign plots = job.o.outdir, '*.png' | *glob1: first = False %}
	{%- if plots %}
##{{hash}} Plots
	{% for png in plots %}
::: {.tab .force-tab}

		{%- case png | stem | ext | [1:] %}
			{%- when 'scree' %}
###{{hash}} Scree plot
A scree plot to show the proportion of explained variance by PC

![Scree plot]({{png}})
			{%- when 'bi' %}
###{{hash}} Bi-plot
![Bi-plot]({{png}})
			{%- when 'pairs' %}
###{{hash}} Pairs plot
![Pairs plot]({{png}})
			{%- when 'loadings' %}
###{{hash}} Loadings plot
Component loadings and label genes most responsible for variation

![Loadings plot]({{png}})
			{%- when 'eigencor' %}
###{{hash}} Eigen correlation plot
![Eigen correlation plot]({{png}})
		{%- endcase %}
:::
	{% endfor  %}
	{%- endif %}

##{{hash}} Selected PCs

```table
file: {{job.o.outfile}}
caption: Selected components
```

::::
{% endfor %}

[1]: Blighe, Kevin. "Haplotype Classification Using Copy Number Variation and Principal Components Analysis." Open Bioinformatics Journal 7 (2013): 19-24.