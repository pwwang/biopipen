# {{report.title}}

One of the aim of RNAseq data analysis is the detection of differentially expressed genes. The analysis takes the normalised read count data and performing statistical analysis to discover quantitative changes in expression levels between experimental groups.

{% case args.tool %}
	{%- when 'deseq2' %}
Here we use DESeq2[[1]] to perform the analysis. DESeq2 provides methods to test for differential expression by use of negative binomial generalized linear models; the estimates of dispersion and logarithmic fold changes incorporate data-driven prior distributions.
	{%- when 'edger' %}
Here we use edgeR[[1]] to perform the analysis. edgeR is a Bioconductor software package for examining differential expression of replicated count data. An overdispersed Poisson model is used to account for both biological and technical variability.
{% endcase %}

{% for job in jobs %}
:::: {.panel .accordion}
{%- if forloop.length > 1 %}
## {{ job.i.efile | stem2 | @append: "-" + stem(job.i.gfile) }}
{%-	endif %}

{% 	assign punc = forloop.length | @minus: 1 | ?=: '#' | $?! ''  %}

{{punc}}## Compairson

```table
file: {{job.i.gfile}}
caption: Compairson {{job.i.gfile | stem}}
```

{{punc}}## Differentially expressed genes
::: {.panel .tab}
{{punc}}### Top DEGs

```table
file: {{job.o.outfile}}
caption: Top DEGs
rows: 100
download: true
```
:::

::: {.panel}
{{punc}}### Top up-regulated genes

```table
file: {{job.o.outfile | prefix | prefix | @append: ".up.xls" }}
caption: Top up-regulated genes
rows: 100
download: true
```
:::

::: {.panel}
{{punc}}### Top down-regulated genes

```table
file: {{job.o.outfile | prefix | prefix | @append: ".down.xls" }}
caption: Top down-regulated genes
rows: 100
download: true
```
:::

{{punc}}## Plots

{% 	from os import path %}

{% 	if `job.o.outdir | @append: "/mdsplot.png" | path.exists` -%}
![MDS Plot]({{job.o.outdir | @append: "/mdsplot.png"}})
{%- endif %}

{%-	if `job.o.outdir | @append: "/volcano.png" | path.exists` -%}
![Volcano plot]({{job.o.outdir | @append: "/volcano.png"}})
{%- endif %}

{%-	if `job.o.outdir | @append: "/maplot.png" | path.exists` -%}
![MA plot]({{job.o.outdir | @append: "/maplot.png"}})
{%- endif %}

{%-	if `job.o.outdir | @append: "/heatmap.png" | path.exists`-%}
![Heatmap]({{job.o.outdir | @append: "/heatmap.png"}})
{%- endif %}

{%-	if `job.o.outdir | @append: "/qqplot.png" | path.exists` -%}
![QQplot on p-values]({{job.o.outdir | @append: "/qqplot.png"}})
{%- endif %}

::::
{% endfor %}


{% case args.tool %}
	{%- when 'deseq2' %}
[[1]]: Love, Michael I., Wolfgang Huber, and Simon Anders. "Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2." Genome biology 15.12 (2014): 550.
	{%- when 'edger' -%}
[[1]]: Robinson, Mark D., Davis J. McCarthy, and Gordon K. Smyth. "edgeR: a Bioconductor package for differential expression analysis of digital gene expression data." Bioinformatics 26.1 (2010): 139-140.
{% endcase %}
