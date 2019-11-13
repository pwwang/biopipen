# {{title}}

LICHeE[1] is a combinatorial method designed to reconstruct multi-sample cell lineage trees and infer the subclonal composition of the given samples based on variant allele frequencies (VAFs) of deep-sequencing somatic single nucleotide variants (SSNVs). The program accepts as input a list of SNVs with specified per-sample VAFs and outputs the inferred cell lineage tree(s) and the sample subclone decomposition. It provides an optional simple GUI to allow users to interact with the trees dynamically.

It:

- Doesn't require mutations to be clustered before analysis
- Doesn't require mutations to be shared in all samples
- Provides sample decomposition
- Reports one tree preferably

{% from os import path %}
{% for job in jobs %}
:::: {.tab}

	{%- if forloop.length > 1: %}
## {{job.i.infile | stem}}
	{%- endif %}

	{%- assign hash = forloop.length | ?:_>1 | :'#' | :'' %}

	{%- if `job.errfile | read | :"Found 0 valid trees after network adjustments" not in _` %}
##{{hash}} Top tree

![Top tree]({{job.o.outdir, '*.png' | *glob1}})

##{{hash}} Cluster details

```table
file: {{job.o.outdir, '*.snvinfo.txt' | *glob1}}
caption: Variant information
```
	{%- else %}
No valid trees generated.
	{%- endif %}

::::
{% endfor %}


[1]: Popic, Victoria, et al. "Fast and scalable inference of multi-sample cancer lineages." Genome biology 16.1 (2015): 91.