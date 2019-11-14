# {{title}}

SCHISM[1] is a computational tool designed to infer subclonal hierarchy and the tumor evolution from somatic mutations. The inference process involves computational assessment of two fundamental properties of tumor evolution: lineage precedence rule and lineage divergence rule.

First, SCHISM combines information about somatic mutation cellularity (aka mutation cancer cell fraction) across all tumor sample(s) available from a patient in a hypothesis testing framework to identify the statistical support for the lineage relationship between each pair of mutations or mutation clusters. The results of the hypothesis test are represented as Cluster Order Precedence Violation (CPOV) matrix which informs the subsequent step in SCHISM and ensures compliance of candidate tree topologies with lineage precedence rule.

Next, an implementation of genetic algorithm (GA) based on a fitness function that incorporates considerations for both lineage precedence (CPOV) rule and lineage divergence rule explores the space of tree topologies and returns a prioritized list of candidate subclonal phylogenetic trees, most compatible with observed cellularity data.

It:

- Requires mutations to be clustered somewhere else
- Requires mutations to be shared in all samples
- Doesn't report sample decomposition
- Reports multiple trees

{% from os import path %}
{% for job in jobs %}
:::: {.tab}

	{%- if forloop.length > 1: %}
## {{job.i.infile | stem}}
	{%- endif %}

{% assign hash = forloop.length | ?:_>1 | :'#' | :'' %}
##{{hash}} Consensus trees

![Consensus trees]({{job.o.outdir, '*.tree.png' | *glob1}})

##{{hash}} Variants in each cluster

```table
file: {{job.o.outdir | str | path.join: 'mutation_to_cluster_assignment.tsv'}}
caption: 'Variants in each cluster'
```

::::
{% endfor %}

[1]: Niknafs, Noushin, et al. "Subclonal hierarchy inference from somatic mutations: automatic reconstruction of cancer evolutionary trees from multi-region next generation sequencing." PLoS computational biology 11.10 (2015): e1004416.