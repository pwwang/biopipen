# {{title}}

This analysis is done using PyClone[1], which is designed to infer the prevalence of point mutations in heterogeneous cancer samples.

{% python from os import path %}
{% for job in jobs %}

:::: {.tab}
{% 	if path.exists(path.join(job.o.outdir, 'tables')) %}

{# tab title #}
{% 		if len(jobs) > 1 %}
## {{job.i.muts | stem}}
{% 		endif %}

{% 	if len(jobs) > 1 %}#{% endif %}## Cluster
::: {.tab}
{% 	if len(jobs) > 1 %}#{% endif %}### Density
![Cluster density]({{job.o.outdir}}/plots/cluster/density.svg)
:::

::: {.tab}
{% 	if len(jobs) > 1 %}#{% endif %}### Parallel coordinates
{% 	if path.exists(path.join(job.o.outdir, 'plots/cluster/parallel_coordinates.svg')) %}
![Parallel coordinates for clusters]({{job.o.outdir}}/plots/cluster/parallel_coordinates.svg)
{% 	endif %}
:::

::: {.tab}
{% 	if len(jobs) > 1 %}#{% endif %}### Scatter
{% 	if path.exists(path.join(job.o.outdir, 'plots/cluster/scatter.svg')) %}
![Scatter for clusters]({{job.o.outdir}}/plots/cluster/scatter.svg)
{% 	endif %}
:::

::: {.tab}
{% 	if len(jobs) > 1 %}#{% endif %}### Table
```table
file: {{job.o.outdir}}/tables/cluster.tsv
caption: Cellular density of clusters
rows: 100
```
:::

{% 	if len(jobs) > 1 %}#{% endif %}## Loci
::: {.tab}
{% 	if len(jobs) > 1 %}#{% endif %}### Density
{% 	if path.exists(path.join(job.o.outdir, 'plots/loci/density.svg')) %}
![cluster-density]({{job.o.outdir}}/plots/loci/density.svg)
{%  endif %}
:::

::: {.tab}
{% 	if len(jobs) > 1 %}#{% endif %}### Parallel coordinates
{% 	if path.exists(path.join(job.o.outdir, 'plots/loci/parallel_coordinates.svg')) %}
![Cellular prevalence]({{job.o.outdir}}/plots/loci/parallel_coordinates.svg)
{% 	endif %}
{% 	if path.exists(path.join(job.o.outdir, 'plots/loci/vaf_parallel_coordinates.svg')) %}
![VAF]({{job.o.outdir}}/plots/loci/vaf_parallel_coordinates.svg)
{% 	endif %}
:::

::: {.tab}
{% 	if len(jobs) > 1 %}#{% endif %}### Scatter
{% 	if path.exists(path.join(job.o.outdir, 'plots/loci/scatter.svg')) %}
![Cellular prevalence]({{job.o.outdir}}/plots/loci/scatter.svg)
{% 	endif %}
{% 	if path.exists(path.join(job.o.outdir, 'plots/loci/vaf_scatter.svg')) %}
![VAF]({{job.o.outdir}}/plots/loci/vaf_scatter.svg)
{% 	endif %}
:::

::: {.tab}
{% 	if len(jobs) > 1 %}#{% endif %}### Similarity
{% 	if path.exists(path.join(job.o.outdir, 'plots/loci/similarity_matrix.svg')) %}
![Loci similarity by cellular prevalence]({{job.o.outdir}}/plots/loci/similarity_matrix.svg)
{% 	endif %}
:::

::: {.tab}
{% 	if len(jobs) > 1 %}#{% endif %}### Table
```table
file: {{job.o.outdir}}/tables/loci.tsv
caption: Cellular prevalence for loci
rows: 100
```
:::

{% else %}{# if path.exists #}

{# Show the last line of error messages #}
{{job.errfile | readlines | lambda x: x[-1] if x else ''}}

{% endif %}{# if path.exists #}
::::
{% endfor %}{# for jobs #}

[1]: Roth, Andrew, et al. "PyClone: statistical inference of clonal population structure in cancer." Nature methods 11.4 (2014): 396.
