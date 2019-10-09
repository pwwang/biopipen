# {{title}}

This analysis is done using PyClone[1], which is designed to infer the prevalence of point mutations in heterogeneous cancer samples.

{% python from os import path %}
{% for job in jobs %}

:::: {.tab}
{% if path.exists(path.join(job.o.outdir, 'tables')) %}

{# tab title #}
{% if len(jobs) > 1 %}
## {{job.i.infile | stem}}
{% endif %}

{% if len(jobs) > 1 %}#{% endif %}## Cluster
::: {.tab}
{% if len(jobs) > 1 %}#{% endif %}### Density
![Cluster density]({{job.o.outdir}}/plots/cluster/density.svg)
:::

::: {.tab}
{% if len(jobs) > 1 %}#{% endif %}### Parallel coordinates
![Parallel coordinates for clusters]({{job.o.outdir}}/plots/cluster/parallel_coordinates.svg)
:::

::: {.tab}
{% if len(jobs) > 1 %}#{% endif %}### Scatter
![Scatter for clusters]({{job.o.outdir}}/plots/cluster/scatter.svg)
:::

::: {.tab}
{% if len(jobs) > 1 %}#{% endif %}### Table
```table
file: {{job.o.outdir}}/tables/cluster.tsv
caption: Cellular density of clusters
```
:::

{% if len(jobs) > 1 %}#{% endif %}## Loci
::: {.tab}
{% if len(jobs) > 1 %}#{% endif %}### Density
![cluster-density]({{job.o.outdir}}/plots/loci/density.svg)
:::

::: {.tab}
{% if len(jobs) > 1 %}#{% endif %}### Parallel coordinates
![Cellular prevalence]({{job.o.outdir}}/plots/loci/parallel_coordinates.svg)
![VAF]({{job.o.outdir}}/plots/loci/vaf_parallel_coordinates.svg)
:::

::: {.tab}
{% if len(jobs) > 1 %}#{% endif %}### Scatter
![Cellular prevalence]({{job.o.outdir}}/plots/loci/scatter.svg)
![VAF]({{job.o.outdir}}/plots/loci/vaf_scatter.svg)
:::

::: {.tab}
{% if len(jobs) > 1 %}#{% endif %}### Similarity
![Loci similarity by cellular prevalence]({{job.o.outdir}}/plots/loci/similarity_matrix.svg)
:::

::: {.tab}
{% if len(jobs) > 1 %}#{% endif %}### Table
```table
file: {{job.o.outdir}}/tables/loci.tsv
caption: Cellular prevalence for loci
```
:::

{% else %}{# if path.exists #}

{# Show the last line of error messages #}
{{job.stderr | readlines | [-1]}}

{% endif %}{# if path.exists #}
::::
{% endfor %}{# for jobs #}

[1]: Roth, Andrew, et al. "PyClone: statistical inference of clonal population structure in cancer." Nature methods 11.4 (2014): 396.
