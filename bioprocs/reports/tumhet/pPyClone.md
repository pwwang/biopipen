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

### Cluster
::: {.tab}
#### Density
![cluster-density]({{job.o.outdir}}/plots/cluster/density.svg)
:::

::: {.tab}
#### Parallel coordinates
![cluster-density]({{job.o.outdir}}/plots/cluster/parallel_coordinates.svg)
:::

::: {.tab}
#### Scatter
![scatter]({{job.o.outdir}}/plots/cluster/scatter.svg)
:::

::: {.tab}
#### Table
```table
file: {{job.o.outdir}}/tables/cluster.csv
```
:::

### Loci
::: {.tab}
#### Density
![cluster-density]({{job.o.outdir}}/plots/loci/density.svg)
:::

::: {.tab}
#### Parallel coordinates
![Cellular prevalence]({{job.o.outdir}}/plots/loci/parallel_coordinates.svg)
![VAF]({{job.o.outdir}}/plots/loci/vaf_parallel_coordinates.svg)
:::

::: {.tab}
#### Scatter
![Cellular prevalence]({{job.o.outdir}}/plots/loci/scatter.svg)
![VAF]({{job.o.outdir}}/plots/loci/vaf_scatter.svg)
:::

::: {.tab}
#### Similarity
![scatter]({{job.o.outdir}}/plots/loci/similarity_matrix.svg)
:::

::: {.tab}
#### Table
```table
file: {{job.o.outdir}}/tables/loci.csv
```
:::

{% else %}{# if path.exists #}

{# Show the last line of error messages #}
{{job.stderr | readlines | [-1]}}

{% endif %}{# if path.exists #}
::::
{% endfor %}{# for jobs #}

[1]: Roth, Andrew, et al. "PyClone: statistical inference of clonal population structure in cancer." Nature methods 11.4 (2014): 396.
