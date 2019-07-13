## Clonanity analysis using PyClone[1]

{% python from os import path %}
{% for job in jobs %}
:::: {.tab}
### Job {{job.index}}

![Similarity matrix]({{path.join(job.o.outdir, "plots/loci/similarity_matrix.svg")}})

```table
caption: Clusters
file: "{{path.join(job.o.outdir, "tables/cluster.tsv")}}"
rows: 10
```

::::
{% endfor %}

[1]: Roth, Andrew, et al. "PyClone: statistical inference of clonal population structure in cancer." Nature methods 11.4 (2014): 396.
