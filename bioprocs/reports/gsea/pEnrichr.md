# {{report.title}}

Enrichr[[1]] is an easy to use intuitive enrichment analysis web-based tool providing various types of visualization summaries of collective functions of gene lists. Here we use its API to perform the analysis of out gene sets.


{% python from pyppl.utils import always_list %}
{% python from os import path %}
{% for job in jobs %}

:::: {.panel}

## {{job.i.infile | stem}}

{% 	for lib in always_list(args.libs) %}

::: {.panel .tab}

### {{lib}}

{% if path.isfile(glob1(job.o.outdir, '*.%s.png' % lib)) %}
![Gene set enrichment analysis against {{lib}}]({{glob1(job.o.outdir, '*.%s.png' % lib)}})

```table
caption: Details of GSEA against {{lib}}
file: {{glob1(job.o.outdir, '*.%s.txt' % lib)}}
```
{% else %}
Nothing enriched.
{% endif %}

:::

{% 	endfor %}

::::

{% endfor %}


[[1]]: Kuleshov, Maxim V., et al. "Enrichr: a comprehensive gene set enrichment analysis web server 2016 update." Nucleic acids research 44.W1 (2016): W90-W97.
