# {{title}}

Enrichr[1] is an easy to use intuitive enrichment analysis web-based tool providing various types of visualization summaries of collective functions of gene lists. Here we use its API to perform the analysis of out gene sets.


{% python from pyppl.utils import alwaysList %}
{% python from os import path %}
{% for job in jobs %}

:::: {.tab .force-collapse}

{% if len(jobs) > 1 %}
## {{job.i.infile | stem}}
{% endif %}

{% for lib in alwaysList(args.libs) %}

::: {.tab .force-tab}

{% if len(jobs) > 1 %}#{% endif %}## {{lib}}

{% if path.isfile(glob1(job.o.outdir, '*.%s.png' % lib)) %}
![Gene set enrichment analysis against {{lib}}]({{glob1(job.o.outdir, '*.%s.png' % lib)}})

```table
caption: Details of GSEA against {{lib}}
file: {{glob1(job.o.outdir, '*.%s.txt' % lib)}}
dtargs:
	order: [[3, 'asc']]
```
{% else %}
Nothing enriched.

{% endif %}

:::

{% endfor %}

::::

{% endfor %}


[1]: Kuleshov, Maxim V., et al. "Enrichr: a comprehensive gene set enrichment analysis web server 2016 update." Nucleic acids research 44.W1 (2016): W90-W97.
