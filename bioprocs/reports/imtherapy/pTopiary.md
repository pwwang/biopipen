# {{title}}

{% python from pathlib import Path %}
{% for job in jobs %}

:::: {.tab}
## {{job.i.infile | stem}}

```table
file: {{Path(job.o.outfile).with_suffix('.neat.txt') | str}}
caption: Prioritized neoantigens by given HLA alleles
rows: 100
```

::::
{% endfor %}