# {{title}}

{% python from os import path %}
{% for job in jobs %}
::: {.tab}

{% if len(jobs) > 1 %}
## {{job.i.infile | stem}}
{% endif %}

{% if path.exists(glob1(job.o.outdir, '*.combined.png')) %}
![Survival plot]({{glob1(job.o.outdir, '*.combined.png')}})\
{% else %}

{% for survplot in glob1(job.o.outdir, '*.surv.png', first = False) %}
![{{survplot.split('.')[0]}}]({{survplot}})
{% endfor %}

{% endif %}{# if path.exists #}

```table
file: {{glob1(job.o.outdir, '*.survival.txt')}}
caption: Details
```

:::
{% endfor %}
