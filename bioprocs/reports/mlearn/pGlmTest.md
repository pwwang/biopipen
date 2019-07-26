# {{title}}

{% for job in jobs %}

::: {.tab}
{% if len(jobs) > 1 %}
## {{job.i.infile | stem}}
{% endif %}


![]({{glob(job.o.outdir, '*.roc.png')}})

{% if report.get('detail') %}
```table
caption: Details of model performance
file: {{glob(job.o.outdir, '*.result.txt')}}
```
{% endif %}

:::
{% endfor %}
