# {{title}}

{% for job in jobs %}

::: {.tab}
{% if len(jobs) > 1 %}
## {{job.i.infile | stem}}
{% endif %}


![]({{glob1(job.o.outdir, '*.roc.png')}})

{% if report.get('detail') %}
```table
caption: Details of model performance
file: {{glob1(job.o.outdir, '*.result.txt')}}
```
{% endif %}

:::
{% endfor %}
