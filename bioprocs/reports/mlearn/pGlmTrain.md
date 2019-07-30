# {{title}}

{% for job in jobs %}

::: {.tab}
{% if len(jobs) > 1 %}
## {{job.i.infile | stem}}
{% endif %}

```table
caption: Model features
file: {{glob1(job.o.outdir, '*.features.txt')}}
```

{% if args.plot %}
![]({{glob1(job.o.outdir, '*.glm.png')}})
{% endif %}

:::
{% endfor %}
