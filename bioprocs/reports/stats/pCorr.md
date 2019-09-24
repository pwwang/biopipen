# {{title}}
Using method: {{args.method | @capitalize}}

{% for job in jobs %}
:::: {.tab}

{% if len(jobs) > 1 %}
## {{job.i.infile | stem}}
{% endif %}

::: {.tab}
{% 	if len(jobs) > 1 %}#{% endif %}## Table
```table
file: {{ job.o.outfile | prefix | @append: ".mat.txt" }}
caption: {{args.method | @capitalize}} Correlation Coefficient between variables
```

:::

{% if args.plot %}
::: {.tab}
{% 	if len(jobs) > 1 %}#{% endif %}## Plot

![Heatmap of {{args.method | @capitalize}} Correlation Coefficient]({{ job.o.outfile | prefix | @append: ".png" }})

:::
{% endif %}

::::

{% endfor %}
