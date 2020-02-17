# {{report.title}}

{% for job in jobs %}
:::: {.panel}

{%	if len(jobs) > 1 %}
## {{ job.i.infiles | [0] | stem }}
{%	endif %}

![Venn/UpSet Plot]({{job.outdir, '*.png' | *glob1}})

{%	if report.get('details') %}
```table
file: {{job.outdir, 'data.txt' | *glob1}}
caption: Details in the plot
```
{%	endif %}

::::
{% endfor %}
