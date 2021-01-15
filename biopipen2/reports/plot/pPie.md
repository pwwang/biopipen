# {{report.title}}

{% for job in jobs %}
::::: {.tab}
{% if len(jobs) > 1 %}
## {{job.i.infile | stem}}
{% else %}
## Pie chart
{% endif %}

![Pie chart]({{job.o.outfile}})
:::::
{% endfor %}
