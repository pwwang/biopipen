# {{title}}

{% for job in jobs %}
::::: {.tab}
{% if len(jobs) > 1 %}
## {{job.i.infile | stem}}
{% else %}
## ROC curve
{% endif %}

![ROC curve]({{job.o.outfile}})
:::::
{% endfor %}