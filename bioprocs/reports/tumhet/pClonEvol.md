# {{title}}

{% python from os import path %}
{% for job in jobs %}

:::: {.tab .force-collapse}
{% if len(jobs) > 1 %}
## {{job.i.mutfile | stem}}
{% endif %}

{% for modelfile in glob1(job.o.outdir, 'model-*.png', first = False) %}

::: {.tab}
{% if len(jobs) > 1 %}#{% endif %}## {{modelfile | stem}}
![{{modelfile | stem}}]({{modelfile}})
:::

{% endfor %}

::::

{% endfor %}