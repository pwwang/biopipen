# {{title}}
{% python from pathlib import Path %}
{% for job in jobs %}

{%  if len(jobs) == 1 %}

{%  else %}
:::: {.tab}

## {{job.i.infile | stem}}

{% for sumfile in sorted(job.o.outdir.glob('summary-*.png')) %}
::: {.tab .force-tab}
### {{sumfile.stem}}
![{{sumfile.stem}}]({{sumfile}})
:::
{% endfor %}

::::
{%  endif %}

{% endfor %}