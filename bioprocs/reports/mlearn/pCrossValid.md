# {{title}}

{% python from os import path %}
{% for job in jobs %}
::: {.tab}
{% if len(jobs) > 1 %}
## {{job.i.infile | stem}}
{% endif %}


{% assign modelfigs = glob1(job.o.outdir, '*.model.png', first = False) %}
{% if report.get('modelfig', True) and modelfigs %}
![The model]({{modelfigs[0]}})
{% endif %}

{% assign rocfigs = glob1(job.o.outdir, '*.roc.png', first = False) %}
{% if report.get('roc', True) and rocfigs %}
![ROC curve]({{rocfigs[0]}})
{% endif %}

{% if report.get('auc', True) %}
```table
caption: AUCs
file: {{glob1(job.o.outdir, '*.aucs.txt')}}
```
{% endif %}

{% assign varimpfigs = glob1(job.o.outdir, '*.varimp.png', first = False) %}
{% if report.get('varimpfig', True) and varimpfigs %}
![Feature importance]({{varimpfigs[0]}})
{% endif %}

{% if report.get('varimp', True) %}
```table
caption: Feature importance
file: {{glob1(job.o.outdir, '*.varimp.txt')}}
```
{% endif %}

:::
{% endfor %}{# for job #}