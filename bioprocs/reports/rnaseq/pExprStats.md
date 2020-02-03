# {{report.title}}

{% for job in jobs %}
:::: {.tab}
	{%- if forloop.length > 1 %}
## {{ job.i.efile | stem | stem | @append: "-" + stem(job.i.gfile) }}
	{%- endif %}

##{% if forloop.length > 1 %}#{%endif%} Stats for all samples

{% assign all_sample_plot = lambda plots: [
	x
	for x in plots
	if '.batch.' not in str(x) and '.group.' not in str(x)
][0] %}
![PCA]({{job.o.outdir, '*.pca.png' | *glob1: first=False | all_sample_plot}})

![Boxplot]({{job.o.outdir, '*.boxplot.png' | *glob1: first=False | all_sample_plot}})

![Violin Plot]({{job.o.outdir, '*.violin.png' | *glob1: first=False | all_sample_plot}})

![Histogram]({{job.o.outdir, '*.histo.png' | *glob1: first=False | all_sample_plot}})

{% 	if `job.o.outdir, '*.group.*.png' | *glob1` %}
##{% if forloop.length > 1 %}#{%endif%} Stats for groups

![Boxplot]({{job.o.outdir, '*.group.boxplot.png' | *glob1}})

![Violin Plot]({{job.o.outdir, '*.group.violin.png' | *glob1}})

![Histogram]({{job.o.outdir, '*.group.histo.png' | *glob1}})

![QQ Plot]({{job.o.outdir, '*.group.qq.png' | *glob1}})
{% 	endif %}

{% 	if `job.o.outdir, '*.batch.*.png' | *glob1` %}
##{% if forloop.length > 1 %}#{%endif%} Stats for batches

![Boxplot]({{job.o.outdir, '*.batch.boxplot.png' | *glob1}})

![Violin Plot]({{job.o.outdir, '*.batch.violin.png' | *glob1}})

![Histogram]({{job.o.outdir, '*.batch.histo.png' | *glob1}})

![QQ Plot]({{job.o.outdir, '*.batch.qq.png' | *glob1}})
{% 	endif %}

::::
{% endfor %}
