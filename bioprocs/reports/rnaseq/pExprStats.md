# {{report.title}}

{% for job in jobs %}
:::: {.panel}
## {{ job.i.infile | stem | stem | @append: "-" + stem(job.i.gfile) }}

::: {.panel}
### Stats for all samples

{% assign all_sample_plot = lambda plots: [
	x
	for x in plots
	if '.batch.' not in str(x) and '.group.' not in str(x)
][0] %}

{% if args.plot %}
![First 2 components from PCA analysis]({{job.o.outdir, '*.pca.png' | *glob1: first=False | all_sample_plot}})
![Boxplot of expressions of all samples]({{job.o.outdir, '*.boxplot.png' | *glob1: first=False | all_sample_plot}})
![Violin plot of expressions of all samples]({{job.o.outdir, '*.violin.png' | *glob1: first=False | all_sample_plot}})
![Distribution of expressions of all samples]({{job.o.outdir, '*.histo.png' | *glob1: first=False | all_sample_plot}})
{% endif %}

:::

{% 	if `job.o.outdir, '*.group.*.png' | *glob1` %}
::: {.panel}
### Stats for groups

![Boxplot of expressions of different groups]({{job.o.outdir, '*.group.boxplot.png' | *glob1}})
![Violin plot of expressions of different groups]({{job.o.outdir, '*.group.violin.png' | *glob1}})
![Distribution of expressions of different groups]({{job.o.outdir, '*.group.histo.png' | *glob1}})
![QQ Plot of expressions of different groups]({{job.o.outdir, '*.group.qq.png' | *glob1}})
:::
{% 	endif %}

{% 	if `job.o.outdir, '*.batch.*.png' | *glob1` %}
::: {.panel}
### Stats for batches

![Boxplot of expressions of different batches]({{job.o.outdir, '*.batch.boxplot.png' | *glob1}})
![Violin plot of expressions of different groups]({{job.o.outdir, '*.batch.violin.png' | *glob1}})
![Distribution of expressions of different groups]({{job.o.outdir, '*.batch.histo.png' | *glob1}})
![QQ Plot of expressions of different groups]({{job.o.outdir, '*.batch.qq.png' | *glob1}})
{% 	endif %}
:::
::::
{% endfor %}
