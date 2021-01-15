# {{report.title}}

{% from pathlib import Path %}

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

{% 	if args.plot %}
{% 		if args.plot.pca %}
![First 2 components from PCA analysis]({{job.o.outdir, '*.pca.png' | *glob1: first=False | all_sample_plot}})
{% 		endif %}
{% 		if args.plot.boxplot %}
![Boxplot of expressions of all samples]({{job.o.outdir, '*.boxplot.png' | *glob1: first=False | all_sample_plot}})
{% 		endif %}
{% 		if args.plot.violin %}
![Violin plot of expressions of all samples]({{job.o.outdir, '*.violin.png' | *glob1: first=False | all_sample_plot}})
{% 		endif %}
{% 		if args.plot.histogram %}
![Distribution of expressions of all samples]({{job.o.outdir, '*.histo.png' | *glob1: first=False | all_sample_plot}})
{% 		endif %}
{% endif %}

:::

{% 	if `job.o.outdir, '*.group.*.png' | *glob1 | isinstance: (Path, list)` %}
::: {.panel}
### Stats for groups

{% 		if args.plot %}
{% 			if args.plot.boxplot %}
![Boxplot of expressions of different groups]({{job.o.outdir, '*.group.boxplot.png' | *glob1}})
{% 			endif %}
{% 			if args.plot.violin %}
![Violin plot of expressions of different groups]({{job.o.outdir, '*.group.violin.png' | *glob1}})
{% 			endif %}
{% 			if args.plot.histogram %}
![Distribution of expressions of different groups]({{job.o.outdir, '*.group.histo.png' | *glob1}})
{% 			endif %}
{% 			if args.plot.qqplot %}
![QQ Plot of expressions of different groups]({{job.o.outdir, '*.group.qq.png' | *glob1}})
{% 			endif %}
:::
{% 		endif %}
{% 	endif %}

{% 	if `job.o.outdir, '*.batch.*.png' | *glob1 | isinstance: (Path, list)` %}
::: {.panel}
### Stats for batches

{% 		if args.plot %}
{% 			if args.plot.boxplot %}
![Boxplot of expressions of different batches]({{job.o.outdir, '*.batch.boxplot.png' | *glob1}})
{% 			endif %}
{% 			if args.plot.violin %}
![Violin plot of expressions of different groups]({{job.o.outdir, '*.batch.violin.png' | *glob1}})
{% 			endif %}
{% 			if args.plot.histogram %}
![Distribution of expressions of different groups]({{job.o.outdir, '*.batch.histo.png' | *glob1}})
{% 			endif %}
{% 			if args.plot.qqplot %}
![QQ Plot of expressions of different groups]({{job.o.outdir, '*.batch.qq.png' | *glob1}})
{% 			endif %}
{% 		endif %}
{% 	endif %}
:::
::::
{% endfor %}
