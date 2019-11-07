# {{title}}

All of the basic anlyses were performed by Maftools[1], which attempts to summarize, analyze, annotate and visualize MAF files in an efficient manner from either TCGA sources or any in-house studies as long as the data is in MAF format.

{% python from os import path %}
{% for job in jobs %}

:::: {.tab .force-collapse}

{% 	if len(jobs) > 1 %}
## {{job.i.indir | stem}}
{% 	endif %}

{% 	if len(jobs) > 1 %}#{% endif %}## Summary
![Summary of mutations]({{glob1(job.o.outdir, 'summary.png')}})

{% 	if path.isfile(glob1(job.o.outdir, 'titv.png')) %}

![Transition-transversion ratio]({{glob1(job.o.outdir, 'titv.png')}})

- Transitions are interchanges of two-ring purines (A &lt;-&gt; G) or of one-ring pyrimidines (C &lt;-&gt; T): they therefore involve bases of similar shape. \
- Transversions are interchanges of purine for pyrimidine bases, which therefore involve exchange of one-ring and two-ring structures.
{% 	endif %}

{% 	if path.isfile(glob1(job.o.outdir, 'genecloud.png')) %}
Gene clound: Size of each gene is proportional to the total number of samples in which it is mutated/altered. \

![Gene Cloud]({{glob1(job.o.outdir, 'genecloud.png')}})

{% 	endif %}

{% 	if path.isfile(glob1(job.o.outdir, 'tcgacomp.png')) %}
TCGA contains over 30 different cancer cohorts and median mutation load across them varies from as low as 7 per exome (Pheochromocytoma and Paraganglioma arising from Adrenal Gland) to as high as 315 per exome (Skin Cutaneoys Melanoma). It is informative to see how mutation load in given maf stands against TCGA cohorts. This draws distribution of variants compiled from over 10,000 WXS samples across 33 TCGA landmark cohorts. Plot generated is similar to the one described in [3]. \

![Tumor mutation burden compared with TCGA cohorts]({{glob1(job.o.outdir, 'tcgacomp.png')}})

{% 	endif %}


{% 	if path.isfile(glob1(job.o.outdir, 'rainfalls', '*.png')) %}
{% 		if len(jobs) > 1 %}#{% endif %}## Rainfall plots
Cancer genomes, especially solid tumors are characterized by genomic loci with localized hyper-mutations[4]. Such hyper mutated genomic regions can be visualized by plotting inter variant distance on a linear genomic scale. \
{% 		for rainfall in glob1(job.o.outdir, 'rainfalls', '*.png', first = False) %}

::: {.tab .force-tab}
{% 			if len(jobs) > 1 %}#{% endif %}### {{path.basename(rainfall)[:-13]}}

![Rainfall plot for sample {{path.basename(rainfall)[:-13]}}]({{rainfall}})

:::

{% 		endfor %}
{% 	endif %}

{% 	if path.isfile(glob1(job.o.outdir, 'oncoplot.png')) %}
{% 		if len(jobs) > 1 %}#{% endif %}## Oncoplot

![Oncoplot of top mutated genes]({{glob1(job.o.outdir, 'oncoplot.png')}})

NOTE: Variants annotated as Multi_Hit are those genes which are mutated more than once in the same sample.
{% 	endif %}

{% 	if path.isfile(glob1(job.o.outdir, 'lollipops', '*.png')) %}
{% 		if len(jobs) > 1 %}#{% endif %}## Mutations on protein domains
{% 		if path.isfile(glob1(job.o.outdir, 'pfam.png')) %}

::: {.tab .force-tab}
{% 			if len(jobs) > 1 %}#{% endif %}### Overall
This serves the purpose of knowing what domain in given cancer cohort, is most frequently affected. \

![Mutated domains]({{glob1(job.o.outdir, 'pfam.png')}})

```table
caption: Details of mutated domains
file: {{glob1(job.o.outdir, 'pfam.csv')}}
rows: 100
csvargs:
	delimiter: ','
```
:::

{% 		endif %}
{% 		for lollipop in glob1(job.o.outdir, 'lollipops', '*.png', first = False) %}

::: {.tab .force-tab}
{% 			if len(jobs) > 1 %}#{% endif %}### {{path.basename(lollipop)[:-13]}}
Lollipop plots are simple and most effective way showing mutation spots on protein structure. Many oncogenes have a preferential sites which are mutated more often than any other locus. These spots are considered to be mutational hot-spots and lollipop plots can be used to display them along with rest of the mutations. \

![Lollipop plot for gene {{path.basename(lollipop)[:-13]}}]({{lollipop}})

:::

{% 		endfor %}
{% 	endif %}

{% 	if path.isfile(glob1(job.o.outdir, 'oncodrive.png')) or path.isfile(glob1(job.o.outdir, 'pancan.png')) %}
{% 		if path.isfile(glob1(job.o.outdir, 'oncodrive.png')) %}
{% 		if len(jobs) > 1 %}#{% endif %}## Detecting cancer driver genes based on positional clustering

Cancer genes (driver) detection from a given MAF is a based on algorithm oncodriveCLUST[2]. Concept is based on the fact that most of the variants in cancer causing genes are enriched at few specific loci (aka hot-spots). This method takes advantage of such positions to identify cancer genes. \

![Driver genes]({{glob1(job.o.outdir, 'oncodrive.png')}})

{%  	endif %}
{% 		if path.isfile(glob1(job.o.outdir, 'pancan.png')) %}
{% 		if len(jobs) > 1 %}#{% endif %}## Comparison with pancancer mutated genes

[5] performed MutSigCV[6] analysis on 21 cancer cohorts and identified over 200 genes to be significantly mutated which consists of previously un-subscribed novel genes. Their results show only few genes are mutated in multiple cohort while many of them are tissue/cohort specific. We can compare mutSig results against this pan-can list of significantly mutated genes to see genes specifically mutated in given cohort. \

![Comparison with pancancer mutated genes]({{glob1(job.o.outdir, 'pancan.png')}})

{%  	endif %}
{%  endif %}

{% 	if path.isfile(glob1(job.o.outdir, 'signature.png')) %}
{% 		if len(jobs) > 1 %}#{% endif %}## Mutation signatures

Every cancer, as it progresses leaves a signature characterized by specific pattern of nucleotide substitutions. Studies[3] have shown such mutational signatures, derived from over 7000 cancer samples 5. Such signatures can be extracted by decomposing matrix of nucleotide substitutions, classified into 96 substitution classes based on immediate bases surrounding the mutated base. Extracted signatures can also be compared to those [validated COSMIC signatures](http://cancer.sanger.ac.uk/cosmic/signatures).

![Identified signatures]({{glob1(job.o.outdir, 'signature.png')}})

![Mutation signature similarities]({{glob1(job.o.outdir, 'signature-sim.png')}})
{%  endif %}

::::

{% endfor %}{# for job #}

[1]: Mayakonda, Anand, et al. "Maftools: efficient and comprehensive analysis of somatic variants in cancer." Genome research 28.11 (2018): 1747-1756.
[2]: Tamborero, David, Abel Gonzalez-Perez, and Nuria Lopez-Bigas. "OncodriveCLUST: exploiting the positional clustering of somatic mutations to identify cancer genes." Bioinformatics 29.18 (2013): 2238-2244.
[3]: Alexandrov, Ludmil B., et al. "Signatures of mutational processes in human cancer." Nature 500.7463 (2013): 415.
[4]: Leiserson, Mark DM, et al. "CoMEt: a statistical approach to identify combinations of mutually exclusive alterations in cancer." Genome biology 16.1 (2015): 160.
[5]: Lawrence, Michael S., et al. "Discovery and saturation analysis of cancer genes across 21 tumour types." Nature 505.7484 (2014): 495.
[6]: Lawrence, Michael S., et al. "Mutational heterogeneity in cancer and the search for new cancer-associated genes." Nature 499.7457 (2013): 214.
