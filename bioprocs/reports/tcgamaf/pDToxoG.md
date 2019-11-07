# {{title}}

The OxoG artifact from SNV calls was removed by D-ToxoG[1]. The goal of D-ToxoG is to limit the output mutation calls to less than {{args.params.artifactThresholdRate}} artifact by following steps:

1. Make an estimate of the number of OxoG artifacts in a given MAF file \
2. Calculate the likelihood that a given SNV is an OxoG artifact \
3. Apply threshold ({{args.params.artifactThresholdRate}}) to false discovery rate \

{% for job in jobs %}
:::: {.tab}

{% if len(jobs) > 1 %}
## {{job.i.infile | stem}}
{% endif %}

{% if report.get('download') %}
[Download MAF file]({{job.o.outfile}} "file-download")
{% endif %}

::: {.tab}
### All cases
```table
file: {{glob1(job.outdir, 'caseTableData.tsv')}}
caption: Summary of all cases
```

![Lego plot before filtering for all cases]({{glob1(job.outdir, 'figures', '*.dtoxog.maf*_lego_before.png')}})

![Lego plot after filtering for all cases]({{glob1(job.outdir, 'figures', '*.dtoxog.maf*_lego_after.png')}})

:::

{% for legoplot in glob1(job.outdir, 'figures', '*_lego_before.png', first = False) %}
{% if '.dtoxog.maf' not in legoplot %}
::: {.tab}
### {{legoplot | bn | [:-16]}}

![Lego plot before filtering for {{legoplot | bn | [:-16]}}]({{legoplot}})


![Lego plot after filtering]({{legoplot | bn | [:-16] | @append: "_lego_after.png" | @prepend: dirname(legoplot) + '/'}})
:::
{% endif %}{# if not legoplot.endswith #}
{% endfor %}{# for legoplot #}

::::
{% endfor %}{# for job #}

[1]: Costello, Maura, et al. "Discovery and characterization of artifactual mutations in deep coverage targeted capture sequencing data due to oxidative DNA damage during sample preparation." Nucleic acids research 41.6 (2013): e67-e67.
