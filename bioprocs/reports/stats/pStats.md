# {{title}}
## Glimpse

{% mode loose %}
{% python from os import path %}
{% for job in jobs %}
:::: {.tab}

{% if len(jobs) > 1 %}
### {{job.i.infile | stem}}
{% endif %}

```table
caption: First 10 samples
file: {{glob1(job.o.outdir, "*.top10.txt")}}
```

::::
{% endfor %}

## Statistics

{% for job in jobs %}
:::: {.tab}

{% if len(jobs) > 1 %}
### {{job.i.infile | stem}}
{% endif %}

{% for i, statfile in enumerate(sorted(glob1(job.o.outdir, "feature-stat.*.txt", first = False))) %}
::: {.tab}
{% assign feature = statfile | stem | [13:] %}
{% if len(jobs) > 1 %}#{%endif%}### {{feature}}

```table
caption: Basic statistics
file: {{statfile}}
```

![{{feature}} distribution]({{path.join(job.o.outdir, "feature-plot." + feature + ".png")}})

{% assign testfile = path.join(job.o.outdir, "feature-test." + feature + ".txt") %}
{% if path.isfile(testfile) %}
```table
caption: Tests
file: {{testfile}}
```

{% if i == 0 %}
{# Explain tests #}
- T test: Used to determine if there is a significant difference between the means of two groups. [Ref](https://en.wikipedia.org/wiki/Student%27s_t-test)
- Wilcoxon rank-sum test (a.k.a: Mann–Whitney U test): Used to determine whether two independent samples were selected from populations having the same distribution. [Ref](https://en.wikipedia.org/wiki/Mann%E2%80%93Whitney_U_test)
- ANOVA: Provides a statistical test of whether two or more population means are equal. [Ref](https://en.wikipedia.org/wiki/Analysis_of_variance)
- Kruskal Wallis test: Tests whether the mean ranks are the same in all the groups. [Ref](https://en.wikipedia.org/wiki/Kruskal%E2%80%93Wallis_one-way_analysis_of_variance)
- Chi-square test: Used to determine whether there is a significant difference between the expected frequencies and the observed frequencies in one or more categories. [Ref](https://en.wikipedia.org/wiki/Chi-squared_test)
- Fisher's exact test: Used to determine if there are nonrandom associations between two categorical variables. [Ref](https://en.wikipedia.org/wiki/Fisher%27s_exact_test)
{% endif %}

{% endif %}
:::
{% endfor %}


::::
{% endfor %}