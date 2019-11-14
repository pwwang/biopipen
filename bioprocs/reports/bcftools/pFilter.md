# {{title}}

The filters were applied using bcftools[1].

{% import yaml %}
{% from os import path %}
{% if args.stat %}
```table
caption: Number of varaints with filters.
csvargs:
	delimiter: ','
---
{% for i, job in enumerate(jobs) -%}
{%- assign fstable = glob1(job.outdir, '*.filterstats.txt') -%}
{%- for j, line in enumerate(readlines(fstable)) -%}
{%- if j > 0 or i == 0 %}
{{line.replace('\t', ',')}}
{%- endif -%}
{%- endfor -%}
{%- endfor %}
```
{%- assign tableannfile = open(glob1(jobs[0].outdir, '*.filterstats.txt.ann')) -%}
{%- assign tableann = yaml.safe_load(tableannfile) -%}
{%- python tableannfile.close() -%}
{%- if report.get('tableann') -%}
{%- python tableann.update(report.tableann) -%}
{%- endif -%}{# if report.get('tableann') #}
{%- for key, value in tableann.items() %}
- {{key}}: {{value}}
{%- endfor -%}
{% endif %}{# if args.stat #}

[1]: Li, Heng, et al. "The sequence alignment/map format and SAMtools." Bioinformatics 25.16 (2009): 2078-2079.