{% from "utils/misc.liq" import report_jobs, table_of_images -%}
{% from "utils/gsea.liq" import fgsea_report, gsea_report -%}

<script>
    import { Image, DataTable } from "@@";
</script>

{% python %}
def python_map(func, iter):
    return [func(elem) for elem in iter]
{% endpython %}

{%- macro report_job(job, h=1) -%}

{%  for subset_dir in job.out.outdir | joinpaths: "*" | glob %}
<h{{h}}>{{subset_dir | basename}}</h{{h}}>
{%      for groupdir in subset_dir | joinpaths: "*" | glob %}
<h{{h+1}}>{{ groupdir | basename }}</h{{h+1}}>
{%          if envs.fgsea %}
{{            fgsea_report(groupdir, h+2, envs, 8) }}
{%          else %}
{{            gsea_report(groupdir, h+2, envs, 8) }}
{%          endif %}
{%      endfor %}
{%  endfor %}

{%- endmacro -%}

{%- macro head_job(job) -%}
{% assign config = job.in.configfile | read | toml_loads %}
{% assign name = config.name or stem(job.out.outdir) %}
<h1>{{name | escape}}</h1>
{%- endmacro -%}

{{ report_jobs(jobs, head_job, report_job) }}
