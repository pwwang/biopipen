{% from "utils/misc.liq" import report_jobs, table_of_images -%}
{% from "utils/gsea.liq" import fgsea_report_script, fgsea_report, gsea_report -%}

<script>
{{ fgsea_report_script() }}
</script>

{%- macro report_job(job, h=2) -%}

{%  for cldir in job.out.outdir | glob: "*" | first | glob: '*' %}
<h{{h}}>{{ cldir | basename }}</h{{h}}>
{%      if envs.fgsea %}
{{        fgsea_report(cldir, h+1, envs, 10) }}
{%      else %}
{{        gsea_report(cldir, h+1, envs, 10) }}
{%      endif %}
{%  endfor %}

{%- endmacro -%}

{%- macro head_job(job) -%}
{%- set name = job.out.outdir | glob: '*' | first | stem -%}
<h1>{{name | escape}}</h1>
{%- endmacro -%}

{{ report_jobs(jobs, head_job, report_job) }}
