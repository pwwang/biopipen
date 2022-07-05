{% from "utils/misc.liq" import report_jobs, table_of_images -%}
{% from "utils/gsea.liq" import fgsea_report_script, fgsea_report, gsea_report -%}

<script>
{{ fgsea_report_script() }}
</script>

{%- macro report_job(job, h=2) -%}
{%- set name = job.in.configfile | config: "toml" | attr: "name" -%}
{%- if name or proc.size > 1 -%}
{%- else -%}
{%- set h = 1 -%}
{%- endif -%}

{%- for ssdir in job.out.outdir | glob: "*" -%}
<h{{h}}>{{ ssdir | stem }}</h{{h}}>

{%  for cldir in ssdir | glob: '*' %}
<h{{h+1}}>{{ cldir | basename }}</h{{h+1}}>
{%      if envs.fgsea %}
{{        fgsea_report(cldir, h+2, envs, 10) }}
{%      else %}
{{        gsea_report(cldir, h+2, envs, 10) }}
{%      endif %}
{%  endfor %}

{%- endfor -%}
{%- endmacro -%}

{%- macro head_job(job) -%}
{%- set name = job.in.configfile | config: "toml" | attr: "name" -%}
{%- if name or proc.size > 1 -%}
<h1>{{name | default: "Job #" + str(job.index+1) | escape}}</h1>
{%- endif -%}
{%- endmacro -%}

{{ report_jobs(jobs, head_job, report_job) }}
