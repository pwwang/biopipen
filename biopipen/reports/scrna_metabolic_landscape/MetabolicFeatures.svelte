{% from "utils/misc.liq" import report_jobs, table_of_images -%}
{% from "utils/gsea.liq" import fgsea_report_script, fgsea_report, gsea_report -%}

<script>
{{ fgsea_report_script() }}
</script>

{%- macro report_job(job, h=2) -%}
  {%- for ssdir in job.out.outdir | glob: "*" -%}
    {%- if basename(ssdir) == "ALL" -%}
      {%- set h = 1 -%}
    {%- else -%}
      <h{{h}}>{{ ssdir | stem }}</h{{h}}>
    {%- endif -%}

    {% for cldir in ssdir | glob: '*' %}
      <h{{h+1}}>{{ cldir | basename }}</h{{h+1}}>
      {% if envs.fgsea %}
          {{ fgsea_report(cldir, h+2, envs, envs.top) }}
      {% else %}
          {{ gsea_report(cldir, h+2, envs, envs.top) }}
      {% endif %}
    {% endfor %}
  {%- endfor -%}
{%- endmacro -%}

{%- macro head_job(job) -%}
  <h1>{{job.in.sobjfile | stem | escape}}</h1>
{%- endmacro -%}

{{ report_jobs(jobs, head_job, report_job) }}
