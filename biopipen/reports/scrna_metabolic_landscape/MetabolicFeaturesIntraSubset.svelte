{% from "utils/misc.liq" import report_jobs, table_of_images -%}
{% from "utils/gsea.liq" import fgsea_report, gsea_report -%}

<script>
  import { Image, DataTable, Descr } from "$libs";
  import { Tabs, Tab, TabContent, Accordion, AccordionItem, InlineNotification } from "$ccs";
</script>

{%- macro report_job(job, h=2) -%}
  {% if envs.fgsea %}
    {{ job | render_job: h=h }}
  {% else %}
    {%  for groupdir in job.out.outdir | glob: "*" %}
      <h{{h}}>{{groupdir | basename}}</h{{h}}>
      {%- set dsdirs = groupdir | glob: "*" -%}
      {% for dsdir in groupdir | glob: "*" %}
          <h{{h+1}}>{{ dsdir | basename }}</h{{h+1}}>
          {{ gsea_report(dsdir, h+2, envs, envs.top) }}
      {% endfor %}
    {% endfor %}
  {% endif %}
{%- endmacro -%}

{%- macro head_job(job) -%}
  <h1>{{job.in.sobjfile | stem | escape}}</h1>
{%- endmacro -%}

{{ report_jobs(jobs, head_job, report_job) }}
