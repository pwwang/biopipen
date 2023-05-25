{% from "utils/misc.liq" import report_jobs, table_of_images -%}
{% from "utils/gsea.liq" import fgsea_report, gsea_report -%}

<script>
  import { Image, DataTable } from "$libs";
</script>

{%- macro report_job(job, h=2) -%}
  {%  for groupdir in job.out.outdir | glob: "*" %}
    <h{{h}}>{{groupdir | basename}}</h{{h}}>
    {%- set dsdirs = groupdir | glob: "*" -%}
    {% for dsdir in groupdir | glob: "*" %}
        <h{{h+1}}>{{ dsdir | basename }}</h{{h+1}}>
        {% if envs.fgsea %}
            {% if dsdir | joinpaths: "fgsea.txt" | isfile %}
                {{ fgsea_report(dsdir, h+2, envs, envs.top) }}
            {% else %}
                <p>Not enough events.</p>
            {% endif %}
        {% else %}
            {{ gsea_report(dsdir, h+2, envs, envs.top) }}
        {% endif %}
    {% endfor %}
  {% endfor %}
{%- endmacro -%}

{%- macro head_job(job) -%}
  <h1>{{job.in.sobjfile | stem | escape}}</h1>
{%- endmacro -%}

{{ report_jobs(jobs, head_job, report_job) }}
