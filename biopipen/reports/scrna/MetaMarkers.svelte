{% from "utils/misc.liq" import report_jobs, table_of_images -%}
{% from "utils/gsea.liq" import enrichr_report -%}
<script>
    import { Image, DataTable, Descr } from "$libs";
    import { Tabs, Tab, TabContent, InlineNotification, Accordion, AccordionItem } from "$ccs";
</script>


{%- macro report_job(job, h=1) -%}
    {{ job | render_job: h=h }}
{%- endmacro -%}


{%- macro head_job(job) -%}
  <h1>{{job.in.srtobj | stem | escape}}</h1>
{%- endmacro -%}

{{ report_jobs(jobs, head_job, report_job) }}