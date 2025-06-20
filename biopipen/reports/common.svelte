{% from "utils/misc.liq" import report_jobs, table_of_images -%}
<script>
    import { Image, DataTable, Descr } from "$libs";
    import { Tabs, Tab, TabContent, UnorderedList, ListItem, InlineNotification } from "$ccs";
</script>

{%- macro report_job(job, h=1) -%}
    {{ job | render_job: h=h }}
{%- endmacro -%}

{%- macro head_job(job) -%}
    <h1>{{job.in | attr: "values" | call | first | stem0 | escape}}</h1>
{%- endmacro -%}

{{ report_jobs(jobs, head_job, report_job) }}
