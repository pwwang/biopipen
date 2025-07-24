{% from "utils/misc.liq" import report_jobs -%}
<script>
    import { Image } from "$libs";
</script>

{%- macro report_job(job, h=1) -%}
    {%- for pngfile in job.out.outdir | glob: '*.png' -%}
    <h{{h+1}}>{{pngfile | stem0 | title}}</h{{h+1}}>
    <Image src="{{pngfile}}" />
    {%- endfor -%}
{%- endmacro -%}

{%- macro head_job(job) -%}
<h1>
    Samples: {{job.in.segfiles | first | stem0 }}-etc
    ({{job.in.segfiles | len}})
</h1>
{%- endmacro -%}

{{ report_jobs(jobs, head_job, report_job) }}
