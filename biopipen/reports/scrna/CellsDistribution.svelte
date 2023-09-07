{% from "utils/misc.liq" import report_jobs -%}
<script>
    import { Image } from "$libs";
</script>

{%- macro report_job(job, h=1) -%}
    {%- for secdir in job.out.outdir | glob: "*" -%}
        {%- set sec = secdir | basename -%}

        {%- if sec != "DEFAULT" -%}
            <h{{h}}>{{sec | escape}}</h{{h}}>
        {%- else -%}
            {%- set h = h - 1 -%}
        {%- endif -%}

        {%- for plotfile in secdir | glob: "*.png" -%}
            <h{{h+1}}>{{ plotfile | stem }}</h{{h+1}}>
            <Image src={{plotfile | quote}} />
        {%- endfor -%}
    {%- endfor -%}
{%- endmacro -%}

{%- macro head_job(job) -%}
<h1>{{job.in.srtobj | stem0 | escape}}</h1>
{%- endmacro -%}

{{ report_jobs(jobs, head_job, report_job) }}
