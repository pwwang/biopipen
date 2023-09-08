{% from "utils/misc.liq" import report_jobs -%}
<script>
    import { Image } from "$libs";
</script>

{%- macro report_job(job, h=1) -%}
    {%- set secdirs = job.out.outdir | glob: "*" -%}
    {%- if len(secdirs) == 1 -%}
        {%- set secname = secdirs | first | basename -%}
        {%- for plotfile in secdirs[0] | glob: "*.png" -%}
            {%- if secname == "DEFAULT" -%}
                <h{{h}}>{{ plotfile | stem | escape }}</h{{h}}>
            {%- else -%}
                <h{{h}}>{{ secname | escape }} - {{ plotfile | stem | escape }}</h{{h}}>
            {%- endif -%}
            <Image src={{plotfile | quote}} />
        {%- endfor -%}
    {%- else -%}
        {%- for secdir in secdirs -%}
            {%- set sec = secdir | basename -%}
            <h{{h}}>{{sec | escape}}</h{{h}}>
            {%- for plotfile in secdir | glob: "*.png" -%}
                <h{{h+1}}>{{ plotfile | stem }}</h{{h+1}}>
                <Image src={{plotfile | quote}} />
            {%- endfor -%}
        {%- endfor -%}
    {%- endif -%}
{%- endmacro -%}

{%- macro head_job(job) -%}
<h1>{{job.in.srtobj | stem0 | escape}}</h1>
{%- endmacro -%}

{{ report_jobs(jobs, head_job, report_job) }}
