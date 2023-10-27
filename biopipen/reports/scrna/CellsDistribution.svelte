{% from "utils/misc.liq" import report_jobs, table_of_images -%}
<script>
    import { Image } from "$libs";
</script>

{%- macro report_job(job, h=1) -%}
    {%- set secdirs = job.out.outdir | glob: "*" -%}
    {%- if len(secdirs) == 1 -%}
        {%- set secname = secdirs | first | basename -%}
        {%- if secdirs[0] | joinpaths: "venn.png" | exists -%}
            {%- if secname == "DEFAULT" -%}
                <h{{h}}>Case overlapping</h{{h}}>
            {%- else -%}
                <h{{h}}>{{ secname | escape }} - Case overlapping</h{{h}}>
            {%- endif -%}
            {{ table_of_images(
                [joinpaths(secdirs[0], "venn.png"), joinpaths(secdirs[0], "upset.png")],
                ["Venn plot", "Upset plot"]) }}
        {%- endif -%}
        {%- for plotfile in secdirs[0] | glob: "case-*.png" -%}
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
            {%- if secdir | joinpaths: "venn.png" | exists -%}
                <h{{h+1}}>Case overlapping</h{{h+1}}>
                {{ table_of_images(
                    [joinpaths(secdir, "venn.png"), joinpaths(secdir, "upset.png")],
                    ["Venn plot", "Upset plot"]) }}
            {%- endif -%}
            {%- for plotfile in secdir | glob: "case-*.png" -%}
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
