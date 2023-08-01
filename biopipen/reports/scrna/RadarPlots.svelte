{% from "utils/misc.liq" import report_jobs, table_of_images -%}
<script>
    import { Image } from "$libs";
</script>

{%- macro report_job(job, h=1) -%}
    {%- set section_file = job.out.outdir | joinpaths: "sections.toml" -%}
    {%- if section_file | exists -%}
        {%- set sections = section_file | config: "toml" -%}
        {%- for section, cases in sections.items() -%}
            <h{{h}}>{{section | escape}}</h{{h}}>
            {%- set imgs = [] -%}
            {%- for case in cases -%}
                {%- set img = job.out.outdir | joinpaths: case + ".png" -%}
                {%- set _ = imgs.append(img) -%}
            {%- endfor -%}
            {{ table_of_images(imgs) }}
        {%- endfor -%}
    {%- else -%}
        {%- for img in job.out.outdir | glob: "*.png" -%}
            <h{{h}}>{{img | stem | escape}}</h{{h}}>
            <Image src={{img | quote}} />
        {%- endfor -%}
    {%- endif -%}
{%- endmacro -%}

{%- macro head_job(job) -%}
<h1>{{job.in.srtobj | stem}}</h1>
{%- endmacro -%}

{{ report_jobs(jobs, head_job, report_job) }}
