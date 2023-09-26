{% from "utils/misc.liq" import report_jobs, table_of_images -%}

<script>
    import { Image, DataTable } from "$libs";
</script>

{%- macro report_job(job, h=1) -%}

    {% if envs.stats %}
    <h{{h}}>Sample Information</h{{h}}>
    {% endif %}
    {% if envs.exclude_cols and isinstance(envs.exclude_cols, str) %}
        {% set excluded_cols = envs.exclude_cols | replace: " ", "" | split: "," %}
    {% else %}
        {% set excluded_cols = envs.exclude_cols %}
    {% endif %}

    <DataTable
        data={ {{ job.out.outfile | datatable: sep="\t", excluded=excluded_cols }} }
        pageSize={50}
        />

    {% if envs.stats %}
        <h{{h}}>Statistics</h{{h}}>
        {%- set stat_imgs = job.outdir | glob: "*.png" -%}
        {{- table_of_images(stat_imgs) -}}
    {% endif %}

{%- endmacro -%}

{%- macro head_job(job) -%}
    <h1>{{job.in.infile | stem | escape }}</h1>
{%- endmacro -%}

{{ report_jobs(jobs, head_job, report_job) }}

