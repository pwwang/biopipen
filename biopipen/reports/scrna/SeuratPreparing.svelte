{% from "utils/misc.liq" import report_jobs, table_of_images -%}
{% from_ os import path %}
<script>
    import { Image, DataTable } from "$libs";
    import { Tile } from "$ccs";
</script>

{%- macro report_job(job, h=1) -%}
    <h{{h}}>Applied filters</h{{h}}>
    <Tile>
        <p>Cell filters: {{envs.cell_qc | str | escape}}</p>
        <p>Gene filters: {{
            proc.envs.gene_qc
            | str
            | replace: "{", "&#123"
            | replace: "}", "&#125"
        }}</p>

        <DataTable
            src={{job.outdir | joinpaths: 'plots', 'dim.txt' | quote}}
            data={ {{job.outdir | joinpaths: 'plots', 'dim.txt' | datatable: sep="\t"}} } />
    </Tile>

    <h{{h}}>Violin plots</h{{h}}>
    {% set qcimgs = job.outdir | glob: "plots", "*.vln.png" %}
    {{ table_of_images(qcimgs) }}

    <h{{h}}>Scatter plots</h{{h}}>
    {% set qcimgs = job.outdir | glob: "plots", "*.scatter.png" %}
    {{ table_of_images(qcimgs) }}

{%- endmacro -%}

{%- macro head_job(job) -%}
    <h1>{{job.in.metafile | stem}}</h1>
{%- endmacro -%}

{{ report_jobs(jobs, head_job, report_job) }}
