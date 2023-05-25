{% from "utils/misc.liq" import report_jobs, table_of_images -%}
{% from_ os import path %}
<script>
    import { Image } from "$libs";
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
        <p>Dimension (genes x cells) before filtering: {{job.outdir | joinpaths: "before-qc", "dim.txt" | read}}</p>
        <p>Dimension (genes x cells) after filtering: {{job.outdir | joinpaths: "after-qc", "dim.txt" | read}}</p>
    </Tile>

    <h{{h}}>QC features before filtering</h{{h}}>
    {% set qcimgs = job.outdir | glob: "before-qc", "*.png" %}
    {{ table_of_images(qcimgs, 3) }}

    <h{{h}}>QC features after filtering</h{{h}}>
    {% set qcimgs = job.outdir | glob: "after-qc", "*.png" %}
    {{ table_of_images(qcimgs, 3) }}

{%- endmacro -%}

{%- macro head_job(job) -%}
    <h1>{{job.in.metafile | stem}}</h1>
{%- endmacro -%}

{{ report_jobs(jobs, head_job, report_job) }}
