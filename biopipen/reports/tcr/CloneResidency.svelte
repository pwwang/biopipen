{% from "utils/misc.liq" import report_jobs, table_of_images -%}
<script>
    import { Image } from "@@";
</script>

{%- macro report_job(job, h=1) -%}
<h{{h}}>Residency plots</h{{h}}>

{% assign scatter_pngs = job.out.outdir | joinpaths: "scatter", "scatter_*.png" | glob %}
{{ table_of_images(scatter_pngs, col=3) }}

<h{{h}}>Clonotype overlapping</h{{h}}>
{% assign venn_pngs = job.out.outdir | joinpaths: "venn", "*.png" | glob %}
{{ table_of_images(venn_pngs, col=3) }}

{%- endmacro -%}

{%- macro head_job(job) -%}
<h1>{{job.out.outdir | stem | replace: ".immunarch", ""}}</h1>
{%- endmacro -%}

{{ report_jobs(jobs, head_job, report_job) }}

