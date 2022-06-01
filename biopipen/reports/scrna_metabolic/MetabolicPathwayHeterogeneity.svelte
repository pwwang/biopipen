{% from "utils/misc.liq" import report_jobs, table_of_images -%}

<script>
    import { Image } from "@@";
</script>


{%- macro report_job(job, h=1) -%}
<h{{h}}>Metabolic pathways enriched in genes with highest contribution to the metabolic heterogeneities</h{{h}}>
<Image src="{{job.out.outdir | glob: '*' | first | joinpaths: 'pathway_heterogeneity.png'}}" />
{%- endmacro -%}

{%- macro head_job(job) -%}
{%- set name = job.out.outdir | glob: "*" | first | stem -%}
<h1>{{name | escape}}</h1>
{%- endmacro -%}

{{ report_jobs(jobs, head_job, report_job) }}
