{% from "utils/misc.liq" import report_jobs, table_of_images -%}

<script>
    import { Image } from "@@";
</script>

{%- macro report_job(job, h=2) -%}
{%- set name = job.in.configfile | config: "toml" | attr: "name" -%}
{%- if name or proc.size > 1 -%}
{%- else -%}
{%- set h = 1 -%}
{%- endif -%}

{%- for ssdir in job.out.outdir | glob: "*" -%}
<h{{h}}>{{ ssdir | stem }}</h{{h}}>

<h{{h+1}}>Metabolic pathways enriched in genes with highest contribution to the metabolic heterogeneities</h{{h+1}}>
<Image src="{{job.out.outdir | glob: '*' | first | joinpaths: 'pathway_heterogeneity.png'}}" />

{%- endfor -%}
{%- endmacro -%}

{%- macro head_job(job) -%}
{%- set name = job.in.configfile | config: "toml" | attr: "name" -%}
{%- if name or proc.size > 1 -%}
<h1>{{name | default: "Job #" + str(job.index+1) | escape}}</h1>
{%- endif -%}
{%- endmacro -%}

{{ report_jobs(jobs, head_job, report_job) }}
