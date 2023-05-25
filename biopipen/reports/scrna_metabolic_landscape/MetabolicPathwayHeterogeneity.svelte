{% from "utils/misc.liq" import report_jobs, table_of_images -%}

<script>
    import { Image } from "$libs";
</script>

{%- macro report_job(job, h=2) -%}
  {%- for ssdir in job.out.outdir | glob: "*" -%}
  <h{{h}}>{{ ssdir | stem }}</h{{h}}>

  <h{{h+1}}>Metabolic pathways enriched in genes with highest contribution to the metabolic heterogeneities</h{{h+1}}>
  <Image src="{{ssdir | joinpaths: 'pathway_heterogeneity.png'}}" />

  {%- endfor -%}
{%- endmacro -%}

{%- macro head_job(job) -%}
  <h1>{{job.in.sobjfile | stem | escape}}</h1>
{%- endmacro -%}

{{ report_jobs(jobs, head_job, report_job) }}
