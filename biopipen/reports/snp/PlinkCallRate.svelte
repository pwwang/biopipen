{% from "utils/misc.liq" import report_jobs -%}
<script>
    import { Image, Descr } from "$libs";
</script>

{%- macro report_job(job, h=1) -%}
    <h{{h+1}}>Sample Call Rate</h{{h+1}}>
    {%- for pngfile in job.out.outdir | glob: '*.samplecr.png' -%}
    <Descr>Cutoff: {{envs.samplecr}}</Descr>
    <Image src="{{pngfile}}" />
    {%- endfor -%}

    <h{{h+1}}>Variant Call Rate</h{{h+1}}>
    {%- for pngfile in job.out.outdir | glob: '*.varcr.png' -%}
    <Descr>Cutoff: {{envs.varcr}}</Descr>
    <Image src="{{pngfile}}" />
    {%- endfor -%}
{%- endmacro -%}

{%- macro head_job(job) -%}
<h1>Sample: {{job.in.cnrfile | stem0 }}</h1>
{%- endmacro -%}

{{ report_jobs(jobs, head_job, report_job) }}
