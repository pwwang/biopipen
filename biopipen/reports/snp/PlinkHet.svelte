{% from "utils/misc.liq" import report_jobs -%}
<script>
    import { Image, Descr } from "$libs";
</script>

{%- macro report_job(job, h=1) -%}
    {%- for pngfile in job.out.outdir | glob: '*.png' -%}
    <h{{h+1}}>Distribution</h{{h+1}}>
    <Descr>Cutoff: [mean - {{envs.cutoff}} x sd, mean + {{envs.cutoff}} x sd]</Descr>
    <Image src="{{pngfile}}" />
    {%- endfor -%}
{%- endmacro -%}

{%- macro head_job(job) -%}
<h1>Sample: {{job.in.cnrfile | stem0 }}</h1>
{%- endmacro -%}

{{ report_jobs(jobs, head_job, report_job) }}
