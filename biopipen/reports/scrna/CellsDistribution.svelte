{% from "utils/misc.liq" import report_jobs -%}
<script>
    import { Image } from "@@";
</script>

{%- macro report_job(job, h=1) -%}
{%- for casefile in job.out.outdir | glob: "*.png" -%}
<h{{h}}>{{casefile | stem | escape}}</h{{h}}>
<Image src={{casefile | quote}} />
{%- endfor -%}
{%- endmacro -%}

{%- macro head_job(job) -%}
{%- if job.in.casefile -%}
{%-     set cases = job.in.casefile | config: "toml" -%}
{%- else -%}
{%-     set cases = envs -%}
{%- endif -%}
<h1>{{cases | attr: "name" | escape}}</h1>
{%- endmacro -%}

{{ report_jobs(jobs, head_job, report_job) }}
