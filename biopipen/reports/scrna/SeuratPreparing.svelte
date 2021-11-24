{% from "utils/misc.liq" import report_jobs -%}
{% from_ os import path %}
<script>
    import { Image } from "@@";
</script>

{%- macro report_job(job, h=1) -%}
{% for plotdir in job.out.outdir | joinpaths: "*" | glob %}
{%    if path.isdir(plotdir) %}
<h{{h}}>{{plotdir | basename}}</h{{h}}>
{%      for img in plotdir | joinpaths: "*.png" | glob %}
<Image src={{img | quote}} />
{%      endfor %}
{%    endif %}
{% endfor %}
{%- endmacro -%}

{%- macro head_job(job) -%}
<h1>{{job.in.metafile | stem}}</h1>
{%- endmacro -%}

{{ report_jobs(jobs, head_job, report_job) }}
