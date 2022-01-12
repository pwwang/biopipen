{% from "utils/misc.liq" import report_jobs, table_of_images -%}

<script>
    import { Image } from "@@";
</script>

{%- macro report_job(job, h=1) -%}
{% for ssdir in job.out.outdir | joinpaths: "*" | glob %}
<h{{h}}>{{ssdir | basename}}</h{{h}}>
{{ table_of_images(sorted(glob(joinpaths(ssdir, "*.png")))) }}
{% endfor %}
{%- endmacro -%}

{%- macro head_job(job) -%}
{% assign config = job.in.configfile | read | toml_loads %}
{% assign name = config.name or stem(job.out.outdir) %}
<h1>{{name | escape}}</h1>
{%- endmacro -%}

{{ report_jobs(jobs, head_job, report_job) }}
