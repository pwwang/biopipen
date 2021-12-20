{% from "utils/misc.liq" import report_jobs, table_of_images -%}

<script>
    import { Image } from "@@";
</script>

{% python %}
def python_map(func, iter):
    return [func(elem) for elem in iter]
{% endpython %}

{%- macro report_job(job, h=1) -%}
{{ table_of_images(
    glob(joinpaths(job.out.outdir, "*", "pathway_heterogeneity.png")),
    list(python_map(basename, glob(joinpaths(job.out.outdir, "*")))),
) }}
{%- endmacro -%}

{%- macro head_job(job) -%}
{% assign config = job.in.configfile | read | toml_loads %}
{% assign name = config.name or stem(job.out.outdir) %}
<h1>{{name | escape}}</h1>
{%- endmacro -%}

{{ report_jobs(jobs, head_job, report_job) }}
