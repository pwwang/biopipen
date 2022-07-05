{% from "utils/misc.liq" import report_jobs, table_of_images -%}
{% from "utils/gsea.liq" import fgsea_report, gsea_report -%}

<script>
    import { Image, DataTable } from "@@";
</script>

{% python %}
def python_map(func, iter):
    return [func(elem) for elem in iter]
{% endpython %}

{%- macro report_job(job, h=2) -%}
{%- set name = job.in.configfile | config: "toml" | attr: "name" -%}
{%- if name or proc.size > 1 -%}
{%- else -%}
{%- set h = 1 -%}
{%- endif -%}

{%  for groupdir in job.out.outdir | glob: "*" %}
<h{{h}}>{{groupdir | basename}}</h{{h}}>
    {%- set dsdirs = groupdir | glob: "*" -%}
    {% for dsdir in groupdir | glob: "*" %}
        <h{{h+1}}>{{ dsdir | basename }}</h{{h+1}}>
        {% if envs.fgsea %}
            {% if dsdir | joinpaths: "fgsea.txt" | as_path | attr: "is_file" | call %}
                {{ fgsea_report(dsdir, h+2, envs, 10) }}
            {% else %}
                <p>Not enough events.</p>
            {% endif %}
        {% else %}
            {{ gsea_report(dsdir, h+2, envs, 10) }}
        {% endif %}
    {% endfor %}
{% endfor %}

{%- endmacro -%}

{%- macro head_job(job) -%}
{%- set name = job.in.configfile | config: "toml" | attr: "name" -%}
{%- if name or proc.size > 1 -%}
<h1>{{name | default: "Job #" + str(job.index+1) | escape}}</h1>
{%- endif -%}
{%- endmacro -%}

{{ report_jobs(jobs, head_job, report_job) }}
