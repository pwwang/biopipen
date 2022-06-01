{% macro report_jobs(jobs, job_head_macro, job_report_macro) %}

{% if len(jobs) == 1 %}

{{  job_report_macro(jobs[0], h=1) }}

{% else %}

{%  for job in jobs %}
{{      job_head_macro(job)}}
{{      job_report_macro(job, h=2) }}
{%  endfor %}

{% endif %}
{% endmacro %}


{%- macro table_of_images_script() -%}
import { Image } from "../components";
{%- endmacro -%}

{%- macro table_of_images(srcs, caps=None, col=2, table_width=100) -%}
<table class="image-table" style="background-color: var(--cds-background-hover); width: 100%">

{% for batch_srcs in srcs | batch: col, "" %}
{% set outer_loop = loop %}
<tr>
    {% for src in batch_srcs %}
        {% set i = col * outer_loop.index0 + loop.index0 %}
        {% if i >= len(srcs) %}
            <td style="width: {{table_width / col}}%"></td>
        {% else %}
        <td style="width: {{table_width / col}}%; vertical-align:top;">
            {% if caps is none %}
            <div
                style="padding-left: 28px; font-weight: bold; padding-top: 10px; margin-bottom: -10px;">
                {{ src | stem }}
            </div>
            {% elif caps %}
            <div
                style="padding-left: 28px; font-weight: bold; padding-top: 10px; margin-bottom: -10px;">
                {{ caps[i] }}
            </div>
            {% endif %}
            <Image
                style="max-width: 90%"
                src={{src | quote}}
            />
        </td>
        {% endif %}
    {% endfor %}
    {% if table_width < 100 %}
    <td style="width: {{100 - table_width}}%"></td>
    {% endif %}
</tr>
{% endfor %}

</table>
{%- endmacro -%}
