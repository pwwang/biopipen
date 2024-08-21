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
import { Image } from "$libs";
{%- endmacro -%}

{%- macro table_of_images(srcs, caps=None, col=2, table_width=100) -%}
<table class="image-table" style="background-color: var(--cds-background-hover); width: 100%">

{% for batch_srcs in srcs | batch: col, "" %}
{% set outer_loop = loop %}
<tr>
    {% for srcinfo in batch_srcs %}
        {% if srcinfo | isinstance: str %}
            {% set src = srcinfo %}
            {% set download = None %}
        {% else %}
            {% set src = srcinfo['src'] %}
            {% set download = srcinfo.get('download', None) %}
        {% endif %}
        {% set i = col * outer_loop.index0 + loop.index0 %}
        {% if i >= len(srcs) %}
            <td style="width: {{table_width / col}}%"></td>
        {% else %}
        <td style="width: {{table_width / col}}%; vertical-align:top;">
            {% if caps is none %}
            <div
                style="padding-left: 28px; font-weight: bold; padding-top: 16px;">
                {{ src | stem }}
            </div>
            {% elif caps is false %}
            {% else %}
            <div
                style="padding-left: 28px; font-weight: bold; padding-top: 16px;">
                {{ caps[i] }}
            </div>
            {% endif %}
            {% if download %}
                <Image
                    style="max-width: 90%"
                    src={{src | quote}}
                    download={ {{download | json}} }
                />
            {% else %}
                <Image
                    style="max-width: 90%"
                    src={{src | quote}}
                />
            {% endif %}
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
