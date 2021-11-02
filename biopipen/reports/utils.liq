{% macro report_jobs(jobs, job_head_macro, job_report_macro) %}
{% if len(jobs) == 1 %}
{{  job_report_macro(jobs[0], h=1) }}
{% else %}
{%  for job in jobs %}
<div class="job_wrapper>
{{      job_head_macro(job)}}
{{      job_report_macro(job, h=2) }}
</div>
{%  endfor %}
{% endif %}
{% endmacro %}


{%- macro table_of_images(srcs, caps=None, col=2, max_width="90%") -%}
<table class="image-table">
{% for i, src in enumerate(srcs) %}
{%  if i % col == 0 %}
<tr>
{%  endif %}
<td style="width: {{100.0/col}}%; background-color: var(--cds-background-hover); vertical-align:top;">
    {% if caps %}
    <div style="text-align: center; padding-top: 10px; margin-bottom: -10px;">{{ caps[i] }}</div>
    {% endif %}
    <Image src={{src | quote}} style="max-width: {{max_width}}" />
</td>
{%  if (i+1) % col == 0 %}
</tr>
{%  endif %}
{% endfor  %}
{% if len(srcs) % col > 0 %}
{%  for i in range(col - len(srcs) % col) %}
<td style="width: {{100.0/col}}%; background-color: var(--cds-background-hover); vertical-align:top;">
</td>
{%  endfor %}
</tr>
{% endif %}
</table>
{%- endmacro -%}
