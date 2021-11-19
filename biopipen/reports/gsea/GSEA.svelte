{% from "utils.liq" import report_jobs, table_of_images -%}
<script>
    import { Image, DataTable } from "@@";
    import { Tile } from "carbon-components-svelte";
</script>

{%- macro report_job(job, h=1) -%}
<h{{h}}>Global view</h{{h}}>

<embed src={{job.out.outdir | joinpaths: "*.global.plots.pdf" | glob | first | quote}}
    width="100%"
    height="1000"
    type="application/pdf" />

<h{{h}}>Summary</h{{h}}>
{% for sumfile in job.out.outdir | joinpaths: "*.SUMMARY.RESULTS.REPORT.*.txt" | glob %}
{%   set klass = stem(sumfile).split(".")[-1] %}
<h{{h+1}}>{{klass}}</h{{h+1}}>
<DataTable data={ {{sumfile | datatable: sep="\t"}} } />
{% endfor %}

<h{{h}}>Enrichment details</h{{h}}>
{% for sumfile in job.out.outdir | joinpaths: "*.SUMMARY.RESULTS.REPORT.*.txt" | glob %}
{%   set klass = stem(sumfile).split(".")[-1] %}
<h{{h+1}}>{{klass}}</h{{h+1}}>
{%   set sumdata = sumfile | datatable: sep="\t" | json_loads %}
{%   set has_signif = [] %}
{%   for row in sumdata %}
{%      if row["FDR_q_val"] < 0.25 %}
{%          set _ = has_signif.append(1) %}
<embed src={{job.out.outdir | joinpaths: "*." + row["GS"] + ".plot." + klass + ".*.pdf" | glob | first | quote}}
    width="100%"
    height="700"
    type="application/pdf" />
{%      endif %}
{%   endfor %}
{%   if len(has_signif) == 0 %}
<Tile>No significantly (FDR_q_val &lt; 0.25) enriched pathways found.</Tile>
{%   endif %}
{% endfor %}

{%- endmacro -%}

{%- macro head_job(job) -%}
<h1>{{job.in.infile | stem}}</h1>
{%- endmacro -%}

{{ report_jobs(jobs, head_job, report_job) }}
