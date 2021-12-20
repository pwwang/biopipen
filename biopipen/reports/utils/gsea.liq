{% from "utils/misc.liq" import table_of_images -%}



{%- macro fgsea_report(fgsea_dir, h, envs, nrows=100) -%}

<h{{h}}>Enrichment table</h{{h}}>
<Image src={{ fgsea_dir | joinpaths: "gsea_table.png" | quote }} />

<h{{h}}>Enrichment pathways</h{{h}}>
<DataTable src={{ fgsea_dir | joinpaths: "fgsea.txt" | quote }}
    data={ {{ fgsea_dir | joinpaths: "fgsea.txt" | datatable: sep="\t", nrows=nrows }} } />

<h{{h}}>Enrichment plot of pathways</h{{h}}>
{{ table_of_images(glob(joinpaths(fgsea_dir, "fgsea_*.png"))[:nrows]) }}

{%- endmacro -%}



{%- macro gsea_report(gsea_dir, h, envs, nrows=100) -%}
<h{{h}}>Global view</h{{h}}>

<embed src={{gsea_dir | joinpaths: "*.global.plots.pdf" | glob | first | quote}}
    width="100%"
    height="1000"
    type="application/pdf" />

<h{{h}}>Summary</h{{h}}>
{% for sumfile in gsea_dir | joinpaths: "*.SUMMARY.RESULTS.REPORT.*.txt" | glob %}
{%   set klass = stem(sumfile).split(".")[-1] %}
<h{{h+1}}>{{klass}}</h{{h+1}}>
<DataTable data={ {{sumfile | datatable: sep="\t", nrows=nrows}} } />
{% endfor %}

<h{{h}}>Enrichment details</h{{h}}>
{% set cutoff = envs.get("fdr.q.val.threshold", envs.get("fdr_q_val_threshold", 0.25)) %}
{% for sumfile in gsea_dir | joinpaths: "*.SUMMARY.RESULTS.REPORT.*.txt" | glob %}
{%   set klass = stem(sumfile).split(".")[-1] %}
<h{{h+1}}>{{klass}}</h{{h+1}}>
{%   set sumdata = sumfile | datatable: sep="\t" | json_loads %}
{%   set has_signif = [] %}
{%   for row in sumdata %}
{%      if row["FDR_q_val"] < cutoff %}
{%          set _ = has_signif.append(1) %}
<embed src={{gsea_dir | joinpaths: "*." + row["GS"] + ".plot." + klass + ".*.pdf" | glob | first | quote}}
    width="100%"
    height="700"
    type="application/pdf" />
{%      endif %}
{%   endfor %}
{%   if len(has_signif) == 0 %}
<Tile>No significantly (FDR_q_val &lt; {{cutoff}}) enriched pathways found.</Tile>
{%   endif %}
{% endfor %}

{%- endmacro -%}
