{% from "utils/misc.liq" import table_of_images -%}

{%- macro fgsea_report_script() -%}
import { Image, DataTable } from "../components";
{%- endmacro -%}

{%- macro fgsea_report(fgsea_dir, h, envs, nrows=100) -%}
{%- addfilter splitgenes -%}
def splitgenes(data):
    for dat in data:
        dat["leadingEdge"] = dat["leadingEdge"].replace(",", " ")
    return json_dumps(data)
{%- endaddfilter -%}

<h{{h}}>Enrichment table</h{{h}}>
<Image src={{ fgsea_dir | joinpaths: "gsea_table.png" | quote }} />

{% set data = fgsea_dir | joinpaths: "fgsea.txt" | datatable: sep="\t", nrows=nrows | json_loads %}

<h{{h}}>Enrichment pathways</h{{h}}>
<DataTable src={{ fgsea_dir | joinpaths: "fgsea.txt" | quote }}
    data={ {{ data | splitgenes: }} }
    pageSize={10} />

<h{{h}}>Enrichment plot of pathways</h{{h}}>
{%- python -%}
import os
def fgsea_plots(pathways, fgsea_dir):
    out = []
    for pathway in pathways:
        pathway = pathway.replace("/", "-")
        pwfig = joinpaths(fgsea_dir, f"fgsea_{pathway}.png")
        if os.path.exists(pwfig):
            out.append(pwfig)
    return out
{%- endpython -%}
{{ table_of_images(
    fgsea_plots(liquid_map(data, "pathway"), fgsea_dir),
    liquid_map(data, "pathway"),
    table_width=75
) }}

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
