{% from "utils/misc.liq" import report_jobs -%}
<script>
    import { Image, DataTable } from "@@";
    import { Tabs, Tab, TabContent } from "carbon-components-svelte";
</script>


{%- macro report_job(job, h=1) -%}
{% for casedir in job.out.outdir | glob: "*" %}
{%  set case = casedir | basename %}
<h{{h}}>{{case}}</h{{h}}>

<h{{h+1}}>Markers</h{{h+1}}>
<DataTable
    src={{ casedir | joinpaths: "markers.txt" | quote }}
    data={ {{ casedir | joinpaths: "markers.txt" | datatable: sep="\t", nrows=100 }} }
    />

<h{{h+1}}>Enrichment analysis</h{{h+1}}>
<Tabs>
    {% for enrtxt in casedir | as_path | attr: "glob" | call: "Enrichr-*.txt"  %}
    {%  set db = enrtxt | stem | replace: "Enrichr-", "" %}
    <Tab label="{{db}}" title="{{db}}" />
    {% endfor %}
    <div slot="content">
        {% for enrtxt in casedir | as_path | attr: "glob" | call: "Enrichr-*.txt" %}
        {%  set db = enrtxt | stem | replace: "Enrichr-", "" %}
        <TabContent>
            <Image src={{casedir | joinpaths: "Enrichr-" + db + ".png" | quote}} />
            <DataTable
                src={{ enrtxt | quote }}
                data={ {{ enrtxt | datatable: sep="\t", nrows=100 }} }
                />
        </TabContent>
        {% endfor %}
    </div>
</Tabs>

{% endfor %}
{%- endmacro -%}


{%- macro head_job(job) -%}
{% if in.casefile %}
{%  set name = in.casefile | config: "toml" | attr: "name" %}
{% else %}
{%  set name = envs | attr: "name" %}
{% endif %}
<h1>{{name | escape}}</h1>
{%- endmacro -%}

{{ report_jobs(jobs, head_job, report_job) }}