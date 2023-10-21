{% from "utils/misc.liq" import report_jobs -%}
<script>
    import { Image, DataTable } from "$libs";
    import { Tabs, Tab, TabContent } from "$ccs";
</script>

{%- macro report_job(job, h=1) -%}
<h{{h}}>TCR Cluster size distribution</h{{h}}>

{% for casedir in job.out.outdir | glob: "ClusterSizeDistribution", "*" %}
    {% set casename = casedir | basename %}
    {% if casename != "DEFAULT" %}
        <h{{h+1}}>{{casename}}</h{{h+1}}>
    {% endif %}
    <Image src={{casedir | joinpaths: "cluster_size_distribution.png" | quote}} />
{% endfor %}

<h{{h}}>Shared TCR clusters among samples</h{{h}}>

{% for casedir in job.out.outdir | glob: "SharedClusters", "*" %}
    {% set casename = casedir | basename %}
    {% if casename != "DEFAULT" %}
        <h{{h+1}}>{{casename}}</h{{h+1}}>
    {% endif %}
    <Image src={{casedir | joinpaths: "shared_clusters.png" | quote}} />
{% endfor %}


<h{{h}}>Sample diversity using TCR clusters</h{{h}}>

{% for casedir in job.out.outdir | glob: "SampleDiversity", "*" %}
    {% set casename = casedir | basename %}
    {% if casename != "DEFAULT" %}
        <h{{h+1}}>{{casename}}</h{{h+1}}>
    {% endif %}

    <Tabs>
        <Tab label="Plot" />
        <Tab label="Table" />
        <svelte:fragment slot="content">
            <TabContent>
                <Image src={{casedir | joinpaths: "diversity.png" | quote}} />
            </TabContent>
            <TabContent>
                <DataTable src={{casedir | joinpaths: "diversity.txt" | quote}}
                    data={ {{ casedir | joinpaths: "diversity.txt" | datatable: sep="\t", index_col=0 }} } />
            </TabContent>
        </svelte:fragment>
    </Tabs>
{% endfor %}
{%- endmacro -%}

{%- macro head_job(job) -%}
<h1>{{job.in.immfile | stem | replace: ".immunarch", ""}}</h1>
{%- endmacro -%}

{{ report_jobs(jobs, head_job, report_job) }}
