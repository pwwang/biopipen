{% from "utils/misc.liq" import report_jobs -%}
<script>
    import { Image, DataTable } from "@@";
    import { Tabs, Tab, TabContent } from "carbon-components-svelte";
</script>

{%- macro report_job(job, h=1) -%}
<h{{h}}>TCR Cluster size distribution</h{{h}}>

<Image src={{job.out.outdir | joinpaths: "ClusterSizeDistribution/cluster_size_distribution.png" | quote}} />

<h{{h}}>Shared TCR clusters among samples</h{{h}}>

<Tabs>
    <Tab label="Plot" />
    <Tab label="Table" />
    <svelte:fragment slot="content">
        <TabContent>
            <Image src={{job.out.outdir | joinpaths: "SharedClusters/shared_clusters.png" | quote}} />
        </TabContent>
        <TabContent>
            <DataTable src={{ job.out.outdir | joinpaths: "SharedClusters/shared_clusters.txt" | quote }}
                data={ {{ job.out.outdir | joinpaths: "SharedClusters/shared_clusters.txt" | datatable: sep="\t", index_col=0 }} } />
        </TabContent>
    </svelte:fragment>
</Tabs>

{%- if job.out.outdir | joinpaths: "SharedClustersByGrouping" | as_path | attr: "is_dir" | call -%}
<h{{h}}>Shared TCR clusters from groups</h{{h}}>
<Image src={{job.out.outdir | joinpaths: "SharedClustersByGrouping/shared_clusters_by_grouping.png" | quote}} />
{%- endif -%}


<h{{h}}>Sample diversity using TCR clusters</h{{h}}>

{%- for divfile in job.out.outdir | glob: "SampleDiversity/diversity_*.txt" -%}
{%-     set method = divfile | stem | replace: "diversity_", "" -%}
{%-     set plotname = divfile | basename | replace: ".txt", ".png" -%}
{%-     set plotfile = divfile | dirname | joinpaths: plotname -%}
<h{{h+1}}>Method: {{method}}</h{{h+1}}>

{%- if method == "chao1" -%}
<p>Chao1 estimator is a nonparameteric asymptotic estimator of species richness (number of species in a population).</p>
{%- elif method == "hill" -%}
<p>Hill numbers are a mathematically unified family of diversity indices (differing only by an exponent q).</p>
{%- elif method == "div" -%}
<p>True diversity, or the effective number of types, refers to the number of equally-abundant types needed for the average proportional abundance of the types to equal that observed in the dataset of interest where all types may not be equally abundant.</p>
{%- elif method == "gini.simp" -%}
<p>The Gini-Simpson index is the probability of interspecific encounter, i.e., probability that two entities represent different types.</p>
{%- elif method == "inv.simp" -%}
<p>Inverse Simpson index is the effective number of types that is obtained when the weighted arithmetic mean is used to quantify average proportional abundance of types in the dataset of interest.</p>
{%- elif method == "gini" -%}
<p>The Gini coefficient measures the inequality among values of a frequency distribution (for example levels of income). A Gini coefficient of zero expresses perfect equality, where all values are the same (for example, where everyone has the same income). A Gini coefficient of one (or 100 percents ) expresses maximal inequality among values (for example where only one person has all the income).</p>
{%- elif method == "raref " -%}
<p>Rarefaction is a technique to assess species richness from the results of sampling through extrapolation.</p>
{%- endif -%}

<Tabs>
    <Tab label="Plot" />
    <Tab label="Table" />
    <svelte:fragment slot="content">
        <TabContent>
            <Image src={{plotfile | quote}} />
        </TabContent>
        <TabContent>
            <DataTable src={{divfile | quote}}
                data={ {{ divfile | datatable: sep="\t", index_col=0 }} } />
        </TabContent>
    </svelte:fragment>
</Tabs>

{%- endfor -%}

{%- endmacro -%}

{%- macro head_job(job) -%}
<h1>{{job.in.immfile | stem | replace: ".immunarch", ""}}</h1>
{%- endmacro -%}

{{ report_jobs(jobs, head_job, report_job) }}

