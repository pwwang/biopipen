{% from "utils/misc.liq" import report_jobs -%}
{% from_ os import path %}
<script>
    import { DataTable, Image } from "@@";
    import { Tabs, Tab, TabContent } from "carbon-components-svelte";
</script>

{%- macro report_job(job, h=1) -%}
<h{{h}}>Number of cells for clusters</h{{h}}>
<Tabs>
    <Tab label="Bar Plot" />
    <Tab label="Table" />
    <svelte:fragment slot="content">
        <TabContent>
            <Image src={{job.out.outdir | joinpaths: "stats/ncells.png" | quote}} />
        </TabContent>
        <TabContent>
            <DataTable src={{job.out.outdir | joinpaths: "stats/ncells.txt" | quote}}
                data={ {{job.out.outdir | joinpaths: "stats/ncells.txt" | datatable: sep="\t", nrows=100 }} } />
        </TabContent>
    </svelte:fragment>
</Tabs>

<h{{h}}>Number of cells from each sample for clusters</h{{h}}>
<Image src={{job.out.outdir | joinpaths: "stats/ncellspersample.png" | quote}} />

<h{{h}}>Fraction of cells from each sample for clusters</h{{h}}>
<Image src={{job.out.outdir | joinpaths: "stats/perccellspersample.png" | quote}} />

{% set ridgeplots = job.out.outdir | joinpaths: "exprs/ridgeplots.png" %}
{%- if ridgeplots | as_path | attr: "is_file" | call -%}
<h{{h}}>Ridge plots of gene expressions</h{{h}}>
<Image src={{job.out.outdir | joinpaths: "exprs/ridgeplots.png" | quote}} />
{%- endif -%}

{% set vlnplots = job.out.outdir | joinpaths: "exprs/vlnplots.png" %}
{%- if vlnplots | as_path | attr: "is_file" | call -%}
<h{{h}}>Violin plots of gene expressions</h{{h}}>
<Image src={{job.out.outdir | joinpaths: "exprs/vlnplots.png" | quote}} />
{%- endif -%}

{% set featureplots = job.out.outdir | joinpaths: "exprs/featureplots.png" %}
{%- if featureplots | as_path | attr: "is_file" | call -%}
<h{{h}}>Feature plots of gene expressions</h{{h}}>
<Image src={{job.out.outdir | joinpaths: "exprs/featureplots.png" | quote}} />
{%- endif -%}

{% set dotplot = job.out.outdir | joinpaths: "exprs/dotplot.png" %}
{%- if dotplot | as_path | attr: "is_file" | call -%}
<h{{h}}>Dot plot of gene expressions</h{{h}}>
<p>Intuitive way of visualizing how feature expression changes across different identity classes (clusters). The size of the dot encodes the percentage of cells within a class, while the color encodes the AverageExpression level across all cells within a class (blue is high).</p>
<Image src={{job.out.outdir | joinpaths: "exprs/dotplot.png" | quote}} />
{%- endif -%}

{% set heatmap = job.out.outdir | joinpaths: "exprs/heatmap.png" %}
{%- if heatmap | as_path | attr: "is_file" | call -%}
<h{{h}}>Heatmap of gene expressions</h{{h}}>
<Image src={{job.out.outdir | joinpaths: "exprs/heatmap.png" | quote}} />
{%- endif -%}

{%- endmacro -%}

{%- macro head_job(job) -%}
<h1>{{job.in.srtobj | stem}}</h1>
{%- endmacro -%}

{{ report_jobs(jobs, head_job, report_job) }}