{% from "utils/misc.liq" import report_jobs, table_of_images -%}
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

{%- if job.out.outdir | glob: "exprs/ridgeplots-*.png" -%}
<h{{h}}>Ridge plots of gene expressions</h{{h}}>
    {%- set figures = job.out.outdir | glob: "exprs/ridgeplots-*.png" -%}
    {%- for figure in figures -%}
        {%- set title = figure | append: ".title" | read | escape -%}
        {%- if not title.startswith("ridgeplots-") or len(figures) > 1 -%}
            <h{{h+1}}>{{ title | escape }}</h{{h+1}}>
        {%- endif -%}
        <p><Image src={{figure | quote}} /></p>
    {%- endfor -%}
{%- endif -%}

{%- if job.out.outdir | glob: "exprs/vlnplots-*.png" -%}
<h{{h}}>Violin plots of gene expressions</h{{h}}>
    {%- set figures = job.out.outdir | glob: "exprs/vlnplots-*.png" -%}
    {%- for figure in figures -%}
        {%- set title = figure | append: ".title" | read | escape -%}
        {%- if not title.startswith("vlnplots-") or len(figures) > 1 -%}
            <h{{h+1}}>{{ title | escape }}</h{{h+1}}>
        {%- endif -%}
        <p><Image src={{figure | quote}} /></p>
    {%- endfor -%}
{%- endif -%}

{%- if job.out.outdir | glob: "exprs/featureplots-*.png" -%}
<h{{h}}>Feature plots of gene expressions</h{{h}}>
    {%- set figures = job.out.outdir | glob: "exprs/featureplots-*.png" -%}
    {%- for figure in figures -%}
        {%- set title = figure | append: ".title" | read | escape -%}
        {%- if not title.startswith("featureplots-") or len(figures) > 1 -%}
            <h{{h+1}}>{{ title | escape }}</h{{h+1}}>
        {%- endif -%}
        <p><Image src={{figure | quote}} /></p>
    {%- endfor -%}
{%- endif -%}

{%- if job.out.outdir | glob: "exprs/dotplot-*.png" -%}
<h{{h}}>Dot plot of gene expressions</h{{h}}>
    {%- set figures = job.out.outdir | glob: "exprs/dotplot-*.png" -%}
    {%- for figure in figures -%}
        {%- set title = figure | append: ".title" | read | escape -%}
        {%- if not title.startswith("dotplot-") or len(figures) > 1 -%}
            <h{{h+1}}>{{ title | escape }}</h{{h+1}}>
        {%- endif -%}
        <p><Image src={{figure | quote}} /></p>
    {%- endfor -%}
{%- endif -%}

{%- if job.out.outdir | glob: "exprs/heatmap-*.png" -%}
    <h{{h}}>Heatmap of gene expressions</h{{h}}>
    {%- set figures = job.out.outdir | glob: "exprs/heatmap-*.png" -%}
    {%- for figure in figures -%}
        {%- set title = figure | append: ".title" | read | escape -%}
        {%- if not title.startswith("heatmap-") or len(figures) > 1 -%}
            <h{{h+1}}>{{ title | escape }}</h{{h+1}}>
        {%- endif -%}
        <p><Image src={{figure | quote}} /></p>
    {%- endfor -%}
{%- endif -%}


{%- if job.out.outdir | joinpaths: "dimplots" | as_path | attr: "is_dir" | call -%}
<h{{h}}>Dimensional reduction plots</h{{h}}>
{%-   set dpimgs = job.out.outdir | glob: "dimplots", "*.png" -%}
{{ table_of_images(dpimgs, caps=[]) }}
{%- endif -%}

{%- endmacro -%}

{%- macro head_job(job) -%}
<h1>{{job.in.srtobj | stem}}</h1>
{%- endmacro -%}

{{ report_jobs(jobs, head_job, report_job) }}