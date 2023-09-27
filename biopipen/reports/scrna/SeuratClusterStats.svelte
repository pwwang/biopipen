{% from "utils/misc.liq" import report_jobs, table_of_images -%}
{% from_ os import path %}
<script>
    import { DataTable, Image } from "$libs";
    import { Tabs, Tab, TabContent } from "$ccs";
</script>

{%- macro report_job(job, h=1) -%}
{%- for statname in proc.envs.stats -%}
    {%- if statname.startswith("nCells") -%}
        {%- set num_or_frac = "Number" -%}
        {%- set rest_title = statname | replace: "nCells_", "" | replace: "_", " " -%}
    {%- else -%}
        {%- set num_or_frac = "Fraction" -%}
        {%- set rest_title = statname | replace: "fracCells_", "" | replace: "_", " " -%}
    {%- endif -%}
    {%- if rest_title == "All" or rest_title == "ALL" -%}
        {%- set rest_title = "" -%}
    {%- else -%}
        {%- set rest_title = "(" + rest_title + ")" -%}
    {%- endif -%}

    {%- set plotfile = job.out.outdir | joinpaths: "stats", statname + ".png" -%}
    {%- set tablefile = plotfile + ".txt" -%}
    <h{{h}}>{{num_or_frac}} of cells {{rest_title}}</h{{h}}>
    <Tabs>
        <Tab label="Plot" />
        <Tab label="Table" />
        <svelte:fragment slot="content">
            <TabContent>
                <Image src={{plotfile | quote}} />
            </TabContent>
            <TabContent>
                <DataTable src={{tablefile | quote}}
                    data={ {{tablefile | datatable: sep="\t", nrows=100 }} } />
            </TabContent>
        </svelte:fragment>
    </Tabs>
{%- endfor -%}

{%- if job.out.outdir | glob: "features/table-*.tsv" -%}
<h{{h}}>Feature matrix</h{{h}}>
    {%- set tabfiles = job.out.outdir | glob: "features/table-*.tsv" -%}
    {%- for tabfile in tabfiles -%}
        {%- set title = tabfile | append: ".title" | read | escape -%}
        {%- if not title.startswith("table-") or len(tabfiles) > 1 -%}
            <h{{h+1}}>{{ title | escape }}</h{{h+1}}>
        {%- endif -%}
        <p><DataTable src={{tabfile | quote}}
            data={ {{tabfile | datatable: sep="\t", nrows=100 }} } /></p>
    {%- endfor -%}
{%- endif -%}

{%- if job.out.outdir | glob: "features/ridgeplots-*.png" -%}
<h{{h}}>Ridge plots of features</h{{h}}>
    {%- set figures = job.out.outdir | glob: "features/ridgeplots-*.png" -%}
    {%- for figure in figures -%}
        {%- set title = figure | append: ".title" | read | escape -%}
        {%- if not title.startswith("ridgeplots-") or len(figures) > 1 -%}
            <h{{h+1}}>{{ title | escape }}</h{{h+1}}>
        {%- endif -%}
        <p><Image src={{figure | quote}} /></p>
    {%- endfor -%}
{%- endif -%}

{%- if job.out.outdir | glob: "features/vlnplots-*.png" -%}
<h{{h}}>Violin plots of features</h{{h}}>
    {%- set figures = job.out.outdir | glob: "features/vlnplots-*.png" -%}
    {%- for figure in figures -%}
        {%- set title = figure | append: ".title" | read | escape -%}
        {%- if not title.startswith("vlnplots-") or len(figures) > 1 -%}
            <h{{h+1}}>{{ title | escape }}</h{{h+1}}>
        {%- endif -%}
        <p><Image src={{figure | quote}} /></p>
    {%- endfor -%}
{%- endif -%}

{%- if job.out.outdir | glob: "features/featureplots-*.png" -%}
<h{{h}}>Feature plots of features</h{{h}}>
    {%- set figures = job.out.outdir | glob: "features/featureplots-*.png" -%}
    {%- for figure in figures -%}
        {%- set title = figure | append: ".title" | read | escape -%}
        {%- if not title.startswith("featureplots-") or len(figures) > 1 -%}
            <h{{h+1}}>{{ title | escape }}</h{{h+1}}>
        {%- endif -%}
        <p><Image src={{figure | quote}} /></p>
    {%- endfor -%}
{%- endif -%}

{%- if job.out.outdir | glob: "features/dotplot-*.png" -%}
<h{{h}}>Dot plot of features</h{{h}}>
    {%- set figures = job.out.outdir | glob: "features/dotplot-*.png" -%}
    {%- for figure in figures -%}
        {%- set title = figure | append: ".title" | read | escape -%}
        {%- if not title.startswith("dotplot-") or len(figures) > 1 -%}
            <h{{h+1}}>{{ title | escape }}</h{{h+1}}>
        {%- endif -%}
        <p><Image src={{figure | quote}} /></p>
    {%- endfor -%}
{%- endif -%}

{%- if job.out.outdir | glob: "features/heatmap-*.png" -%}
    <h{{h}}>Heatmap of features</h{{h}}>
    {%- set figures = job.out.outdir | glob: "features/heatmap-*.png" -%}
    {%- for figure in figures -%}
        {%- set title = figure | append: ".title" | read | escape -%}
        {%- if not title.startswith("heatmap-") or len(figures) > 1 -%}
            <h{{h+1}}>{{ title | escape }}</h{{h+1}}>
        {%- endif -%}
        <p><Image src={{figure | quote}} /></p>
    {%- endfor -%}
{%- endif -%}


{%- if job.out.outdir | joinpaths: "dimplots" | isdir -%}
<h{{h}}>Dimensional reduction plots</h{{h}}>
{%-   set dpimgs = job.out.outdir | glob: "dimplots", "*.png" -%}
{{ table_of_images(dpimgs, caps=[]) }}
{%- endif -%}

{%- endmacro -%}

{%- macro head_job(job) -%}
<h1>{{job.in.srtobj | stem}}</h1>
{%- endmacro -%}

{{ report_jobs(jobs, head_job, report_job) }}