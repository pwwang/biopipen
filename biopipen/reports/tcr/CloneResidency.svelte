{% from "utils/misc.liq" import report_jobs, table_of_images -%}
<script>
    import { Image, DataTable } from "@@";
    import { Dropdown } from "carbon-components-svelte";

    let count_sample;

</script>


{%- macro report_job(job, h=1) -%}

    <h{{h}}>Count table</h{{h}}>
    <Dropdown
        selectedId="-1"
        items={[
            { id: "-1", text: "Select a sample" },
            {% for i, countfile in job.out.outdir | glob: "counts", "*.txt" | enumerate -%}
                { id: "{{i}}", text: "{{countfile | stem}}" },
            {% endfor %}
        ]}
        on:select={({ detail }) => {
            count_sample = detail.selectedItem.text;
        }}
    />
    {% for countfile in job.out.outdir | glob: "counts", "*.txt" -%}
        {#if count_sample == {{ countfile | stem | quote }} }
        <DataTable src={{
            job.out.outdir
            | joinpaths: "counts", stem(countfile) + ".txt"
            | quote
        }} data={ {{
            job.out.outdir
            | joinpaths: "counts", stem(countfile) + ".txt"
            | datatable: sep="\t", nrows=20, index_col=None
        }} } />
        {/if}
    {%- endfor %}

    <h{{h}}>Residency plots</h{{h}}>

    {% if job.out.outdir | joinpaths: "sample_groups" | as_path | attr: "exists" | call %}
        {% for groupfile in job.out.outdir | glob: "sample_groups", "*.txt" %}
            {% set group = groupfile | stem %}
            <h{{h+1}}>{{group}}</h{{h+1}}>
            {% set scatter_pngs = [] %}
            {% for sample in groupfile | readlines %}
                {% set spngs = job.out.outdir | glob: "scatter", "scatter_" + sample + "_*.png" %}
                {% set _ = scatter_pngs.extend(spngs) %}
            {% endfor %}
            {{ table_of_images(scatter_pngs, col=3, caps=false) }}
        {% endfor %}
    {% else %}
        {% assign scatter_pngs = job.out.outdir | joinpaths: "scatter", "scatter_*.png" | glob %}
        {{ table_of_images(scatter_pngs, col=3, caps=false) }}
    {% endif %}

    <h{{h}}>Clonotype overlapping</h{{h}}>

    {% if job.out.outdir | joinpaths: "sample_groups" | as_path | attr: "exists" | call %}
        {% for groupfile in job.out.outdir | glob: "sample_groups", "*.txt" %}
            {% set group = groupfile | stem %}
            <h{{h+1}}>{{group}}</h{{h+1}}>
            {% set venn_pngs = [] %}
            {% for sample in groupfile | readlines %}
                {% set spng = job.out.outdir | glob0: "venn", "venn_" + sample + ".png" %}
                {% set _ = venn_pngs.append(spng) %}
            {% endfor %}
            {{ table_of_images(venn_pngs, col=3) }}
        {% endfor %}
    {% else %}
        {% assign venn_pngs = job.out.outdir | joinpaths: "venn", "*.png" | glob %}
        {{ table_of_images(venn_pngs, col=3) }}
    {% endif %}
{%- endmacro -%}

{%- macro head_job(job) -%}
    <h1>{{job.out.outdir | stem | replace: ".immunarch", ""}}</h1>
{%- endmacro -%}

{{ report_jobs(jobs, head_job, report_job) }}
