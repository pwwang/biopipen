{% from "utils/misc.liq" import report_jobs, table_of_images -%}
<script>
    import { Image, DataTable } from "$libs";
    import { Dropdown } from "$ccs";

    let count_sample;

</script>

{%- macro report_case(casename, h=2) -%}
    <h{{h}}>Count table</h{{h}}>
    <Dropdown
        selectedId="-1"
        items={[
            { id: "-1", text: "Select a sample" },
            {% for i, countfile in job.out.outdir | glob: casename, "counts", "*.txt" | enumerate -%}
                { id: "{{i}}", text: "{{countfile | stem}}" },
            {% endfor %}
        ]}
        on:select={({ detail }) => {
            count_sample = `{{casename}}-${detail.selectedItem.text}`;
        }}
    />
    {% for countfile in job.out.outdir | glob: casename, "counts", "*.txt" -%}
        {#if count_sample == {{ countfile | stem | prepend: casename + '-' | quote }} }
        <DataTable src={{
            job.out.outdir
            | joinpaths: casename, "counts", stem(countfile) + ".txt"
            | quote
        }} data={ {{
            job.out.outdir
            | joinpaths: casename, "counts", stem(countfile) + ".txt"
            | datatable: sep="\t", nrows=20, index_col=None
        }} } />
        {/if}
    {%- endfor %}

    <h{{h}}>Residency plots</h{{h}}>

    {% if job.out.outdir | joinpaths: casename, "sample_groups" | exists %}
        {% for groupfile in job.out.outdir | glob: casename, "sample_groups", "*.txt" %}
            {% set group = groupfile | stem %}
            <h{{h+1}}>{{group}}</h{{h+1}}>
            {% set scatter_pngs = [] %}
            {% for sample in groupfile | readlines %}
                {% set spngs = job.out.outdir | glob: casename, "scatter", "scatter_" + sample + "_*.png" %}
                {% set _ = scatter_pngs.extend(spngs) %}
            {% endfor %}
            {{ table_of_images(scatter_pngs, col=3, caps=false) }}
        {% endfor %}
    {% else %}
        {% assign scatter_pngs = job.out.outdir | joinpaths: casename, "scatter", "scatter_*.png" | glob %}
        {{ table_of_images(scatter_pngs, col=3, caps=false) }}
    {% endif %}

    <h{{h}}>Clonotype overlapping</h{{h}}>

    {% if job.out.outdir | joinpaths: casename, "sample_groups" | exists %}
        {% for groupfile in job.out.outdir | glob: casename, "sample_groups", "*.txt" %}
            {% set group = groupfile | stem %}
            <h{{h+1}}>{{group}}</h{{h+1}}>
            {% set venn_pngs = [] %}
            {% for sample in groupfile | readlines %}
                {% set spng = job.out.outdir | glob0: casename, "venn", "venn_" + sample + ".png" %}
                {% set _ = venn_pngs.append(spng) %}
            {% endfor %}
            {{ table_of_images(venn_pngs, col=3) }}
        {% endfor %}
    {% else %}
        {% assign venn_pngs = job.out.outdir | joinpaths: casename, "venn", "*.png" | glob %}
        {{ table_of_images(venn_pngs, col=3) }}
    {% endif %}
{%- endmacro -%}

{%- macro report_job(job, h=1) -%}
    {% set casedirs = job.out.outdir | glob: "*" %}
    {% if len(casedirs) == 1 %}
        {{ report_case(stem(casedirs[0]), h) }}
    {% else %}
        {% for casedir in casedirs %}
            <h{{h}}>{{casedir | stem}}</h{{h}}>
            {{ report_case(stem(casedir), h+1) }}
        {% endfor %}
    {% endif %}
{%- endmacro -%}

{%- macro head_job(job) -%}
    <h1>{{job.out.outdir | stem | replace: ".immunarch", ""}}</h1>
{%- endmacro -%}

{{ report_jobs(jobs, head_job, report_job) }}
