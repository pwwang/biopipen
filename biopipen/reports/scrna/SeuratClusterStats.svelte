{% from "utils/misc.liq" import report_jobs, table_of_images -%}
{% from_ os import path %}
<script>
    import { DataTable, Image } from "$libs";
    import { Tabs, Tab, TabContent } from "$ccs";
</script>

{%- macro report_job(job, h=1) -%}
    {%- set stats_reports_file = job.out.outdir | joinpaths: "stats", "report_toc.json" -%}
    {%- set features_reports_file = job.out.outdir | joinpaths: "features", "report_toc.json" -%}
    {%- set dimplots_reports_file = job.out.outdir | joinpaths: "dimplots", "report_toc.json" -%}

    {%- if stats_reports_file | exists -%}
        {%- set stats = stats_reports_file | config: "json" -%}
        {% for key, value in stats.items() -%}
            <h{{h}}>{{key | escape}}</h{{h}}>
            <Tabs>
                {% if 'bar' in value -%}
                    <Tab label="Bar plot" />
                {% endif -%}
                {% if 'pie' in value -%}
                    <Tab label="Pie chart" />
                {% endif -%}
                {% if 'table' in value -%}
                    <Tab label="Table" />
                {% endif -%}
                <svelte:fragment slot="content">
                    {% if 'bar' in value -%}
                        <TabContent>
                            <Image src="{{job.out.outdir}}/stats/{{value.bar}}" />
                        </TabContent>
                    {% endif -%}
                    {% if 'pie' in value -%}
                        <TabContent>
                            <Image src="{{job.out.outdir}}/stats/{{value.pie}}" />
                        </TabContent>
                    {% endif -%}
                    {% if 'table' in value -%}
                        <TabContent>
                            <DataTable src="{{job.out.outdir}}/stats/{{value.table}}"
                                data={ {{job.out.outdir | joinpaths: "stats", value.table | datatable: sep="\t", nrows=100 }} }
                                />
                        </TabContent>
                    {% endif -%}
                </svelte:fragment>
            </Tabs>
        {%- endfor -%}
    {%- endif -%}

    {%- if features_reports_file | exists -%}
        {%- set features = features_reports_file | config: "json" %}
        {% for key, value in features.items() -%}
            <h{{h}}>{{key | escape}}</h{{h}}>
            {% for val in value -%}
                {% if "name" in val -%}
                    <h{{h+1}}>{{val.name | escape}}</h{{h+1}}>
                {%- endif -%}
                {% if val.kind == "table" -%}
                    <DataTable src="{{job.out.outdir}}/features/{{val.file}}"
                        data={ {{job.out.outdir | joinpaths: "features", val.file | datatable: sep="\t", nrows=100 }} }
                        />
                {% else -%}
                    <Image src="{{job.out.outdir}}/features/{{val.file}}" />
                {% endif -%}
            {%- endfor -%}
        {%- endfor -%}
    {%- endif -%}

    {%- if dimplots_reports_file | exists -%}
        {%- set dimplots = dimplots_reports_file | config: "json" %}
        {% for key, value in dimplots.items() -%}
            <h{{h}}>{{key | escape}}</h{{h}}>
            <Image src="{{job.out.outdir}}/dimplots/{{value}}" />
        {%- endfor -%}
    {%- endif -%}
{%- endmacro -%}

{%- macro head_job(job) -%}
<h1>{{job.in.srtobj | stem}}</h1>
{%- endmacro -%}

{{ report_jobs(jobs, head_job, report_job) }}