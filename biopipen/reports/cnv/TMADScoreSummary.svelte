{% from "utils/misc.liq" import report_jobs -%}

<script>
    import { Image, DataTable } from "$libs";
    import { Tabs, Tab, TabContent } from "$ccs";
</script>

{%- macro report_job(job, h=1) -%}
<h{{h}}>tMAD score for each sample</h{{h}}>
<Tabs>
    <Tab label="Barplot" />
    <Tab label="Table" />
    <svelte:fragment slot="content">
        <TabContent>
            <Image src="{{job.out.outdir | joinpaths: 'tMAD.png'}}" />
        </TabContent>
        <TabContent>
            <DataTable src="{{job.out.outdir | joinpaths: 'tMAD.txt'}}"
                data={ {{ job.out.outdir | joinpaths: 'tMAD.txt' | datatable: sep="\t" }} } />
        </TabContent>
    </svelte:fragment>
</Tabs>

{% if envs.group_cols | isinstance: str %}
    <h{{h}}>By {{envs.group_cols | replace: ",", ", then "}}</h{{h}}>

    <Tabs>
        <Tab label="Barplot" />
        <Tab label="Voilin Plot" />
        <svelte:fragment slot="content">
            <TabContent>
                <Image src="{{job.out.outdir | joinpaths: 'tMAD_' + envs.group_cols + '_bar.png'}}" />
            </TabContent>
            <TabContent>
                <Image src="{{job.out.outdir | joinpaths: 'tMAD_' + envs.group_cols + '_box_violin.png'}}" />
            </TabContent>
        </svelte:fragment>
    </Tabs>
{% elif envs.group_cols %}
    {% for group_col in envs.group_cols %}
    <h{{h}}>By {{group_col | replace: ",", ", then "}}</h{{h}}>

    <Tabs>
        <Tab label="Barplot" />
        <Tab label="Voilin Plot" />
        <svelte:fragment slot="content">
            <TabContent>
                <Image src="{{job.out.outdir | joinpaths: 'tMAD_' + group_col + '_bar.png'}}" />
            </TabContent>
            <TabContent>
                <Image src="{{job.out.outdir | joinpaths: 'tMAD_' + group_col + '_box_violin.png'}}" />
            </TabContent>
        </svelte:fragment>
    </Tabs>
    {% endfor %}
{% endif %}
{%- endmacro -%}

{%- macro head_job(job) -%}
<h1>{{job.in.metafile | stem0 }}</h1>
{%- endmacro -%}

{{ report_jobs(jobs, head_job, report_job) }}
