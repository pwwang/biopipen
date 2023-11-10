{% from "utils/misc.liq" import report_jobs -%}

<script>
    import { Image, DataTable } from "$libs";
    import { Tabs, Tab, TabContent } from "$ccs";
</script>

{%- macro report_job(job, h=1) -%}
<h{{h}}>Aneuploidy Score for each sample</h{{h}}>
<Tabs>
    <Tab label="Barplot (arm)" />
    <Tab label="Table (arm)" />
    <Tab label="Barplot (seg)" />
    <Tab label="Table (seg)" />
    <svelte:fragment slot="content">
        <TabContent>
            <Image src="{{job.out.outdir | joinpaths: 'AS_arm.png'}}" />
        </TabContent>
        <TabContent>
            <DataTable src="{{job.out.outdir | joinpaths: 'AS_arm.txt'}}"
                data={ {{ job.out.outdir | joinpaths: 'AS_arm.txt' | datatable: sep="\t" }} } />
        </TabContent>
        <TabContent>
            <Image src="{{job.out.outdir | joinpaths: 'AS_seg.png'}}" />
        </TabContent>
        <TabContent>
            <DataTable src="{{job.out.outdir | joinpaths: 'AS_seg.txt'}}"
                data={ {{ job.out.outdir | joinpaths: 'AS_seg.txt' | datatable: sep="\t" }} } />
        </TabContent>
    </svelte:fragment>
</Tabs>

{% if envs.group_cols | isinstance: str %}
    <h{{h}}>Aneuploidy Score By {{envs.group_cols | replace: ",", ", then "}}</h{{h}}>

    <Tabs>
        <Tab label="Barplot (arm)" />
        <Tab label="Voilin Plot (arm)" />
        <Tab label="Barplot (seg)" />
        <Tab label="Voilin Plot (seg)" />
        <svelte:fragment slot="content">
            <TabContent>
                <Image src="{{job.out.outdir | joinpaths: 'AS_arm_bar_' + envs.group_cols + '.png'}}" />
            </TabContent>
            <TabContent>
                <Image src="{{job.out.outdir | joinpaths: 'AS_arm_violin_' + envs.group_cols + '.png'}}" />
            </TabContent>
            <TabContent>
                <Image src="{{job.out.outdir | joinpaths: 'AS_seg_bar_' + envs.group_cols + '.png'}}" />
            </TabContent>
            <TabContent>
                <Image src="{{job.out.outdir | joinpaths: 'AS_seg_violin_' + envs.group_cols + '.png'}}" />
            </TabContent>
        </svelte:fragment>
    </Tabs>
{% elif envs.group_cols %}
    {% for group_col in envs.group_cols %}
    <h{{h}}>Aneuploidy Score By {{group_col | replace: ",", ", then "}}</h{{h}}>

    <Tabs>
        <Tab label="Barplot (arm)" />
        <Tab label="Voilin Plot (arm)" />
        <Tab label="Barplot (seg)" />
        <Tab label="Voilin Plot (seg)" />
        <svelte:fragment slot="content">
            <TabContent>
                <Image src="{{job.out.outdir | joinpaths: 'AS_arm_bar_' + group_col + '.png'}}" />
            </TabContent>
            <TabContent>
                <Image src="{{job.out.outdir | joinpaths: 'AS_arm_violin_' + group_col + '.png'}}" />
            </TabContent>
            <TabContent>
                <Image src="{{job.out.outdir | joinpaths: 'AS_seg_bar_' + group_col + '.png'}}" />
            </TabContent>
            <TabContent>
                <Image src="{{job.out.outdir | joinpaths: 'AS_seg_violin_' + group_col + '.png'}}" />
            </TabContent>
        </svelte:fragment>
    </Tabs>
    {% endfor %}
{% endif %}

<h{{h}}>Arm Aneuploidy for each chromosome</h{{h}}>
<Tabs>
    <Tab label="arm" />
    <Tab label="seg" />
    <svelte:fragment slot="content">
        <TabContent>
            <DataTable src="{{job.out.outdir | joinpaths: 'CAA_arm.txt'}}"
                data={ {{ job.out.outdir | joinpaths: 'CAA_arm.txt' | datatable: sep="\t" }} } />
        </TabContent>
        <TabContent>
            <DataTable src="{{job.out.outdir | joinpaths: 'CAA_seg.txt'}}"
                data={ {{ job.out.outdir | joinpaths: 'CAA_seg.txt' | datatable: sep="\t" }} } />
        </TabContent>
    </svelte:fragment>
</Tabs>

{% set hms = job.out.outdir | glob: "Heatmap_*.png" %}
{% if hms %}
    <h{{h+1}}>Heatmaps</h{{h+1}}>
    <Tabs>
        {% for hmfile in hms %}
        <Tab label={{hmfile | stem | replace: "Heatmap_", "" | quote}} />
        {% endfor %}
        <svelte:fragment slot="content">
            {% for hmfile in hms %}
            <TabContent>
                <Image src="{{hmfile}}" />
            </TabContent>
            {% endfor %}
        </svelte:fragment>
    </Tabs>
{% endif %}

{%- endmacro -%}

{%- macro head_job(job) -%}
<h1>{{job.in.metafile | stem0 }}</h1>
{%- endmacro -%}

{{ report_jobs(jobs, head_job, report_job) }}

