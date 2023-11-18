{% from "utils/misc.liq" import report_jobs -%}
{% from "utils/gsea.liq" import enrichr_report -%}
<script>
    import { Image, DataTable } from "$libs";
    import { Tabs, Tab, TabContent, InlineNotification } from "$ccs";
</script>


{%- macro report_job(job, h=1) -%}
    {%- set alldirs = job.out.outdir | glob: "*" -%}
    {%- set ovdir = job.out.outdir | joinpaths: "OVERLAPS" -%}
    {%- set secdirs = [] -%}
    {%- for adir in alldirs -%}
        {%- if basename(adir) != "OVERLAPS" -%}
            {%- set _ = secdirs.append(adir) -%}
        {%- endif -%}
    {%- endfor -%}

    {%- if len(secdirs) == 1 -%}
        {%- set secname = secdirs | first | basename -%}
        {%- for casedir in secdirs[0] | glob: "*" -%}
            {%- if secname == "DEFAULT" -%}
                <h{{h}}>{{casedir | basename | escape}}</h{{h}}>
            {%- else -%}
                <h{{h}}>{{secname | escape}} - {{casedir | basename | escape}}</h{{h}}>
            {%- endif -%}
            {%- if casedir | joinpaths: "error.txt" | exists -%}
                <InlineNotification
                    hideCloseButton
                    lowContrast
                    kind="warning"
                    subtitle={{ casedir | joinpaths: "error.txt" | read | quote }}
                    />
            {%- else -%}
                <h{{h+1}}>Markers</h{{h+1}}>
                <Tabs>
                    <Tab label="Markers" />
                    <Tab label="Volcano Plot" />
                    <Tab label="Dot Plot" />
                    <svelte:fragment slot="content">
                        <TabContent>
                            <DataTable
                                src={{ casedir | joinpaths: "markers.txt" | quote }}
                                data={ {{ casedir | joinpaths: "markers.txt" | datatable: sep="\t", nrows=100 }} }
                                />
                        </TabContent>
                        <TabContent>
                            <Image src={{ casedir | joinpaths: "volcano.png" | quote }} />
                        </TabContent>
                        <TabContent>
                            <Image src={{ casedir | joinpaths: "dotplot.png" | quote }} />
                        </TabContent>
                    </svelte:fragment>
                </Tabs>

                <h{{h+1}}>Enrichment analysis</h{{h+1}}>
                {{ enrichr_report(casedir) }}
            {%- endif -%}
        {%- endfor -%}
    {%- else -%}
        {%- for secdir in secdirs -%}
            {%- set sec = secdir | basename -%}
            <h{{h}}>{{sec | escape}}</h{{h}}>
            {%- for casedir in secdir | glob: "*" -%}
                <h{{h+1}}>{{casedir | basename | escape}}</h{{h+1}}>
                {%- if casedir | joinpaths: "error.txt" | exists -%}
                    <InlineNotification
                        hideCloseButton
                        lowContrast
                        kind="warning"
                        subtitle={{ casedir | joinpaths: "error.txt" | read | quote }}
                        />
                {%- else -%}
                    <h{{h+2}}>Markers</h{{h+2}}>
                    <Tabs>
                        <Tab label="Markers" />
                        <Tab label="Volcano Plot" />
                        <Tab label="Dot Plot" />
                        <svelte:fragment slot="content">
                            <TabContent>
                                <DataTable
                                    src={{ casedir | joinpaths: "markers.txt" | quote }}
                                    data={ {{ casedir | joinpaths: "markers.txt" | datatable: sep="\t", nrows=100 }} }
                                    />
                            </TabContent>
                            <TabContent>
                                <Image src={{ casedir | joinpaths: "volcano.png" | quote }} />
                            </TabContent>
                            <TabContent>
                                <Image src={{ casedir | joinpaths: "dotplot.png" | quote }} />
                            </TabContent>
                        </svelte:fragment>
                    </Tabs>

                    <h{{h+2}}>Enrichment analysis</h{{h+2}}>
                    {{ enrichr_report(casedir) }}
                {%- endif -%}
            {%- endfor -%}
        {%- endfor -%}
    {%- endif -%}

    {%- if ovdir | exists -%}
        <h{{h}}>Overlapping Markers</h{{h}}>
        {%- for casedir in ovdir | glob: "*" -%}
            <h{{h+1}}>{{casedir | basename | escape}}</h{{h+1}}>
            <Tabs>
                {%- if casedir | joinpaths: "venn.png" | exists -%}
                <Tab label="Venn Diagram" />
                {%- endif -%}
                <Tab label="UpSet Plot" />
                <Tab label="Marks" />
                <svelte:fragment slot="content">
                    {%- if casedir | joinpaths: "venn.png" | exists -%}
                    <TabContent>
                        <Image src={{ casedir | joinpaths: "venn.png" | quote }} />
                    </TabContent>
                    {%- endif -%}
                    <TabContent>
                        <Image src={{ casedir | joinpaths: "upset.png" | quote }} />
                    </TabContent>
                    <TabContent>
                        <DataTable
                            src={{ casedir | joinpaths: "markers.txt" | quote }}
                            data={ {{ casedir | joinpaths: "markers.txt" | datatable: sep="\t" }} }
                            />
                    </TabContent>
                </svelte:fragment>
            </Tabs>
        {%- endfor -%}
    {%- endif -%}
{%- endmacro -%}


{%- macro head_job(job) -%}
  <h1>{{job.in.srtobj | stem | escape}}</h1>
{%- endmacro -%}

{{ report_jobs(jobs, head_job, report_job) }}