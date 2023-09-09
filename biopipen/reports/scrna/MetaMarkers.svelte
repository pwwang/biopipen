{% from "utils/misc.liq" import report_jobs, table_of_images -%}
{% from "utils/gsea.liq" import enrichr_report -%}
<script>
    import { Image, DataTable } from "$libs";
    import { Tabs, Tab, TabContent, InlineNotification } from "$ccs";
</script>


{%- macro report_job(job, h=1) -%}
    {%- set secdirs = job.out.outdir | glob: "*" -%}
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
                    <Tab label="Top 10 genes" />
                    <Tab label="Table of top 100 genes" />
                    <div slot="content">
                    <TabContent>
                        {{ table_of_images(glob(casedir, "plots", "*.png"), col = 4) }}
                    </TabContent>
                    <TabContent>
                        <DataTable
                            src={{ casedir | joinpaths: "markers.txt" | quote }}
                            data={ {{ casedir | joinpaths: "markers.txt" | datatable: sep="\t", nrows=100 }} }
                            />
                    </TabContent>
                    </div>
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
                        <Tab label="Top 10 genes" />
                        <Tab label="Table of top 100 genes" />
                        <TabContent>
                            {{ table_of_images(glob(casedir, "plots", "*.png"), col = 4) }}
                        </TabContent>
                        <TabContent>
                            <DataTable
                                src={{ casedir | joinpaths: "markers.txt" | quote }}
                                data={ {{ casedir | joinpaths: "markers.txt" | datatable: sep="\t", nrows=100 }} }
                                />
                        </TabContent>
                    </Tabs>

                    <h{{h+2}}>Enrichment analysis</h{{h+2}}>
                    {{ enrichr_report(casedir) }}
                {%- endif -%}
            {%- endfor -%}
        {%- endfor -%}
    {%- endif -%}
{%- endmacro -%}


{%- macro head_job(job) -%}
  <h1>{{job.in.srtobj | stem | escape}}</h1>
{%- endmacro -%}

{{ report_jobs(jobs, head_job, report_job) }}