{% from "utils/misc.liq" import report_jobs -%}

<script>
    import { Image, DataTable } from "$libs";
    import { Tabs, Tab, TabContent, Tile } from "$ccs";
</script>

{%- macro report_job(job, h=1) -%}
<Tabs>
    <Tab label="Arm-level gains/losses" />
    <Tab label="Table" />
    <Tab label="Total number of arm gains/losses" />
    <div slot="content">
        <TabContent>
            <Image src="{{job.out.outdir}}/AneuploidyScore.png" />
        </TabContent>
        <TabContent>
            <DataTable src="{{ job.out.outdir | joinpaths: 'CAA.txt' }}"
                data={ {{ job.out.outdir | joinpaths: 'CAA.txt' | datatable: sep="\t" }} } />
        </TabContent>
        <TabContent>
            <Tile>
                <pre>{{job.out.outdir | joinpaths: "AS.txt" | read}}
                </pre>
            </Tile>
        </TabContent>
    </div>
</Tabs>
{%- endmacro -%}

{%- macro head_job(job) -%}
<h1>{{job.in.segfile | stem0 }}</h1>
{%- endmacro -%}

{{ report_jobs(jobs, head_job, report_job) }}
