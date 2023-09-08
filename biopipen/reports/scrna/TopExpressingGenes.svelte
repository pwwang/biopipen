{% from "utils/misc.liq" import report_jobs -%}
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

            <h{{h+1}}>Markers</h{{h+1}}>
            <DataTable
                src={{ casedir | joinpaths: "exprn.txt" | quote }}
                data={ {{ casedir | joinpaths: "exprn.txt" | datatable: sep="\t", nrows=100 }} }
                />

            <h{{h+1}}>Enrichment analysis</h{{h+1}}>
            {{ enrichr_report(casedir) }}

        {%- endfor -%}
    {%- else -%}
        {%- for secdir in secdirs -%}
            {%- set sec = secdir | basename -%}
            <h{{h}}>{{sec | escape}}</h{{h}}>
            {%- for casedir in secdir | glob: "*" -%}
                <h{{h+1}}>{{casedir | basename | escape}}</h{{h+1}}>

                <h{{h+2}}>Markers</h{{h+2}}>
                <DataTable
                    src={{ casedir | joinpaths: "exprn.txt" | quote }}
                    data={ {{ casedir | joinpaths: "exprn.txt" | datatable: sep="\t", nrows=100 }} }
                    />

                <h{{h+2}}>Enrichment analysis</h{{h+2}}>
                {{ enrichr_report(casedir) }}

            {%- endfor -%}
        {%- endfor -%}
    {%- endif -%}
{%- endmacro -%}

{%- macro head_job(job) -%}
  <h1>{{job.in.srtobj | stem0 | escape}}</h1>
{%- endmacro -%}

{{ report_jobs(jobs, head_job, report_job) }}