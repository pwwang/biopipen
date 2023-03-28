{% from "utils/misc.liq" import report_jobs -%}
<script>
    import { Image, DataTable } from "$libs";
    import { Tabs, Tab, TabContent, InlineNotification } from "$ccs";
</script>


{%- macro report_job(job, h=1) -%}
  {% for casedir in job.out.outdir | glob: "*" %}
    {% set case = casedir | basename %}
    <h{{h}}>{{case}}</h{{h}}>

    {% if casedir | joinpaths: "error.txt" | exists %}
      <InlineNotification
        hideCloseButton
        lowContrast
        kind="warning"
        subtitle={{ casedir | joinpaths: "error.txt" | read | quote }}
        />
    {% else %}
      <h{{h+1}}>Markers</h{{h+1}}>
      <DataTable
          src={{ casedir | joinpaths: "markers.txt" | quote }}
          data={ {{ casedir | joinpaths: "markers.txt" | datatable: sep="\t", nrows=100 }} }
          />

      <h{{h+1}}>Enrichment analysis</h{{h+1}}>
      <Tabs>
          {% for enrtxt in casedir | glob: "Enrichr-*.txt"  %}
          {%  set db = enrtxt | stem | replace: "Enrichr-", "" %}
          <Tab label="{{db}}" title="{{db}}" />
          {% endfor %}
          <div slot="content">
              {% for enrtxt in casedir | glob: "Enrichr-*.txt" %}
              {%  set db = enrtxt | stem | replace: "Enrichr-", "" %}
              <TabContent>
                  <Image src={{casedir | joinpaths: "Enrichr-" + db + ".png" | quote}} />
                  <DataTable
                      src={{ enrtxt | quote }}
                      data={ {{ enrtxt | datatable: sep="\t", nrows=100 }} }
                      />
              </TabContent>
              {% endfor %}
          </div>
      </Tabs>
    {% endif %}
  {% endfor %}
{%- endmacro -%}


{%- macro head_job(job) -%}
  <h1>{{job.in.srtobj | stem | escape}}</h1>
{%- endmacro -%}

{{ report_jobs(jobs, head_job, report_job) }}