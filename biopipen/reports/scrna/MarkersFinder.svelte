{% from "utils/misc.liq" import report_jobs -%}
<script>
    import { Image, DataTable } from "$libs";
    import { Tabs, Tab, TabContent, InlineNotification } from "$ccs";
</script>


{%- macro report_job(job, h=1) -%}
  {% set cases = {} %}
  {% for casedir in job.out.outdir | glob: "*:*" %}
    {% set case = casedir | basename %}
    {% set section = case | split: ":" | first %}
    {% set name = case | split: ":" | last %}
    {% set _ = cases.setdefault(section, []) %}
    {% set _ = cases[section].append(name) %}
  {% endfor %}

  {% for section, names in cases.items() %}
    {% if section != "DEFAULT" or len(cases) > 1 %}
    <h{{h}}>{{section}}</h{{h}}>
    {% else %}
    {% set h = h - 1 %}
    {% endif %}

    {% for name in names %}
      {% set casedir = job.out.outdir | joinpaths: section + ":" + name %}
      <h{{h+1}}>{{name}}</h{{h+1}}>
      {% if casedir | joinpaths: "error.txt" | exists %}
        <InlineNotification
          hideCloseButton
          lowContrast
          kind="warning"
          subtitle={{ casedir | joinpaths: "error.txt" | read | quote }}
          />
      {% else %}
        <h{{h+2}}>Markers</h{{h+2}}>
        <DataTable
            src={{ casedir | joinpaths: "markers.txt" | quote }}
            data={ {{ casedir | joinpaths: "markers.txt" | datatable: sep="\t", nrows=100 }} }
            />

        <h{{h+2}}>Enrichment analysis</h{{h+2}}>
        <Tabs>
            {% for enrtxt in casedir | glob: "Enrichr-*.txt"  %}
              {% set db = enrtxt | stem | replace: "Enrichr-", "" %}
              <Tab label="{{db}}" title="{{db}}" />
            {% endfor %}
            <div slot="content">
                {% for enrtxt in casedir | glob: "Enrichr-*.txt" %}
                  {% set db = enrtxt | stem | replace: "Enrichr-", "" %}
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
  {% endfor %}
{%- endmacro -%}


{%- macro head_job(job) -%}
  <h1>{{job.in.srtobj | stem | escape}}</h1>
{%- endmacro -%}

{{ report_jobs(jobs, head_job, report_job) }}