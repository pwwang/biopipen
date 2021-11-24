{% from "utils/misc.liq" import table_of_images -%}
<script>
    import { Image } from "@@";
    import { Tabs, Tab, TabContent } from "carbon-components-svelte";
</script>

{% for case in envs.cases %}
<h1>{{case}}</h1>

{%  for binsize in envs.cases[case].binsizes %}
<h2>Binsize: {{binsize}}</h2>

{% from_ os.path import join, basename %}
{% assign manplots = [] %}
{% assign circplots = [] %}
{% assign samples = [] %}
{% for job in jobs %}
{%  set manplot = job.out.outdir | joinpaths: case, "manhattan."+str(binsize)+".*.png" | glob %}
{%  set circplot = job.out.outdir | joinpaths: case, "circular."+str(binsize)+".*.png" | glob %}
{%  set _ = manplots.append(manplot[0]) %}
{%  if len(circplot) > 0 %}
{%      set _ = circplots.append(circplot[0]) %}
{%  endif %}
{%  set _ = samples.append(basename(job.out.outdir).replace(".cnvpytor", "")) %}
{% endfor %}

{% if len(circplots) != len(samples) %}
{%  set circplots = [] %}
{% endif %}

<Tabs>
    <Tab label="Manhattan plot" />
    {% if circplots %}
    <Tab label="Circular plot" />
    {% endif %}
    <div slot="content">
        <TabContent>
            {{ table_of_images(manplots, samples, 1) }}
        </TabContent>
        {% if circplots %}
        <TabContent>
            {{ table_of_images(circplots, samples, 2) }}
        </TabContent>
        {% endif %}
    </div>
</Tabs>

{%  endfor %}

{% endfor %}
