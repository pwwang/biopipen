{% from "utils/misc.liq" import table_of_images -%}
<script>
    import { Image } from "$libs";
    import { Tabs, Tab, TabContent } from "$ccs";
</script>

{% for binsize in envs.binsizes %}
<h1>Binsize: {{binsize}}</h1>

{% from_ os.path import join, basename %}
{% assign manplots = [] %}
{% assign circplots = [] %}
{% assign samples = [] %}
{% for job in jobs %}
{%  set manplot = job.out.outdir | glob: "manhattan."+str(binsize)+".*.png" %}
{%  set circplot = job.out.outdir | glob: "circular."+str(binsize)+".*.png" %}
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

{% endfor %}
