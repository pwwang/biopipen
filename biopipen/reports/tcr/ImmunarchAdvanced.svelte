{% from "utils.liq" import report_jobs, table_of_images -%}
<script>
    import { Image } from "@@";
    import { Tabs, Tab, TabContent, Tile, UnorderedList, ListItem } from "carbon-components-svelte";
</script>

{%- macro report_job(job, h=1) -%}

<h{{h}}>Gene usage</h{{h}}>
<h{{h+1}}>
    Top {{envs.gu_top}} genes
</h{{h+1}}>
{% for gupng in job.out.outdir | joinpaths: "gene_usage", "gene_usage*.png" | glob %}
<Image src={{gupng | quote}} />
{% endfor %}

<h{{h+1}}>
    Gene usage analysis
</h{{h+1}}>
{% assign guapngs = job.out.outdir | joinpaths: "gene_usage_analysis", "gene_usage_analysis*.png" | glob %}
{{ table_of_images(guapngs) }}

<h{{h}}>Spectratyping</h{{h}}>
{% for spect_sam_dir in job.out.outdir | joinpaths: "spectratyping", "*" | glob %}
<h{{h+1}}>
    {{ spect_sam_dir | basename }}
</h{{h+1}}>
{% assign spectpngs = spect_sam_dir | joinpaths: "spectratyping*.png" | glob %}
{{ table_of_images(spectpngs) }}
{% endfor %}

<h{{h}}>Diversity estimation</h{{h}}>
<h{{h+1}}>
    Sample diversity
</h{{h+1}}>
{% assign div_met_dirs = job.out.outdir | joinpaths: "diversity", "*" | glob %}
<Tabs>
    {% for dm_dir in div_met_dirs %}
    <Tab label="Method: {{ dm_dir | basename }}" />
    {% endfor %}
    <div slot="content">
        {% for dm_dir in div_met_dirs %}
        <TabContent>
            {% assign divpngs = dm_dir| joinpaths: "diversity*.png" | glob %}
            {{ table_of_images(divpngs) }}
        </TabContent>
        {% endfor %}
    </div>
</Tabs>

<h{{h+1}}>
    Rarefaction analysis
</h{{h+1}}>
{% assign rfpngs = job.out.outdir | joinpaths: "raref", "raref*.png" | glob %}
{{ table_of_images(rfpngs) }}

<h{{h}}>Tracking of clonotypes</h{{h}}>
{% for name, value in envs.tracking_target.items() %}
<h{{h+1}}>Clonotypes: {{name}}</h{{h+1}}>
<Image src={{ job.out.outdir | joinpaths: "tracking", "tracking_" + name + ".png" | quote }} />
{% endfor %}

<h{{h}}>Kmer and sequence motif analysis</h{{h}}>
{% for kmerdir in job.out.outdir | joinpaths: "kmer", "kmer_*" | glob %}
{%  assign k = kmerdir | stem | replace: "kmer_", "" | float | int %}
<h{{h+1}}>K = {{k}}</h{{h+1}}>

{%  for kmerpng in kmerdir | joinpaths: "head_*.png" | glob %}
        {% assign head = kmerpng | stem0 | replace: "head_", "" | int %}
        <h{{h+2}}>The most {{head}} abundant kmers</h{{h+2}}>
        <Image src={{kmerpng | quote}} />
{%  endfor %}

<h{{h+2}}>Motif analysis</h{{h+2}}>

{% if job.index == 0 %}
<Tile>
<p>
    The method determines which matrix to compute:
</p>

<UnorderedList style="padding-left: 1rem">
    <ListItem>
        `freq` - position frequency matrix (PFM);
    </ListItem>
    <ListItem>
        `prob` - position probability matrix (PPM);
    </ListItem>
    <ListItem>
        `wei` - position weight matrix (PWM)
    </ListItem>
    <ListItem>
        `self` - self-information matrix.
    </ListItem>
</UnorderedList>
</Tile>
{% endif %}

{%  for motdir in kmerdir | joinpaths: "motif_*" | glob %}
{%      assign motmethod = motdir | stem | replace: "motif_", "" %}
{%      assign motpngs = motdir | joinpaths: "*.png" | glob %}
<h{{h+3}}>Method: {{ motmethod }}</h{{h+3}}>
{{ table_of_images(motpngs) }}
{%  endfor %}

{% endfor %}

{%- endmacro -%}

{%- macro head_job(job) -%}
<h1>{{job.out.outdir | stem | replace: ".immunarch", ""}}</h1>
{%- endmacro -%}

{{ report_jobs(jobs, head_job, report_job) }}
