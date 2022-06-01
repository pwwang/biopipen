{% from "utils/misc.liq" import report_jobs, table_of_images -%}
<script>
    import { Image } from "@@";
    import { Tabs, Tab, TabContent, Tile, UnorderedList, ListItem } from "carbon-components-svelte";
</script>

{%- macro report_job(job, h=1) -%}

<h{{h}}>Exploratory analysis</h{{h}}>

<h{{h+1}}>CDR3 length distribution</h{{h+1}}>

<Tabs>
    <Tab label="By aa" />
    <Tab label="By nt" />
    <div slot="content">
        <TabContent>
            {% assign lenpngs = job.out.outdir | joinpaths: "len-aa", "len*.png" | glob %}
            {{ table_of_images(lenpngs) }}
        </TabContent>
        <TabContent>
            {% assign lenpngs = job.out.outdir | joinpaths: "len-nt", "len*.png" | glob %}
            {{ table_of_images(lenpngs) }}
        </TabContent>
    </div>
</Tabs>


<h{{h+1}}>Clonotype volume (# clonotypes)</h{{h+1}}>

<p>Number of clonotypes in (groups of) samples.

{% assign volpngs = job.out.outdir | joinpaths: "volume", "volume*.png" | glob %}
{{ table_of_images(volpngs) }}


<h{{h+1}}>Clonotype abundances</h{{h+1}}>

{% assign cntpngs = job.out.outdir | joinpaths: "count", "count*.png" | glob %}
{{ table_of_images(cntpngs) }}

<h{{h}}>Clonality</h{{h}}>

<h{{h+1}}>
    Top clones
</h{{h+1}}>

{% assign tcpngs = job.out.outdir | joinpaths: "top_clones", "top_clones*.png" | glob %}
{{ table_of_images(tcpngs) }}

<h{{h+1}}>Rare clones</h{{h+1}}>

{% assign rcpngs = job.out.outdir | joinpaths: "rare_clones", "rare_clones*.png" | glob %}
{{ table_of_images(rcpngs) }}

<h{{h+1}}>
    Clonal space homeostasis
</h{{h+1}}>

<p>The proportion of the repertoire occupied by the clones of a given size</p>

{% assign hcpngs = job.out.outdir | joinpaths: "homeo_clones", "hom_clones*.png" | glob %}
{{ table_of_images(hcpngs) }}

<h{{h}}>Repertoire overlaps</h{{h}}>

{% if job.index == 0 %}
<Tile>
<p>
    Repertoire overlap is the most common approach to measure repertoire similarity.
    It is achieved by computation of specific statistics on clonotypes shared between
    given repertoires, also called “public” clonotypes.
    immunarch provides several indices: - number of public clonotypes (.method = "public")
    - a classic measure of overlap similarity.
</p>

<UnorderedList style="padding-left: 1rem">
    <ListItem>
        overlap coefficient (.method = "overlap") - a normalised measure of overlap
        similarity. It is defined as the size of the intersection divided by the smaller of the size of the two sets.
    </ListItem>
    <ListItem>
        Jaccard index (.method = "jaccard") - it measures similarity between finite
        sample sets, and is defined as the size of the intersection divided by the size of the union of the sample sets.
    </ListItem>
    <ListItem>
        Tversky index (.method = "tversky") - an asymmetric similarity measure on sets
        that compares a variant to a prototype. If using default arguments, it’s similar to Dice’s coefficient.
    </ListItem>
    <ListItem>
        cosine similarity (.method = "cosine") - a measure of similarity between two non-zero vectors
    </ListItem>
    <ListItem>
        Morisita’s overlap index (.method = "morisita") - a statistical measure of dispersion of individuals in a population.
        It is used to compare overlap among samples.
    </ListItem>
    <ListItem>
        incremental overlap - overlaps of the N most abundant clonotypes with incrementally growing N
        (.method = "inc+METHOD", e.g., "inc+public" or "inc+morisita").
    </ListItem>
</UnorderedList>
</Tile>
{% endif %}

{% assign ovpngs = job.out.outdir | joinpaths: "overlap", "overlap-*.png" | glob %}
<Tabs>
    {% for ovpng in ovpngs %}
    {%  assign ovmethod = ovpng | stem | replace: "overlap-", "" %}
    <Tab label="Method: {{ovmethod}}" />
    {% endfor %}

    <div slot="content">
        {% for ovpng in ovpngs %}
        {%  assign ovmethod = ovpng | stem | replace: "overlap-", "" %}
        {%  assign ovredpngs = job.out.outdir | joinpaths: "overlap", "overlapanalysis-"+ ovmethod +"-*.png" | glob %}
        <TabContent>
            <Image src={{ovpng | quote}} />
            {{ table_of_images(ovredpngs) }}
        </TabContent>
        {% endfor %}
    </div>
</Tabs>


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
{% assign div_met_dirs = job.out.outdir | joinpaths: "diversity", "*" | glob %}
<h{{h+1}}>
    Sample diversity (all_clones)
</h{{h+1}}>
<Tabs>
    {% for dm_dir in div_met_dirs %}
    <Tab label="Method: {{ dm_dir | basename }}" />
    {% endfor %}
    <div slot="content">
        {% for dm_dir in div_met_dirs %}
        <TabContent>
            {% assign divpngs = dm_dir| joinpaths: "diversity-1-*.png" | glob %}
            {{ table_of_images(divpngs) }}
        </TabContent>
        {% endfor %}
    </div>
</Tabs>

<h{{h+1}}>
    Sample diversity (clone size >= 2)
</h{{h+1}}>
<Tabs>
    {% for dm_dir in div_met_dirs %}
    <Tab label="Method: {{ dm_dir | basename }}" />
    {% endfor %}
    <div slot="content">
        {% for dm_dir in div_met_dirs %}
        <TabContent>
            {% assign divpngs = dm_dir| joinpaths: "diversity-2-*.png" | glob %}
            {{ table_of_images(divpngs) }}
        </TabContent>
        {% endfor %}
    </div>
</Tabs>

<h{{h+1}}>
    Sample diversity (clone size >= 3)
</h{{h+1}}>
<Tabs>
    {% for dm_dir in div_met_dirs %}
    <Tab label="Method: {{ dm_dir | basename }}" />
    {% endfor %}
    <div slot="content">
        {% for dm_dir in div_met_dirs %}
        <TabContent>
            {% assign divpngs = dm_dir| joinpaths: "diversity-3-*.png" | glob %}
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
