{% from "utils.liq" import report_jobs, table_of_images -%}
<script>
    import { Image } from "@@";
    import { Tabs, Tab, TabContent, Tile, UnorderedList, ListItem } from "carbon-components-svelte";
</script>

{%- macro report_job(job, h=1) -%}

<h{{h}}>Exploratory analysis</h{{h}}>

<h{{h+1}}>
    CDR3 length distribution
</h{{h+1}}>

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


<h{{h+1}}>
    Clonotype volume (# clonotypes)
</h{{h+1}}>

<p>Number of clonotypes in (groups of) samples.

{% assign volpngs = job.out.outdir | joinpaths: "volume", "volume*.png" | glob %}
{{ table_of_images(volpngs) }}


<h{{h+1}}>
    Clonotype abundances
</h{{h+1}}>

{% assign cntpngs = job.out.outdir | joinpaths: "count", "count*.png" | glob %}
{{ table_of_images(cntpngs) }}

<h{{h}}>Clonality</h{{h}}>

<h{{h+1}}>
    Top clones
</h{{h+1}}>

{% assign tcpngs = job.out.outdir | joinpaths: "top_clones", "top_clones*.png" | glob %}
{{ table_of_images(tcpngs) }}

<h{{h+1}}>
    Rare clones
</h{{h+1}}>

{% assign rcpngs = job.out.outdir | joinpaths: "rare_clones", "rare_clones*.png" | glob %}
{{ table_of_images(rcpngs) }}

<h{{h+1}}>
    Clonal space homeostasis
</h{{h+1}}>

<p>The proportion of the repertoire occupied by the clones of a given size</p>

{% assign hcpngs = job.out.outdir | joinpaths: "homeo_clones", "hom_clones*.png" | glob %}
{{ table_of_images(hcpngs) }}

<h{{h}}>
    Repertoire overlaps
</h{{h}}>

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

{%- endmacro -%}


{%- macro head_job(job) -%}
<h1>{{job.out.outdir | stem | replace: ".immunarch", ""}}</h1>
{%- endmacro -%}

{{ report_jobs(jobs, head_job, report_job) }}
