{% from "utils/misc.liq" import report_jobs, table_of_images -%}
<script>
    import { Image, DataTable } from "$libs";
    import { Tabs, Tab, TabContent, Accordion, AccordionItem } from "$ccs";
</script>

{%- macro report_job(job, h=1) -%}

<h{{h}}>Exploratory analysis</h{{h}}>

    <h{{h+1}}>CDR3 length distribution</h{{h+1}}>

    {% assign lenpngs = job.out.outdir | glob: "len", "*.png" %}
    {{ table_of_images(lenpngs) }}

    <h{{h+1}}>Clonotype volume (# clonotypes)</h{{h+1}}>

    {% assign volpngs = job.out.outdir | glob: "volume", "*.png" %}
    {{ table_of_images(volpngs) }}

    <h{{h+1}}>Clonotype abundances</h{{h+1}}>

    {% assign cntpngs = job.out.outdir | glob: "count", "*.png" %}
    {{ table_of_images(cntpngs) }}

<h{{h}}>Clonality</h{{h}}>

    <h{{h+1}}>Top clones</h{{h+1}}>

    {% assign tcpngs = job.out.outdir | glob: "top_clones", "*.png" %}
    {{ table_of_images(tcpngs) }}

    <h{{h+1}}>Rare clones</h{{h+1}}>

    {% assign rcpngs = job.out.outdir | glob: "rare_clones", "*.png" %}
    {{ table_of_images(rcpngs) }}

    <h{{h+1}}>Clonal space homeostasis</h{{h+1}}>

    <p>The proportion of the repertoire occupied by the clones of a given size</p>

    {% assign hcpngs = job.out.outdir | glob: "homeo_clones", "*.png" %}
{{ table_of_images(hcpngs) }}

<h{{h}}>Repertoire overlaps</h{{h}}>

    {% if job.index == 0 %}
    <Accordion>
        <AccordionItem title="Overlapping methods">
            <p>
                Repertoire overlap is the most common approach to measure repertoire similarity.
                Immunarch provides several indices:
            </p>
            <p>
                - number of public clonotypes (.method = "public") - a classic measure of overlap similarity.
            </p>
            <p>
                - overlap coefficient (.method = "overlap") - a normalised measure of overlap
                similarity. It is defined as the size of the intersection divided by the smaller of the size of the two sets.
            </p>
            <p>
                - Jaccard index (.method = "jaccard") - it measures similarity between finite
                sample sets, and is defined as the size of the intersection divided by the size of the union of the sample sets.
            </p>
            <p>
                - Tversky index (.method = "tversky") - an asymmetric similarity measure on sets
                that compares a variant to a prototype. If using default arguments, it’s similar to Dice’s coefficient.
            </p>
            <p>
                - cosine similarity (.method = "cosine") - a measure of similarity between two non-zero vectors
            </p>
            <p>
                - Morisita’s overlap index (.method = "morisita") - a statistical measure of dispersion of individuals in a population.
                It is used to compare overlap among samples.
            </p>
            <p>
                - incremental overlap - overlaps of the N most abundant clonotypes with incrementally growing N
                (.method = "inc+METHOD", e.g., "inc+public" or "inc+morisita").
            </p>
        </AccordionItem>
    </Accordion>
    {% endif %}

    {% for ovdir in job.out.outdir | glob: "overlap", "*" | sort %}
        {% set ovname = ovdir | basename %}
        <h{{h+1}}>{{ovname}}</h{{h+1}}>
        {% assign ovpngs = ovdir | glob: "*.png" | sort %}
        {{ table_of_images(ovpngs) }}
    {% endfor %}

<h{{h}}>Gene usage</h{{h}}>
    {% for gu_dir in job.out.outdir | glob: "gene_usage", "*" | sort %}
        {% set gu_name = gu_dir | basename %}
        {% if gu_name != "DEFAULT" %}
            <h{{h+1}}>{{gu_name}}</h{{h+1}}>
        {% endif %}
        {% assign gupngs = gu_dir | glob: "*.png" | sort %}
        {{ table_of_images(gupngs) }}
    {% endfor %}

<h{{h}}>Spectratyping</h{{h}}>
    {% for spect_sam_dir in job.out.outdir | glob: "spectratyping", "*" | sort %}
        <h{{h+1}}>{{ spect_sam_dir | basename }}</h{{h+1}}>
        {% assign spectpngs = spect_sam_dir | glob: "*.png" | sort %}
        {{ table_of_images(spectpngs) }}
    {% endfor %}

<h{{h}}>Diversity estimation</h{{h}}>
    {% assign div_met_dirs = job.out.outdir | glob: "diversity", "*" | sort %}
    {% for dm_dir in div_met_dirs %}
        <h{{h+1}}>{{dm_dir | basename}}</h{{h+1}}>
        {% if dm_dir | glob: "diversity.test.*.txt" %}
            {% assign dm_test_file = dm_dir | glob0: "diversity.test.*.txt" %}
            {% assign dm_test_method = dm_test_file | stem | replace: "diversity.test.", "" %}
            <Tabs>
                <Tab label="Plot" />
                <Tab label="Test: {{dm_test_method}}" />
                <div slot="content">
                    <TabContent>
                        <Image src={{ dm_dir | joinpaths: "diversity.png" | quote }} />
                    </TabContent>
                    <TabContent>
                        <DataTable src={{ dm_test_file | quote }} data={ {{ dm_test_file | datatable: sep="\t" }} } />
                    </TabContent>
                </div>
            </Tabs>
        {% else %}
            <Image src={{ dm_dir | joinpaths: "diversity.png" | quote }} />
        {% endif %}
    {% endfor %}

{% if job.out.outdir | glob: "rarefraction", "*" %}
<h{{h}}>Rarefaction analysis</h{{h}}>
    {% for rfdir in job.out.outdir | glob: "rarefraction", "*" | sort %}
        {% assign rfname = rfdir | basename %}
        {% if rfname != "DEFAULT" %}
            <h{{h+1}}>{{rfname}}</h{{h+1}}>
        {% endif %}
        {% assign rfpngs = rfdir | glob: "*.png" | sort %}
        {{ table_of_images(rfpngs) }}
    {% endfor %}
{% endif %}

{% if job.out.outdir | glob: "tracking", "*.png" %}
<h{{h}}>Tracking of clonotypes</h{{h}}>
    {% assign trackpngs = job.out.outdir | glob: "tracking", "*.png" | sort %}
    {{ table_of_images(trackpngs) }}
{% endif %}

<h{{h}}>Kmer and sequence motif analysis</h{{h}}>
    {% for kmerdir in job.out.outdir | glob: "kmer", "*" | sort %}
        {% assign kmercase = kmerdir | basename %}
        {% if kmercase != "DEFAULT" %}
            <h{{h+1}}>{{kmercase}}</h{{h+1}}>
        {% endif %}
        {% assign kmerpngs = kmerdir | glob: "*.png" | sort %}
        {{ table_of_images(kmerpngs) }}
    {% endfor %}

{%- endmacro -%}


{%- macro head_job(job) -%}
<h1>{{job.out.outdir | stem | replace: ".immunarch", ""}}</h1>
{%- endmacro -%}

{{ report_jobs(jobs, head_job, report_job) }}
