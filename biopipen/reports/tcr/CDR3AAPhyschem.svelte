{% from "utils/misc.liq" import report_jobs -%}
<script>
    import { Image, DataTable } from "$libs";
    import { Tile, UnorderedList, ListItem, Link, Tabs, Tab, TabContent } from "$ccs";
</script>


{%- macro report_job(job, h=1) -%}

    <h{{h}}>Summary</h{{h}}>
    <Tile>
        The idea is to perform a regression between two groups of cells
        (e.g. Treg vs Tconv) at different length of CDR3 AA sequences.
        The regression will be performed for each physicochemical feature of the
        AA (hydrophobicity, volume and isolectric point).
    </Tile>

    <h{{h+1}}>Reference</h{{h+1}}>
    <UnorderedList style="padding-left: 1rem">
        <ListItem class="ccs-li"><Link href="https://www.nature.com/articles/ni.3491" target="_blank">https://www.nature.com/articles/ni.3491</Link></ListItem>
        <ListItem class="ccs-li"><Link href="https://www.nature.com/articles/s41590-022-01129-x" target="_blank">https://www.nature.com/articles/s41590-022-01129-x</Link></ListItem>
        <ListItem class="ccs-li">Wimley, W. C. &amp; White, S. H. Experimentally determined hydrophobicity scale for proteins at membrane - interfaces. Nat. Struct. Biol. 3, 842-848 (1996).</ListItem>
        <ListItem class="ccs-li">Hdbk of chemistry &amp; physics 72nd edition. (CRC Press, 1991).</ListItem>
        <ListItem class="ccs-li">Zamyatnin, A. A. Protein volume in solution. Prog. Biophys. Mol. Biol. 24, 107-123 (1972).</ListItem>
    </UnorderedList>

    <h{{h}}>Available Cells</h{{h}}>
    <DataTable
        src={{job.out.outdir | joinpaths: "stats.txt" | quote}}
        data={ {{job.out.outdir | joinpaths: "stats.txt" | datatable: sep="\t", index_col=None}} }
        />

    {% for subsetdir in job.out.outdir | glob: "*" %}
        {% if not subsetdir | isdir %}
            {% continue %}
        {% endif %}
        {% if basename(subsetdir) == "ALL" %}
            {% set h = h - 1 %}
        {% else %}
            <h{{h}}>Subset: {{subsetdir | stem}}</h{{h}}>
        {% endif %}

        <h{{h+1}}>Estimated OR (per s.d.) for each physicochemical feature</h{{h+1}}>
        <Tabs>
            <Tab label="Plot" />
            <Tab label="Table" />
            <svelte:fragment slot="content">
                <TabContent>
                    <Image src={{subsetdir | joinpaths: "estimated_coefficients.png" | quote}} />
                </TabContent>
                <TabContent>
                    <DataTable
                        src={{subsetdir | joinpaths: "estimates.txt" | quote}}
                        data={ {{subsetdir | joinpaths: "estimates.txt" | datatable: sep="\t", index_col=None}} }
                        />
                </TabContent>
            </svelte:fragment>
        </Tabs>

        <h{{h+1}}>Hydrophobicity Distribution</h{{h+1}}>
        <Image src={{subsetdir | joinpaths: "distribution.png" | quote}} />

    {% endfor %}

{%- endmacro -%}

{%- macro head_job(job) -%}
    <h1>{{job.out.outdir | stem | replace: ".immunarch", ""}}</h1>
{%- endmacro -%}

{{ report_jobs(jobs, head_job, report_job) }}

<style>
:global(.ccs-li) {
    /* default will have messy bullet points */
    list-style-type: square;
    position: unset !important;
}
</style>
