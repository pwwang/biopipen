{% from "utils/misc.liq" import report_jobs -%}
<script>
    import { Image, DataTable, Descr } from "$libs";
    import { Tile, UnorderedList, ListItem, Link, Tabs, Tab, TabContent } from "$ccs";
</script>

<h1>Summary</h1>
<Descr>
    The idea is to perform a regression between two groups of cells
    (e.g. Treg vs Tconv) at different length of CDR3 AA sequences.
    The regression will be performed for each physicochemical feature of the
    AA (hydrophobicity, volume and isolectric point).

    <h2>Reference</h2>
    <UnorderedList>
        <ListItem><Link href="https://www.nature.com/articles/ni.3491" target="_blank">https://www.nature.com/articles/ni.3491</Link></ListItem>
        <ListItem><Link href="https://www.nature.com/articles/s41590-022-01129-x" target="_blank">https://www.nature.com/articles/s41590-022-01129-x</Link></ListItem>
        <ListItem>Wimley, W. C. &amp; White, S. H. Experimentally determined hydrophobicity scale for proteins at membrane - interfaces. Nat. Struct. Biol. 3, 842-848 (1996).</ListItem>
        <ListItem>Hdbk of chemistry &amp; physics 72nd edition. (CRC Press, 1991).</ListItem>
        <ListItem>Zamyatnin, A. A. Protein volume in solution. Prog. Biophys. Mol. Biol. 24, 107-123 (1972).</ListItem>
    </UnorderedList>
</Descr>

{%- macro report_job(job, h=1) -%}
    {{ job | render_job: h=h }}
{%- endmacro -%}

{%- macro head_job(job) -%}
    <h1>{{job.out.outdir | stem | replace: ".scRep", ""}}</h1>
{%- endmacro -%}

{{ report_jobs(jobs, head_job, report_job) }}
