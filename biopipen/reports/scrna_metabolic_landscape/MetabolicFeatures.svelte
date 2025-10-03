{% from "utils/misc.liq" import report_jobs -%}
<script>
    import { Image, DataTable, Descr } from "$libs";
    import { Tabs, Tab, TabContent, UnorderedList, ListItem, InlineNotification, Tile } from "$ccs";
</script>

<h1>Introduction</h1>

<Descr>
    Metabolic landscape of single cells in the tumor microenvironment.
</Descr>

<h2>Workflow of the original analysis</h2>
<Image src="https://raw.githubusercontent.com/LocasaleLab/Single-Cell-Metabolic-Landscape/master/pipeline.png" />

<h2>Reference</h2>
<UnorderedList>
    <ListItem><a href="https://www.nature.com/articles/s41467-019-11738-0" target="_blank">
        Zhengtao, Ziwei Dai, and Jason W. Locasale.
        "Metabolic landscape of the tumor microenvironment at single cell resolution."
        Nature communications 10.1 (2019): 1-12.
    </a></ListItem>
    <ListItem><a href="https://github.com/LocasaleLab/Single-Cell-Metabolic-Landscape" target="_blank">
        Orginal pipeline
    </a></ListItem>
</UnorderedList>

<h2>Analyses with this pipeline</h2>

<Descr>
The cells are grouped at 2 dimensions: `subset_by`, usually the clinic groups that bring biological meaning
(i.e. different timepoints or sample types (tumor/normal)), and `group_by`, usually the cell types.
</Descr>

<UnorderedList>
<ListItem>
    <a href="?proc=MetabolicPathwayActivity" class="listitem">MetabolicPathwayActivity</a>
    <Tile><p>Investigating the metabolic pathways of the cells in different subsets and groups.</p></Tile>
</ListItem>
<ListItem>
    <a href="?proc=MetabolicPathwayHeterogeneity" class="listitem">MetabolicPathwayHeterogeneity</a>
    <Tile><p>Showing metabolic pathways enriched in genes with highest contribution to the metabolic heterogeneities</p></Tile>
</ListItem>
<ListItem>
    <span class="listitem">MetabolicFeatures (this page)</span>
    <Tile>
    <p>Gene set enrichment analysis against the metabolic pathways for comparisons by different groups in different subsets.</p>
    <p>The metabolic features are actual gene set enrichment analysis (GSEA) results for the metabolic pathways with given comparisons.</p>
    </Tile>
</ListItem>
</UnorderedList>

<style>
.listitem {
    font-size: large;
    font-weight: bold;
    margin: 1rem 0 0.5rem 0;
    display: inline-block;
}
</style>

{%- macro report_job(job, h=1) -%}
    {{ job | render_job: h=h }}
{%- endmacro -%}

{%- macro head_job(job) -%}
    <h1>{{job.in | attr: "values" | call | first | stem0 | escape}}</h1>
{%- endmacro -%}

{{ report_jobs(jobs, head_job, report_job) }}
