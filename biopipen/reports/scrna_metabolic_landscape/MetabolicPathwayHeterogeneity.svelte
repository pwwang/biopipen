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
    <a href="../MetabolicPathwayActivity/index.html">MetabolicPathwayActivity</a>
    <Tile>
    <p>Investigating the metabolic pathways of the cells in different subsets and groups.</p>
    </Tile>
</ListItem>
<ListItem>
    MetabolicPathwayHeterogeneity (this page)
    <Tile>
    <p>Showing metabolic pathways enriched in genes with highest contribution to the metabolic heterogeneities</p>
    <p>
        The PCA analysis was applied on normalized expression values.
        The function prcomp in R was used to perform the PCA analysis.
        For each metabolic gene, we computed its PCA score defined as the sum of absolute values of the loadings of this gene in the top PCs
        that in total account for certain variance to measure variability of gene expression across cells.
        We then sorted the PCA scores of the genes in descending order and applied GSEA analysis to the ranked list of genes to identify metabolic pathways
        enriched in genes with highest variability.
    </p>
    </Tile>
</ListItem>
<ListItem>
    <a href="../MetabolicFeatures/index.html">MetabolicFeatures</a>
    <Tile>
    <p>Gene set enrichment analysis against the metabolic pathways for comparisons by different groups in different subsets.</p>
    </Tile>
</ListItem>
</UnorderedList>

{%- macro report_job(job, h=1) -%}
    {{ job | render_job: h=h }}
{%- endmacro -%}

{%- macro head_job(job) -%}
    <h1>{{job.in | attr: "values" | call | first | stem0 | escape}}</h1>
{%- endmacro -%}

{{ report_jobs(jobs, head_job, report_job) }}
