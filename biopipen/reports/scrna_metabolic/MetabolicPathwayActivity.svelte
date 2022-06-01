{% from "utils/misc.liq" import report_jobs -%}

<script>
    import { Image } from "@@";
</script>

<h1>Introduction</h1>

<div class="markdown-wrap">
<Markdown>

Metabolic landscape of single cells in the tumor microenvironment.

## Workflow of the original analysis
<Image src="https://github.com/LocasaleLab/Single-Cell-Metabolic-Landscape/raw/master/pipeline.png" />

## Reference
- [Xiao, Zhengtao, Ziwei Dai, and Jason W. Locasale. "Metabolic landscape of the tumor microenvironment at single cell resolution." Nature communications 10.1 (2019): 1-12.](https://www.nature.com/articles/s41467-019-11738-0)
- [Original pipeline](https://github.com/LocasaleLab/Single-Cell-Metabolic-Landscape)

## Analyses with this pipeline

The cells are grouped at 2 dimensions: `grouping`, usually the cell types, and `subsetting`, usually
the groups that bring biological meaning (i.e. different timepoints or sample types (tumor/normal)).

- **MetabolicPathwayActivity**

    Investigating the metabolic pathways of the cells in different groups and subsets.

    The cells are first grouped by subsets and then the metabolic activities are examined for each groups in different subsets.

- **MetabolicPathwayHeterogeneity**

    Showing metabolic pathways enriched in genes with highest contribution to the metabolic heterogeneities

- **MetabolicFeatures**

    Gene set enrichment analysis against the metabolic pathways for groups in different subsets.

- **MetabolicFeaturesIntraSubsets**

    Gene set enrichment analysis against the metabolic pathways for subsets based on the designed comparison in different groups.

</Markdown>
</div>

{%- macro report_job(job, h=2) -%}
{%- set ssdir = job.out.outdir | glob: "*" | first -%}

<h{{ h }}>Metabolic pathway activities in cell types</h{{ h }}>
<Image src="{{ssdir | joinpaths: 'KEGGpathway_activity_heatmap.png'}}" />

<h{{ h }}>Distributions of pathway activities in different cell types</h{{ h }}>
<Image src="{{ssdir | joinpaths: 'pathway_activity_violinplot.png'}}" />

{%- endmacro -%}

{%- macro head_job(job) -%}
{% set name = job.out.outdir | glob: "*" | first | stem %}
<h1>{{name | escape}}</h1>
{%- endmacro -%}

{{ report_jobs(jobs, head_job, report_job) }}

<style>

.markdown-wrap ul {
    list-style-position: inside;
    list-style-type: disc;
    margin-top: 0.2rem;
}

.markdown-wrap ul li p {
    padding-left: 1.2rem;
}

.markdown-wrap ul li > p:first-of-type {
    padding-left: 0;
    display: inline;
}
</style>
