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

- **MetabolicFeaturesIntraSubsets** (optional)

    Gene set enrichment analysis against the metabolic pathways for subsets based on the designed comparison in different groups.

</Markdown>
</div>

{%- macro report_job(job, h=2) -%}
{%- set name = job.in.configfile | config: "toml" | attr: "name" -%}
{%- if name or proc.size > 1 -%}
{%- else -%}
{%- set h = 1 -%}
{%- endif -%}

{%- for ssdir in job.out.outdir | glob: "*" -%}
<h{{h}}>{{ ssdir | stem }}</h{{h}}>

<h{{ h+1 }}>Metabolic pathway activities in cell types</h{{ h+1 }}>
<Image src="{{ssdir | joinpaths: 'KEGGpathway_activity_heatmap.png'}}" />

<h{{ h+1 }}>Distributions of pathway activities in different cell types</h{{ h+1 }}>
<Image src="{{ssdir | joinpaths: 'pathway_activity_violinplot.png'}}" />

{%- endfor -%}
{%- endmacro -%}

{%- macro head_job(job) -%}
{%- set name = job.in.configfile | config: "toml" | attr: "name" -%}
{%- if name or proc.size > 1 -%}
<h1>{{name | default: "Job #" + str(job.index+1) | escape}}</h1>
{%- endif -%}
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
