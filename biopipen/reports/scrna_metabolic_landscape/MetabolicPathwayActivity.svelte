{% from "utils/misc.liq" import report_jobs -%}

<script>
    import { Image, Descr } from "$libs";
    import { ListItem, UnorderedList } from "$ccs";
</script>

<h1>Introduction</h1>

<Descr>
    Metabolic landscape of single cells in the tumor microenvironment.
</Descr>

<h2>Workflow of the original analysis</h2>
<Image src="https://github.com/LocasaleLab/Single-Cell-Metabolic-Landscape/raw/master/pipeline.png" />

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
The cells are grouped at 2 dimensions: `grouping`, usually the cell types, and `subsetting`, usually
the groups that bring biological meaning (i.e. different timepoints or sample types (tumor/normal)).
</Descr>

<UnorderedList>
<ListItem>
    MetabolicPathwayActivity (this page)
    <p>Investigating the metabolic pathways of the cells in different groups and subsets.</p>
    <p>The cells are first grouped by subsets and then the metabolic activities are examined for each groups in different subsets.</p>
</ListItem>
<ListItem>
    <a href="../MetabolicPathwayHeterogeneity/index.html">MetabolicPathwayHeterogeneity</a>
    <p>Showing metabolic pathways enriched in genes with highest contribution to the metabolic heterogeneities</p>
</ListItem>
<ListItem>
    <a href="../MetabolicFeatures/index.html">MetabolicFeatures</a>
    <p>Gene set enrichment analysis against the metabolic pathways for groups in different subsets.</p>
</ListItem>
<ListItem>
    <a href="../MetabolicFeaturesIntraSubsets/index.html">MetabolicFeaturesIntraSubsets</a>
    <p>Gene set enrichment analysis against the metabolic pathways for subsets in different groups.</p>
</ListItem>
</UnorderedList>


{%- macro report_job(job, h=2) -%}
  {%- for ssdir in job.out.outdir | glob: "*" -%}
  {%- if not isdir(ssdir) -%}
    {%- continue -%}
  {%- endif -%}
  <h{{h}}>{{ ssdir | stem }}</h{{h}}>

  <h{{ h+1 }}>Metabolic pathway activities by <code>{{envs.grouping}}</code></h{{ h+1 }}>
  <Image src="{{ssdir | joinpaths: 'KEGGpathway_activity_heatmap.png'}}" />

  <h{{ h+1 }}>Distributions of pathway activities by <code>{{envs.grouping}}</code></h{{ h+1 }}>
  <Image src="{{ssdir | joinpaths: 'pathway_activity_violinplot.png'}}" />
  {%- endfor -%}

  {% if job.out.outdir | glob: "*.group-*.png" -%}
  <h{{h}}>Merged heatmaps</h{{h}}>
  {% for group_hm in job.out.outdir | glob: "*.group-*.png" -%}
    {%- if group_hm.endswith(".group-unclustered.png") -%}
        <h{{h+1}}>{{group_hm | stem | replace: ".group-unclustered", " (Group Unclustered)"}}</h{{h+1}}>
        <Image src="{{group_hm}}" />
    {%- else -%}
        <h{{h+1}}>{{group_hm | stem | replace: ".group-clustered", " (Group Clustered)"}}</h{{h+1}}>
        <Image src="{{group_hm}}" />
    {%- endif -%}
  {%- endfor -%}
  {%- endif -%}
{%- endmacro -%}

{%- macro head_job(job) -%}
  <h1>{{job.in.sobjfile | stem | escape}}</h1>
{%- endmacro -%}

{{ report_jobs(jobs, head_job, report_job) }}
