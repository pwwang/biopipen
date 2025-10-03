{% from "utils/misc.liq" import report_jobs -%}
<script>
    import { Image, DataTable, Descr, Math } from "$libs";
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
    <span class="listitem">MetabolicPathwayActivity (this page)</span>
    <Tile>
        <p>Investigating the metabolic pathways of the cells in different subsets and groups.</p>
        <p>The cells are first subset by subsets and then the metabolic activities are examined for each groups in different subsets.</p>
        <p> </p>
        <p>A pathway activity score defined as the relative gene expression value averaged over all genes in the pathway and all cells of the group.</p>
        <p> </p>
        <p>For the i-th metabolic gene, we first calculated its mean expression level across cells of the j-th cell group:
        <Math displayMode>E_{i,j} = \frac{ {\mathop {\sum }\nolimits_{k = 1}^{n_j} g_{i,k}}}{ {n_j}},\,i \in 1 \ldots M,j \in 1 \ldots N</Math>
        <p>
            In which n<sub>j</sub> is the number of cells in the j-th cell group, g<sub>i,k</sub> is the expression level of the i-th gene in the k-th cell in this cell group,
            M is the number of metabolic genes, and N is the number of cell groups. The relative expression level of the i-th gene in the j-th cell group was then
            defined as the ratio of E<sub>i,j</sub> to its average over all cell groups:
        </p>
        <Math displayMode>r_{i,j} = \frac{ {E_{i,j}}}{ {\frac{1}{N}\mathop {\sum }\nolimits_j^N E_{i,j}}}</Math>
        <p>
            Here r<sub>i,j</sub> quantifies the relative expression level of gene i in cell group j comparing to the average expression level of this gene in all cell groups.
            A r<sub>i,j</sub> value &gt;1 means that expression level of gene i is higher in cell group j compared to its average expression level over all cell groups.
            The pathway activity score for the t-th pathway and the j-th cell group was then defined as the weighted average of r<sub>i,j</sub> over all genes included
            in this pathway:
        </p>
        <Math displayMode>p_{t,j} = \frac{ {\mathop {\sum }\nolimits_{i = 1}^{m_t} w_i \times r_{i,j}}}{ {\mathop {\sum }\nolimits_{i = 1}^{m_t} w_i}}</Math>
        <p>Where p<sub>t,j</sub> represents the activity of the t-th pathway in the j-th cell group, m<sub>t</sub> is the number of genes in the pathway t, w<sub>i</sub>
        is the weighting factor equal to the reciprocal of number of pathways that include the i-th gene.
        To avoid the possibility that pathway activity scores were affected by genes with low expression level or high drop-out rates,
        we excluded the outliers in each pathway defined by genes with relative expression levels greater than three times 75th percentile or below 1/3 times 25th percentile.
        Statistical significance of higher or lower pathway activity in a specific cell group was then evaluated by a random permutation test,
        in which the cell group labels were randomly shuffled for 5000 (for the scRNA datasets) to simulate a null distribution of the pathway activity scores
        and compare to the pathway activity scores in the original, non-shuffled dataset.
        For the pathway activity score p<sub>t,j</sub>, we then calculated a p-value defined as the fraction of random pathway activity scores larger than pt,j
        (if p<sub>t,j</sub> is &gt;1) or smaller than p<sub>t,j</sub> (if p<sub>t,j</sub> is &lt;1) to assess if activity of this pathway is significantly
        higher or lower in this cell group than average.</p>
    </Tile>
</ListItem>
<ListItem>
    <a href="?proc=MetabolicPathwayHeterogeneity" class="listitem">MetabolicPathwayHeterogeneity</a>
    <Tile>
        <p>Showing metabolic pathways enriched in genes with highest contribution to the metabolic heterogeneities</p>
    </Tile>
</ListItem>
<ListItem>
    <a href="?proc=MetabolicFeatures" class="listitem">MetabolicFeatures</a>
    <Tile>
        <p>Gene set enrichment analysis against the metabolic pathways for comparisons by different groups in different subsets.</p>
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

