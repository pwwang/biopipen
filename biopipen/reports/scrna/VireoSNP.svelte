{% from "utils/misc.liq" import report_jobs -%}

<script>
    import { Image, DataTable, Descr } from "$libs";
</script>

{%- macro report_job(job, h=1) -%}
    {% if job.outdir | joinpath: "best_n_clones.txt" | exists  %}
    <h{{h}}>Identification for a suitable number of clones</h{{h}}>
    <Descr>
        It is generally difficult to identify the number of clones, which is a balance between subclone resolution and analysis reliability. More clones maybe preferred, but there could be higher risk that the subclones are not genuine but rather technical noise.

        Here, we could use ELBO for different number of clones as an indictor for model selection. However, this is still imperfect. One empirical suggestion is to choose the n_clones when ELBO stops increasing dramatically.

        The best n_clones identified based on the ELBO elbow method is {{ job.outdir | joinpath: 'best_n_clones.txt' | read }}.
    </Descr>
    <Image src="{{ job.outdir | joinpath: 'ELBO_n_clones.png' }}" />
    {% endif %}

    <h{{h}}>Model Fitting</h{{h}}>
    <Image src="{{ job.outdir | joinpath: 'model_fitting.png' }}" />

    <h{{h}}>Clone Assignment and  allele frequency</h{{h}}>
    <Image src="{{ job.outdir | joinpath: 'assignment_allele_freq.png' }}" />

    <h{{h}}>Clone Visualization with Allele Frequency</h{{h}}>
    <Image src="{{ job.outdir | joinpath: 'clone_allele_heatmap.png' }}" />
{%- endmacro -%}

{%- macro head_job(job) -%}
    <h1>{{job.in.cellsnpout | stem}}</h1>
{%- endmacro -%}

{{ report_jobs(jobs, head_job, report_job) }}
