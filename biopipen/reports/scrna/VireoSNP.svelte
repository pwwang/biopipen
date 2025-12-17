{% from "utils/misc.liq" import report_jobs -%}

<script>
    import { Image, DataTable, Descr } from "$libs";
</script>

{%- macro report_job(job, h=1) -%}
    {% if job.out.outdir | joinpath: "best_n_clones.txt" | exists  %}
    <h{{h}}>Identification for a suitable number of clones</h{{h}}>
    <Descr>
        <p>
        It is generally difficult to identify the number of clones, which is a balance between subclone resolution and analysis reliability. More clones maybe preferred, but there could be higher risk that the subclones are not genuine but rather technical noise.
        </p>

        <p>
        Here, we could use ELBO for different number of clones as an indictor for model selection. However, this is still imperfect. One empirical solution is to choose the n_clones when ELBO stops increasing dramatically.
        </p>

        <p>
        The best n_clones identified based on the ELBO elbow method is
        <code>{{ job.out.outdir | joinpath: 'best_n_clones.txt' | read }}</code>.
        </p>
    </Descr>
    <Image src="{{ job.out.outdir | joinpath: 'ELBO_n_clones.png' }}" />
    {% endif %}

    <h{{h}}>Model Fitting</h{{h}}>
    <Image src="{{ job.out.outdir | joinpath: 'model_fitting.png' }}" />

    <h{{h}}>Clone Assignment and allele frequency</h{{h}}>
    <Image src="{{ job.out.outdir | joinpath: 'assignment_allele_freq.png' }}" />

    <h{{h}}>Clone Visualization with Allele Frequency</h{{h}}>
    <Image src="{{ job.out.outdir | joinpath: 'clone_allele_heatmap.png' }}" />

    <h{{h}}>Cell-Clone Assignment Table</h{{h}}>
    <DataTable src="{{ job.out.outdir | joinpath: 'cell_clone_assignment.tsv' }}"
        data={ {{ job.out.outdir | joinpath: 'cell_clone_assignment.tsv' | datatable: sep="\t", nrows=100, header=0, index_col=None }} }
        pageSize={50}
     />
{%- endmacro -%}

{%- macro head_job(job) -%}
    <h1>{{job.in.cellsnpout | stem}}</h1>
{%- endmacro -%}

{{ report_jobs(jobs, head_job, report_job) }}
