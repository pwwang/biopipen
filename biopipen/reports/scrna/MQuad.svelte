{% from "utils/misc.liq" import report_jobs -%}

<script>
    import { DataTable, Iframe, Descr } from "$libs";
    import { Link, Tile, InlineNotification } from "$ccs";
</script>

{%- macro report_job(job, h=1) -%}
    <h{{h}}>Heatmap of the allele frequency of qualified variants</h{{h}}>
    <Iframe src="{{ job.out.outdir | joinpath: 'top variants heatmap.pdf' }}" />

    <h{{h}}>deltaBIC cutoff</h{{h}}>
    {% if job.out.outdir | joinpath: "deltaBIC_cdf.pdf" | exists  %}
    <Descr>
    A cdf plot of deltaBIC distribution of all variants, including the cutoff determined by MQuad.
    </Descr>
    <Iframe src="{{ job.out.outdir | joinpath: 'deltaBIC_cdf.pdf' }}" />
    {% else %}
    <Descr>
    <p>No deltaBIC distribution plot available.</p>
    <p>The cutoff is specified as <code>{{ envs.t | default: envs.tenx }}</code> in the parameters.</p>
    </Descr>
    {% endif %}

    <h{{h}}>Passed Variants</h{{h}}>
    {% if job.out.outdir | joinpath: "passed_variant_names.txt" | read | strip == ""  %}
    <InlineNotification
        hideCloseButton
        lowContrast
        kind="warning">
        No variants passed the filtering criteria.
        </InlineNotification>
    {% else %}
    <DataTable src="{{ job.out.outdir | joinpath: 'passed_variant_names.txt' }}"
        data={ {{ job.out.outdir | joinpath: 'passed_variant_names.txt' | datatable: sep="\t", nrows=100, header=None, index_col=None }} }
     />
    {% endif %}

{%- endmacro -%}

{%- macro head_job(job) -%}
    <h1>{{job.in.cellsnpout | stem}}</h1>
{%- endmacro -%}

<h1>Reference</h1>

<Descr>
    <Link href="https://www.nature.com/articles/s41467-022-28845-0#Sec9" target="_blank">
    Kwok, A. W. C., Qiao, C., Huang, R., Sham, M. H., Ho, J. W., & Huang, Y. (2022).
    MQuad enables clonal substructure discovery using single cell mitochondrial variants.
    Nature communications, 13(1), 1205.
    </Link>
</Descr>
<Tile>
    <h4>Abstract</h4>
    <br />
    <p>
    Mitochondrial mutations are increasingly recognised as informative endogenous genetic markers that
    can be used to reconstruct cellular clonal structure using single-cell RNA or DNA sequencing data.
    However, identifying informative mtDNA variants in noisy and sparse single-cell sequencing data is
    still challenging with few computation methods available. Here we present an open source computational
    tool MQuad that accurately calls clonally informative mtDNA variants in a population of single cells,
    and an analysis suite for complete clonality inference, based on single cell RNA, DNA or ATAC
    sequencing data. Through a variety of simulated and experimental single cell sequencing data,
    we showed that MQuad can identify mitochondrial variants with both high sensitivity and specificity,
    outperforming existing methods by a large extent. Furthermore, we demonstrate its wide applicability
    in different single cell sequencing protocols, particularly in complementing single-nucleotide and
    copy-number variations to extract finer clonal resolution.
    </p>
</Tile>

{{ report_jobs(jobs, head_job, report_job) }}
