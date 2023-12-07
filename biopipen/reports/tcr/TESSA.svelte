{% from "utils/misc.liq" import report_jobs, table_of_images -%}
<script>
    import { Image, DataTable, Descr } from "$libs";
    import { UnorderedList, ListItem } from "$ccs";
</script>

<Descr>
    <h1>Introduction</h1>
    <p><a href="https://github.com/jcao89757/TESSA" target="_blank">Tessa</a> is a Bayesian model to integrate T cell receptor (TCR) sequence profiling with transcriptomes of T cells. Enabled by the recently developed single cell sequencing techniques, which provide both TCR sequences and RNA sequences of each T cell concurrently, Tessa maps the functional landscape of the TCR repertoire, and generates insights into understanding human immune response to diseases. As the first part of tessa, BriseisEncoder is employed prior to the Bayesian algorithm to capture the TCR sequence features and create numerical embeddings. We showed that the reconstructed Atchley Factor matrices and CDR3 sequences, generated through the numerical embeddings, are highly similar to their original counterparts. The CDR3 peptide sequences are constructed via a RandomForest model applied on the reconstructed Atchley Factor matrices.</p>
    <p></p>
    <p>For more information, please refer to the following papers:</p>
    <UnorderedList>
        <ListItem><a href="https://www.nature.com/articles/s41592-020-01020-3" target="_blank">Mapping the Functional Landscape of TCR Repertoire</a>, Zhang, Z., Xiong, D., Wang, X. et al. 2021.</ListItem>
        <ListItem><a href="https://www.nature.com/articles/s42256-021-00383-2" target="_blank">Deep learning-based prediction of the T cell receptor–antigen binding specificity</a>, Lu, T., Zhang, Z., Zhu, J. et al. 2021.</ListItem>
    </UnorderedList>
</Descr>

{%- macro report_job(job, h=1) -%}
    {{ job | render_job: h=h }}
{%- endmacro -%}

{%- macro head_job(job) -%}
    <h1>{{job.in.immdata | stem | escape}}</h1>
{%- endmacro -%}

{{ report_jobs(jobs, head_job, report_job) }}
