"""Cellranger pipeline module for BioPipen"""
from ..core.proc import Proc
from ..core.config import config


class CellRangerCount(Proc):
    """Run cellranger count

    to count gene expression and/or feature barcode reads
    requires cellranger v7+.

    Input:
        fastqs: The input fastq files
            Either a list of fastq files or a directory containing fastq files
            If a directory is provided, it should be passed as a list with one
            element.
        id: The id defining output directory. If not provided, it is inferred
            from the fastq files.
            Note that, unlike the `--id` argument of cellranger, this will not select
            the samples from `in.fastqs`. In stead, it will symlink the fastq files
            to a temporary directory with this `id` as prefix and pass that to
            cellranger.

    Output:
        outdir: The output directory

    Envs:
        ncores: Number of cores to use
        cellranger: Path to cellranger
        ref: Path of folder containing 10x-compatible transcriptome reference
        tmpdir: Path to temporary directory, used to save the soft-lined fastq files
            to pass to cellranger
        outdir_is_mounted (flag): A flag indicating whether the output directory is
            on a mounted filesystem. As of `cellranger` v9.0.1, `cellranger vdj` will
            fail when trying to copy/operate files to a mounted filesystem.
            See <https://github.com/10XGenomics/cellranger/issues/210> and
            <https://github.com/10XGenomics/cellranger/issues/250> for similar issues.
            If that is the case, set this flag to `True` to use `envs.tmpdir` as
            the output directory for `cellranger vdj`, and then move the results
            to the final output directory after `cellranger vdj` finishes.
            In this case, make sure that `envs.tmpdir` must have enough space and
            it must be a local filesystem.
        copy_outs_only (flag): If `outdir_is_mounted` is `True`, set this flag to `True`
            to only copy the `outs` folder from the temporary output directory
            to the final output directory, instead of the whole output directory.
        include_introns (flag): Set to false to exclude intronic reads in count.
        create_bam (flag): Enable or disable BAM file generation.
            This is required by cellrange v8+. When using cellrange v8-, it will be
            transformed to `--no-bam`.
        <more>: Other environment variables required by `cellranger count`
            See `cellranger count --help` for more details or
            <https://www.10xgenomics.com/support/software/cell-ranger/advanced/cr-command-line-arguments#count>
    """  # noqa: E501
    input = "fastqs:files, id"
    output = """outdir:dir:
        {%- set fastqs = in.fastqs -%}
        {%- if len(fastqs) == 1 and isdir(fastqs[0]) -%}
            {%- set fastqs = fastqs[0] | glob: "*.fastq.gz" -%}
        {%- endif -%}
        {%- if in.id -%}
            {{in.id}}
        {%- else -%}
            {%- set id = commonprefix(*fastqs) |
                regex_replace: "_L\\d+(:?_.*)?$", "" |
                regex_replace: "_S\\d+$", "" -%}
            {{- id -}}
        {%- endif -%}
    """
    lang = config.lang.python
    envs = {
        "ncores": config.misc.ncores,
        "cellranger": config.exe.cellranger,
        "ref": config.ref.ref_cellranger_gex,
        "tmpdir": config.path.tmpdir,
        "outdir_is_mounted": False,
        "copy_outs_only": True,
        "include_introns": True,
        "create_bam": False,
    }
    script = "file://../scripts/cellranger/CellRangerCount.py"
    plugin_opts = {
        "report": "file://../reports/cellranger/CellRangerCount.svelte",
        "report_paging": 5,
    }


class CellRangerVdj(Proc):
    """Run cellranger vdj

    to perform sequence assembly and paired clonotype calling.
    requires cellranger v7+.

    Input:
        fastqs: The input fastq files
            Either a list of fastq files or a directory containing fastq files
            If a directory is provided, it should be passed as a list with one
            element.
        id: The id determining the output directory. If not provided, it is inferred
            from the fastq files.

    Output:
        outdir: The output directory

    Envs:
        ncores: Number of cores to use
        cellranger: Path to cellranger
        ref: Path of folder containing 10x-compatible transcriptome reference
        tmpdir: Path to temporary directory, used to save the soft-lined fastq files
            to pass to cellranger.
        outdir_is_mounted (flag): A flag indicating whether the output directory is
            on a mounted filesystem. As of `cellranger` v9.0.1, `cellranger vdj` will
            fail when trying to copy the VDJ reference files to a mounted filesystem.
            See <https://github.com/10XGenomics/cellranger/issues/210> and
            <https://github.com/10XGenomics/cellranger/issues/250> for similar issues.
            If that is the case, set this flag to `True` to use `envs.tmpdir` as
            the output directory for `cellranger vdj`, and then move the results
            to the final output directory after `cellranger vdj` finishes.
            In this case, make sure that `envs.tmpdir` must have enough space and
            it must be a local filesystem.
        copy_outs_only (flag): If `outdir_is_mounted` is `True`, set this flag to `True`
            to only copy the `outs` folder from the temporary output directory
            to the final output directory, instead of the whole output directory.
        <more>: Other environment variables required by `cellranger vdj`
            See `cellranger vdj --help` for more details or
            <https://www.10xgenomics.com/support/software/cell-ranger/advanced/cr-command-line-arguments#vdj>
    """  # noqa: E501
    input = "fastqs:files, id"
    output = """outdir:dir:
        {%- set fastqs = in.fastqs -%}
        {%- if len(fastqs) == 1 and isdir(fastqs[0]) -%}
            {%- set fastqs = fastqs[0] | glob: "*.fastq.gz" -%}
        {%- endif -%}
        {%- if in.id -%}
            {{in.id}}
        {%- else -%}
            {%- set id = commonprefix(*fastqs) |
                regex_replace: "_L\\d+(:?_.*)?$", "" |
                regex_replace: "_S\\d+$", "" -%}
            {{- id -}}
        {%- endif -%}
    """
    lang = config.lang.python
    envs = {
        "ncores": config.misc.ncores,
        "cellranger": config.exe.cellranger,
        "ref": config.ref.ref_cellranger_vdj,
        "outdir_is_mounted": False,
        "copy_outs_only": True,
        "tmpdir": config.path.tmpdir,
    }
    script = "file://../scripts/cellranger/CellRangerVdj.py"
    plugin_opts = {
        "report": "file://../reports/cellranger/CellRangerVdj.svelte",
        "report_paging": 5,
    }


class CellRangerSummary(Proc):
    """Summarize cellranger metrics

    Input:
        indirs: The directories containing cellranger results
            from `CellRangerCount`/`CellRangerVdj`.

    Output:
        outdir: The output directory

    Envs:
        group (type=auto): The group of the samples for boxplots.
            If `None`, don't do boxplots.
            It can be a dict of group names and sample names, e.g.
            `{"group1": ["sample1", "sample2"], "group2": ["sample3"]}`
            or a file containing the group information, with the first column
            being the sample names and the second column being the group names.
            The file should be tab-delimited with no header.
    """
    input = "indirs:dirs"
    input_data = lambda ch: [list(ch.iloc[:, 0])]
    output = "outdir:dir:{{in.indirs | first | stem | append: '-etc.summary'}}"
    lang = config.lang.rscript
    script = "file://../scripts/cellranger/CellRangerSummary.R"
    envs = {"group": None}
    plugin_opts = {
        "report": "file://../reports/cellranger/CellRangerSummary.svelte",
        "report_paging": 8,
    }
