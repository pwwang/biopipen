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

    Output:
        outdir: The output directory

    Envs:
        ncores: Number of cores to use
        cellranger: Path to cellranger
        ref: Path of folder containing 10x-compatible transcriptome reference
        tmpdir: Path to temporary directory, used to save the soft-lined fastq files
            to pass to cellranger
        include_introns (flag): Set to false to exclude intronic reads in count.
        create_bam (flag): Enable or disable BAM file generation.
            This is required by cellrange v8+. When using cellrange v8-, it will be
            transformed to `--no-bam`.
        <more>: Other environment variables required by `cellranger count`
            See `cellranger count --help` for more details or
            https://www.10xgenomics.com/support/software/cell-ranger/advanced/cr-command-line-arguments#count
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
            to pass to cellranger
        <more>: Other environment variables required by `cellranger vdj`
            See `cellranger vdj --help` for more details or
            https://www.10xgenomics.com/support/software/cell-ranger/advanced/cr-command-line-arguments#vdj
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
    output = "outdir:dir:{{in.indirs | first | stem | append: '-etc.summary'}}"
    lang = config.lang.rscript
    script = "file://../scripts/cellranger/CellRangerSummary.R"
    envs = {"group": None}
    plugin_opts = {
        "report": "file://../reports/cellranger/CellRangerSummary.svelte",
        "report_paging": 8,
    }
