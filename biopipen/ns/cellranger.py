"""Cellranger pipeline module for BioPipen"""
from ..core.proc import Proc
from ..core.config import config


class CellRangerCount(Proc):
    """Run cellranger count

    to count gene expression and/or feature barcode reads

    Input:
        fastqs: The input fastq files
            Either a list of fastq files or a directory containing fastq files
            If a directory is provided, it should be passed as a list with one
            element.

    Output:
        outdir: The output directory

    Envs:
        ncores: Number of cores to use
        cellranger: Path to cellranger
        ref: Path of folder containing 10x-compatible transcriptome reference
        tmpdir: Path to temporary directory, used to save the soft-lined fastq files
            to pass to cellranger
        include_introns: Set to false to exclude intronic reads in count.
        <more>: Other environment variables required by `cellranger count`
            See `cellranger count --help` for more details or
            https://www.10xgenomics.com/support/software/cell-ranger/advanced/cr-command-line-arguments#count
    """  # noqa: E501
    input = "fastqs:files"
    output = """outdir:dir:
        {%- set fastqs = in.fastqs -%}
        {%- if len(fastqs) == 1 and isdir(fastqs[0]) -%}
            {%- set fastqs = fastqs[0] | glob: "*.fastq.gz" -%}
        {%- endif -%}
        {%- set sample = commonprefix(*fastqs) |
            regex_replace: "_L\\d+_$", "" |
            regex_replace: "_S\\d+$", "" -%}
        {{- sample -}}
    """
    lang = config.lang.python
    envs = {
        "ncores": config.misc.ncores,
        "cellranger": config.exe.cellranger,
        "ref": config.ref.ref_cellranger_gex,
        "tmpdir": config.path.tmpdir,
        "include_introns": "true",
    }
    script = "file://../scripts/cellranger/CellRangerCount.py"
    plugin_opts = {
        "report": "file://../reports/cellranger/CellRangerCount.svelte",
    }


class CellRangerVdj(Proc):
    """Run cellranger vdj

    to perform sequence assembly and paired clonotype calling

    Input:
        fastqs: The input fastq files
            Either a list of fastq files or a directory containing fastq files
            If a directory is provided, it should be passed as a list with one
            element.

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
    input = "fastqs:files"
    output = """outdir:dir:
        {%- set fastqs = in.fastqs -%}
        {%- if len(fastqs) == 1 and isdir(fastqs[0]) -%}
            {%- set fastqs = fastqs[0] | glob: "*.fastq.gz" -%}
        {%- endif -%}
        {%- set sample = commonprefix(*fastqs) |
            regex_replace: "_L\\d+_$", "" |
            regex_replace: "_S\\d+$", "" -%}
        {{- sample -}}
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
    }
