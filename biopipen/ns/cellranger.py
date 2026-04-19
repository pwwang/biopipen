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
    output_flatten = True
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
    output_flatten = True
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


class CellRangerMulti(Proc):
    """Run cellranger multi

    to analyze Gene Expression, V(D)J, Feature Barcode (CRISPR, Antibody,
    Antigen), and Flex data from a single GEM well.
    Requires cellranger v7+.

    The multi config CSV is generated automatically from the `envs.gex`,
    `envs.vdj`, `envs.feature`, and `envs.libraries` options.

    Input:
        fastqs: The input FASTQ files or directories for all libraries in the
            GEM well. All FASTQs are symlinked into a temporary directory and
            used as the `fastqs` path in the `[libraries]` section of the
            generated multi config CSV (unless an explicit `fastqs` path is
            given per library in `envs.libraries`).
        id: The id defining the output directory. If not provided, it is
            inferred from the common prefix of the FASTQ filenames.

    Output:
        outdir: The output directory

    Envs:
        ncores: Number of cores to use
        cellranger: Path to cellranger
        tmpdir: Path to temporary directory, used to store the symlinked FASTQ
            files and the generated multi config CSV, and also used as the
            output directory when `outdir_is_mounted` is `True`.
        outdir_is_mounted (flag): A flag indicating whether the output directory
            is on a mounted filesystem. If `True`, cellranger multi will run in
            a local tmpdir and results are moved back to the final output
            directory. Make sure `envs.tmpdir` has enough space and is on a
            local filesystem.
        copy_outs_only (flag): If `outdir_is_mounted` is `True`, set this flag
            to `True` to only copy the `outs` folder from the temporary output
            directory to the final output directory.
        gex (ns): Options for the `[gene-expression]` section of the multi
            config CSV. Set to `False` or `null` to omit this section entirely.
            - reference: Path to the 10x Genomics-compatible transcriptome
                reference. Required for Gene Expression libraries (not needed
                for Flex with `probe_set`). Defaults to
                `config.ref.ref_cellranger_gex`.
            - create_bam (flag): Enable or disable BAM file generation.
                Default: false.
            - <more>: Any other `[gene-expression]` CSV options (underscore
                keys are converted to kebab-case). See
                <https://www.10xgenomics.com/support/software/cell-ranger/latest/advanced/cr-multi-config-csv-opts>
        vdj (ns): Options for the `[vdj]` section. Set to `False`/`null` to
            omit (default). Provide a dict with at least `reference` when V(D)J
            libraries are present.
            - reference: Path to the 10x-compatible V(D)J reference.
            - <more>: Any other `[vdj]` CSV options.
        feature (ns): Options for the `[feature]` section. Set to
            `False`/`null` to omit (default). Provide a dict with at least
            `reference` when Antibody Capture, CRISPR Guide Capture, or Antigen
            Capture libraries are present.
            - reference: Path to the Feature reference CSV.
            - <more>: Any other `[feature]` CSV options.
        libraries (list): A list of library definitions for the `[libraries]`
            section of the multi config CSV. Each entry is a dict with:
            - fastq_id (required): Sample name prefix used in the FASTQ
                filenames (the part before `_S<N>_`).
            - feature_types (required): Library type, e.g. `Gene Expression`,
                `VDJ-T`, `VDJ-B`, `Antibody Capture`, `CRISPR Guide Capture`.
            - fastqs: Optional explicit path to the FASTQ directory. If
                omitted, the temporary symlink directory created from
                `in.fastqs` is used.
            - lanes: Optional pipe-separated lane numbers, e.g. `1|2`.
    """  # noqa: E501
    input = "fastqs:files, id"
    output = """outdir:dir:
        {%- if in.id -%}
            {{in.id}}
        {%- else -%}
            {%- set fastqs = in.fastqs -%}
            {%- if len(fastqs) == 1 and isdir(fastqs[0]) -%}
                {%- set fastqs = fastqs[0] | glob: "*.fastq.gz" -%}
            {%- endif -%}
            {%- set id = commonprefix(*fastqs) |
                regex_replace: "_S\\d+.*$", "" |
                regex_replace: "_+$", "" -%}
            {{- id if id else "cellranger_multi" -}}
        {%- endif -%}
    """
    output_flatten = True
    lang = config.lang.python
    envs = {
        "ncores": config.misc.ncores,
        "cellranger": config.exe.cellranger,
        "tmpdir": config.path.tmpdir,
        "outdir_is_mounted": False,
        "copy_outs_only": True,
        "gex": {
            "reference": config.ref.ref_cellranger_gex,
            "create_bam": False,
        },
        "vdj": {},
        "feature": {},
        "libraries": [],
    }
    script = "file://../scripts/cellranger/CellRangerMulti.py"
    plugin_opts = {
        "report": "file://../reports/cellranger/CellRangerMulti.svelte",
        "report_paging": 5,
    }


class CellRangerMultiSummary(Proc):
    """Summarize cellranger multi metrics

    Reads per-sample metrics from the `outs/per_sample_outs/<sample>/metrics_summary.csv`
    of each `CellRangerMulti` output directory and generates a merged metrics table
    and QC plots.

    Input:
        indirs: The directories containing cellranger multi results
            from `CellRangerMulti`.

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
    output = "outdir:dir:{{in.indirs | first | stem | append: '-etc.multi_summary'}}"
    lang = config.lang.rscript
    script = "file://../scripts/cellranger/CellRangerMultiSummary.R"
    envs = {"group": None}
    plugin_opts = {
        "report": "file://../reports/cellranger/CellRangerMultiSummary.svelte",
        "report_paging": 8,
    }
