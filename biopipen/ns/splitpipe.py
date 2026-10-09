"""Split-pipe processes for Parse Biosciences split-pool data"""
from ..core.proc import Proc
from ..core.config import config


class SplitPipeRef(Proc):
    """Build an indexed reference genome for split-pipe

    Runs `split-pipe --mode mkref`, which is required before any data can be
    processed. The output of this process is used by `SplitPipeCount` with
    `envs.genome_dir`.
    Requires the Parse Biosciences Pipeline (split-pipe) v1.9.0+.

    Input:
        fasta: The input reference sequence file(s) in fasta format
            Multiple files are allowed, e.g. for a species-mix reference. They
            should correspond to `envs.genome_name` and `in.genes` in the same
            order.
        genes: The input reference gene annotation file(s) in GTF format
            Should correspond to `envs.genome_name` and `in.fasta` in the same
            order.

    Output:
        outdir: The output directory, named by `envs.genome_name`

    Envs:
        splitpipe: Path to the split-pipe executable
        nthreads: Number of threads to use
        ncores: Alias for `nthreads`. Takes precedence when both are given.
        genome_name (type=list): The name(s) of the genome(s), e.g. `hg38`, or
            `["hg38", "GRCm39"]` for a species-mix reference. When multiple
            names are given, the output directory is named after them joined
            by `_`.
        parfile: Path to the parameter file with one `NAME VALUE` entry per
            line. Required for small genomes, for which STAR rejects its
            default `genomeSAindexNbases` (14): STAR's options are passed
            through the `mkref_star_args` parameter as a comma-separated
            list, e.g.
            `mkref_star_args  --genomeSAindexNbases,8,--genomeSAsparseD,2`.
        <more>: Other environment variables required by `split-pipe --mode mkref`
            See <https://support.parsebiosciences.com/hc/en-us/articles/23060102930580>
            for more details.
    """  # noqa: E501
    input = "fasta:files, genes:files"
    output = """outdir:dir:
        {%- if envs.genome_name is string -%}
            {{- envs.genome_name -}}
        {%- else -%}
            {{- envs.genome_name | join: "_" -}}
        {%- endif -%}
    """
    output_flatten = True
    lang = config.lang.python
    envs = {
        "splitpipe": config.exe.splitpipe,
        "nthreads": config.misc.ncores,
        "ncores": None,
        "genome_name": None,
        "parfile": None,
    }
    script = "file://../scripts/splitpipe/SplitPipeRef.py"

class SplitPipeCount(Proc):
    """Run split-pipe on one sublibrary

    Runs `split-pipe --mode all` to pre-process, align and post-process the
    sequencing data of a single sublibrary, and to generate the expression
    matrices and the report for it.
    Requires the Parse Biosciences Pipeline (split-pipe) v1.9.0+.

    Terminology: a "sublibrary" is one well-pool of a kit processed by
    `--mode all` (i.e. the fastq files of one sample well-pool of one
    sequencing run). Sublibraries are always processed independently and must
    never be merged at this step: fastq files from multiple sequencing lanes
    of the SAME sublibrary are passed together as a list to `in.fq1`, while
    fastq files of DIFFERENT sublibraries belong in different jobs of this
    process. Sublibraries are merged into a single dataset afterwards, using
    `SplitPipeCombine` (`--mode comb`).

    Input:
        fq1: The input fastq (R1) file(s) of one sublibrary
            Multiple files are allowed when the sublibrary was sequenced over
            multiple lanes, and they are passed together to `--fq1`.
        fq2: The input fastq (R2) file(s) of the sublibrary
            Multiple files are allowed when the sublibrary was sequenced over
            multiple lanes, and they are passed together to `--fq2`. If not
            given, the R2 file names are derived from the R1 ones by
            split-pipe, which requires the R1 and R2 files to differ only by
            `R1`/`R2` in their names.
        id: The id defining the output directory. If not provided, it is
            inferred from the fq1 file names, e.g. `s1_S1_R1_001.fastq.gz`
            gives `s1`.

    Output:
        outdir: The output directory

    Envs:
        splitpipe: Path to the split-pipe executable
        nthreads: Number of threads to use
        ncores: Alias for `nthreads`. Takes precedence when both are given.
        chemistry: The chemistry version of the kit, e.g. `v1`, `v2`, `v3` or
            `v4`. Required by split-pipe.
        kit: The kit name, e.g. `WT`, `WT_mini`, `WT_mega`, `WT_mega_384`,
            `WT_penta`, `WT_penta_384`, `MAX_WT_S`, `MAX_WT_M`, `MAX_WT_L` or
            `MAX_WT_L_384`. Required by split-pipe.
        genome_dir: Path to the indexed genome directory built by
            `SplitPipeRef` (`--mode mkref`). It can also be passed by the
            input channel, see `SplitPipePipeline`.
        samples (type=list): The sample definitions, as a list of `NAME WELLS`
            entries, e.g. `["xcond_1 A1-A6", "xcond_2 C1-C6"]`. Each entry is
            passed to split-pipe as one `--sample NAME WELLS` argument.
            If omitted, the pipeline processes all samples (`all-sample`).
        samp_sltab: Path to the SampleLoadingTable (excel) file from the
            Support Suite.
        samp_slt: Alias of `samp_sltab` in Parse's own documentation, where
            both spellings are used for the SampleLoadingTable.
            Don't set both of them.
        samp_list: Path to the sample specification list (text) file, with one
            `NAME WELLS` entry per line.
        parfile: Path to the parameter file with one `NAME VALUE` entry per
            line, e.g. `post_min_map_frac 0.01`.
        targeted_list: Path to the target gene list (csv) for Gene Select data.
        crispr (flag): Whether this is a CRISPR Detect sublibrary. It must be
            used together with `envs.crsp_guides` and `envs.parent_dir`.
        crsp_guides: Path to the csv file with the guide RNA sequences for
            CRISPR Detect data.
        crsp_read_thresh (type=int): The `guide_read_cutoff` for CRISPR Detect
            data. The default is 3.
        crsp_tscp_thresh (type=int): The `guide_tscp_cutoff` for CRISPR Detect
            data. The default is 5.
        ffpe (flag): Whether this is an Evercode WT FFPE sublibrary.
            Only v4 chemistry is supported for FFPE data.
        tcr_analysis (flag): Whether to run TCR analysis. Immune databases must
            have been installed with the pipeline, see `split-pipe --immune_check`.
        bcr_analysis (flag): Whether to run BCR analysis.
        immune_genome: The immune genome to align to for TCR/BCR analysis.
            `human` and `mouse` are supported for TCR, and `human`, `mouse`
            and `transgenic_mouse` for BCR analysis.
        parent_dir: Path to the whole transcriptome (parent) output directory
            of the sublibrary. If provided, the parent results are used to
            filter the cells of the TCR/BCR/CRISPR analysis.
        use_imgt_db (flag): Whether to use the IMGT reference for the
            TCR/BCR analysis. The IMGT database must have been installed with
            `bash install_immune_dbs.sh -I -y` and followed by
            `pip install --no-cache-dir .` in the pipeline directory.
        min_map_fraction (type=float): The minimum mapping fraction of the
            reads to be retained as transcripts. Parse documents this
            parameter as `post_min_map_frac` in the `--parfile` parameter
            file, so a parameter file is generated for this run; it cannot be
            used together with `envs.parfile`.
        <more>: Other environment variables required by `split-pipe --mode all`
            See <https://support.parsebiosciences.com/hc/en-us/articles/23060220494228>
            for more details.
    """  # noqa: E501
    input = "fq1:files, fq2:files, id"
    output = """outdir:dir:
        {%- if in.id -%}
            {{in.id}}
        {%- else -%}
            {%- set id = in.fq1 | first | basename |
                regex_replace: "(_L\\d+)?(_R\\d+)?(_\\d+)?\\.f(ast)?q(\\.gz)?$", "" |
                regex_replace: "_S\\d+$", "" -%}
            {{- id -}}
        {%- endif -%}
    """
    output_flatten = True
    lang = config.lang.python
    envs = {
        "splitpipe": config.exe.splitpipe,
        "nthreads": config.misc.ncores,
        "ncores": None,
        "chemistry": None,
        "kit": None,
        "genome_dir": None,
        "samples": None,
        "samp_sltab": None,
        "samp_slt": None,
        "samp_list": None,
        "parfile": None,
        "targeted_list": None,
        "crispr": False,
        "crsp_guides": None,
        "crsp_read_thresh": None,
        "crsp_tscp_thresh": None,
        "ffpe": False,
        "tcr_analysis": False,
        "bcr_analysis": False,
        "immune_genome": None,
        "parent_dir": None,
        "use_imgt_db": False,
        "min_map_fraction": None,
    }
    script = "file://../scripts/splitpipe/SplitPipeCount.py"
    plugin_opts = {
        "report": "file://../reports/splitpipe/SplitPipeCount.svelte",
        "report_paging": 5,
    }


class SplitPipeCombine(Proc):
    """Combine sublibraries into a single dataset

    Runs `split-pipe --mode comb` to merge the outputs of the sublibrary
    (one well-pool processed by `--mode all`) directories of one sequencing
    run into a single dataset with a new report.
    Requires the Parse Biosciences Pipeline (split-pipe) v1.9.0+.

    Terminology: each sublibrary is one well-pool processed by
    `SplitPipeCount` (`--mode all`). Only sublibraries of the same library
    preparation can be combined (for example multiple WT sublibraries), and
    sublibraries must have been processed with the same version of the
    pipeline as the one used here. Chemistry, kit and sample definitions are
    parsed from the sublibrary directories, so they are not passed to
    `--mode comb`.

    Input:
        sublibs: The output directories of the sublibraries to combine, from
            `SplitPipeCount`.
        id: The id defining the output directory. If not provided, it is
            inferred from the sublibrary directories (see `SplitPipePipeline`).

    Output:
        outdir: The output directory

    Envs:
        splitpipe: Path to the split-pipe executable
        nthreads: Number of threads to use
        ncores: Alias for `nthreads`. Takes precedence when both are given.
        sublib_list: Path to a file with one sublibrary directory per line.
            If provided, it is passed as `--sublib_list` instead of passing
            `in.sublibs` as `--sublibraries`.
        immune_genome: The immune genome to align to for the TCR/BCR analysis
            combined here, either `human` or `mouse`.
        parent_dir: Path to the combined whole transcriptome (parent) output
            directory. The parent outputs must be combined before the
            corresponding TCR/BCR/CRISPR ones, and both must combine the same
            set of sublibraries.
        use_imgt_db (flag): Whether to use the IMGT reference for the
            TCR/BCR analysis. The IMGT database must have been installed with
            `bash install_immune_dbs.sh -I -y` and followed by
            `pip install --no-cache-dir .` in the pipeline directory.
        <more>: Other environment variables required by `split-pipe --mode comb`.
            Note that `--tcr_analysis`/`--bcr_analysis`/`--crispr` are not
            documented for `--mode comb`; any other option can be passed
            through here if the installed version of split-pipe supports it.
            See <https://support.parsebiosciences.com/hc/en-us/articles/23060220494228>
            for more details.
    """  # noqa: E501
    input = "sublibs:dirs, id"
    output = """outdir:dir:
        {%- if in.id -%}
            {{in.id}}
        {%- else -%}
            {{- in.sublibs | first | basename | append: "-comb" -}}
        {%- endif -%}
    """
    output_flatten = True
    lang = config.lang.python
    envs = {
        "splitpipe": config.exe.splitpipe,
        "nthreads": config.misc.ncores,
        "ncores": None,
        "sublib_list": None,
        "immune_genome": None,
        "parent_dir": None,
        "use_imgt_db": False,
    }
    script = "file://../scripts/splitpipe/SplitPipeCombine.py"
    plugin_opts = {
        "report": "file://../reports/splitpipe/SplitPipeCombine.svelte",
        "report_paging": 5,
    }


class SplitPipeSummary(Proc):
    """Summarize split-pipe metrics

    Merges the metrics files of the sublibrary (or combined) output
    directories of `SplitPipeCount`/`SplitPipeCombine` into one table. The
    names of the metrics files are version-dependent, so they are discovered
    with `envs.metrics_glob` and the process fails if nothing matches.

    Input:
        indirs: The directories containing split-pipe results, from
            `SplitPipeCount`/`SplitPipeCombine`.

    Output:
        outdir: The output directory

    Envs:
        metrics_glob (type=str): The glob pattern of the metrics file names.
            split-pipe 1.9.0 writes `agg_sample_summary.csv` at the top level
            of each output directory (the older `metrics_summary.csv` name is
            not written by this version),
            relative to each input directory, e.g. `*metrics_summary.csv`.
            The directories are searched recursively, so metrics files nested
            in per-sample subdirectories are found as well.
    """
    input = "indirs:dirs"
    input_data = lambda ch: [list(ch.iloc[:, 0])]
    output = "outdir:dir:{{in.indirs | first | stem | append: '-etc.splitpipe_summary'}}"
    lang = config.lang.rscript
    script = "file://../scripts/splitpipe/SplitPipeSummary.R"
    envs = {"metrics_glob": "*agg_sample_summary.csv"}
    plugin_opts = {
        "report": "file://../reports/splitpipe/SplitPipeSummary.svelte",
        "report_paging": 8,
    }
