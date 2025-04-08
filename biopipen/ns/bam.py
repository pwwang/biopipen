"""Tools to process sam/bam/cram files"""

from ..core.proc import Proc
from ..core.config import config


# +-------------------------------------------------------------------+
# | CNV callers                                                       |
# +-------------------------------------------------------------------+
class CNVpytor(Proc):
    """Detect CNV using CNVpytor

    Input:
        bamfile: The bam file
            Will try to index it if it's not indexed.
        snpfile: The snp file

    Output:
        outdir: The output directory

    Envs:
        cnvpytor: Path to cnvpytor
        samtools: Path to samtools, used to index bam file in case it's not
        ncores: Number of cores to use (`-j` for cnvpytor)
        refdir: The directory containing the fasta file for each chromosome
        genome: The genome assembly to put in the VCF file
        chrsize: The geome size file to fix missing contigs in VCF header
        chrom: The chromosomes to run on
        binsizes: The binsizes
        snp: How to read snp data
        filters: The filters to filter the result
            See - https://github.com/abyzovlab/CNVpytor/blob/master/GettingStarted.md#predicting-cnv-regions
        mask_snps: Whether mask 1000 Genome snps
        baf_nomask: Do not use P mask in BAF histograms

    Requires:
        cnvpytor:
           - check: {{proc.envs.cnvpytor}} --version
    """  # noqa: E501
    input = "bamfile:file, snpfile:file"
    output = "outdir:dir:{{in.bamfile | stem}}.cnvpytor"
    lang = config.lang.python
    envs = {
        "cnvpytor": config.exe.cnvpytor,
        "samtools": config.exe.samtools,
        "ncores": config.misc.ncores,
        "refdir": config.ref.refdir,
        "genome": config.ref.genome,
        "chrsize": config.ref.chrsize,
        "chrom": [],
        "binsizes": [10000, 100000],
        # set False to disable snp data importing
        "snp": {
            "sample": "",
            "name1": [],
            "ad": "AD",
            "gt": "GT",
            "noAD": False,
        },
        "filters": {
            # https://github.com/abyzovlab/CNVpytor/blob/master/
            # GettingStarted.md#predicting-cnv-regions
            "CNVsize": [50000, "inf"],  # CNV size
            "eval1": [0, 0.0001],  # filter on e-val1
            "eval2": [],  # filter on e-val2
            "eval3": [],  # filter on e-val3
            "eval4": [],  # filter on e-val4
            "q0": [-1, 0.5],  # filter on Q0
            "pN": [0, 0.5],  # filter on pN
            "dG": [100000, "inf"],  # filter on dG
        },
        "mask_snps": True,
        "baf_nomask": False,
        # other arguments for -rd
    }
    script = "file://../scripts/bam/CNVpytor.py"
    plugin_opts = {"report": "file://../reports/bam/CNVpytor.svelte"}


class ControlFREEC(Proc):
    """Detect CNVs using Control-FREEC

    Input:
        bamfile: The bam file
        snpfile: The snp file

    Output:
        outdir: The output directory

    Envs:
        freec: Path to Control-FREEC executable
        ncores: Number of cores to use
        arggs: Other arguments for Control-FREEC
    """
    input = "bamfile:file, snpfile:file"
    output = "outdir:dir:{{in.bamfile | stem}}.freec"
    lang = config.lang.python
    envs = {
        "freec": config.exe.freec,
        "ncores": config.misc.ncores,
        "tabix": config.exe.tabix,
        "bedtools": config.exe.bedtools,
        "sambamba": config.exe.sambamba,
        "samtools": config.exe.samtools,
        # The "<ref>.fai" file will be used as chrLenFile
        "ref": config.ref.reffa,
        "refdir": config.ref.refdir,
        "rscript": config.lang.rscript,
        "binsize": 50_000,  # shortcut for args.general.window
        "args": {
            "general": {
                "ploidy": 2,
                "breakPointThreshold": 0.8,
            },
            "sample": {"mateOrientation": "FR"},
            "control": {},
            "BAF": {},
            "target": {},
        },
    }
    script = "file://../scripts/bam/ControlFREEC.py"
    plugin_opts = {"report": "file://../reports/bam/ControlFREEC.svelte"}


class CNAClinic(Proc):
    """Detect CNVs using CNAClinic

    Input:
        metafile: The meta file, header included, tab-delimited, including
            following columns:
            - Bam: The path to bam file
            - Sample: Optional. The sample names,
                if you don't want filename of bam file to be used
            - Group: Optional. The group names, either "Case" or "Control"
            - Patient: Optional. The patient names. Since CNAClinic only
                supports paired samples, you need to provide the patient names
                for each sample. Required if "Group" is provided.
            - Binsizer: Optional. Samples used to estimate the bin size
                "Y", "Yes", "T", "True", will be treated as True
                If not provided, will use `envs.binsizer` to get the samples
                to use. Either this column or `envs.binsizer` should be
                provided.

    Output:
        outdir: The output directory

    Envs:
        ncores: Number of cores to use
        seed: The seed for random number generator for choosing samples
            for estimating bin size
        binsizer: The samples used to estimate the bin size, it could be:
            A list of sample names
            A float number (0 < x <= 1), the fraction of samples to use
            A integer number (x > 1), the number of samples to use
        binsize: Directly use this binsize for CNAClinic, in bp.
        genome: The genome assembly
        run_args: The arguments for CNAClinic::runSegmentation
        plot_args: The arguments for CNAClinic::plotSampleData
        plot_multi_args: The arguments for CNAClinic::plotMultiSampleData
    """
    input = "metafile:file"
    output = "outdir:dir:{{in.metafile | stem}}.cnaclinic"
    lang = config.lang.rscript
    envs = {
        "ncores": config.misc.ncores,
        "binsizer": None,
        "binsize": None,
        "seed": 123,
        "genome": config.ref.genome,
        "run_args": {
            # HMM is errored
            "segmentType": ["CBS", "LACBS", "PLS"],
            "segmentsToSummarise": ["CBS", "LACBS", "PLS"],
            "summaryMethod": "mean",
        },
        "plot_args": {},
        "plot_multi_args": False,
    }
    script = "file://../scripts/bam/CNAClinic.R"
    plugin_opts = {
        "report": "file://../reports/bam/CNAClinic.svelte",
        "report_paging": 20,
    }


# +-------------------------------------------------------------------+
# | Bam processing tools                                              |
# +-------------------------------------------------------------------+
class BamSplitChroms(Proc):
    """Split bam file by chromosomes

    Input:
        bamfile: The bam file

    Output:
        outdir: The output directory with bam files for each chromosome

    Envs:
        ncores: Number of cores to use
        samtools: Path to samtools executable
        sambamba: Path to sambamba executable
        tool: The tool to use, either "samtools" or "sambamba"
        keep_other_sq: Keep other chromosomes in "@SQ" field in header
        chroms: The chromosomes to keep, if not provided, will use all
        index: Whether to index the output bam files. Requires the input bam
            file to be sorted.
    """
    input = "bamfile:file"
    output = "outdir:dir:{{in.bamfile | stem}}.split"
    lang = config.lang.python
    envs = {
        "ncores": config.misc.ncores,
        "samtools": config.exe.samtools,
        "sambamba": config.exe.sambamba,
        "tool": "samtools",
        "keep_other_sq": False,
        "chroms": [],
        "index": True,
    }
    script = "file://../scripts/bam/BamSplitChroms.py"


class BamMerge(Proc):
    """Merge bam files

    Input:
        bamfiles: The bam files

    Output:
        outfile: The output bam file

    Envs:
        ncores: Number of cores to use
        tool: The tool to use, either "samtools" or "sambamba"
        samtools: Path to samtools executable
        sambamba: Path to sambamba executable
        sort: Whether to sort the output bam file
        index: Whether to index the output bam file
            Requires envs.sort to be True
        merge_args: The arguments for merging bam files
            `samtools merge` or `sambamba merge`, depending on `tool`
            For `samtools`, these keys are not allowed: `-o`, `-O`,
            `--output-fmt`, `-@`, and `--threads`, as they are managed by
            the script
            For `sambamba`, these keys are not allowed: `-t`, and `--nthreads`,
            as they are managed by the script
        sort_args: The arguments for sorting bam files
            `samtools sort` or `sambamba sort`, depending on `tool`
            For `samtools`, these keys are not allowed: `-o`, `-@`,
            and `--threads`, as they are managed by the script
            For `sambamba`, these keys are not allowed: `-t`, `--nthreads`,
            `-o` and `--out`, as they are managed by the script
    """
    input = "bamfiles:files"
    output = "outfile:file:{{in.bamfiles | first | stem}}.merged.bam"
    lang = config.lang.python
    envs = {
        "ncores": config.misc.ncores,
        "samtools": config.exe.samtools,
        "sambamba": config.exe.sambamba,
        "tool": "samtools",
        "sort": True,
        "index": True,
        "merge_args": [],
        "sort_args": [],
    }
    script = "file://../scripts/bam/BamMerge.py"


class BamSampling(Proc):
    """Keeping only a fraction of read pairs from a bam file

    Input:
        bamfile: The bam file

    Output:
        outfile: The output bam file

    Envs:
        ncores: Number of cores to use
        samtools: Path to samtools executable
        tool: The tool to use, currently only "samtools" is supported
        fraction (type=float): The fraction of reads to keep.
            If `0 < fraction <= 1`, it's the fraction of reads to keep.
            If `fraction > 1`, it's the number of reads to keep.
            Note that when fraction > 1, you may not get the exact number
            of reads specified but a close number.
        seed: The seed for random number generator
        index: Whether to index the output bam file
        sort: Whether to sort the output bam file
        sort_args: The arguments for sorting bam file using `samtools sort`.
            These keys are not allowed: `-o`, `-@`,
            and `--threads`, as they are managed by the script.
    """
    input = "bamfile:file"
    output = "outfile:file:{{in.bamfile | stem}}.sampled{{envs.fraction}}.bam"
    lang = config.lang.python
    envs = {
        "ncores": config.misc.ncores,
        "samtools": config.exe.samtools,
        "tool": "samtools",
        "fraction": None,
        "seed": 8525,
        "index": True,
        "sort": True,
        "sort_args": [],
    }
    script = "file://../scripts/bam/BamSampling.py"


class BamSubsetByBed(Proc):
    """Subset bam file by the regions in a bed file

    Input:
        bamfile: The bam file
        bedfile: The bed file

    Output:
        outfile: The output bam file

    Envs:
        ncores: Number of cores to use
        samtools: Path to samtools executable
        tool: The tool to use, currently only "samtools" is supported
        index: Whether to index the output bam file
    """
    input = "bamfile:file, bedfile:file"
    output = "outfile:file:{{in.bamfile | stem}}-subset.bam"
    lang = config.lang.python
    envs = {
        "ncores": config.misc.ncores,
        "samtools": config.exe.samtools,
        "tool": "samtools",
        "index": True,
    }
    script = "file://../scripts/bam/BamSubsetByBed.py"


class BamSort(Proc):
    """Sort bam file

    Input:
        bamfile: The bam file

    Output:
        outfile: The output bam file

    Envs:
        tool (choice): The tool to use.
            - samtools: Use `samtools`
            - sambamba: Use `sambamba`
        ncores (type=int): Number of cores to use
        samtools: Path to samtools executable
        sambamba: Path to sambamba executable
        tmpdir: The temporary directory to use
        byname (flag): Whether to sort by read name
        index (flag): Whether to index the output bam file
            The index file will be created in the same directory as the output
            bam file
        <more>: Other arguments passed to the sorting tool
            See `samtools sort` or `sambamba sort`
    """
    input = "bamfile:file"
    output = "outfile:file:{{in.bamfile | stem}}.sorted.bam"
    lang = config.lang.python
    envs = {
        "tool": "samtools",
        "ncores": config.misc.ncores,
        "samtools": config.exe.samtools,
        "sambamba": config.exe.sambamba,
        "tmpdir": config.path.tmpdir,
        "byname": False,
        "index": True,
    }
    script = "file://../scripts/bam/BamSort.py"


class SamtoolsView(Proc):
    """View bam file using samtools, mostly used for filtering

    This is a wrapper for `samtools view` command.
    It will create a new bam file with the same name as the input bam file.

    Input:
        bamfile: The bam file

    Output:
        outfile: The output bam file

    Envs:
        ncores: Number of cores to use
        samtools: Path to samtools executable
        index: Whether to index the output bam file
            Requires the input bam file to be sorted.
        <more>: Other arguments passed to the view tool
            See `samtools view` or `sambamba view`.
    """
    input = "bamfile:file"
    output = "outfile:file:{{in.bamfile | stem}}.bam"
    lang = config.lang.python
    envs = {
        "ncores": config.misc.ncores,
        "samtools": config.exe.samtools,
        "index": True,
    }
    script = "file://../scripts/bam/SamtoolsView.py"
