"""CNVkit commnads"""

from ..core.proc import Proc
from ..core.config import config


class CNVkitAccess(Proc):
    """Calculate the sequence-accessible coordinates in chromosomes from the
    given reference genome using `cnvkit.py access`

    Input:
        excfiles: Additional regions to exclude, in BED format

    Output:
        outfile: The output file

    Envs:
        cnvkit: Path to `cnvkit.py`
        min_gap_size (type=int): Minimum gap size between accessible sequence
            regions
        ref: The reference genome fasta file

    Requires:
        cnvkit:
            - check: {{proc.envs.cnvkit}} version
    """
    input = "excfiles:files"
    output = (
        "outfile:file:{{envs.ref | stem0}}.access.{{envs.min_gap_size}}.bed"
    )
    lang = config.lang.python
    envs = {
        "cnvkit": config.exe.cnvkit,
        "min_gap_size": 5000,
        "ref": config.ref.reffa,
    }
    script = "file://../scripts/cnvkit/CNVkitAccess.py"


class CNVkitAutobin(Proc):
    """Quickly estimate read counts or depths in a BAM file to estimate
    reasonable on- and (if relevant) off-target bin sizes.

    Using `cnvkit.py autobin`.

    If multiple BAMs are given, use the BAM with median file size.

    Input:
        bamfiles: The bamfiles
        accfile: The access file
        baitfile: Potentially targeted genomic regions.
            E.g. all possible exons for the reference genome.
            Format - BED, interval list, etc.

    Output:
        target_file: The target BED output
        antitarget_file: The antitarget BED output

    Envs:
        cnvkit: Path to `cnvkit.py`
        method (choice): Sequencing protocol. Determines whether and how to use
            antitarget bins.
            - hybrid: Hybridization capture
            - amplicon: Targeted amplicon sequencing
            - wgs: Whole genome sequencing
        bp_per_bin (type=int): Desired average number of sequencing read bases
            mapped to each bin.
        target_max_size (type=int): Maximum size of target bins.
        target_min_size (type=int): Minimum size of target bins.
        antitarget_max_size (type=int): Maximum size of antitarget bins.
        antitarget_min_size (type=int): Minimum size of antitarget bins.
        annotate: Use gene models from this file to assign names to the target
            regions. Format: UCSC refFlat.txt or ensFlat.txt file (preferred),
            or BED, interval list, GFF, or similar.
        short_names (flag): Reduce multi-accession bait labels to
            be short and consistent.
        ref: The reference genome fasta file

    Requires:
        cnvkit:
            - check: {{proc.envs.cnvkit}} version
    """
    input = "bamfiles:files, accfile:file, baitfile:file"
    output = [
        "target_file:file:{{in.bamfiles | first | stem0}}-etc.target.bed",
        "antitarget_file:file:{{in.bamfiles | first | stem0}}"
        "-etc.antitarget.bed",
    ]
    lang = config.lang.python
    envs = {
        "cnvkit": config.exe.cnvkit,
        "method": "hybrid",
        "bp_per_bin": 100000,
        "target_max_size": 20000,
        "target_min_size": 20,
        "antitarget_max_size": 500000,
        "antitarget_min_size": 500,
        "annotate": None,
        "short_names": False,
        "ref": config.ref.reffa,
    }
    script = "file://../scripts/cnvkit/CNVkitAutobin.py"


class CNVkitCoverage(Proc):
    """Run cnvkit coverage

    Input:
        bamfile: The bamfile
        target_file: The target file or anti-target file

    Output:
        outfile: The output coverage file

    Envs:
        cnvkit: Path to cnvkit.py
        count (flag): Get read depths by counting read midpoints
            within each bin. (An alternative algorithm).
        min_mapq (type=int): Minimum mapping quality to include a read.
        ncores (type=int): Number of subprocesses to calculate coverage
            in parallel
        ref: The reference genome fasta file

    Requires:
        cnvkit:
            - check: {{proc.envs.cnvkit}} version
    """
    input = "bamfile:file, target_file:file"
    output = """outfile:file:
        {%- if "antitarget" in basename(in.target_file) -%}
            {{in.bamfile | stem0}}.antitargetcoverage.cnn
        {%- else -%}
            {{in.bamfile | stem0}}.targetcoverage.cnn
        {%- endif -%}
    """
    lang = config.lang.python
    envs = {
        "cnvkit": config.exe.cnvkit,
        "count": False,
        "min_mapq": 0,
        "ncores": config.misc.ncores,
        "ref": config.ref.reffa,
    }
    script = "file://../scripts/cnvkit/CNVkitCoverage.py"


class CNVkitReference(Proc):
    """Run cnvkit reference

    To genearte a reference file from normal samples, provide the cnn coverage
    files from the normal samples. To generate a flat reference file, provide
    the target/antitarget file.

    Input:
        covfiles: The coverage files from normal samples
        target_file: Target intervals (.bed or .list)
        antitarget_file: Antitarget intervals (.bed or .list)
        sample_sex: Specify the chromosomal sex of all given samples as male or
            female. Guess each sample from coverage of X and Y chromosomes if
            not given.

    Output:
        outfile: The reference cnn file

    Envs:
        cnvkit: Path to cnvkit.py
        cluster (flag): Calculate and store summary stats for
            clustered subsets of the normal samples with similar coverage
            profiles.
        min_cluster_size (type=int): Minimum cluster size to keep in reference
            profiles.
        male_reference (flag): Create a male reference: shift
            female samples chrX log-coverage by -1, so the reference chrX
            average is -1. Otherwise, shift male samples chrX by +1, so the
            reference chrX average is 0.
        no_gc (flag): Skip GC correction.
        no_edge (flag): Skip edge-effect correction.
        no_rmask (flag): Skip RepeatMasker correction.
        ref: The reference genome fasta file

    Requires:
        cnvkit:
            - check: {{proc.envs.cnvkit}} version
    """
    input = [
        "covfiles:files",
        "target_file:file",
        "antitarget_file:file",
        "sample_sex:var",
    ]
    output = """outfile:file:
        {%- if not in.covfiles -%}
            flat.reference.cnn
        {%- else -%}
            {{in.covfiles | first | stem0 }}-etc.reference.cnn
        {%- endif -%}
    """
    lang = config.lang.python
    envs = {
        "cnvkit": config.exe.cnvkit,
        "cluster": False,
        "min_cluster_size": 4,
        "male_reference": False,
        "no_gc": False,
        "no_edge": False,
        "no_rmask": False,
        "ref": config.ref.reffa,
    }
    script = "file://../scripts/cnvkit/CNVkitReference.py"


class CNVkitFix(Proc):
    """Run cnvkit.py fix

    Input:
        target_file: The target file
        antitarget_file: The antitarget file
        reference: The refence cnn file
        sample_id: Sample ID for target/antitarget files.
            Otherwise inferred from file names.

    Output:
        outfile: The fixed coverage files (.cnr)

    Envs:
        cnvkit: Path to cnvkit.py
        cluster (flag): Compare and use cluster-specific values
            present in the reference profile.
            (requires `envs.cluster=True` for `CNVkitReference`).
        no_gc (flag): Skip GC correction.
        no_edge (flag): Skip edge-effect correction.
        no_rmask (flag): Skip RepeatMasker correction.

    Requires:
        cnvkit:
            - check: {{proc.envs.cnvkit}} version
    """
    input = (
        "target_file:file, antitarget_file:file, reference:file, sample_id:var"
    )
    output = (
        "outfile:file:{{in.sample_id | default: stem0(in.target_file)}}.cnr"
    )
    lang = config.lang.python
    envs = {
        "cnvkit": config.exe.cnvkit,
        "cluster": False,
        "no_gc": False,
        "no_edge": False,
        "no_rmask": False,
    }
    script = "file://../scripts/cnvkit/CNVkitFix.py"


class CNVkitSegment(Proc):
    """Run cnvkit.py segment

    For segmentation methods, see
    https://cnvkit.readthedocs.io/en/stable/pipeline.html#segmentation-methods

    Input:
        cnrfile: The fixed coverage files (.cnr)
        vcf: VCF file name containing variants for segmentation
            by allele frequencies (optional).
        sample_id: Specify the name of the sample in the VCF to use for b-allele
            frequency extraction and as the default plot title.
        normal_id: Corresponding normal sample ID in the input VCF.
            This sample is used to select only germline SNVs to plot
            b-allele frequencies.

    Output:
        outfile: The segmentation file (.cns)

    Envs:
        cnvkit: Path to cnvkit.py
        method: Method to use for segmentation.
            Candidates - cbs, flasso, haar, none, hmm, hmm-tumor, hmm-germline
        threshold: Significance threshold (p-value or FDR, depending on method)
            to accept breakpoints during segmentation. For HMM methods,
            this is the smoothing window size.
        drop_low_coverage (flag): Drop very-low-coverage bins
            before segmentation to avoid false-positive deletions in
            poor-quality tumor samples.
        drop_outliers (type=int): Drop outlier bins more than this many
            multiples of the 95th quantile away from the average within a
            rolling window. Set to 0 for no outlier filtering.
        rscript: Path to Rscript
        ncores (type=int): Number of subprocesses to segment in parallel.
            0 or negative for all available cores
        smooth_cbs (flag): Perform an additional smoothing before
            CBS segmentation, which in some cases may increase the sensitivity.
            Used only for CBS method.
        min_variant_depth (type=int): Minimum read depth for a SNV to be
            displayed in the b-allele frequency plot.
        zygosity_freq (type=float): Ignore VCF's genotypes (GT field) and
            instead infer zygosity from allele frequencies.

    Requires:
        cnvkit:
            - check: {{proc.envs.cnvkit}} version
        r-DNAcopy:
            - check: {{proc.envs.rscript}} <(echo "library(DNAcopy)")
    """
    input = "cnrfile:file, vcf:file, sample_id:var, normal_id:var"
    output = "outfile:file:{{in.cnrfile | stem0}}.cns"
    lang = config.lang.python
    envs = {
        "cnvkit": config.exe.cnvkit,
        "method": "cbs",
        "threshold": None,
        "drop_low_coverage": False,
        "drop_outliers": 10,
        "rscript": config.lang.rscript,
        "ncores": config.misc.ncores,
        "smooth_cbs": False,
        "min_variant_depth": 20,
        "zygosity_freq": 0.25,
    }
    script = "file://../scripts/cnvkit/CNVkitSegment.py"


class CNVkitScatter(Proc):
    """Run cnvkit.py scatter

    Input:
        cnrfile: The fixed cnr file (.cnr)
        cnsfile: The segmentation file (.cns)
        vcf: VCF file name containing variants for segmentation
            by allele frequencies (optional).
        sample_id: Specify the name of the sample in the VCF to use for b-allele
            frequency extraction and as the default plot title.
        normal_id: Corresponding normal sample ID in the input VCF.
            This sample is used to select only germline SNVs to plot
            b-allele frequencies.

    Output:
        outdir: Output directory with plots for multiple cases

    Envs:
        cnvkit: Path to cnvkit.py
        convert: Path to `convert` to convert pdf to png file
        convert_args (ns): The arguments for `convert`
            - density (type=int): Horizontal and vertical density of the image
            - quality (type=int): JPEG/MIFF/PNG compression level
            - background: Background color
            - alpha: Activate, deactivate, reset, or set the alpha channel
            - <more>: See `convert -help` and also:
                https://linux.die.net/man/1/convert
        chromosome: Chromosome or chromosomal range,
            e.g. 'chr1' or 'chr1:2333000-2444000', to display.
            If a range is given, all targeted genes in this range will be
            shown, unless -g/--gene is also given.
        gene: Name of gene or genes (comma-separated) to display.
        width (type=int): Width of margin to show around the selected gene(s)
            (-g/--gene) or small chromosomal region (-c/--chromosome).
        antitarget_marker (flag): Plot antitargets using this
            symbol when plotting in a selected chromosomal region
            (-g/--gene or -c/--chromosome).
        by_bin (flag): Plot data x-coordinates by bin indices
            instead of genomic coordinates. All bins will be shown with equal
            width, no blank regions will be shown, and x-axis values indicate
            bin number (within chromosome) instead of genomic position.
        segment_color: Plot segment lines in this color. Value can be
            any string accepted by matplotlib, e.g. 'red' or '#CC0000'.
        trend (flag): Draw a smoothed local trendline on the
            scatter plot.
        y_max (type=int): y-axis upper limit.
        y_min (tyoe=int): y-axis lower limit.
        min_variant_depth (type=int): Minimum read depth for a SNV to be
            displayed in the b-allele frequency plot.
        zygosity_freq (typ=float): Ignore VCF's genotypes (GT field) and
            instead infer zygosity from allele frequencies.
        title: Plot title. Sample ID if not provided.
        cases (type=json): The cases for different plots with keys as case names
            and values to overwrite the default args given by `envs.<args>`,
            including  `convert_args`, `by_bin`, `chromosome`, `gene`, `width`
            `antitarget_marker`, `segment_color`, `trend`, `y_max`, `y_min`,
            `min_variant_depth`, `zygosity_freq` and `title.
            By default, an `all` case will be created with default arguments
            if no case specified

    Requires:
        cnvkit:
            - check: {{proc.envs.cnvkit}} version
        convert:
            - check: {{proc.envs.convert}} -version
    """
    input = (
        "cnrfile:file, cnsfile:file, config:var, "
        "vcf:file, sample_id:var, normal_id:var"
    )
    output = "outdir:dir:{{in.cnrfile | stem0}}.scatter"
    lang = config.lang.python
    envs = {
        "cnvkit": config.exe.cnvkit,
        "convert": config.exe.convert,
        "convert_args": {
            "density": 150,
            "quality": 90,
            "background": "white",
            "alpha": "remove",
        },
        "chromosome": None,
        "gene": None,
        "width": 1000000,
        "antitarget_marker": False,
        "by_bin": False,
        "segment_color": None,
        "trend": False,
        "y_max": None,
        "y_min": None,
        "min_variant_depth": 20,
        "zygosity_freq": 0.25,
        "title": None,
        "cases": {},
    }
    script = "file://../scripts/cnvkit/CNVkitScatter.py"
    plugin_opts = {
        "report": "file://../reports/cnvkit/CNVkitScatter.svelte",
        "report_paging": 10,
    }


class CNVkitDiagram(Proc):
    """Run cnvkit.py diagram

    Input:
        cnrfile: The fixed cnr file (.cnr)
        cnsfile: The segmentation file (.cns)
        sample_sex: Specify the sample's chromosomal sex as male or female.
            (Otherwise guessed from X and Y coverage).

    Output:
        outdir: Output directory with the scatter plots

    Envs:
        cnvkit: Path to cnvkit.py
        convert: Path to `convert` to convert pdf to png file
        convert_args (ns): The arguments for `convert`
            - density (type=int): Horizontal and vertical density of the image
            - quality (type=int): JPEG/MIFF/PNG compression level
            - background: Background color
            - alpha: Activate, deactivate, reset, or set the alpha channel
            - <more>: See `convert -help` and also:
                https://linux.die.net/man/1/convert
        threshold (type=float): Copy number change threshold to label genes.
        min_probes (type=int): Minimum number of covered probes to label a gene.
        male_reference (flag): Assume inputs were normalized to a
            male reference (i.e. female samples will have +1 log-CNR of chrX;
            otherwise male samples would have -1 chrX).
        no_shift_xy (flag): Don't adjust the X and Y chromosomes
            according to sample sex.
        title: Plot title. Sample ID if not provided.
        cases (type=json): The cases with keys as names and values as different
            configs, including `threshold`, `min_probes`, `male_reference`,
            `no_shift_xy` and `title`

    Requires:
        cnvkit:
            - check: {{proc.envs.cnvkit}} version
        convert:
            - check: {{proc.envs.convert}} -version
    """
    input = "cnrfile:file, cnsfile:file, sample_sex:var"
    output = "outdir:dir:{{in.cnrfile | stem0}}.diagram"
    lang = config.lang.python
    envs = {
        "cnvkit": config.exe.cnvkit,
        "convert": config.exe.convert,
        "convert_args": {
            "density": 150,
            "quality": 90,
            "background": "white",
            "alpha": "remove",
        },
        "threshold": 0.5,
        "min_probes": 3,
        "male_reference": False,
        "no_shift_xy": False,
        "title": None,
        "cases": {},
    }
    script = "file://../scripts/cnvkit/CNVkitDiagram.py"
    plugin_opts = {
        "report": "file://../reports/cnvkit/CNVkitDiagram.svelte",
        "report_paging": 10,
    }


class CNVkitHeatmap(Proc):
    """Run cnvkit.py heatmap for multiple cases

    Input:
        segfiles: Sample coverages as raw probes (.cnr) or segments (.cns).
        sample_sex: Specify the chromosomal sex of all given samples as male
            or female. Separated by comma. (Default: guess each sample from
            coverage of X and Y chromosomes).

    Output:
        outdir: Output directory with heatmaps of multiple cases

    Envs:
        cnvkit: Path to cnvkit.py
        convert: Path to `convert` to convert pdf to png file
        convert_args (ns): The arguments for `convert`
            - density (type=int): Horizontal and vertical density of the image
            - quality (type=int): JPEG/MIFF/PNG compression level
            - background: Background color
            - alpha: Activate, deactivate, reset, or set the alpha channel
            - <more>: See `convert -help` and also:
                https://linux.die.net/man/1/convert
        by_bin (flag): Plot data x-coordinates by bin indices
            instead of genomic coordinates. All bins will be shown with equal
            width, no blank regions will be shown, and x-axis values indicate
            bin number (within chromosome) instead of genomic position.
        chromosome: Chromosome (e.g. 'chr1') or chromosomal range
            (e.g. 'chr1:2333000-2444000') to display.
        desaturate (flag): Tweak color saturation to focus on
            significant changes.
        male_reference (flag): Assume inputs were normalized to
            a male reference. (i.e. female samples will have +1 log-CNR of chrX;
            otherwise male samples would have -1 chrX).
        no_shift_xy (flag): Don't adjust the X and Y chromosomes
            according to sample sex.
        order: A file with sample names in the desired order.
        cases (type=json): The cases for different plots with keys as case names
            and values to overwrite the default args given by `envs.<args>`,
            including `convert_args`, `by_bin`, `chromosome`, `desaturate`,
            `male_reference`, and, `no_shift_xy`.
            By default, an `all` case will be created with default arguments
            if no case specified

    Requires:
        cnvkit:
            - check: {{proc.envs.cnvkit}} version
        convert:
            - check: {{proc.envs.convert}} -version
    """
    input = "segfiles:files, sample_sex: var"
    output = "outdir:dir:{{in.segfiles | first | stem0}}-etc.heatmap"
    lang = config.lang.python
    envs = {
        "cnvkit": config.exe.cnvkit,
        "convert": config.exe.config,
        "convert_args": {
            "density": 150,
            "quality": 90,
            "background": "white",
            "alpha": "remove",
        },
        "by_bin": False,
        "chromosome": False,
        "desaturate": False,
        "male_reference": False,
        "no_shift_xy": False,
        "order": None,
        "cases": {},
    }
    script = "file://../scripts/cnvkit/CNVkitHeatmap.py"
    plugin_opts = {"report": "file://../reports/cnvkit/CNVkitHeatmap.svelte"}


class CNVkitCall(Proc):
    """Run cnvkit.py call

    Input:
        cnrfile: The fixed cnr file (.cnr), used to generate VCF file
        cnsfile: The segmentation file (.cns)
        vcf: VCF file name containing variants for segmentation
            by allele frequencies (optional).
        sample_id: Specify the name of the sample in the VCF to use for b-allele
            frequency extraction and as the default plot title.
        normal_id: Corresponding normal sample ID in the input VCF.
            This sample is used to select only germline SNVs to plot
            b-allele frequencies.
        sample_sex: Specify the sample's chromosomal sex as male or female.
            (Otherwise guessed from X and Y coverage).
        purity: Estimated tumor cell fraction, a.k.a. purity or cellularity.

    Output:
        outdir: The output directory including the call file (.call.cns)
            bed file, and the vcf file

    Envs:
        cnvkit: Path to cnvkit.py
        center: Re-center the log2 ratio values using this estimator of
            the center or average value.
        center_at (type=float): Subtract a constant number from all log2 ratios.
            For "manual" re-centering, in case the --center option gives
            unsatisfactory results.)
        filter: Merge segments flagged by the specified
            filter(s) with the adjacent segment(s).
        method (choice): Calling method (threshold, clonal or none).
            - threshold: Using hard thresholds for calling each integer copy
                number.
                Use `thresholds` to set a list of threshold log2 values for
                each copy number state
            - clonal: Rescaling and rounding.
                For a given known tumor cell fraction and normal ploidy,
                then simple rounding to the nearest integer copy number
            - none: Do not add a “cn” column or allele copy numbers.
                But still performs rescaling, re-centering, and extracting
                b-allele frequencies from a VCF (if requested).
        thresholds: Hard thresholds for calling each integer copy number,
            separated by commas.
        ploidy (type=float): Ploidy of the sample cells.
        drop_low_coverage (flag): Drop very-low-coverage bins
            before segmentation to avoid false-positive deletions in
            poor-quality tumor samples.
        male_reference (flag): Assume inputs were normalized to a
            male reference.
            (i.e. female samples will have +1 log-CNR of chrX; otherwise
            male samples would have -1 chrX).
        min_variant_depth (type=int): Minimum read depth for a SNV to be
            displayed in the b-allele frequency plot.
        zygosity_freq (type=float): Ignore VCF's genotypes (GT field) and
            instead infer zygosity from allele frequencies.

    Requires:
        cnvkit:
            - check: {{proc.envs.cnvkit}} version
    """
    input = [
        "cnrfile:file",
        "cnsfile:file",
        "vcf:file",
        "sample_id:var",
        "normal_id:var",
        "sample_sex:var",
        "purity:var",
    ]
    output = "outdir:dir:{{in.cnsfile | stem0}}.cnvkit"
    lang = config.lang.python
    envs = {
        "cnvkit": config.exe.cnvkit,
        "center": "median",
        "center_at": None,
        "filter": None,
        "method": "threshold",
        "thresholds": "-1.1,-0.25,0.2,0.7",
        "ploidy": 2,
        "drop_low_coverage": False,
        "male_reference": False,
        "min_variant_depth": 20,
        "zygosity_freq": 0.25,
    }
    script = "file://../scripts/cnvkit/CNVkitCall.py"


class CNVkitBatch(Proc):
    """Run cnvkit batch

    If you need in-depth control of the parameters, for example, multiple
    scatter plots in different regions, or you need to specify sample-sex for
    different samples, take a look at `biopipen.ns.cnvkit_pipeline`

    Input:
        metafile: The meta data file containing the sample information
            Two columns BamFile and `envs.type_col` are required.
            The tumor samples should be labeled as `envs.type_tumor` and the
            normal samples should be labeled as `envs.type_normal` in the
            `envs.type_col` column. If normal samples are not found, a
            flat reference will be used.
            The could be other columns in the meta file, but they could be
            used in `biopipen.ns.cnvkit_pipeline`.

    Output:
        outdir: The output directory

    Envs:
        cnvkit: Path to cnvkit.py
        method: Sequencing assay type: hybridization capture ('hybrid'),
            targeted amplicon sequencing ('amplicon'), or whole genome
            sequencing ('wgs'). Determines whether and how to use antitarget
            bins.
        segment_method: cbs,flasso,haar,none,hmm,hmm-tumor,hmm-germline
            Method used in the 'segment' step.
        male_reference: Use or assume a male reference (i.e. female samples
            will have +1 log-CNR of chrX; otherwise male samples would have
            -1 chrX).
        count_reads: Get read depths by counting read midpoints within each bin.
            (An alternative algorithm).
        drop_low_coverage: Drop very-low-coverage bins before segmentation to
            avoid false-positive deletions in poor-quality tumor samples.
        ncores: Number of subprocesses used to running each of the BAM files
            in parallel
        rscript: Path to the Rscript excecutable to use for running R code.
            Use this option to specify a non-default R installation.
        ref: Path to a FASTA file containing the reference genome.
        targets: Target intervals (.bed or .list) (optional for wgs)
        antitargets: Anti-target intervals (.bed or .list) (optional for wgs)
        annotate: Use gene models from this file to assign names to the
            target regions. Format: UCSC refFlat.txt or ensFlat.txt file
            (preferred), or BED, interval list, GFF, or similar.
        short_names: Reduce multi-accession bait labels to be short
            and consistent.
        target_avg_size: Average size of split target bins
            (results are approximate).
        access: Regions of accessible sequence on chromosomes (.bed),
            as output by the 'access' command.
        access_min_gap_size: Minimum gap size between accessible
            sequence regions if `envs.access` is not specified.
        access_excludes: Exclude these regions from the accessible genome
            Used when `envs.access` is not specified.
        antitarget_avg_size: Average size of antitarget bins
            (results are approximate).
        antitarget_min_size: Minimum size of antitarget bins
            (smaller regions are dropped).
        cluster: Calculate and use cluster-specific summary stats in the
            reference pool to normalize samples.
        reference: Copy number reference file (.cnn) to reuse
        scatter: Create a whole-genome copy ratio profile as a PDF scatter plot.
        diagram: Create an ideogram of copy ratios on chromosomes as a PDF.
        type_col: type_col: The column name in the metafile that
            indicates the sample type.
        type_tumor: The type of tumor samples in `envs.type_col` column of
            `in.metafile`
        type_normal: The type of normal samples in `envs.type_col` column of
            `in.metafile`

    Requires:
        cnvkit:
            - check: {{proc.envs.cnvkit}} version
        r-DNAcopy:
            - check: {{proc.envs.rscript}} <(echo "library(DNAcopy)")
    """
    input = "metafile:file"
    output = "outdir:dir:{{in.metafile | stem0}}.cnvkit"
    lang = config.lang.python
    envs = {
        "cnvkit": config.exe.cnvkit,
        "method": "hybrid",
        "segment_method": "cbs",
        "male_reference": False,
        "count_reads": False,
        "drop_low_coverage": False,
        "ncores": config.misc.ncores,
        "rscript": config.lang.rscript,
        "ref": config.ref.reffa,
        "targets": False,
        "antitargets": False,
        "annotate": False,
        "short_names": False,
        "target_avg_size": False,
        "access": False,
        "access_min_gap_size": 5000,
        "access_excludes": False,
        "antitarget_avg_size": False,
        "antitarget_min_size": False,
        "cluster": False,
        "reference": False,
        "scatter": True,
        "diagram": True,
        "type_tumor": "Tumor",
        "type_normal": "Normal",
        "type_col": "SampleType",
    }
    script = "file://../scripts/cnvkit/CNVkitBatch.py"


class CNVkitGuessBaits(Proc):
    """Guess the bait intervals from the bam files

    It runs scripts/guess_baits.py from the cnvkit repo.

    Input:
        bamfiles: The bam files
        atfile: The potential target file or access file
            e.g. all known exons in the reference genome or
            from `cnvkit.py access`

    Output:
        targetfile: The target file

    Envs:
        cnvkit: Path to cnvkit.py
        guided (flag): `in.atfile` is a potential target file when
            `True`, otherwise it is an access file.
        samtools: Path to samtools executable
        ncores (type=int): Number of subprocesses to segment in parallel
            `0` to use the maximum number of available CPUs.
        ref: Path to a FASTA file containing the reference genome.
        min_depth (type=int): Minimum sequencing read depth to accept as
            captured. For guided only.
        min_gap (type=int): Merge regions separated by gaps smaller than this.
        min_length (type=int): Minimum region length to accept as captured.
            `min_gap` and `min_length` are for unguided only.
    """
    input = "bamfiles:files, atfile:file"
    output = "targetfile:file:{{in.bamfiles | first | stem}}_etc.baits.bed"
    lang = config.lang.python
    envs = {
        "cnvkit": config.exe.cnvkit,
        "samtools": config.exe.samtools,
        "ncores": config.misc.ncores,
        "ref": config.ref.reffa,
        "guided": None,
        "min_depth": 5,
        "min_gap": 25,
        "min_length": 50,
    }
    script = "file://../scripts/cnvkit/CNVkitGuessBaits.py"
