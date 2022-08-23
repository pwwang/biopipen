"""CNVkit commnads"""

from ..core.proc import Proc
from ..core.config import config


class CNVkitAccess(Proc):
    """Run cnvkit access

    Input:
        excfiles: Additional regions to exclude, in BED format

    Output:
        outfile: The output file

    Envs:
        cnvkit: Path to cnvkit.py
        min_gap_size: Minimum gap size between accessible sequence regions
        ref: The reference genome fasta file

    Requires:
        - name: cnvkit
          check: |
            {{proc.envs.cnvkit}} version
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
    """Run cnvkit autobin

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
        cnvkit: Path to cnvkit.py
        method: Sequencing protocol: hybridization capture ('hybrid'),
            targeted amplicon sequencing ('amplicon'),
            or whole genome sequencing ('wgs'). Determines
            whether and how to use antitarget bins.
        bp_per_bin: Desired average number of sequencing read bases mapped to
            each bin.
        target_max_size: Maximum size of target bins.
        target_min_size: Minimum size of target bins.
        antitarget_max_size: Maximum size of antitarget bins.
        antitarget_min_size: Minimum size of antitarget bins.
        annotate: Use gene models from this file to assign names to the target
            regions. Format: UCSC refFlat.txt or ensFlat.txt file (preferred),
            or BED, interval list, GFF, or similar.
        short_names: Reduce multi-accession bait labels to be short and
            consistent.
        ref: The reference genome fasta file

    Requires:
        - name: cnvkit
          check: |
            {{proc.envs.cnvkit}} version
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
        count: Get read depths by counting read midpoints within each bin. (An
            alternative algorithm).
        min_mapq: Minimum mapping quality to include a read.
        ncores: Number of subprocesses to calculate coverage in parallel
        ref: The reference genome fasta file

    Requires:
        - name: cnvkit
          check: |
            {{proc.envs.cnvkit}} version
    """

    input = "bamfile:file, target_file:file"
    output = """
        {%- if "antitarget" in basename(in.target_file) -%}
            outfile:file:{{in.bamfile | stem0}}.antitargetcoverage.cnn
        {%- else -%}
            outfile:file:{{in.bamfile | stem0}}.targetcoverage.cnn
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
        cluster: Calculate and store summary stats for clustered subsets of
            the normal samples with similar coverage profiles.
        min_cluster_size: Minimum cluster size to keep in reference profiles.
        male_reference: Create a male reference: shift female samples
            chrX log-coverage by -1, so the reference chrX average is -1.
            Otherwise, shift male samples chrX by +1, so the reference
            chrX average is 0.
        no_gc: Skip GC correction.
        no_edge: Skip edge-effect correction.
        no_rmask: Skip RepeatMasker correction.
        ref: The reference genome fasta file

    Requires:
        - name: cnvkit
          check: |
            {{proc.envs.cnvkit}} version
    """

    input = (
        "covfiles:files, target_file:file, antitarget_file:file, sample_sex:var"
    )
    output = "outfile:file:{{in.covfiles | first | stem0 }}-etc.reference.cnn"
    lang = config.lang.python
    envs = {
        "cnvkit": config.exe.cnvkit,
        "cluster": False,
        "min_cluster_size": False,
        "sample_sex": False,
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
        cluster: Compare and use cluster-specific values present in the
            reference profile
            (requires `envs.cluster=True` for `CNVkitReference`).
        no_gc: Skip GC correction.
        no_edge: Skip edge-effect correction.
        no_rmask: Skip RepeatMasker correction.

    Requires:
        - name: cnvkit
          check: |
            {{proc.envs.cnvkit}} version
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
        drop_low_coverage: Drop very-low-coverage bins before segmentation to
            avoid false-positive deletions in poor-quality tumor samples.
        drop_outliers: Drop outlier bins more than this many multiples of
            the 95th quantile away from the average within a rolling window.
            Set to 0 for no outlier filtering.
        rscript: Path to Rscript
        ncores: Number of subprocesses to segment in parallel
            0 or negative for all available cores
        smooth_cbs: Perform an additional smoothing before CBS segmentation,
            which in some cases may increase the sensitivity. Used only for
            CBS method.
        min_variant_depth: Minimum read depth for a SNV to be displayed
            in the b-allele frequency plot.
        zygosity_freq: Ignore VCF's genotypes (GT field) and instead infer
            zygosity from allele frequencies.

    Requires:
        - name: cnvkit
          check: |
            {{proc.envs.cnvkit}} version
        - name: r-DNAcopy
          check: |
            {{proc.envs.rscript}} <(echo "library(DNAcopy)")
    """

    input = "cnrfile:file, vcf:file, sample_id:var, normal_id:var"
    output = "outfile:file:{{in.cnrfile | stem0}}.cns"
    lang = config.lang.python
    envs = {
        "cnvkit": config.exe.cnvkit,
        "method": "cbs",
        "threshold": False,
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
        convert_args: The arguments for `convert`
        chromosome: Chromosome or chromosomal range,
            e.g. 'chr1' or 'chr1:2333000-2444000', to display.
            If a range is given, all targeted genes in this range will be
            shown, unless -g/--gene is also given.
        gene: Name of gene or genes (comma-separated) to display.
        width: Width of margin to show around the selected gene(s) (-g/--gene)
            or small chromosomal region (-c/--chromosome).
        antitarget_marker: Plot antitargets using this symbol when plotting
            in a selected chromosomal region (-g/--gene or -c/--chromosome).
            Same as targets if not given
        by_bin: Plot data x-coordinates by bin indices instead of genomic
            coordinates. All bins will be shown with equal width, no blank
            regions will be shown, and x-axis values indicate bin number
            (within chromosome) instead of genomic position.
        segment_color: Plot segment lines in this color. Value can be
            any string accepted by matplotlib, e.g. 'red' or '#CC0000'.
        trend: Draw a smoothed local trendline on the scatter plot.
        y_max: y-axis upper limit.
        y_min: y-axis lower limit.
        min_variant_depth: Minimum read depth for a SNV to be displayed
            in the b-allele frequency plot.
        zygosity_freq: Ignore VCF's genotypes (GT field) and instead infer
            zygosity from allele frequencies.
        title: Plot title. Sample ID if not provided.
        cases: The cases for different plots with keys as case names and values
            to overwrite the default args given by `envs.<args>`, including -
            `convert_args`, `by_bin`, `chromosome`, `gene`, `width`
            `antitarget_marker`, `segment_color`, `trend`, `y_max`, `y_min`,
            `min_variant_depth`, `zygosity_freq` and `title.
            By default, an `all` case will be created with default arguments
            if no case specified

    Requires:
        - name: cnvkit
          check: |
            {{proc.envs.cnvkit}} version
        - name: convert
          check: |
            {{proc.envs.convert}} -version
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
        "chromosome": False,
        "gene": False,
        "width": 1000000,
        "antitarget_marker": False,
        "by_bin": False,
        "segment_color": False,
        "trend": False,
        "y_max": False,
        "y_min": False,
        "min_variant_depth": 20,
        "zygosity_freq": 0.25,
        "title": False,
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
        convert: Path to `convert` to convert PDF to png
        convert_args: The default arguments for convert
        threshold: Copy number change threshold to label genes.
        min_probes: Minimum number of covered probes to label a gene.
        male_reference: Assume inputs were normalized to a male reference
            (i.e. female samples will have +1 log-CNR of chrX; otherwise
            male samples would have -1 chrX).
        no_shift_xy: Don't adjust the X and Y chromosomes according
            to sample sex.
        title: Plot title. Sample ID if not provided.
        cases: The cases with keys as names and values as different configs,
            including `threshold`, `min_probes`, `male_reference`, `no_shift_xy`
            and `title`

    Requires:
        - name: cnvkit
          check: |
            {{proc.envs.cnvkit}} version
        - name: convert
          check: |
            {{proc.envs.convert}} -version
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
        "title": False,
        "cases": {},
    }
    script = "file://../scripts/cnvkit/CNVkitDiagram.py"
    plugin_opts = {
        "report": "file://../reports/cnvkit/CNVkitScatter.svelte",
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
        convert: Path to `convert` to convert PDF to png
        convert_args: The default arguments for convert
        by_bin: Plot data x-coordinates by bin indices instead of genomic
            coordinates. All bins will be shown with equal width,
            no blank regions will be shown, and x-axis values indicate
            bin number (within chromosome) instead of genomic position.
        chromosome: Chromosome (e.g. 'chr1') or chromosomal range
            (e.g. 'chr1:2333000-2444000') to display.
        desaturate: Tweak color saturation to focus on significant changes.
        male_reference: Assume inputs were normalized to a male reference
            (i.e. female samples will have +1 log-CNR of chrX; otherwise
            male samples would have -1 chrX).
        no_shift_xy: Don't adjust the X and Y chromosomes according to
            sample sex.
        cases: The cases for different plots with keys as case names and values
            to overwrite the default args given by `envs.<args>`, including -
            `convert_args`, `by_bin`, `chromosome`, `desaturate`,
            `male_reference`, and, `no_shift_xy`.
            By default, an `all` case will be created with default arguments
            if no case specified

    Requires:
        - name: cnvkit
          check: |
            {{proc.envs.cnvkit}} version
        - name: convert
          check: |
            {{proc.envs.convert}} -version
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

    Output:
        outdir: The output directory including the call file (.call.cns)
            bed file, and the vcf file

    Envs:
        cnvkit: Path to cnvkit.py
        center: Re-center the log2 ratio values using this estimator of
            the center or average value.
        center_at: Subtract a constant number from all log2 ratios.
            For "manual" re-centering, in case the --center option gives
            unsatisfactory results.)
        filter: Merge segments flagged by the specified filter(s) with
            the adjacent segment(s).
        method: Calling method (threshold, clonal or none).
        thresholds: Hard thresholds for calling each integer copy number,
            separated by commas.
        ploidy: Ploidy of the sample cells.
        purity: Estimated tumor cell fraction, a.k.a. purity or cellularity.
        drop_low_coverage: Drop very-low-coverage bins before segmentation
            to avoid false-positive deletions in poor-quality tumor samples.
        male_reference: Assume inputs were normalized to a male reference
            (i.e. female samples will have +1 log-CNR of chrX; otherwise
            male samples would have -1 chrX).
        min_variant_depth: Minimum read depth for a SNV to be displayed
            in the b-allele frequency plot.
        zygosity_freq: Ignore VCF's genotypes (GT field) and instead infer
            zygosity from allele frequencies.

    Requires:
        - name: cnvkit
          check: |
            {{proc.envs.cnvkit}} version
    """
    input = (
        "cnrfile:file, cnsfile:file, "
        "vcf:file, sample_id:var, normal_id:var, sample_sex:var"
    )
    output = "outdir:dir:{{in.cnsfile | stem0}}.cnvkit"
    lang = config.lang.python
    envs = {
        "cnvkit": config.exe.cnvkit,
        "center": "median",
        "center_at": False,
        "filter": False,
        "method": "threshold",
        "thresholds": "-1.1,-0.25,0.2,0.7",
        "ploidy": 2,
        "purity": False,
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
        - name: cnvkit
          check: |
            {{proc.envs.cnvkit}} version
        - name: r-DNAcopy
          check: |
            {{proc.envs.rscript}} <(echo "library(DNAcopy)")
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
