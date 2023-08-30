# CNVkit pipeline

The pipeline decouples `cnvkit.py batch` so that we get detailed control over each step.

## Pipeline overview

The pipeline consists of the following steps:

1. `cnvkit.py access` to generate a BED file of accessible regions if not given
2. Guess baits from bam files if baitfile is not given
3. `cnvkit.py autobin` to generate target and antitarget files
4. `cnvkit.py coverage` to generate coverage files for target region
5. `cnvkit.py coverage` to generate coverage files for antitarget region
6. `cnvkit.py reference` to generate a reference.cnn file using normal samples (or a "flat" reference file if no normal samples are given)
7. `cnvkit.py fix` to combine the uncorrected target and antitarget coverage tables (.cnn) and correct for biases in regional coverage and GC content, according to the given reference.
8. `cnvkit.py segment` to infer discrete copy number segments from the given coverage table:
9. `cnvkit.py call` to call copy number alterations from the given segments file
10. `cnvkit.py scatter` to generate scatter plots of log2 ratios
11. `cnvkit.py diagram` to generate a diagram of copy number alterations on all chromosomes
12. `cnvkit.py heatmap` to generate a heatmap of segment-level log2 ratios
13. `cnvkit.py heatmap` to generate a heatmap of bin-level log2 ratios

See also the flowchart below:

![CNVkit pipeline](./CNVkit-pipeline.png)

## Input files

- `metafile`: a tab-separated file (see the next section) containing sample information
- `baitfile`: Potentially targeted genomic regions.
    E.g. all possible exons for the reference genome.
    This is optional when `method` is `wgs`.
- `accfile`: The accessible genomic regions.
    If not given, use `cnvkit.py access` to generate one. You can control the details by configuration items `[CNVkitAccess.envs]`

## Configurations

### Special configurations

- `access_excludes`: File(s) with regions to be excluded for
    `cnvkit.py access`.
- `guessbaits_guided`: Whether to use guided mode for guessing baits.
- `metacols`: The column names for each type of information in metafile
  - `group`: The column name in the metafile that indicates the sample group
      Default: `Group`
  - `purity`: The column name in the metafile that indicates the sample
      purity. Default: `Purity`
  - `snpvcf`: The column name in the metafile that indicates the path to
      the SNP VCF file. Default: `SnpVcf`
  - `bam`: The column name in the metafile that indicates the path to the
      BAM file. Default: `Bam`
  - `vcf_sample_id`: The column name in the metafile that indicates the
      sample ID in the VCF file. Default: `VcfSampleId`
  - `vcf_normal_id`: The column name in the metafile that indicates the
      normal sample ID in the VCF file. Default: `VcfNormalId`
  - `sex`: The column name in the metafile that indicates the sample
      sex. Default: `Sex`
  - `guess_baits`: The column name in the metafile that indicates whether
      to guess the bait file from the bam files. Default: `GuessBaits`
- `guessbaits`: Guess the bait file from the bam files, either guided or
    unguided.
    If `False`, `baitfile` is used. Otherwise, if `baitfile` is given, use it
    (guided), otherwise use `accfile` (unguided).
    The bam files with `metacols.guess_baits` column set to `True`, `TRUE`,
    `true`, `1`, `Yes`, `YES`, or `yes` will be used to guess the bait file.
- `case`: The group name of samples in `metacols.group` to call CNVs for.
    If not specified, use all samples. In such a case, `control` must not be
    specified, as we are using a flat reference.
- `control`: The group name of samples in `metacols.group` to use as reference
    if not specified, use a flat reference.

### Global configurations

The options that are used by multiple processes (can be overriden individually by `[<proc>.envs.xxx]`):

- `cnvkit`: the path to the `cnvkit.py` executable, defaults to
    `config.exe.cnvkit` from `./.biopipen.toml` or `~/.biopipen.toml`.
- `rscript`: Path to the Rscript excecutable to use for running R code.
    Requires `DNAcopy` to be installed in R, defaults to
    `config.lang.rscript`
- `samtools`: Path to samtools, used for guessing bait file.
- `convert`: Linux `convert` command to convert pdf to png
    So that they can be embedded in the HTML report.
- `ncores`: number of cores to use, defaults to `config.misc.ncores`
- `reffa`: the reference genome (e.g. hg19.fa), defaults to `config.ref.reffa`
    Used by `CNVkitAccess`, `CNVkitAutobin` and `CNVkitReference`
- `annotate`: Use gene models from this file to assign names to the
    target regions. Format: UCSC refFlat.txt or ensFlat.txt file
    (preferred), or BED, interval list, GFF, or similar. Defaults to `config.ref.refflat`
- `short_names`: Reduce multi-accession bait labels to be short and consistent
- `method`: Sequencing protocol: hybridization capture ('hybrid'),
    targeted amplicon sequencing ('amplicon'),
    or whole genome sequencing ('wgs'). Determines
    whether and how to use antitarget bins.
- `male_reference`: Use or assume a male reference (i.e. female samples
    will have +1 log-CNR of chrX; otherwise male samples would have
    -1 chrX).
    Used by `CNVkitReference`, `CNVkitCall`, `CNVkitHeatmapCns` and
    `CNVkitHeatmapCnr`.
- `drop_low_coverage`: Drop very-low-coverage bins before segmentation to
    avoid false-positive deletions in poor-quality tumor samples.
    Used by `CNVkitSegment` and `CNVkitCall`
- `no_gc`: Skip GC correction for `cnvkit.py reference/fix`.
- `no_edge`: Skip edge-effect correction for `cnvkit.py reference/fix`.
- `no_rmask`: Skip RepeatMasker correction for `cnvkit.py reference/fix`.
    no_* options are used by `CNVkitReference` and `CNVkitFix`
- `min_variant_depth`: Minimum read depth for a SNV to be displayed
    in the b-allele frequency plot.
    Used by `CNVkitSegment` and `CNVkitCall`
- `zygosity_freq`: Ignore VCF's genotypes (GT field) and instead infer
    zygosity from allele frequencies.
    Used by `CNVkitSegment` and `CNVkitCall`

### Process-specific configurations

The options that are used by a single process. See the process-specific documentation for details. You can configure them by `[<proc>.envs.xxx]` in the config file.

- [`CNVkitAccess`][1]
- [`CNVkitGuessBaits`][2]
- [`CNVkitAutobin`][3]
- [`CNVkitCoverageTarget`][4]
- [`CNVkitCoverageAntitarget`][4]
- [`CNVkitReference`][5]
- [`CNVkitFix`][6]
- [`CNVkitSegment`][7]
- [`CNVkitScatter`][8]
- [`CNVkitDiagram`][9]
- [`CNVkitHeatmapCns`][10]
- [`CNVkitHeatmapCnr`][10]
- [`CNVkitHeatmapCall`][11]

## The metafile

A metafile should be with the following columns:

- Sample: The sample_id used for target/antitarget files. If not provided,
    the sample_id will be the first part of basename of the bam file.
    For exapmle: `D123.tumor.bam -> D123`
- `<bam>`: The path to the bam file, better using absolute path.
- `<group>`: The type of the sample, defining the tumor/normal samples.
- `<sex>`: Guess each sample from coverage of X and Y chromosomes if
    not given.
- `<purity>`: Estimated tumor cell fraction, a.k.a. purity or cellularity.
- `<snpvcf>`: file name containing variants for segmentation by allele
    frequencies.
- `<vcf_sample_id>`: Sample ID in the VCF file.
- `<vcf_normal_id>`: Normal sample ID in the VCF file.
- `<guess_baits>`: Guess the bait file from the bam file

[1]: https://pwwang.github.io/biopipen/api/biopipen.ns.cnvkit/#biopipen.ns.cnvkit.CNVkitAccess
[2]: https://pwwang.github.io/biopipen/api/biopipen.ns.cnvkit/#biopipen.ns.cnvkit.CNVkitGuessBaits
[3]: https://pwwang.github.io/biopipen/api/biopipen.ns.cnvkit/#biopipen.ns.cnvkit.CNVkitAutobin
[4]: https://pwwang.github.io/biopipen/api/biopipen.ns.cnvkit/#biopipen.ns.cnvkit.CNVkitCoverage
[5]: https://pwwang.github.io/biopipen/api/biopipen.ns.cnvkit/#biopipen.ns.cnvkit.CNVkitReference
[6]: https://pwwang.github.io/biopipen/api/biopipen.ns.cnvkit/#biopipen.ns.cnvkit.CNVkitFix
[7]: https://pwwang.github.io/biopipen/api/biopipen.ns.cnvkit/#biopipen.ns.cnvkit.CNVkitSegment
[8]: https://pwwang.github.io/biopipen/api/biopipen.ns.cnvkit/#biopipen.ns.cnvkit.CNVkitScatter
[9]: https://pwwang.github.io/biopipen/api/biopipen.ns.cnvkit/#biopipen.ns.cnvkit.CNVkitDiagram
[10]: https://pwwang.github.io/biopipen/api/biopipen.ns.cnvkit/#biopipen.ns.cnvkit.CNVkitHeatmap
[11]: https://pwwang.github.io/biopipen/api/biopipen.ns.cnvkit/#biopipen.ns.cnvkit.CNVkitCall
