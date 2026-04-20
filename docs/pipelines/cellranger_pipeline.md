# CellRanger pipeline

Including three pipelines: `CellRangerCountPipeline`, `CellRangerVdjPipeline`, and `CellRangerMultiPipeline`.

## Pipeline overview

Each pipeline contains two processes.

- `CellRangerCount`/`CellRangerVdj`/`CellRangerMulti`: Run cellranger count/vdj/multi on each sample or GEM well.
- `CellRangerSummary`/`CellRangerMultiSummary`: Summarize the results from each sample.

## Input files

Each sample should have a set of fastq files, in the format of:

```python
[
    # fastq files for sample 1
    [
        "sample1_S1_L001_R1_001.fastq.gz",
        "sample1_S1_L001_R2_001.fastq.gz",
        "sample1_S1_L002_R1_001.fastq.gz",
        "sample1_S1_L002_R2_001.fastq.gz",
    ],
    # fastq files for sample 2
    [
        "sample2_S1_L001_R1_001.fastq.gz",
        "sample2_S1_L001_R2_001.fastq.gz",
        "sample2_S1_L002_R1_001.fastq.gz",
        "sample2_S1_L002_R2_001.fastq.gz",
    ],
    ...
]
```

If the ids cannot be inferred from the fastq file names, or you want to use a different
id than the inferred one, you can specify the ids in the input:

```python
input = [
    # fastq files for sample 1
    [
        "sample1_S1_L001_R1_001.fastq.gz",
        "sample1_S1_L001_R2_001.fastq.gz",
        "sample1_S1_L002_R1_001.fastq.gz",
        "sample1_S1_L002_R2_001.fastq.gz",
    ],
    # fastq files for sample 2
    [
        "sample2_S1_L001_R1_001.fastq.gz",
        "sample2_S1_L001_R2_001.fastq.gz",
        "sample2_S1_L002_R1_001.fastq.gz",
        "sample2_S1_L002_R2_001.fastq.gz",
    ],
    ...
]
ids = [
    "sampleA",
    "sampleB",
    ...
]
```

## Configurations

- `input`: The input fastq files for each sample/GEM well. See `Input files` for details.
- `ids`: The ids for each sample/GEM well. See `Input files` for details.

Other than the input, you should provide other configurations to the processes to each individual
process. Check the documentation of each process for more details.

- [`CellRangerCount`](https://pwwang.github.io/biopipen/api/biopipen.ns.cellranger/#biopipen.ns.cellranger.CellRangerCount)
- [`CellRangerVdj`](https://pwwang.github.io/biopipen/api/biopipen.ns.cellranger/#biopipen.ns.cellranger.CellRangerVdj)
- [`CellRangerMulti`](https://pwwang.github.io/biopipen/api/biopipen.ns.cellranger/#biopipen.ns.cellranger.CellRangerMulti)
- [`CellRangerSummary`](https://pwwang.github.io/biopipen/api/biopipen.ns.cellranger/#biopipen.ns.cellranger.CellRangerSummary)
- [`CellRangerMultiSummary`](https://pwwang.github.io/biopipen/api/biopipen.ns.cellranger/#biopipen.ns.cellranger.CellRangerMultiSummary)

### CellRangerMultiPipeline library configuration

`CellRangerMultiPipeline` generates the `cellranger multi` config CSV internally from the `envs`
of the `CellRangerMulti` process. Each element of `input` corresponds to one GEM well (one
`cellranger multi` run). You must specify the library metadata in `envs.libraries`:

```toml
[CellRangerMultiPipeline]
input = [
    [
        "sample_A_S1_L001_R1_001.fastq.gz",
        "sample_A_S1_L001_R2_001.fastq.gz",
        "sample_B_S1_L001_R1_001.fastq.gz",
        "sample_B_S1_L001_R2_001.fastq.gz",
    ],
]
ids = ["GEMwell1"]

[CellRangerMulti.envs.gex]
reference = "/path/to/refdata-gex-GRCh38"
create_bam = false

# Optional: VDJ libraries
# [CellRangerMulti.envs.vdj]
# reference = "/path/to/refdata-cellranger-vdj-GRCh38"

# Optional: Feature Barcode libraries (antibody capture, CRISPR)
# [CellRangerMulti.envs.feature]
# reference = "/path/to/feature_ref.csv"

[[CellRangerMulti.envs.libraries]]
fastq_id = "sample_A"
feature_types = "Gene Expression"

[[CellRangerMulti.envs.libraries]]
fastq_id = "sample_B"
feature_types = "Gene Expression"
```

Each library entry must have at minimum `fastq_id` (matching the FASTQ filename prefix) and
`feature_types` (e.g. `Gene Expression`, `VDJ-T`, `VDJ-B`, `Antibody Capture`, `CRISPR Guide Capture`).
Optionally, `fastqs` (path to the directory containing the FASTQs) and `lanes` can be specified.

### Reference

To run the pipeline, you need to provide the reference genome for the cellranger pipeline. You can
provide the reference genome in the configuration:

```toml
[CellRangerCount.envs]
ref = "/path/to/reference"
```

For the multi pipeline, the reference and library configuration is specified under `CellRangerMulti.envs`:

```toml
[CellRangerMulti.envs.gex]
reference = "/path/to/refdata-gex-GRCh38"

# Optional: VDJ reference
# [CellRangerMulti.envs.vdj]
# reference = "/path/to/refdata-cellranger-vdj-GRCh38"

[[CellRangerMulti.envs.libraries]]
fastq_id = "Sample_A"
feature_types = "Gene Expression"
```

To obtain the reference genome, please refer to:

- <https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/inputs/cr-inputs-overview#count>
- <https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/inputs/cr-inputs-overview#vdj>

You may also make your own reference by `cellranger mkref` for gene expression. See the cellranger documentation:

- <https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/inputs/cr-3p-references>

Also check out docker/cellranger_pipeline/make-examples.sh to see how the references are prepared.

## Docker image

You can use the docker image [`biopipen/cellranger-pipeline`](https://hub.docker.com/r/biopipen/cellranger-pipeline) to run the pipeline. The image contains
the cellranger software and the biopipen package. It is also built with an example dataset for you to
test the pipeline:

```
/example/example-data/Sample_X_S1_L001_R1_001.fastq.gz
/example/example-data/Sample_X_S1_L001_R2_001.fastq.gz
/example/example-data/Sample_Y_S1_L001_R1_001.fastq.gz
/example/example-data/Sample_Y_S1_L001_R2_001.fastq.gz
/example/example-data/Sample_Y_S1_L002_R1_001.fastq.gz
/example/example-data/Sample_Y_S1_L002_R2_001.fastq.gz
```

A sample configuration file is also provided at `/biopipen/docker/cellranger_pipeline/CellrangerCountPipeline.config.toml` for the count pipeline and `/biopipen/docker/cellranger_pipeline/CellrangerMultiPipeline.config.toml` for the multi pipeline.

Note that the docker image was not built with the reference genome. You need to provide the reference genome
by yourself.

Also note that from biopipen v0.34.26, the docker image has been updated to use CellRanger 10.0.0. If you need to use CellRanger 9.0.1, you can use the previous docker image with the tag `<0.34.26`, e.g., `biopipen/cellranger-pipeline:0.34.20`.
