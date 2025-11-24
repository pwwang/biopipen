# CellRanger pipeline

Including two pipelines: `CellRangerCountPipeline` and `CellRangerVdjPipeline`.

## Pipeline overview

Each pipeline contains two processes.

- `CellRangerCount`/`CellRangerVdj`: Run cellranger count/vdj on each sample.
- `CellRangerSummary`: Summarize the results from each sample.

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

- `input`: The input fastq files for each sample. See `Input files` for details.
- `ids`: The ids for each sample. See `Input files` for details.

Other than the input, you should provide other configurations to the processes to each individual
process. Check the documentation of each process for more details.

- [`CellRangerCount`](https://pwwang.github.io/biopipen/api/biopipen.ns.cellranger/#biopipen.ns.cellranger.CellRangerCount)
- [`CellRangerVdj`](https://pwwang.github.io/biopipen/api/biopipen.ns.cellranger/#biopipen.ns.cellranger.CellRangerVdj)
- [`CellRangerSummary`](https://pwwang.github.io/biopipen/api/biopipen.ns.cellranger/#biopipen.ns.cellranger.CellRangerSummary)

### Reference

To run the pipeline, you need to provide the reference genome for the cellranger pipeline. You can
provide the reference genome in the configuration:

```toml
[CellRangerCount.envs]
ref = "/path/to/reference"
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

A sample configuration file is also provided at `/biopipen/docker/cellranger_pipeline/CellrangerCountPipeline.config.toml`.

Note that the docker image was not built with the reference genome. You need to provide the reference genome
by yourself.

Also note that from biopipen v0.34.26, the docker image has been updated to use CellRanger 10.0.0. If you need to use CellRanger 9.0.1, you can use the previous docker image with the tag `<0.34.26`, e.g., `biopipen/cellranger-pipeline:0.34.20`.
