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

## Configurations

- `input`: The input fastq files for each sample.

