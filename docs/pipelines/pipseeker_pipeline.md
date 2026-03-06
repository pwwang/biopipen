# Pipseeker pipeline

## Pipeline overview

Each pipeline contains two processes.

- `PipseekerFull`: Run pipseeker full on each sample.
- `PipseekerSummary`: Summarize the results from each sample.

## Input files

Each sample should have a set of fastq files, in the format of:

```python
[
    # fastq files for sample 1
    [
        "sample1_R1.fastq.gz",
        "sample1_R2.fastq.gz",
    ],
    # fastq files for sample 2
    [
        "sample2_R1.fastq.gz",
        "sample2_R2.fastq.gz",
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
        "sample1_R1.fastq.gz",
        "sample1_R2.fastq.gz",
    ],
    # fastq files for sample 2
    [
        "sample2_R1.fastq.gz",
        "sample2_R2.fastq.gz",
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

- [`PipseekerFull`](https://pwwang.github.io/biopipen/api/biopipen.ns.pipseeker/#biopipen.ns.pipseeker.PipseekerFull)
- [`PipseekerSummary`](https://pwwang.github.io/biopipen/api/biopipen.ns.pipseeker/#biopipen.ns.pipseeker.PipseekerSummary)

### Reference

To run the pipeline, you need to provide the reference genome for the pipseeker pipeline. You can
provide the reference genome in the configuration:

```toml
[PipseekerFull.envs]
ref = "/path/to/reference"
```

To obtain the reference genome, please refer to:

- <https://www.fluentbio.com/resources/pipseeker-downloads/#-mapping-references>

You can also provide cell type annotation reference:

- <https://www.fluentbio.com/resources/pipseeker-downloads/#annotation-references>

## Docker image

You can use the docker image [`biopipen/pipseeker-pipeline`](https://hub.docker.com/r/biopipen/pipseeker-pipeline) to run the pipeline. The image contains
the pipseeker software and the biopipen package. It is also built with an example dataset for you to
test the pipeline:

```
/example/example-data/Sample1_Sub_R1.fastq.gz
/example/example-data/Sample1_Sub_R2.fastq.gz
/example/example-data/Sample2_Sub_R1.fastq.gz
/example/example-data/Sample2_Sub_R2.fastq.gz
```

A reference for test is also provided at `/example/example-data/STAR_test_index`.

A sample configuration file is also provided at `/biopipen/docker/pipseeker_pipeline/PipseekerPipeline.config.toml`.

Note that the docker image was not built with the reference genome. You need to provide the reference genome
by yourself.
