# SplitPipe pipeline

Including the `SplitPipePipeline` and the processes to run
[Parse Biosciences Pipeline (`split-pipe`)](https://support.parsebiosciences.com/hc/en-us/articles/23060220494228)
on Evercode split-pool data.

The Parse Biosciences Pipeline is **proprietary software**. It is
**not** distributed with biopipen, is **not** contained in any published docker
image, and must not be redistributed. You need to obtain it from Parse
Biosciences and install it yourself (see [Installing split-pipe](#installing-split-pipe)).
The processes here only *execute* a `split-pipe` binary that you install, in
the same way other biopipen processes execute their wrapped tools.

## Pipeline overview

- `SplitPipeRef`: Build an indexed reference genome (`split-pipe --mode mkref`).
- `SplitPipeCount`: Run `split-pipe --mode all` on **one sublibrary**, one job per sublibrary.
- `SplitPipeCombine`: Combine the sublibraries of one sequencing run into a single dataset
    (`split-pipe --mode comb`).
- `SplitPipeSummary`: Merge the metrics of the sublibrary (or combined) output directories
    into one table.

`SplitPipePipeline` runs `SplitPipeRef` (unless `genome_dir` is given), then
`SplitPipeCount` for each sublibrary, then `SplitPipeCombine` and `SplitPipeSummary`.
`SplitPipeCount`, `SplitPipeCombine` and `SplitPipeSummary` can also be run individually
with `pipen run splitpipe <Process>`.

## Sublibraries and combining

A **sublibrary** is one well-pool of one sequencing run (the reads of one sample
well-pool). Each sublibrary is processed independently by `SplitPipeCount`
(`--mode all`): fastq files of multiple sequencing *lanes* of the same sublibrary
are passed together in the same list of `fq1`, while fastq files of different
sublibraries belong in different lists.

Sublibraries with the same chemistry/kit are merged into one dataset afterwards
with `SplitPipeCombine` (`--mode comb`), which generates a combined report for
the whole run. `SplitPipePipeline` does this automatically with `combine = true`
(the default); set it to `false` to only process the sublibraries.

## Input files

Each sublibrary should have a set of fastq (R1) files, passed as `fq1`, in the format of:

```python
[
    # fastq (R1) files of sublibrary 1 (multiple lanes)
    [
        "s1_S1_L001_R1_001.fastq.gz",
        "s1_S1_L002_R1_001.fastq.gz",
    ],
    # fastq (R1) files of sublibrary 2
    [
        "s2_S1_L001_R1_001.fastq.gz",
    ],
    ...
]
```

The ids of the sublibraries are inferred from the fastq file names this way
(`s1_S1_R1_001.fastq.gz` -> `s1`). If they cannot be inferred, or you want to
use different ids, you can specify them:

```python
fq1 = [
    ["s1_S1_L001_R1_001.fastq.gz", "s1_S1_L002_R1_001.fastq.gz"],
    ["s2_S1_L001_R1_001.fastq.gz"],
]
ids = ["sublibrary1", "sublibrary2"]
```

The R2 files are passed in the same format with `fq2`:

```python
fq2 = [
    # fastq (R2) files of sublibrary 1 (multiple lanes)
    [
        "s1_S1_L001_R2_001.fastq.gz",
        "s1_S1_L002_R2_001.fastq.gz",
    ],
    # fastq (R2) files of sublibrary 2
    [
        "s2_S1_L001_R2_001.fastq.gz",
    ],
    ...
]
```

If `fq2` is not given, the R2 file names are derived from the R1 ones by
split-pipe, which requires the R1 and R2 files of a sublibrary to differ only
by `R1`/`R2` in their names.

## Configurations

- `fq1`: The fastq (R1) files for each sublibrary. See [Input files](#input-files).
- `fq2`: The fastq (R2) files for each sublibrary, in the same format as `fq1`.
    Optional: if not given, the R2 file names are derived from the R1 ones by
    split-pipe.
- `ids`: The ids for each sublibrary. See [Input files](#input-files).
- `fasta`, `genes`, `genome_name`: The reference sequence/annotation and the genome name(s),
    used to build the genome with `SplitPipeRef`.
- `genome_dir`: The directory of an already built genome. If given, `SplitPipeRef`
    is not run and this directory is passed to `SplitPipeCount` directly.
    If not given, the genome is built by `SplitPipeRef` from `fasta`/`genes`/`genome_name`.
- `combine`: Whether to combine the sublibraries into a single dataset with
    `SplitPipeCombine`.
- `summary`: Whether to summarize the metrics with `SplitPipeSummary`.

Other than the input, you should provide other configurations to each individual
process, e.g. the chemistry and kit of the sublibrary:

```toml
[SplitPipePipeline]
fq1 = [
    ["/data/s1_S1_L001_R1_001.fastq.gz", "/data/s1_S1_L002_R1_001.fastq.gz"],
    ["/data/s2_S1_L001_R1_001.fastq.gz"],
]
fq2 = [
    ["/data/s1_S1_L001_R2_001.fastq.gz", "/data/s1_S1_L002_R2_001.fastq.gz"],
    ["/data/s2_S1_L001_R2_001.fastq.gz"],
]
fasta = ["/genomes/ref.fa"]
genes = ["/genomes/ref.gtf"]
genome_name = ["GRCh38"]

[SplitPipeCount.envs]
chemistry = "v3"
kit = "WT"
nthreads = 16
```

Check the documentation of each process for all the environment variables:

- [`SplitPipeRef`](https://pwwang.github.io/biopipen/api/biopipen.ns.splitpipe/#biopipen.ns.splitpipe.SplitPipeRef)
- [`SplitPipeCount`](https://pwwang.github.io/biopipen/api/biopipen.ns.splitpipe/#biopipen.ns.splitpipe.SplitPipeCount)
- [`SplitPipeCombine`](https://pwwang.github.io/biopipen/api/biopipen.ns.splitpipe/#biopipen.ns.splitpipe.SplitPipeCombine)
- [`SplitPipeSummary`](https://pwwang.github.io/biopipen/api/biopipen.ns.splitpipe/#biopipen.ns.splitpipe.SplitPipeSummary)

### Reference genome

You can either build a genome with `SplitPipeRef` (used automatically by
`SplitPipePipeline` when `genome_dir` is not given):

```toml
[SplitPipeRef.envs]
genome_name = "GRCh38"
```

or pass the directory of a genome already built by `split-pipe --mode mkref`:

```toml
[SplitPipePipeline]
genome_dir = "/genomes/GRCh38"
```

For a species-mix reference, pass one genome name per fasta/genes pair,
e.g. `genome_name = ["hg38", "GRCm39"]`.

## Installing split-pipe

The Parse Biosciences Pipeline is installed from a zip package downloaded from
Parse Biosciences (see the [Pipeline Download](https://support.parsebiosciences.com/hc/en-us/articles/17200056667924)
page). The instructions below follow Parse's installation guide; check it for
the details and for the hardware requirements:

- <https://support.parsebiosciences.com/hc/en-us/articles/23060102930580>

```bash
# Create a new environment with Python 3.12.8
conda create -n spipe conda-forge::python==3.12.8
conda activate spipe

# Unzip the package and install the dependencies
unzip ParseBiosciences-Pipeline.1.9.0.zip -d /opt
cd /opt/ParseBiosciences-Pipeline.1.9.0
bash ./install_dependencies_conda.sh -i -y

# Install the pipeline itself
pip install --no-cache-dir ./

# Check the installation
split-pipe --help
```

To install the TCR/BCR dependencies (and optionally the IMGT reference), run
`bash ./install_immune_dbs.sh -i -y` (`-I -y` for IMGT) in the package
directory, followed by `pip install --no-cache-dir .` again, and pass
`envs.tcr_analysis`/`envs.bcr_analysis` and `envs.use_imgt_db` to
`SplitPipeCount`/`SplitPipeCombine`.

### Locating the binary

The processes run the binary given by their `splitpipe` environment variable,
which defaults to `split-pipe` (resolved from `PATH`). Point it to the
`split-pipe` of your environment, either in your configuration file
(`~/.biopipen.toml`):

```toml
[exe]
splitpipe = "/opt/conda/envs/spipe/bin/split-pipe"
```

or per run:

```bash
pipen run splitpipe SplitPipeCount --envs.splitpipe /opt/conda/envs/spipe/bin/split-pipe ...
```

Note that the pipeline version requirements of split-pipe apply: all
sublibraries combined together must have been processed with the same version of
split-pipe as the one running the combination.

## Docker image

There is no published docker image for this pipeline, since the Parse
Biosciences Pipeline is proprietary and must not be redistributed. Instead,
`docker/splitpipe_pipeline/` contains a `Dockerfile` recipe that you build
locally from **your own** licensed download of the package. The recipe does not
download the Parse software; the extracted package directory is the build
context. See `docker/splitpipe_pipeline/README.md` for details.

The image built from this recipe contains the Parse Biosciences Pipeline and is
for research use only. Do not push it to any public registry and do not
redistribute it.

## Tests

The tests under `tests/test_splitpipe/` are not run in CI (`LOCAL_ONLY`), since
they require `split-pipe` and your own data. Put a small reference
(`ref.fa`, `ref.gtf`) and small sublibraries (`s1_R1.fastq.gz`, `s2_R1.fastq.gz`,
... see the docstrings of the tests) in `tests/test_splitpipe/data/` and run
e.g.:

```bash
python tests/test_splitpipe/SplitPipeCount/test.py
```
