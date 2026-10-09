# splitpipe_pipeline

A docker image for the [`SplitPipePipeline`](../../docs/pipelines/splitpipe_pipeline.md):
biopipen + the Parse Biosciences Pipeline (`split-pipe`).

The image embeds the Parse Biosciences Pipeline, research-use-only software
licensed by Parse Biosciences.

## Build

Download the Parse Biosciences Pipeline package from Parse Biosciences
(`ParseBiosciences-Pipeline.<version>.zip`) and extract it into the biopipen
checkout. The build context is the checkout itself (the extracted package
directory inside it is the default `PARSE_SRC`).

The build pushes nothing by itself, but the recipe is wslc-first (the flags are
the same with plain `docker build`). Create the private repository
`biopipen/splitpipe` on Docker Hub (a private repository works — Docker
Personal includes one, paid tiers more) and log in:

```powershell
wslc registry login
```

Then build with the biopipen checkout as the context:

```powershell
wslc build -f C:\path\to\biopipen\docker\splitpipe_pipeline\Dockerfile `
    -t biopipen/splitpipe:1.9.0 `
    C:\path\to\biopipen
```

If the extracted package is not at the checkout root, point `PARSE_SRC` at it
inside the context:

```powershell
wslc build -f C:\path\to\biopipen\docker\splitpipe_pipeline\Dockerfile `
    --build-arg PARSE_SRC=path/to/ParseBiosciences-Pipeline.1.9.0 `
    -t biopipen/splitpipe:1.9.0 C:\path\to\biopipen
```

The recipe follows Parse's installation instructions: it creates a `spipe`
conda environment with Python 3.12.8, runs Parse's own
`install_dependencies_conda.sh`, and installs the pipeline into it. The biopipen
checkout (with its R dependencies) is copied to `/biopipen` and installed
editable into the base environment, and `split-pipe` is linked into `PATH`, so
the processes work with their default configuration. Since the image installs
this checkout, it embeds the current development version of biopipen, including
the `splitpipe` namespace (not released on PyPI yet).

## Push

```powershell
wslc push biopipen/splitpipe:1.9.0
```

## Run

Mount a working directory with a configuration file
(`SplitPipePipeline.config.toml`) and your data, and run:

```powershell
wslc run --rm -w /workdir `
    -v .:/workdir -v D:\data:/data -v D:\genomes:/genomes `
    biopipen/splitpipe:1.9.0 `
    SplitPipePipeline @SplitPipePipeline.config.toml
```

with `SplitPipePipeline.config.toml`:

```toml
[SplitPipePipeline]
# fastq (R1) files, one list per sublibrary (lanes of the same sublibrary
# belong in the same list)
input = [
    ["/data/s1_S1_L001_R1_001.fastq.gz", "/data/s1_S1_L002_R1_001.fastq.gz"],
    ["/data/s2_S1_L001_R1_001.fastq.gz"],
]
# optional: if not given, the R2 file names are derived from the R1 ones by
# split-pipe (they must differ only by `R1`/`R2`)
fq2 = [
    ["/data/s1_S1_L001_R2_001.fastq.gz", "/data/s1_S1_L002_R2_001.fastq.gz"],
    ["/data/s2_S1_L001_R2_001.fastq.gz"],
]
# use an already built genome ...
genome_dir = "/genomes/GRCh38"
# ... or build one from the reference files:
# fasta = ["/genomes/ref.fa"]
# genes = ["/genomes/ref.gtf"]
# genome_name = ["GRCh38"]

combine = true   # combine the sublibraries into a single dataset
summary = true   # summarize the metrics into one table

[SplitPipeCount.envs]
chemistry = "v3"
kit = "WT"
nthreads = 16
# `ncores` is an alias for `nthreads`, taking precedence when both are given:
# ncores = 16
```

Alternatively, run the pipeline with pipen-board:

```powershell
wslc run --rm -p 18521:18521 -w /workdir `
    -v .:/workdir -v D:\data:/data -v D:\genomes:/genomes `
    biopipen/splitpipe:1.9.0 `
    board biopipen.ns.splitpipe_pipeline:SplitPipePipeline -a /biopipen/docker/splitpipe_pipeline/board.toml
```

Everything (`split-pipe`, `pipen` and `biopipen`) is in the image, and
`split-pipe` is on `PATH`, so no extra configuration is needed. Run
`wslc run -it --rm -w /workdir -v .:/workdir biopipen/splitpipe:1.9.0 bash`
to get a shell in the image.

TCR/BCR analysis is not enabled by default. To enable it, uncomment the
`install_immune_dbs.sh` step in the `Dockerfile` (it downloads the immune
references, a few GB) and rebuild.

See the [SplitPipe pipeline docs](../../docs/pipelines/splitpipe_pipeline.md)
for all the configurations.
