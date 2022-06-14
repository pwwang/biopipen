# biopipen - A set of processes/pipelines for bioinformatics

## Installation

```shell
pip install -U biopipen
```

## Usage

### Use as APIs

```python
from pipen import Proc, Pipen
from biopipen.ns.bed import BedLiftOver

MyBedLiftOver = Proc.from_proc(BedLiftOver)

if __name__ == "__main__":
    Pipen().set_start(MyBedLiftOver).run()
```

### Use as pipen-cli-run plugin

```shell
❯ pipen run bed BedLiftOver

DESCRIPTION:
  Liftover a BED file using liftOver

USAGE:
  pipen [OPTIONS]

OPTIONS FOR <BedLiftOver>:
  --in.inbed <list>               - The input BED file Default: \[]
  --out.outbed <auto>             - The output BED file Default: <awaiting compiling>
  --envs.liftover <str>           - The path to liftOver Default: liftOver
  --envs.chain <str>              - The map chain file for liftover
                                    Default: ~/reference/hg38ToHg19.over.chain.gz

OPTIONAL OPTIONS:
  --config <path>                 - Read options from a configuration file in TOML. Default: None
  -h, --help                      - Print help information for this command
  --full                          - Show full options for this command

PIPELINE OPTIONS:
  --profile <str>                 - The default profile from the configuration to run the pipeline.
                                    This profile will be used unless a profile is specified in the
                                    process or in the .run method of pipen. Default: default
  --outdir <path>                 - The output directory of the pipeline
                                    Default: ~/bedliftover_pipeline_results/
  --workdir <str>                 - The workdir for the pipeline. Default: ./.pipen
  --scheduler <str>               - The scheduler to run the jobs. Default: local
```
