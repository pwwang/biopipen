<p align="center">
  <img height="120" style="height: 120px" src="https://github.com/pwwang/biopipen/blob/dev/docs/img/logo.png?raw=true">
</p>
<p align="center">
  A set of processes/pipelines for bioinformatics based on
  <a href="https://github.com/pwwang/pipen" target="_blank">pipen</a>
</p>
<hr />

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
❯ pipen run bed BedLiftOver --help
Usage: pipen [-h | -h+] [options]

Liftover a BED file using liftOver
Use `@configfile` to load default values for the options.

Pipeline Options:
  --name NAME           The name for the pipeline, will affect the default workdir and
                        outdir. [default: BedLiftOver_pipeline]
  --profile PROFILE     The default profile from the configuration to run the pipeline.
                        This profile will be used unless a profile is specified in the
                        process or in the .run method of pipen. You can check the available
                        profiles by running `pipen profile`
  --outdir OUTDIR       The output directory of the pipeline [default: ./<name>_results]
  --forks FORKS         How many jobs to run simultaneously by the scheduler
  --scheduler SCHEDULER
                        The scheduler to run the jobs

Namespace <envs>:
  --envs ENVS           Environment variables for the process [default: {'liftover':
                        'liftOver', 'chain': ''}]
  --envs.liftover LIFTOVER
                        The path to liftOver [default: liftOver]
  --envs.chain CHAIN    The map chain file for liftover [default: ]

Namespace <in>:
  --in.inbed INBED [INBED ...]
                        The input BED file

Namespace <out>:
  --out.outbed OUTBED   The output BED file [default: {{in.inbed | basename}}]

Options:
  -h, --help, -h+, --help+
                        show help message (with + to show more options) and exit
```
