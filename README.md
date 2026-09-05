<p align="center">
  <img width="240px" src="https://github.com/pwwang/biopipen/blob/dev/docs/img/logo.png?raw=true" />
</p>
<p align="center">
  A set of processes/pipelines for bioinformatics based on
  <a href="https://github.com/pwwang/pipen" target="_blank">pipen</a>
</p>

<p align="center">
  <a href="https://pypi.org/project/biopipen" target="_blank"><img alt="PyPI" src="https://img.shields.io/pypi/v/biopipen.svg"></a>
  <a href="https://pypi.org/project/biopipen" target="_blank"><img alt="Python" src="https://img.shields.io/pypi/pyversions/biopipen.svg"></a>
  <a href="https://pypi.org/project/biopipen" target="_blank"><img alt="License" src="https://img.shields.io/pypi/l/biopipen.svg"></a>
  <a href="https://github.com/pwwang/biopipen/actions/workflows/ci.yml" target="_blank"><img alt="CI" src="https://github.com/pwwang/biopipen/actions/workflows/ci.yml/badge.svg"></a>
  <a href="https://pwwang.github.io/biopipen" target="_blank"><img alt="Docs" src="https://img.shields.io/badge/docs-pwwang.github.io%2Fbiopipen-blue"></a>
</p>

<hr />

`biopipen` is a collection of **170+ ready-made bioinformatics processes and
pipelines**, built on top of the [pipen](https://github.com/pwwang/pipen) workflow
framework. Each process wraps a well-known bioinformatics tool (Seurat, CNVkit,
bcftools, immunarch, ...) into a reusable, composable unit with typed input/output
options, environment handling, and optional interactive reports powered by
[pipen-board](https://github.com/pwwang/pipen-board).

Use processes individually as Python classes, chain them into custom pipelines,
or run any of the built-in pipelines directly from the command line.

## Installation

```shell
pip install -U biopipen
```

Optional extras:

```shell
pip install -U "biopipen[runinfo]"     # write run information to JSON files
pip install -U "biopipen[log2file]"    # log to files
```

Requires Python >= 3.9. Process-level dependencies (e.g. Seurat, CNVkit) are
usually installed per-process in conda environments (see the
[documentation](https://pwwang.github.io/biopipen)).

## Usage

### As Python APIs

Import any process, customize it with `Proc.from_proc()`, and run it with pipen:

```python
from pipen import Proc, Pipen
from biopipen.ns.bed import BedLiftOver

MyBedLiftOver = Proc.from_proc(BedLiftOver)

if __name__ == "__main__":
    Pipen().set_start(MyBedLiftOver).run()
```

Processes are also composable — chain multiple processes together as a
workflow by declaring dependencies with `requires`:

```python
import pipen
from biopipen.ns.vcf import VcfLiftOver, VcfIndex

# Index the VCF lifted over by VcfLiftOver
class MyVcfIndex(VcfIndex):
    requires = VcfLiftOver

if __name__ == "__main__":
    pipen.run(
        name="liftover_pipeline",
        starts=[VcfLiftOver],
        data=[{"invcf": "input.vcf"}],
    )
```

The output of `VcfLiftOver` is automatically routed into `MyVcfIndex` as its
input. Any chain of `requires` builds a full dependency graph that pipen
schedules and runs.

### From the command line

`biopipen` registers each of its namespaces as a `pipen run` plugin. Every
process is runnable with full argument parsing (use `--help` to explore options):

```shell
pipen run bed BedLiftOver --in.inbed input.bed --envs.chain hg19ToHg38.over.chain
```

Available namespaces and their processes:

| Namespace | Description |
| --------- | ----------- |
| `bam` | BAM file processing: sorting, merging, sampling, splitting, subsetting; CNV detection (CNVpytor, ControlFREEC, CNAClinic); Samplot visualization |
| `bed` | BED operations: liftover, merge, intersect, make windows, consensus, BED→VCF |
| `cellranger` | 10x Genomics Cell Ranger: `count`, `vdj`, `multi` and result summaries |
| `cnv` | CNV scores and summaries: aneuploidy score, TMAD score |
| `cnvkit` | Each step of CNVkit as a process: access, autobin, coverage, reference, fix, segment, call, scatter, diagram, heatmap, batch, guess-baits |
| `delim` | Delimited-file utilities: bind rows, generate sample info files |
| `gene` | Gene name conversion, gene promoter regions |
| `gsea` | Gene set enrichment analysis: GSEA, fgsea, pre-ranked, Enrichr |
| `misc` | General utilities: run arbitrary shell commands, plots, file helpers |
| `pipseeker` | Pipseeker-based T-cell repertoire analysis |
| `plot` | Publication-ready plots: Venn diagram, heatmap, ROC, Manhattan, QQ, scatter, density |
| `protein` | Protein structure analysis: Prodigy, RMSD, MMCIF→PDB, PDB→FASTA |
| `regulatory` | Regulatory genomics: motif scanning, motif affinity test, variant-motif plots |
| `rnaseq` | RNA-seq utilities: unit conversion, transcriptome simulation |
| `scrna` | Single-cell RNA-seq (Seurat-centric): loading, QC/filtering, clustering, marker finding, cell-type annotation (SingleR, CelliD, SCINA, cellassign, scSorter, GPTCelltype, scBERT, scHDeepInsight, ...), trajectory (Slingshot), velocity (scVelo), cell-cell communication, demultiplexing (CellSNP-lite, Vireo), doublet detection (MQuad), imputation, reference mapping, 10x conversion |
| `scrna_metabolic_landscape` | Single-cell metabolic landscape: pathway activity, features, heterogeneity |
| `snp` | GWAS utilities on PLINK: simulation, QC (HWE, het, call rate), filtering, frequency, IBD, Matrix eQTL |
| `stats` | Statistical tests: Chow test, mediation analysis, liquid association, differential coexpression, meta-p-value |
| `tcgamaf` | TCGA MAF utilities: MAF→VCF, add chromosome prefix |
| `tcr` | T-cell receptor repertoire: immunarch loading/filtering, diversity, clonality, CDR3 clustering, TESSA, TCRdock, single-cell repertoire |
| `vcf` | VCF processing: liftover, filtering, indexing, splitting, intersecting, annotation, bcftools wrappers, truvari benchmarking |
| `web` | Data download: regular files and Google Cloud Storage buckets |

## Built-in pipelines

| Pipeline | Description |
| -------- | ----------- |
| [CNVkit pipeline](https://pwwang.github.io/biopipen/pipelines/cnvkit_pipeline/) | Full CNV detection workflow, decoupling `cnvkit.py batch` into 13 individually-controllable steps with flowcharts and interactive reports |
| [Cell Ranger pipelines](https://pwwang.github.io/biopipen/pipelines/cellranger_pipeline/) | `CellRangerCountPipeline`, `CellRangerVdjPipeline`, `CellRangerMultiPipeline` — run and summarize Cell Ranger across samples/GEM wells |
| [Pipseeker pipeline](https://pwwang.github.io/biopipen/pipelines/pipseeker_pipeline/) | Run `pipseeker full` per sample and summarize the results |
| Single-cell metabolic landscape | Metabolic pathway activity, features, and heterogeneity for single-cell data |

## Documentation

Full API reference, pipeline guides, and per-process configuration details are
available at [https://pwwang.github.io/biopipen](https://pwwang.github.io/biopipen).
