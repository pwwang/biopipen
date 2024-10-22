from __future__ import annotations

import sys
import hashlib
from pathlib import Path
from typing import Callable

from biopipen.utils.misc import run_command, dict_to_cli_args

LOCAL_FLAG = "--local"

if len(sys.argv) > 1 and sys.argv[1] != LOCAL_FLAG:
    print(f"Usage: python prep_data.py [{LOCAL_FLAG}]")
    sys.exit(1)

LOCAL = len(sys.argv) > 1 and sys.argv[1] == LOCAL_FLAG

REFDIR = Path(__file__).parent / "reference"
DATADIR = Path(__file__).parent / "data"
REFFA_URL = (
    "http://hgdownload.cse.ucsc.edu/"
    "goldenPath/%(genome)s/bigZips/%(genome)s.fa.gz"
)
CHRSIZE_URL = (
    "http://hgdownload.cse.ucsc.edu/"
    "goldenPath/%(genome)s/bigZips/%(genome)s.chrom.sizes"
)
REFGENE_URL = (
    "http://hgdownload.cse.ucsc.edu/"
    "goldenPath/%(genome)s/bigZips/genes/%(genome)s.refGene.gtf.gz"
)
SCTYPE_DB_URL = (
    "https://github.com/IanevskiAleksandr/sc-type/raw/master/ScTypeDB_full.xlsx"
)
PBMC_MULTIMODEL_URL = (
    "https://zenodo.org/records/7779017/files/pbmc_multimodal_2023.rds?download=1"
)
CHROMS = [
    "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9",
    "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17",
    "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY", "chrM",
]
REFFA_HG19_MD5SUM = "0f680f92acef0774052a4a5392817860"
REFFA_HG38_MD5SUM = "31da86e5ed328dec9fb528799f6658e6"
PBMC_MULTIMODEL_MD5SUM = "45f0cc1c4e977bb49698a697e9b71543"
ARIA2C_OPTS = [
    "--file-allocation=falloc",
    "--auto-file-renaming=false",
    "--allow-overwrite=true",
]
IFNB_SUB_URL = "https://github.com/zhanghao-njmu/SCP/raw/main/data/ifnb_sub.rda"


def decor(msg: str, local_only: bool | Callable = False):
    def _decor(func):
        def wrapper(*args, **kwargs):
            print(f"::group::{msg.format(*args, **kwargs)}")
            if callable(local_only):
                is_local_only = local_only(*args, **kwargs)
            else:
                is_local_only = local_only

            if is_local_only and not LOCAL:
                print("Skipping local-only task, as --local flag is not set.")
            elif func(*args, **kwargs):
                print("Cached")
            print("::endgroup::")
            print()
            print()
        return wrapper
    return _decor


def md5sum(file):
    if not Path(file).exists():
        return None
    with open(file, "rb") as f:
        return hashlib.md5(f.read()).hexdigest()


def download_reffa(genome):
    """Download genome reference sequences"""
    outdir = REFDIR / genome
    outdir.mkdir(exist_ok=True, parents=True)
    url = REFFA_URL % {"genome": genome}
    outfile = outdir / "allchrs.fa.gz"
    reffa = outdir / "chrs.fa"

    aria2c_args = dict(
        s=2,
        x=2,
        o=outfile.name,
        d=outdir,
        _=url,
    )
    aria2c_args[""] = ["aria2c", *ARIA2C_OPTS]
    run_command(dict_to_cli_args(aria2c_args, dashify=True), fg=True)

    seqkit_args = {"": ["seqkit", "grep"]}
    seqkit_args["p"] = CHROMS
    seqkit_args["_"] = outfile
    run_command(
        dict_to_cli_args(seqkit_args, dashify=True, dup_key=True),
        stdout=reffa,
    )
    run_command(["samtools", "faidx", reffa])
    run_command(["rm", "-f", outfile])


@decor("Downloading {0} chromosome sizes", local_only=lambda genome: genome == "hg38")
def download_chrsize(genome):
    """Download genome size file"""
    outdir = REFDIR / genome
    outdir.mkdir(exist_ok=True, parents=True)
    url = CHRSIZE_URL % {"genome": genome}
    outfile = outdir / "chrom.sizes"

    aria2c_args = dict(
        s=2,
        x=2,
        o=outfile.name,
        d=outdir,
        _=url,
    )
    aria2c_args[""] = ["aria2c", *ARIA2C_OPTS]
    run_command(dict_to_cli_args(aria2c_args, dashify=True), fg=True)


@decor("Downloading {0} refgenes", local_only=lambda genome: genome == "hg38")
def download_refgene(genome):
    """Download genome size file"""
    outdir = REFDIR / genome
    outdir.mkdir(exist_ok=True, parents=True)
    url = REFGENE_URL % {"genome": genome}
    outfile = outdir / "refgene.gtf.gz"

    aria2c_args = dict(
        s=2,
        x=2,
        o=outfile.name,
        d=outdir,
        _=url,
    )
    aria2c_args[""] = ["aria2c", *ARIA2C_OPTS]
    run_command(dict_to_cli_args(aria2c_args, dashify=True), fg=True)

    refgene_file = outdir / "refgene.gtf"
    refexon_file = outdir / "refexon.gtf"
    run_command(["gunzip", "-f", outfile])
    run_command(["awk", '$3 == "exon"', refgene_file], stdout=refexon_file)


# @decor("Downloading KEGG_metabolism.gmt")
# def download_kegg_metabolism():
#     """Download KEGG_metabolism.gmt"""
#     outfile = REFDIR / "KEGG_metabolism.gmt"
#     url = (
#         "https://raw.githubusercontent.com/"
#         "LocasaleLab/Single-Cell-Metabolic-Landscape/"
#         "master/Data/KEGG_metabolism.gmt"
#     )
#     aria2c_args = dict(
#         s=2,
#         x=2,
#         o=outfile.name,
#         d=outfile.parent,
#         _=url,
#     )
#     aria2c_args[""] = ["aria2c", *ARIA2C_OPTS]
#     run_command(dict_to_cli_args(aria2c_args, dashify=True), fg=True)


@decor("Downloading hg19 reference genome sequences")
def download_reffa_hg19():
    """Download hg19 reference sequences if not cached"""
    genome = "hg19"
    outdir = REFDIR / genome
    outdir.mkdir(exist_ok=True, parents=True)
    reffa = outdir / "chrs.fa"
    if md5sum(reffa) == REFFA_HG19_MD5SUM:
        return True

    return download_reffa(genome)


@decor("Downloading hg38 reference genome sequences", local_only=True)
def download_reffa_hg38():
    """Download hg38 reference sequences if not cached"""
    genome = "hg38"
    outdir = REFDIR / genome
    outdir.mkdir(exist_ok=True, parents=True)
    reffa = outdir / "chrs.fa"
    if md5sum(reffa) == REFFA_HG38_MD5SUM:
        return True

    return download_reffa(genome)


@decor("Downloading scType database")
def download_sctype_db():
    """Download scType database"""
    name = "ScTypeDB_full.xlsx"
    aria2c_args = dict(
        o=name,
        d=REFDIR,
        _=SCTYPE_DB_URL,
    )
    aria2c_args[""] = ["aria2c", *ARIA2C_OPTS]
    run_command(dict_to_cli_args(aria2c_args, dashify=True), fg=True)


@decor("Downloading ifnb_sub.rda and convert to .rds")
def download_ifnb_sub():
    """Download ifnb_sub.rda"""
    name = "ifnb_sub.rda"
    destname = "ifnb_sub.rds"
    file = DATADIR / name
    destfile = DATADIR / destname
    if file.is_file() and destfile.is_file():
        return True

    aria2c_args = dict(
        o=name,
        d=DATADIR,
        _=IFNB_SUB_URL,
    )
    aria2c_args[""] = ["aria2c", *ARIA2C_OPTS]
    run_command(dict_to_cli_args(aria2c_args, dashify=True), fg=True)
    run_command(
        [
            "Rscript",
            "-e",
            f"""
                library(Seurat)
                load('{file}');
                ifnb_sub <- NormalizeData(ifnb_sub)
                ifnb_sub <- FindVariableFeatures(ifnb_sub)
                ifnb_sub <- ScaleData(ifnb_sub)
                ifnb_sub <- RunPCA(ifnb_sub)
                ifnb_sub <- FindNeighbors(ifnb_sub)
                ifnb_sub <- FindClusters(
                    ifnb_sub,
                    resolution = c(setdiff(seq(0.2, 1.2, 0.2), 0.4), 0.4)
                )
                ifnb_sub <- RunUMAP(ifnb_sub, dims = 1:30)
                saveRDS(ifnb_sub, '{destfile}')
            """,
        ],
        fg=True,
    )


@decor("Downloading pbmc_multimodal_2023.rds")
def download_pbmc_multimodal():
    """Download pbmc_multimodal_2023.rds"""
    name = "pbmc_multimodal_2023.rds"
    destfile = REFDIR / name
    if md5sum(destfile) == PBMC_MULTIMODEL_MD5SUM:
        return True

    aria2c_args = dict(
        o=name,
        d=REFDIR,
        _=PBMC_MULTIMODEL_URL,
    )
    aria2c_args[""] = ["aria2c", *ARIA2C_OPTS]
    run_command(dict_to_cli_args(aria2c_args, dashify=True), fg=True)


if __name__ == "__main__":
    download_reffa_hg19()
    download_reffa_hg38()
    download_chrsize("hg19")
    download_chrsize("hg38")
    download_refgene("hg19")
    download_refgene("hg38")
    # download_kegg_metabolism()
    download_sctype_db()
    download_ifnb_sub()
    download_pbmc_multimodal()
