import hashlib
from pathlib import Path

from biopipen.utils.misc import run_command, dict_to_cli_args

DESTDIR = Path(__file__).parent / "reference"
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
ARIA2C_OPTS = [
    "--file-allocation=falloc",
    "--auto-file-renaming=false",
    "--allow-overwrite=true",
]


def echo(msg):
    def _echo(func):
        def wrapper(*args, **kwargs):
            print(f"::group::{msg.format(*args, **kwargs)}")
            if func(*args, **kwargs):
                print("Cached")
            print("::endgroup::")
        return wrapper
    return _echo


def md5sum(file):
    if not Path(file).exists():
        return None
    with open(file, "rb") as f:
        return hashlib.md5(f.read()).hexdigest()


def download_reffa(genome):
    """Download genome reference sequences"""
    outdir = DESTDIR / genome
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


@echo("Downloading {0} chromosome sizes")
def download_chrsize(genome):
    """Download genome size file"""
    outdir = DESTDIR / genome
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


@echo("Downloading {0} refgenes")
def download_refgene(genome):
    """Download genome size file"""
    outdir = DESTDIR / genome
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


# @echo("Downloading KEGG_metabolism.gmt")
# def download_kegg_metabolism():
#     """Download KEGG_metabolism.gmt"""
#     outfile = DESTDIR / "KEGG_metabolism.gmt"
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


@echo("Downloading hg19 reference genome sequences")
def download_reffa_hg19():
    """Download hg19 reference sequences if not cached"""
    genome = "hg19"
    outdir = DESTDIR / genome
    outdir.mkdir(exist_ok=True, parents=True)
    reffa = outdir / "chrs.fa"
    if md5sum(reffa) == REFFA_HG19_MD5SUM:
        return True

    return download_reffa(genome)


@echo("Downloading hg38 reference genome sequences")
def download_reffa_hg38():
    """Download hg38 reference sequences if not cached"""
    genome = "hg38"
    outdir = DESTDIR / genome
    outdir.mkdir(exist_ok=True, parents=True)
    reffa = outdir / "chrs.fa"
    if md5sum(reffa) == REFFA_HG38_MD5SUM:
        return True

    return download_reffa(genome)


@echo("Downloading scType database")
def download_sctype_db():
    """Download scType database"""
    name = "ScTypeDB_full.xlsx"
    aria2c_args = dict(
        o=name,
        d=DESTDIR,
        _=SCTYPE_DB_URL,
    )
    aria2c_args[""] = ["aria2c", *ARIA2C_OPTS]
    run_command(dict_to_cli_args(aria2c_args, dashify=True), fg=True)


@echo("Downloading pbmc_multimodal_2023.rds")
def download_pbmc_multimodal():
    """Download pbmc_multimodal_2023.rds"""
    name = "pbmc_multimodal_2023.rds"
    aria2c_args = dict(
        o=name,
        d=DESTDIR,
        _=PBMC_MULTIMODEL_URL,
    )
    aria2c_args[""] = ["aria2c", *ARIA2C_OPTS]
    run_command(dict_to_cli_args(aria2c_args, dashify=True), fg=True)


if __name__ == "__main__":
    download_reffa_hg19()
    # download_reffa_hg38()
    download_chrsize("hg19")
    # download_chrsize("hg38")
    download_refgene("hg19")
    # download_refgene("hg38")
    # download_kegg_metabolism()
    download_sctype_db()
    download_pbmc_multimodal()
