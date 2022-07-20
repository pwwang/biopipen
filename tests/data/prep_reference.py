import cmdy
from pathlib import Path

DESTDIR = Path(__file__).parent / "reference"
REFFA_URL = (
    "http://hgdownload.cse.ucsc.edu/"
    "goldenpath/%(genome)s/bigZips/%(genome)s.fa.gz"
)
CHRSIZE_URL = (
    "http://hgdownload.cse.ucsc.edu/"
    "goldenpath/%(genome)s/bigZips/%(genome)s.chrom.sizes"
)
CHROMS = [
    "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9",
    "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17",
    "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY", "chrM",
]

def echo(msg):
    def _echo(func):
        def wrapper(*args, **kwargs):
            print(f"::group::{msg.format(*args, **kwargs)}")
            func(*args, **kwargs)
            print("::endgroup::")
        return wrapper
    return _echo


@echo("Downloading {0} reference genome sequences")
def download_reffa(genome):
    """Download genome reference sequences"""
    outdir = DESTDIR / genome
    outdir.mkdir(exist_ok=True, parents=True)
    url = REFFA_URL % {"genome": genome}
    outfile = outdir / "allchrs.fa.gz"
    reffa = outdir / "chrs.fa"
    cmdy.aria2c(
        s=2,
        x=2,
        o=outfile.name,
        d=outdir,
        _=url,
        **{"file-allocation": "falloc"},
    )
    cmdy.seqkit.grep(p=CHROMS, _=outfile, _dupkey=True).r() > reffa
    cmdy.samtools.faidx(reffa)
    cmdy.rm(f=True, _=outfile)


@echo("Downloading {0} chromosome sizes")
def download_chrsize(genome):
    """Download genome size file"""
    outdir = DESTDIR / genome
    outdir.mkdir(exist_ok=True, parents=True)
    url = CHRSIZE_URL % {"genome": genome}
    outfile = outdir / "chrom.sizes"
    cmdy.aria2c(
        s=2,
        x=2,
        o=outfile.name,
        d=outdir,
        _=url,
        **{"file-allocation": "falloc"},
    )


@echo("Downloading KEGG_metabolism.gmt")
def download_kegg_metabolism():
    """Download KEGG_metabolism.gmt"""
    outfile = DESTDIR / "KEGG_metabolism.gmt"
    url = (
        "https://raw.githubusercontent.com/"
        "LocasaleLab/Single-Cell-Metabolic-Landscape/"
        "master/Data/KEGG_metabolism.gmt"
    )
    cmdy.aria2c(
        s=2,
        x=2,
        o=outfile.name,
        d=outfile.parent,
        _=url,
        **{"file-allocation": "falloc"},
    )


if __name__ == "__main__":
    download_reffa("hg19")
    download_reffa("hg38")
    download_chrsize("hg19")
    download_chrsize("hg38")
    download_kegg_metabolism()
