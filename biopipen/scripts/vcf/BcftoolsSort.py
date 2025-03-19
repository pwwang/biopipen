from typing import Literal
from pathlib import Path, PosixPath  # noqa: F401

from biopipen.utils.misc import run_command, logger
from biopipen.scripts.vcf.bcftools_utils import run_bcftools

infile: str = {{in.infile | quote}}  # pyright: ignore # noqa: E999
outfile: str = {{out.outfile | quote}}  # pyright: ignore
envs: dict = {{envs | dict | repr}}  # pyright: ignore

outdir = Path(outfile).parent
bcftools = envs.pop("bcftools")
tabix = envs.pop("tabix")
ncores = envs.pop("ncores")
gz = envs.pop("gz")
index = envs.pop("index")
chrsize = envs.pop("chrsize")
notfound = envs.pop("notfound")

if chrsize:
    class Contig:
        def __init__(self, name: str, length: str):
            self.name = name
            self.length = length

        def __str__(self) -> str:
            return f"##contig=<ID={self.name},length={self.length}>"

    def parse_header(header_file: Path) -> tuple[list[str], dict[str, Contig]]:
        hlines = []
        ctgs = {}
        with open(header_file) as fh:
            for line in fh:
                if line.startswith("##contig"):
                    ctg = line.strip().split("##contig=<ID=")[1].split(",length=")
                    ctgs[ctg[0]] = Contig(ctg[0], ctg[1].replace(">", ""))
                else:
                    hlines.append(line.strip())
        return hlines, ctgs

    def match_contigs(
        ctgs: dict[str, Contig],
        chroms: list[str],
        notfound: Literal["error", "remove", "start", "end"],
    ) -> list[str]:
        if (
            ctgs
            and chroms
            and all(chrom.startswith("chr") for chrom in chroms)
            and not any(chrom.startswith("chr") for chrom in ctgs)
        ):
            logger.warning(
                "Removing 'chr' prefix from chromosomes in envs.chrsize file, "
                "because the input VCF file does not have 'chr' prefix."
            )
            chroms = [chrom[3:] for chrom in chroms]

        new_ctgs = []
        for chrom in chroms:
            if chrom in ctgs:
                new_ctgs.append(str(ctgs[chrom]))
                del ctgs[chrom]

        if ctgs:
            if notfound == "error":
                raise ValueError(
                    "Chromosomes not found in envs.chrsize file: "
                    f"{', '.join(ctgs.keys())}"
                )
            elif notfound == "start":
                new_ctgs = [str(ctg) for ctg in ctgs.values()] + new_ctgs
            elif notfound == "end":
                new_ctgs = new_ctgs + [str(ctg) for ctg in ctgs.values()]

        return new_ctgs

    chroms = []
    with Path(chrsize).expanduser().open() as fh:
        for line in fh:
            chrom = line.strip().split()[0]
            chroms.append(chrom)

    header_file = outdir / "header.txt"
    run_command(f'{bcftools} view -h {infile} > {header_file}', fg=True)
    header_lines, contigs = parse_header(header_file)
    new_contigs = match_contigs(contigs, chroms, notfound=notfound)
    header_lines = [header_lines[0], *new_contigs, *header_lines[1:]]
    reheader_file = outdir / "reheader.txt"
    with open(reheader_file, "w") as fh:
        fh.writelines([f"{line}\n" for line in header_lines])

    reheader_vcf = outdir / f"{Path(infile).stem}_reheader.vcf"
    run_command([
        bcftools, "reheader",
        "--header", reheader_file,
        "-o", reheader_vcf,
        infile
    ], fg=True)

    infile = str(reheader_vcf)

envs[""] = [bcftools, "sort"]
envs["_"] = infile
envs["o"] = outfile

if index and not gz:
    logger.warning("Forcing envs.gz to True because envs.index is True.")
    gz = True

if "O" not in envs and "output-type" not in envs and "output_type" not in envs:
    envs["O"] = "z" if gz else "v"

run_bcftools(envs, bcftools=bcftools, index=index, tabix=tabix)
