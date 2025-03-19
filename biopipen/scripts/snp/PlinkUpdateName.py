from pathlib import Path
from biopipen.utils.misc import run_command, dict_to_cli_args, logger

indir: str = {{in.indir | quote}}  # pyright: ignore # noqa: #999
namefile: str = {{in.namefile | quote}}  # pyright: ignore
outdir: str = {{out.outdir | quote}}  # pyright: ignore
plink = {{envs.plink | repr}}  # pyright: ignore
bcftools = {{envs.bcftools | repr}}  # pyright: ignore
ncores = {{envs.ncores | repr}}  # pyright: ignore
match_alt = {{envs.match_alt | repr}}  # pyright: ignore

bedfile = list(Path(indir).glob("*.bed"))
if len(bedfile) == 0:
    raise FileNotFoundError(f"No .bed file found in `in.indir`")
elif len(bedfile) > 1:
    logger.warning(f"Multiple .bed files found in `in.indir`, using the first one.")

bedfile = bedfile[0]
input = bedfile.with_suffix("")
output = Path(outdir) / bedfile.stem

if namefile.endswith(".vcf") or namefile.endswith(".vcf.gz"):
    logger.info("VCF file received, extracting names")
    def alt_matched(bim_alt, vcf_alt, match_alt):
        if match_alt == "none":
            return True
        if match_alt == "exact":
            return bim_alt == vcf_alt

        bim_alts = bim_alt.split(",")
        vcf_alts = vcf_alt.split(",")
        if match_alt == "all":
            return set(bim_alts) == set(vcf_alts)
        if match_alt == "any":
            return bool(set(bim_alts) & set(vcf_alts))
        if match_alt == "first_included":
            return bim_alts[0] in vcf_alts
        if match_alt == "first":
            return bim_alts[0] == vcf_alts[0]

        raise ValueError(f"Unknown match_alt: {match_alt}")

    def readline(f):
        line = f.readline().strip()
        return line.split("\t") if line else None

    namefile_tmp = Path(outdir) / "_namefile_from_vcf.txt"
    infofile = Path(outdir) / "_information_from_vcf_unsorted.txt"
    sorted_infofile = Path(outdir) / "_information_from_vcf_sorted.txt"
    sorted_bim = Path(outdir) / "_sorted_bim.txt"
    bt_cmd = [
        bcftools, "query",
        "-f", "%CHROM\\t%ID\\t0\\t%POS\\t%ALT\\t%REF\\n",
        "-o", infofile,
        namefile,
    ]
    ## infofile
    # 1	rs10492	0	10492	T	C
    logger.info("- Extracting information from VCF file ...")
    run_command(bt_cmd, fg=True)
    # sort infofile
    logger.info("- Sorting the information from VCF file ...")
    run_command(
        [
            "sort",
            "-k1,1", "-k4,4n", "-k6,6",
            infofile,
            "--parallel", ncores,
            "-o", sorted_infofile
        ],
        env={"LC_ALL": "C"},
        fg=True,
    )

    ## .bim file
    # 1	1_10492	0	10492	T	C
    # sort .bim file
    logger.info("- Sorting the .bim file ...")
    run_command(
        [
            "sort",
            "-k1,1", "-k4,4n", "-k6,6",
            input.with_suffix(".bim"),
            "--parallel", ncores,
            "-o", sorted_bim
        ],
        env={"LC_ALL": "C"},
        fg=True,
    )
    # query namefile for records in sorted bim file
    logger.info("- Matching and generating the name file ...")
    with sorted_bim.open() as fbim, sorted_infofile.open() as finfo, namefile_tmp.open("w") as fout:  # noqa: E501
        bim = readline(fbim)
        info = readline(finfo)
        while bim and info:
            if (
                bim[0] == info[0]
                and bim[3] == info[3]
                and bim[5] == info[5]
                and alt_matched(bim[4], info[4], match_alt)
            ):
                fout.write(f"{bim[1]}\t{info[1]}\n")
                bim = readline(fbim)
                info = readline(finfo)
            elif (
                bim[0] < info[0]
                or (bim[0] == info[0] and bim[3] < info[3])
                or (bim[0] == info[0] and bim[3] == info[3] and bim[5] < info[5])
            ):
                bim = readline(fbim)
            else:
                info = readline(finfo)

    namefile = str(namefile_tmp)

args = {
    "": plink,
    "bfile": input,
    "out": output,
    "make_bed": True,
    "update_name": namefile,
}

run_command(dict_to_cli_args(args, dashify=True), fg=True)
